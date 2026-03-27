#' Bayesian Estimation and Variable Selection for Group Differences in Markov Random Fields
#'
#' @description
#' The \code{bgmCompare} function estimates group differences in category
#' threshold parameters (main effects) and pairwise interactions (pairwise
#' effects) of a Markov Random Field (MRF) for binary and ordinal variables.
#' Groups can be defined either by supplying two separate datasets (\code{x} and
#' \code{y}) or by a group membership vector. Optionally, Bayesian variable
#' selection can be applied to identify differences across groups.
#'
#' @details
#' This function extends the ordinal MRF framework
#' \insertCite{MarsmanVandenBerghHaslbeck_2025;textual}{bgms} to multiple
#' groups. The basic idea of modeling, analyzing, and testing group
#' differences in MRFs was introduced in
#' \insertCite{MarsmanWaldorpSekulovskiHaslbeck_2024;textual}{bgms}, where
#' two–group comparisons were conducted using adaptive Metropolis sampling.
#' The present implementation generalizes that approach to more than two
#' groups and supports additional samplers (NUTS) with staged warmup
#' adaptation.
#'
#' Key components of the model:
#'
#' @seealso \code{vignette("comparison", package = "bgms")} for a worked example.
#' @family model-fitting
#'
#' @section Pairwise Interactions:
#' For variables \eqn{i} and \eqn{j}, the group-specific interaction is
#' represented as:
#' \deqn{\theta_{ij}^{(g)} = \phi_{ij} + \delta_{ij}^{(g)},}
#' where \eqn{\phi_{ij}} is the baseline effect and
#' \eqn{\delta_{ij}^{(g)}} are group differences constrained to sum to zero.
#'
#' @section Ordinal Variables:
#' \strong{Regular ordinal variables}: category thresholds are decomposed into a
#' baseline plus group differences for each category.
#'
#' \strong{Blume–Capel variables}: category thresholds are quadratic in the
#' category index, with both the linear and quadratic terms split into a
#' baseline plus group differences.
#'
#' @section Variable Selection:
#' When \code{difference_selection = TRUE}, spike-and-slab priors are
#' applied to difference parameters:
#' \itemize{
#'   \item \strong{Bernoulli}: fixed prior inclusion probability.
#'   \item \strong{Beta–Bernoulli}: inclusion probability given a Beta prior.
#' }
#'
#' @section Sampling Algorithms and Warmup:
#' Parameters are updated within a Gibbs framework, using the same
#' sampling algorithms and staged warmup scheme described in
#' \code{\link{bgm}}:
#' \itemize{
#'   \item \strong{Adaptive Metropolis–Hastings}: componentwise random–walk
#'     proposals with Robbins–Monro adaptation of proposal SDs.
#'   \item \strong{Hamiltonian Monte Carlo (HMC)} (\emph{deprecated}): joint
#'     updates with fixed leapfrog trajectories. This method is deprecated;
#'     use NUTS instead.
#'   \item \strong{No–U–Turn Sampler (NUTS)}: an adaptive HMC variant with
#'     dynamic trajectory lengths; warmup uses a staged adaptation schedule.
#' }
#'
#' For details on the staged adaptation schedule (fast–slow–fast phases),
#' see \code{\link{bgm}}. In addition, when
#' \code{difference_selection = TRUE}, updates of inclusion indicators are
#' delayed until late warmup. In NUTS, this appends two extra phases
#' (Stage-3b and Stage-3c), so that the total number of warmup iterations
#' exceeds the user-specified \code{warmup}.
#'
#' After warmup, adaptation is disabled: step size and mass matrix are fixed
#' at their learned values, and proposal SDs remain constant.
#'
#' @param x A data frame or matrix of binary and ordinal responses for
#'   Group 1. Variables should be coded as nonnegative integers starting at
#'   0. For ordinal variables, unused categories are collapsed; for
#'   Blume–Capel variables, all categories are retained.
#' @param y Optional data frame or matrix for Group 2 (two-group designs).
#'   Must have the same variables (columns) as \code{x}.
#' @param group_indicator Optional integer vector of group memberships for
#'   rows of \code{x} (multi-group designs). Ignored if \code{y} is supplied.
#' @param difference_selection Logical. If \code{TRUE}, spike-and-slab priors
#'   are applied to difference parameters. Default: \code{TRUE}.
#' @param main_difference_selection Logical. If \code{TRUE}, apply spike-and-slab
#'   selection to main effect (threshold) differences. If \code{FALSE}, main
#'   effect differences are always included (no selection). Since main effects
#'   are often nuisance parameters and their selection can interfere with
#'   pairwise selection under the Beta-Bernoulli prior, the default is
#'   \code{FALSE}. Only used when \code{difference_selection = TRUE}.
#' @param variable_type Character vector specifying type of each variable:
#'   \code{"ordinal"} (default) or \code{"blume-capel"}.
#' @param baseline_category Integer or vector giving the baseline category
#'   for Blume–Capel variables.
#' @param difference_scale Double. Scale of the Cauchy prior for difference
#'   parameters. Default: \code{1}.
#' @param difference_prior Character. Prior for difference inclusion:
#'   \code{"Bernoulli"} or \code{"Beta-Bernoulli"}. Default: \code{"Bernoulli"}.
#' @param difference_probability Numeric. Prior inclusion probability for
#'   differences (Bernoulli prior). Default: \code{0.5}.
#' @param beta_bernoulli_alpha,beta_bernoulli_beta Doubles. Shape parameters
#'   of the Beta prior for inclusion probabilities in the Beta–Bernoulli
#'   model. Defaults: \code{1}.
#' @param pairwise_scale Double. Scale of the Cauchy prior for baseline
#'   pairwise interactions. Default: \code{1}.
#' @param standardize Logical. If \code{TRUE}, the Cauchy prior scale for each
#'   pairwise interaction (both baseline and difference) is adjusted based on
#'   the range of response scores. Without standardization, pairs with more
#'   response categories experience less shrinkage because their naturally
#'   smaller interaction effects make a fixed prior relatively wide.
#'   Standardization equalizes relative shrinkage across all pairs, with
#'   \code{pairwise_scale} itself applying to the unit interval (binary) case.
#'   See \code{\link{bgm}} for details on the adjustment. Default: \code{FALSE}.
#' @param main_alpha,main_beta Doubles. Shape parameters of the beta-prime
#'   prior for baseline threshold parameters. Defaults: \code{0.5}.
#' @param iter Integer. Number of post–warmup iterations per chain.
#'   Default: \code{1e3}.
#' @param warmup Integer. Number of warmup iterations before sampling.
#'   Default: \code{1e3}.
#' @param na_action Character. How to handle missing data:
#'   \code{"listwise"} (drop rows) or \code{"impute"} (impute within Gibbs).
#'   Default: \code{"listwise"}.
#' @param display_progress Character. Controls progress reporting:
#'   \code{"per-chain"}, \code{"total"}, or \code{"none"}.
#'   Default: \code{"per-chain"}.
#' @param progress_callback An optional R function with signature
#'   \code{function(completed, total)} that is called at regular intervals
#'   during sampling, where \code{completed} is the number of iterations
#'   completed across all chains and \code{total} is the total number of
#'   iterations. Useful for external front-ends (e.g., JASP) that supply
#'   their own progress reporting.
#'   When \code{NULL} (the default), no callback is invoked.
#' @param verbose Logical. If \code{TRUE}, prints informational messages
#'   during data processing (e.g., missing data handling, variable recoding).
#'   Defaults to \code{getOption("bgms.verbose", TRUE)}. Set
#'   \code{options(bgms.verbose = FALSE)} to suppress messages globally.
#' @param update_method Character. Sampling algorithm:
#'   \code{"adaptive-metropolis"} or \code{"nuts"}.
#'   \code{"hamiltonian-mc"} is accepted but deprecated; use \code{"nuts"}
#'   instead. Default: \code{"nuts"}.
#' @param target_accept Numeric between 0 and 1. Target acceptance rate.
#'   Defaults: 0.44 (Metropolis), 0.80 (NUTS).
#' @param hmc_num_leapfrogs `r lifecycle::badge("deprecated")` Integer.
#'   Leapfrog steps for HMC (deprecated). Default: \code{100}.
#' @param nuts_max_depth Integer. Maximum tree depth for NUTS. Default: \code{10}.
#' @param learn_mass_matrix Logical. If \code{TRUE}, adapts a diagonal mass
#' matrix during warmup (NUTS only). Default: \code{TRUE}.
#' @param chains Integer. Number of parallel chains. Default: \code{4}.
#' @param cores Integer. Number of CPU cores. Default:
#'   \code{parallel::detectCores()}.
#' @param seed Optional integer. Random seed for reproducibility.
#' @param main_difference_model,reference_category,pairwise_difference_scale,main_difference_scale,pairwise_difference_prior,main_difference_prior,pairwise_difference_probability,main_difference_probability,pairwise_beta_bernoulli_alpha,pairwise_beta_bernoulli_beta,main_beta_bernoulli_alpha,main_beta_bernoulli_beta,interaction_scale,threshold_alpha,threshold_beta,burnin,save
#'   `r lifecycle::badge("deprecated")`
#'   Deprecated arguments as of \strong{bgms 0.1.6.0}.
#'   Use `difference_scale`, `difference_prior`, `difference_probability`,
#'   `beta_bernoulli_alpha`, `beta_bernoulli_beta`, `baseline_category`,
#'   `pairwise_scale`, and `warmup` instead.
#' @return
#' A list of class \code{"bgmCompare"} containing posterior summaries,
#' posterior mean matrices, and raw MCMC samples:
#' \itemize{
#'   \item \code{posterior_summary_main_baseline},
#'     \code{posterior_summary_pairwise_baseline}: summaries of baseline
#'     thresholds and pairwise interactions.
#'   \item \code{posterior_summary_main_differences},
#'     \code{posterior_summary_pairwise_differences}: summaries of group
#'     differences in thresholds and pairwise interactions.
#'   \item \code{posterior_summary_indicator}: summaries of inclusion
#'     indicators (if \code{difference_selection = TRUE}).
#'   \item \code{posterior_mean_main_baseline},
#'     \code{posterior_mean_associations_baseline}: posterior mean matrices
#'     (legacy style).
#'   \item \code{raw_samples}: list of raw draws per chain for main,
#'     pairwise, and indicator parameters.
#'   \item \code{arguments}: list of function call arguments and metadata.
#' }
#'
#' The \code{summary()} method prints formatted summaries, and
#' \code{coef()} extracts posterior means.
#'
#' NUTS diagnostics (tree depth, divergences, energy, E-BFMI) are included
#' in \code{fit$nuts_diag} if \code{update_method = "nuts"}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontrun{
#' # Run bgmCompare on subset of the Boredom dataset
#' x = Boredom[Boredom$language == "fr", 2:6]
#' y = Boredom[Boredom$language != "fr", 2:6]
#'
#' fit = bgmCompare(x, y, chains = 2)
#'
#' # Posterior inclusion probabilities
#' summary(fit)$indicator
#'
#' # Bayesian model averaged main effects for the groups
#' coef(fit)$main_effects_groups
#'
#' # Bayesian model averaged pairwise effects for the groups
#' coef(fit)$pairwise_effects_groups
#' }
#'
#' @export
bgmCompare = function(
  x,
  y,
  group_indicator,
  difference_selection = TRUE,
  main_difference_selection = FALSE,
  variable_type = "ordinal",
  baseline_category,
  difference_scale = 1,
  difference_prior = c("Bernoulli", "Beta-Bernoulli"),
  difference_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  pairwise_scale = 1,
  main_alpha = 0.5,
  main_beta = 0.5,
  iter = 1e3,
  warmup = 1e3,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis", "hamiltonian-mc"),
  target_accept,
  hmc_num_leapfrogs = 100,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  standardize = FALSE,
  verbose = getOption("bgms.verbose", TRUE),
  progress_callback = NULL,
  main_difference_model,
  reference_category,
  main_difference_scale,
  pairwise_difference_scale,
  pairwise_difference_prior,
  main_difference_prior,
  pairwise_difference_probability,
  main_difference_probability,
  pairwise_beta_bernoulli_alpha,
  pairwise_beta_bernoulli_beta,
  main_beta_bernoulli_alpha,
  main_beta_bernoulli_beta,
  interaction_scale,
  threshold_alpha,
  threshold_beta,
  burnin,
  save
) {
  # Set verbose option for internal functions, restore on exit
  old_verbose = getOption("bgms.verbose")
  options(bgms.verbose = verbose)
  on.exit(options(bgms.verbose = old_verbose), add = TRUE)

  if(hasArg(main_difference_model)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(main_difference_model =)")
  }

  if(hasArg(reference_category)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(reference_category =)", "bgmCompare(baseline_category =)")
    if(!hasArg(baseline_category)) baseline_category = reference_category
  }

  if(hasArg(pairwise_difference_scale) || hasArg(main_difference_scale)) {
    lifecycle::deprecate_warn(
      "0.1.6.0", "bgmCompare(pairwise_difference_scale =, main_difference_scale =)",
      "bgmCompare(difference_scale =)"
    )
    if(!hasArg(difference_scale)) {
      difference_scale = if(!missing(pairwise_difference_scale)) pairwise_difference_scale else main_difference_scale
    }
  }

  if(hasArg(pairwise_difference_prior) || hasArg(main_difference_prior)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_difference_prior =, main_difference_prior =)",
      "bgmCompare(difference_prior =)"
    )
    if(!hasArg(difference_prior)) {
      difference_prior = if(!missing(pairwise_difference_prior)) pairwise_difference_prior else main_difference_prior
    }
  }

  if(hasArg(pairwise_difference_probability) || hasArg(main_difference_probability)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_difference_probability =, main_difference_probability =)",
      "bgmCompare(difference_probability =)"
    )
    if(!hasArg(difference_probability)) {
      difference_probability = if(!missing(pairwise_difference_probability)) pairwise_difference_probability else main_difference_probability
    }
  }

  if(hasArg(pairwise_beta_bernoulli_alpha) || hasArg(main_beta_bernoulli_alpha)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_beta_bernoulli_alpha =, main_beta_bernoulli_alpha =)",
      "bgmCompare(beta_bernoulli_alpha =)"
    )
    if(!hasArg(beta_bernoulli_alpha)) {
      beta_bernoulli_alpha = if(!missing(pairwise_beta_bernoulli_alpha)) pairwise_beta_bernoulli_alpha else main_beta_bernoulli_alpha
    }
  }

  if(hasArg(pairwise_beta_bernoulli_beta) || hasArg(main_beta_bernoulli_beta)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(pairwise_beta_bernoulli_beta =, main_beta_bernoulli_beta =)",
      "bgmCompare(beta_bernoulli_beta =)"
    )
    if(!hasArg(beta_bernoulli_beta)) {
      beta_bernoulli_beta = if(!missing(pairwise_beta_bernoulli_beta)) pairwise_beta_bernoulli_beta else main_beta_bernoulli_beta
    }
  }

  if(hasArg(interaction_scale)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(interaction_scale =)", "bgmCompare(pairwise_scale =)")
    if(!hasArg(pairwise_scale)) pairwise_scale = interaction_scale
  }

  if(hasArg(threshold_alpha) || hasArg(threshold_beta)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgmCompare(threshold_alpha =, threshold_beta =)",
      "bgmCompare(main_alpha =, main_beta =)" # = double-check if these are still part of bgmCompare
    )
    if(!hasArg(main_alpha)) main_alpha = threshold_alpha
    if(!hasArg(main_beta)) main_beta = threshold_beta
  }

  if(hasArg(burnin)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(burnin =)", "bgmCompare(warmup =)")
    if(!hasArg(warmup)) warmup = burnin
  }

  if(hasArg(save)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgmCompare(save =)")
  }

  # --- Build spec, sample, build output ----------------------------------------
  spec = bgm_spec(
    x = x,
    model_type = "compare",
    variable_type = variable_type,
    baseline_category = if(hasArg(baseline_category)) baseline_category else 0L,
    y = if(hasArg(y)) y else NULL,
    group_indicator = if(hasArg(group_indicator)) group_indicator else NULL,
    na_action = na_action,
    pairwise_scale = pairwise_scale,
    main_alpha = main_alpha,
    main_beta = main_beta,
    standardize = standardize,
    difference_selection = difference_selection,
    main_difference_selection = main_difference_selection,
    difference_prior = difference_prior,
    difference_scale = difference_scale,
    difference_probability = difference_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    update_method = update_method,
    target_accept = if(hasArg(target_accept)) target_accept else NULL,
    iter = iter,
    warmup = warmup,
    hmc_num_leapfrogs = hmc_num_leapfrogs,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    chains = chains,
    cores = cores,
    seed = seed,
    display_progress = display_progress,
    verbose = verbose,
    progress_callback = progress_callback
  )

  raw = run_sampler(spec)
  output = build_output(spec, raw)

  output$.bgm_spec = spec
  return(output)
}
