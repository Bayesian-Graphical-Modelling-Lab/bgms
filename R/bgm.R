#' Bayesian Estimation or Edge Selection for Markov Random Fields
#'
#' @description
#' The \code{bgm} function estimates the pseudoposterior distribution of the
#' parameters of a Markov Random Field (MRF) for binary, ordinal, continuous,
#' or mixed (discrete and continuous) variables. Depending on the variable
#' types, the model is an ordinal MRF, a Gaussian graphical model (GGM), or a
#' mixed MRF. Optionally, it performs Bayesian edge selection using
#' spike-and-slab priors to infer the network structure.
#'
# ’ @details
# ’ Depending on the variable types, the model is an ordinal MRF, a Gaussian
# ’ graphical model (GGM), or a mixed MRF. Both regular ordinal variables and
# ’ Blume--Capel ordinal variables (with a baseline category) are supported.
# ’
# ’ Edge selection uses spike-and-slab priors with Bernoulli, Beta-Bernoulli,
# ’ or Stochastic-Block inclusion priors. Parameters are sampled with NUTS
# ’ (default) or adaptive Metropolis--Hastings, with a multi-stage warmup
# ’ schedule. Missing data can be handled via listwise deletion or Gibbs
# ’ imputation.
# ’
# ’ For full details on model specification, prior choices, warmup, and
# ’ output interpretation, see the package website at
# ’ \url{https://bayesian-graphical-modelling-lab.github.io/bgms-docs/}.
# ’
# ’ @seealso \code{vignette(“intro”, package = “bgms”)} for a worked example.
# ’ @family model-fitting
#'
#' @param x A data frame or matrix with \code{n} rows and \code{p} columns.
#'   Columns may contain binary, ordinal, or continuous variables (see
#'   \code{variable_type}). Discrete variables are automatically recoded to
#'   non-negative integers (\code{0, 1, ..., m}); for regular ordinal
#'   variables, unobserved categories are collapsed, while Blume–Capel
#'   variables retain all categories. Continuous variables are column-centered
#'   internally so that the GGM likelihood is formulated with a zero-mean
#'   assumption.
#'
#' @param variable_type Character or character vector. Specifies the type of
#'   each variable in \code{x}. Allowed values: \code{"ordinal"},
#'   \code{"blume-capel"}, or \code{"continuous"}. A single string applies
#'   to all variables. A per-variable vector that mixes discrete
#'   (\code{"ordinal"} / \code{"blume-capel"}) and \code{"continuous"}
#'   types fits a mixed MRF. Binary variables are automatically treated as
#'   \code{"ordinal"}. Default: \code{"ordinal"}.
#'
#' @param baseline_category Integer or vector. Baseline category used in
#'   Blume–Capel variables. Can be a single integer (applied to all) or a
#'   vector of length \code{p}. Required if at least one variable is of type
#'   \code{"blume-capel"}.
#'
#' @param iter Integer. Number of post–burn-in iterations (per chain).
#'   Default: \code{1e3}.
#'
#' @param warmup Integer. Number of warmup iterations before collecting
#'   samples. A minimum of 1000 iterations is enforced, with a warning if a
#'   smaller value is requested. Default: \code{1e3}.
#'
#' @param pairwise_scale Double. Scale of the Cauchy prior for pairwise
#'   interaction parameters. Default: \code{1}.
#'
#' @param standardize Logical. If \code{TRUE}, the Cauchy prior scale for each
#'   pairwise interaction is adjusted based on the range of response scores.
#'   Variables with more response categories have larger score products
#'   \eqn{x_i \cdot x_j}, which typically correspond to smaller interaction
#'   effects \eqn{\sigma_{ij}}. Without standardization, a fixed prior scale
#'   is relatively wide for these smaller effects, resulting in less shrinkage
#'   for high-category pairs and more shrinkage for low-category pairs.
#'   Standardization scales the prior proportionally to the maximum score
#'   product, ensuring equivalent relative shrinkage across all pairs.
#'   After internal recoding, regular ordinal variables have scores
#'   \eqn{0, 1, \ldots, m}. The adjusted scale for the interaction between
#'   variables \eqn{i} and \eqn{j} is \code{pairwise_scale * m_i * m_j},
#'   so that \code{pairwise_scale} itself applies to the unit interval case
#'   (binary variables where \eqn{m_i = m_j = 1}). For Blume-Capel variables
#'   with reference category \eqn{b}, scores are centered as
#'   \eqn{-b, \ldots, m-b}, and the adjustment uses the maximum absolute
#'   product of the score endpoints. For mixed pairs, ordinal variables use
#'   raw score endpoints \eqn{(0, m)} and Blume-Capel variables use centered
#'   score endpoints \eqn{(-b, m-b)}.
#'   Default: \code{FALSE}.
#'
#' @param pseudolikelihood Character. Specifies the pseudo-likelihood
#'   approximation used for mixed MRF models (ignored for pure ordinal or
#'   pure continuous data). Options:
#'   \describe{
#'     \item{\code{"conditional"}}{Conditions on the observed continuous
#'       variables when computing the discrete full conditionals. Faster
#'       because the discrete pseudo-likelihood does not depend on the
#'       continuous precision matrix.}
#'     \item{\code{"marginal"}}{Integrates out the continuous variables,
#'       giving discrete full conditionals that account for induced
#'       interactions through the continuous block. More expensive per
#'       iteration.}
#'   }
#'   Default: \code{"conditional"}.
#'
#' @param main_alpha,main_beta Double. Shape parameters of the
#'   beta-prime prior for threshold parameters. Must be positive. If equal,
#'   the prior is symmetric. Defaults: \code{main_alpha = 0.5} and
#'   \code{main_beta = 0.5}.
#'
#' @param edge_selection Logical. Whether to perform Bayesian edge selection.
#'   If \code{FALSE}, the model estimates all edges. Default: \code{TRUE}.
#'
#' @param edge_prior Character. Specifies the prior for edge inclusion.
#'   Options: \code{"Bernoulli"}, \code{"Beta-Bernoulli"}, or
#'   \code{"Stochastic-Block"}. Default: \code{"Bernoulli"}.
#'
#' @param inclusion_probability Numeric scalar. Prior inclusion probability
#'   of each edge (used with the Bernoulli prior). Default: \code{0.5}.
#'
#' @param beta_bernoulli_alpha,beta_bernoulli_beta Double. Shape parameters
#'   for the beta distribution in the Beta–Bernoulli and the Stochastic-Block
#'   priors. Must be positive. For the Stochastic-Block prior these are the shape
#'   parameters for the within-cluster edge inclusion probabilities.
#'   Defaults: \code{beta_bernoulli_alpha = 1} and \code{beta_bernoulli_beta = 1}.
#'
#' @param beta_bernoulli_alpha_between,beta_bernoulli_beta_between Double.
#' Shape parameters for the between-cluster edge inclusion probabilities in the
#' Stochastic-Block prior. Must be positive.
#' Default: \code{beta_bernoulli_alpha_between = 1} and \code{beta_bernoulli_beta_between = 1}
#'
#' @param dirichlet_alpha Double. Concentration parameter of the Dirichlet
#'   prior on block assignments (used with the Stochastic Block model).
#'   Default: \code{1}.
#'
#' @param lambda Double. Rate of the zero-truncated Poisson prior on the
#'   number of clusters in the Stochastic Block Model. Default: \code{1}.
#'
#' @param na_action Character. Specifies missing data handling. Either
#'   \code{"listwise"} (drop rows with missing values) or \code{"impute"}
#'   (perform single imputation during sampling). Default: \code{"listwise"}.
#'
#' @param display_progress Character. Controls progress reporting during
#'   sampling. Options: \code{"per-chain"} (separate bar per chain),
#'   \code{"total"} (single combined bar), or \code{"none"} (no progress).
#'   Default: \code{"per-chain"}.
#'
#' @param progress_callback An optional R function with signature
#'   \code{function(completed, total)} that is called at regular intervals
#'   during sampling, where \code{completed} is the number of iterations
#'   completed across all chains and \code{total} is the total number of
#'   iterations. Useful for external front-ends (e.g., JASP) that supply
#'   their own progress reporting.
#'   When \code{NULL} (the default), no callback is invoked.
#'
#' @param verbose Logical. If \code{TRUE}, prints informational messages
#'   during data processing (e.g., missing data handling, variable recoding).
#'   Defaults to \code{getOption("bgms.verbose", TRUE)}. Set
#'   \code{options(bgms.verbose = FALSE)} to suppress messages globally.
#'
#' @param update_method Character. Specifies how the MCMC sampler updates
#'   the model parameters:
#'   \describe{
#'     \item{"adaptive-metropolis"}{Componentwise adaptive Metropolis–Hastings
#'       with Robbins–Monro proposal adaptation.}
#'     \item{"nuts"}{The No-U-Turn Sampler with RATTLE constrained integration
#'       for Gaussian models with edge selection.}
#'   }
#'   Default: \code{"nuts"}.
#'
#' @param target_accept Numeric between 0 and 1. Target acceptance rate for
#'   the sampler. Defaults are set automatically if not supplied:
#'   \code{0.44} for adaptive Metropolis and \code{0.80} for NUTS.
#'
#' @param nuts_max_depth Integer. Maximum tree depth in NUTS. Must be positive.
#'   Default: \code{10}.
#'
#' @param learn_mass_matrix Logical. If \code{TRUE}, adapt a diagonal mass
#'   matrix during warmup (NUTS only). If \code{FALSE}, use the identity
#'   matrix. Default: \code{TRUE}.
#'
#' @param chains Integer. Number of parallel chains to run. Default: \code{4}.
#'
#' @param cores Integer. Number of CPU cores for parallel execution.
#'   Default: \code{parallel::detectCores()}.
#'
#' @param seed Optional integer. Random seed for reproducibility. Must be a
#'   single non-negative integer.
#'
#' @param interaction_scale,burnin,save,threshold_alpha,threshold_beta
#'   `r lifecycle::badge("deprecated")`
#'   Deprecated arguments as of \strong{bgms 0.1.6.0}.
#'   Use `pairwise_scale`, `warmup`, `main_alpha`, and `main_beta` instead.
#'
#' @return
#' A list of class \code{"bgms"} with posterior summaries, posterior mean
#' matrices, and access to raw MCMC draws. The object can be passed to
#' \code{print()}, \code{summary()}, and \code{coef()}.
#'
#' Main components include:
#' \itemize{
#'   \item \code{posterior_summary_main}: Data frame with posterior summaries
#'     (mean, sd, MCSE, ESS, Rhat) for main-effect parameters.
#'     For OMRF models these are category thresholds;
#'     for mixed MRF models these are discrete thresholds and
#'     continuous means. \code{NULL} for GGM models (no main effects).
#'   \item \code{posterior_summary_quadratic}: Data frame with posterior
#'     summaries for the residual variance parameters (GGM and mixed MRF).
#'     \code{NULL} for OMRF models.
#'   \item \code{posterior_summary_pairwise}: Data frame with posterior
#'     summaries for partial association parameters.
#'   \item \code{posterior_summary_indicator}: Data frame with posterior
#'     summaries for edge inclusion indicators (if \code{edge_selection = TRUE}).
#'
#'   \item \code{posterior_mean_main}: Posterior mean of main-effect
#'     parameters. \code{NULL} for GGM models. For OMRF: a matrix
#'     (p x max_categories) of category thresholds. For mixed MRF: a list
#'     with \code{$discrete} (threshold matrix) and \code{$continuous}
#'     (q x 1 matrix of means).
#'   \item \code{posterior_mean_pairwise}: Symmetric matrix of posterior
#'     mean partial associations (zero diagonal). For continuous variables
#'     these are unstandardized partial correlations; for discrete variables
#'     these are half the log adjacent-category odds ratio. Use
#'     [extract_precision()], [extract_partial_correlations()], or
#'     [extract_log_odds()] to convert to interpretable scales.
#'   \item \code{posterior_mean_residual_variance}: Named numeric vector of
#'     posterior mean residual variances \eqn{1/\Theta_{ii}}{1/Theta_ii}.
#'     Present for GGM and mixed MRF models; \code{NULL} for OMRF.
#'   \item \code{posterior_mean_indicator}: Symmetric matrix of posterior mean
#'     inclusion probabilities (if edge selection was enabled).
#'
#'   \item  Additional summaries returned when
#'     \code{edge_prior = "Stochastic-Block"}. For more details about this prior
#'     see \insertCite{SekulovskiEtAl_2025;textual}{bgms}.
#'    \itemize{
#'       \item \code{posterior_summary_pairwise_allocations}: Data frame with
#'       posterior summaries (mean, sd, MCSE, ESS, Rhat) for the pairwise
#'       cluster co-occurrence of the nodes. This serves to indicate
#'       whether the estimated posterior allocations,co-clustering matrix
#'       and posterior cluster probabilities (see blow) have converged.
#'       \item \code{posterior_coclustering_matrix}: a symmetric matrix of
#'       pairwise proportions of occurrence of every variable. This matrix
#'       can be plotted to visually inspect the estimated number of clusters
#'       and visually inspect nodes that tend to switch clusters.
#'       \item \code{posterior_mean_allocations}: A vector with the posterior mean
#'       of the cluster allocations of the nodes. This is calculated using the method
#'       proposed in \insertCite{Dahl2009;textual}{bgms}.
#'       \item \code{posterior_mode_allocations}: A vector with the posterior
#'        mode of the cluster allocations of the nodes.
#'       \item \code{posterior_num_blocks}: A data frame with the estimated
#'       posterior inclusion probabilities for all the possible number of clusters.
#'       }
#'   \item \code{raw_samples}: A list of raw MCMC draws per chain:
#'     \describe{
#'       \item{\code{main}}{List of main effect samples.}
#'       \item{\code{pairwise}}{List of pairwise effect samples.}
#'       \item{\code{indicator}}{List of indicator samples
#'         (if edge selection enabled).}
#'       \item{\code{allocations}}{List of cluster allocations
#'         (if SBM prior used).}
#'       \item{\code{nchains}}{Number of chains.}
#'       \item{\code{niter}}{Number of post–warmup iterations per chain.}
#'       \item{\code{parameter_names}}{Named lists of parameter labels.}
#'     }
#'
#'   \item \code{arguments}: A list of function call arguments and metadata
#'     (e.g., number of variables, warmup, sampler settings, package version).
#' }
#'
#' The \code{summary()} method prints formatted posterior summaries, and
#' \code{coef()} extracts posterior mean matrices.
#'
#' NUTS diagnostics (tree depth, divergences, energy, E-BFMI) are included
#' in \code{fit$nuts_diag} if \code{update_method = "nuts"}.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' \donttest{
#' # Run bgm on subset of the Wenchuan dataset
#' fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#'
#' # Posterior inclusion probabilities
#' summary(fit)$indicator
#'
#' # Posterior pairwise effects
#' summary(fit)$pairwise
#' }
#'
#' @export
bgm = function(
  x,
  variable_type = "ordinal",
  baseline_category,
  iter = 1e3,
  warmup = 1e3,
  pairwise_scale = 1,
  main_alpha = 0.5,
  main_beta = 0.5,
  edge_selection = TRUE,
  edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
  inclusion_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  beta_bernoulli_alpha_between = 1,
  beta_bernoulli_beta_between = 1,
  dirichlet_alpha = 1,
  lambda = 1,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis"),
  target_accept,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  standardize = FALSE,
  pseudolikelihood = c("conditional", "marginal"),
  verbose = getOption("bgms.verbose", TRUE),
  progress_callback = NULL,
  interaction_scale,
  burnin,
  save,
  threshold_alpha,
  threshold_beta
) {
  # Set verbose option for internal functions, restore on exit

  old_verbose = getOption("bgms.verbose")
  options(bgms.verbose = verbose)
  on.exit(options(bgms.verbose = old_verbose), add = TRUE)

  if(hasArg(interaction_scale)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(interaction_scale =)", "bgm(pairwise_scale =)")
    if(!hasArg(pairwise_scale)) {
      pairwise_scale = interaction_scale
    }
  }

  if(hasArg(burnin)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(burnin =)", "bgm(warmup =)")
    if(!hasArg(warmup)) {
      warmup = burnin
    }
  }

  if(hasArg(save)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(save =)")
  }

  if(hasArg(threshold_alpha) || hasArg(threshold_beta)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      "bgm(threshold_alpha =, threshold_beta =)",
      "bgm(main_alpha =, main_beta =)"
    )
    if(!hasArg(main_alpha)) main_alpha = threshold_alpha
    if(!hasArg(main_beta)) main_beta = threshold_beta
  }

  # --- Build spec, sample, build output ----------------------------------------
  spec = bgm_spec(
    x = x,
    model_type = "omrf",
    variable_type = variable_type,
    baseline_category = if(hasArg(baseline_category)) baseline_category else 0L,
    na_action = na_action,
    pairwise_scale = pairwise_scale,
    main_alpha = main_alpha,
    main_beta = main_beta,
    standardize = standardize,
    edge_selection = edge_selection,
    edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    beta_bernoulli_alpha_between = beta_bernoulli_alpha_between,
    beta_bernoulli_beta_between = beta_bernoulli_beta_between,
    dirichlet_alpha = dirichlet_alpha,
    lambda = lambda,
    update_method = update_method,
    target_accept = if(hasArg(target_accept)) target_accept else NULL,
    iter = iter,
    warmup = warmup,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    chains = chains,
    cores = cores,
    seed = seed,
    display_progress = display_progress,
    verbose = verbose,
    pseudolikelihood = pseudolikelihood,
    progress_callback = progress_callback
  )

  raw = run_sampler(spec)
  output = build_output(spec, raw)
  return(output)
}
