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
#' \insertCite{MarsmanVandenBerghHaslbeck_2024;textual}{bgms} to multiple
#' groups. The basic idea of modeling, analyzing, and testing group
#' differences in MRFs was introduced in
#' \insertCite{MarsmanWaldorpSekulovskiHaslbeck_2024;textual}{bgms}, where
#' two–group comparisons were conducted using adaptive Metropolis sampling.
#' The present implementation generalizes that approach to more than two
#' groups and supports additional samplers (HMC and NUTS) with staged warmup
#' adaptation.
#'
#' Key components of the model:
#'
#' @seealso \code{vignette("comparison", package = "bgms")} for a worked example.
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
#'   \item \strong{Hamiltonian Monte Carlo (HMC)}: joint updates with fixed
#'     leapfrog trajectories; step size and optionally the mass matrix are
#'     adapted during warmup.
#'   \item \strong{No–U–Turn Sampler (NUTS)}: an adaptive HMC variant with
#'     dynamic trajectory lengths; warmup uses the same staged adaptation
#'     schedule as HMC.
#' }
#'
#' For details on the staged adaptation schedule (fast–slow–fast phases),
#' see \code{\link{bgm}}. In addition, when
#' \code{difference_selection = TRUE}, updates of inclusion indicators are
#' delayed until late warmup. In HMC/NUTS, this appends two extra phases
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
#'   pairwise interactions. Default: \code{2.5}.
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
#' @param update_method Character. Sampling algorithm:
#'   \code{"adaptive-metropolis"}, \code{"hamiltonian-mc"}, or \code{"nuts"}.
#'   Default: \code{"nuts"}.
#' @param target_accept Numeric between 0 and 1. Target acceptance rate.
#'   Defaults: 0.44 (Metropolis), 0.65 (HMC), 0.80 (NUTS).
#' @param hmc_num_leapfrogs Integer. Leapfrog steps for HMC. Default: \code{100}.
#' @param nuts_max_depth Integer. Maximum tree depth for NUTS. Default: \code{10}.
#' @param learn_mass_matrix Logical. If \code{TRUE}, adapt the mass matrix
#'   during warmup (HMC/NUTS only). Default: \code{FALSE}.
#' @param chains Integer. Number of parallel chains. Default: \code{4}.
#' @param cores Integer. Number of CPU cores. Default:
#'   \code{parallel::detectCores()}.
#' @param seed Optional integer. Random seed for reproducibility.
#' @param main_difference_model,reference_category,pairwise_difference_scale,main_difference_scale,pairwise_difference_prior,main_difference_prior,pairwise_difference_probability,main_difference_probability,pairwise_beta_bernoulli_alpha,pairwise_beta_bernoulli_beta,main_beta_bernoulli_alpha,main_beta_bernoulli_beta,interaction_scale,threshold_alpha,threshold_beta,burnin,save
#'   `r lifecycle::badge("deprecated")`
#'   Deprecated arguments as of **bgms 0.1.6.0**.
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
#'     \code{posterior_mean_pairwise_baseline}: posterior mean matrices
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
#' fit <- bgmCompare(x, y)
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
  variable_type = "ordinal",
  baseline_category,
  difference_scale = 1,
  difference_prior = c("Bernoulli", "Beta-Bernoulli"),
  difference_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  pairwise_scale = 2.5,
  main_alpha = 0.5,
  main_beta = 0.5,
  iter = 1e3,
  warmup = 1e3,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis", "hamiltonian-mc"),
  target_accept,
  hmc_num_leapfrogs = 100,
  nuts_max_depth = 10,
  learn_mass_matrix = FALSE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
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

  # Check update method
  update_method_input = update_method
  update_method = match.arg(update_method)

  # Check target acceptance rate
  if(hasArg(target_accept)) {
    target_accept = min(target_accept, 1 - sqrt(.Machine$double.eps))
    target_accept = max(target_accept, 0 + sqrt(.Machine$double.eps))
  } else {
    if(update_method == "adaptive-metropolis") {
      target_accept = 0.44
    } else if(update_method == "hamiltonian-mc") {
      target_accept = 0.65
    } else if(update_method == "nuts") {
      target_accept = 0.80
    }
  }

  # Check and preprocess data
  x = data_check(x, "x")
  if(hasArg(y)) {
    y = data_check(y, "y")
    if(ncol(x) != ncol(y)) stop("x and y must have the same number of columns.")
  }

  if(!hasArg(y) & !hasArg(group_indicator)) {
    stop(paste0(
      "For multi-group designs, the bgmCompare function requires input for\n",
      "either y (group 2 data) or group_indicator (group indicator)."
    ))
  }

  # Validate group indicators
  if(!hasArg(y) && hasArg(group_indicator)) {
    group_indicator = as.vector(group_indicator)
    if(anyNA(group_indicator)) stop("group_indicator cannot contain missing values.")
    if(length(group_indicator) != nrow(x)) stop("Length of group_indicator must match number of rows in x.")
  }

  # Model and preprocessing
  if(!hasArg(y)) {
    y = NULL
  }
  if(!hasArg(group_indicator)) {
    group_indicator = NULL
  }

  model = check_compare_model(
    x = x, y = y, group_indicator = group_indicator, difference_selection = difference_selection,
    variable_type = variable_type, baseline_category = baseline_category,
    difference_scale = difference_scale, difference_prior = difference_prior,
    difference_probability = difference_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    pairwise_scale = pairwise_scale, main_alpha = main_alpha,
    main_beta = main_beta
  )

  x = model$x
  group = model$group
  ordinal_variable = model$variable_bool
  baseline_category = model$baseline_category
  difference_prior = model$difference_prior

  # Check Gibbs input
  check_positive_integer(iter, "iter")
  check_non_negative_integer(warmup, "warmup")

  # Check na_action
  na_action_input = na_action
  na_action = try(match.arg(na_action), silent = TRUE)
  if(inherits(na_action, "try-error")) {
    stop(sprintf("Invalid value for `na_action`. Expected 'listwise' or 'impute', got: %s", na_action_input))
  }

  # Check display_progress
  progress_type = progress_type_from_display_progress(display_progress)


  ## Format data
  data = compare_reformat_data(
    x = x, group = group,
    na_action = na_action,
    variable_bool = ordinal_variable,
    baseline_category = baseline_category
  )

  x = data$x
  group = data$group
  num_obs_groups = tabulate(group)
  missing_index = data$missing_index
  num_categories = data$num_categories

  na_impute = data$na_impute
  baseline_category = data$baseline_category
  num_variables = ncol(x)
  num_interactions = num_variables * (num_variables - 1) / 2

  # Compute `counts_per_category`
  counts_per_category = compute_counts_per_category(
    x, num_categories, group
  )

  # Compute sufficient statistics for Blume-Capel variables
  blume_capel_stats = compute_blume_capel_stats(
    x, baseline_category, ordinal_variable, group
  )
  for (i in which(!ordinal_variable)) {
    x[, i] = sum(x[, i] - baseline_category[i])
  }

  # Compute sufficient statistics for pairwise interactions
  pairwise_stats = compute_pairwise_stats(
    x, group
  )


  # Index vector used to sample interactions in a random order -----------------
  Index = matrix(0, nrow = num_interactions, ncol = 3)
  counter = 0
  for(variable1 in 1:(num_variables - 1)) {
    for(variable2 in (variable1 + 1):num_variables) {
      counter = counter + 1
      Index[counter, ] = c(counter, variable1 - 1, variable2 - 1)
    }
  }

  # Gibbs sampling
  # Prepare indices for main and pairwise effects
  main_effect_indices = matrix(NA, nrow = num_variables, ncol = 2)
  for(variable in seq_len(num_variables)) {
    if(variable > 1) {
      main_effect_indices[variable, 1] = 1 + main_effect_indices[variable - 1, 2]
    } else {
      main_effect_indices[variable, 1] = 0 # C++ starts at zero
    }
    if(ordinal_variable[variable]) {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + num_categories[variable] - 1
    } else {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + 1
    }
  }

  pairwise_effect_indices = matrix(NA, nrow = num_variables, ncol = num_variables)
  tel = 0
  for(v1 in seq_len(num_variables - 1)) {
    for(v2 in seq((v1 + 1), num_variables)) {
      pairwise_effect_indices[v1, v2] = tel
      pairwise_effect_indices[v2, v1] = tel
      tel = tel + 1 # C++ starts at zero
    }
  }

  # Compute group-level data
  num_groups = length(unique(group))
  group_indices = matrix(NA, nrow = num_groups, ncol = 2)

  # Align observations with sorted group
  observations = x
  sorted_group = sort(group)
  for(g in unique(group)) {
    observations[which(sorted_group == g), ] = x[which(group == g), ]
    group_indices[g, 1] = min(which(sorted_group == g)) - 1 # C++ starts at zero
    group_indices[g, 2] = max(which(sorted_group == g)) - 1 # C++ starts at zero
  }

  # Compute projection matrix for group differences
  one = matrix(1, nrow = num_groups, ncol = num_groups)
  V = diag(num_groups) - one / num_groups
  projection = eigen(V)$vectors[, -num_groups]
  if(num_groups == 2) {
    projection = matrix(projection, ncol = 1) / sqrt(2)
  }

  # Setting the seed
  if(missing(seed) || is.null(seed)) {
    # Draw a random seed if none provided
    seed = sample.int(.Machine$integer.max, 1)
  }

  if(!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
    stop("Argument 'seed' must be a single non-negative integer.")
  }

  seed <- as.integer(seed)


  # Call the Rcpp function
  out = run_bgmCompare_parallel(
    observations = observations,
    num_groups = num_groups,
    counts_per_category = counts_per_category,
    blume_capel_stats = blume_capel_stats,
    pairwise_stats = pairwise_stats,
    num_categories = num_categories,
    main_alpha = main_alpha,
    main_beta = main_beta,
    pairwise_scale = pairwise_scale,
    difference_scale = difference_scale,
    difference_selection_alpha = beta_bernoulli_alpha,
    difference_selection_beta = beta_bernoulli_beta,
    difference_prior = model$difference_prior, iter = iter, warmup = warmup,
    na_impute = na_impute, missing_data_indices = missing_index,
    is_ordinal_variable = ordinal_variable,
    baseline_category = baseline_category,
    difference_selection = difference_selection,
    main_effect_indices = main_effect_indices,
    pairwise_effect_indices = pairwise_effect_indices,
    target_accept = target_accept,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    projection = projection,
    group_membership = sorted_group - 1,
    group_indices = group_indices,
    interaction_index_matrix = Index,
    inclusion_probability = model$inclusion_probability_difference,
    num_chains = chains, nThreads = cores,
    seed = seed,
    update_method = update_method, hmc_num_leapfrogs = hmc_num_leapfrogs,
    progress_type = progress_type
  )

  userInterrupt = any(vapply(out, FUN = `[[`, FUN.VALUE = logical(1L), "userInterrupt"))
  if(userInterrupt) {
    warning("Stopped sampling after user interrupt, results are likely uninterpretable.")
    output <- tryCatch(
      prepare_output_bgmCompare(
        out = out,
        observations = observations,
        num_categories = num_categories,
        is_ordinal_variable = ordinal_variable,
        num_groups = num_groups,
        group = sorted_group,
        iter = iter,
        warmup = warmup,
        main_effect_indices = main_effect_indices,
        pairwise_effect_indices = pairwise_effect_indices,
        data_columnnames = if(is.null(colnames(x))) paste0("Variable ", seq_len(ncol(x))) else colnames(x),
        difference_selection = difference_selection,
        difference_prior = difference_prior,
        difference_selection_alpha = beta_bernoulli_alpha,
        difference_selection_beta = beta_bernoulli_beta,
        inclusion_probability = model$inclusion_probability_difference,
        pairwise_scale = pairwise_scale,
        difference_scale = difference_scale,
        update_method = update_method,
        target_accept = target_accept,
        nuts_max_depth = nuts_max_depth,
        hmc_num_leapfrogs = hmc_num_leapfrogs,
        learn_mass_matrix = learn_mass_matrix,
        num_chains = chains,
        projection = projection
      ),
      error = function(e) {
        list(partial = out, error = conditionMessage(e))
      },
      warning = function(w) {
        list(partial = out, warning = conditionMessage(w))
      }
    )
    return(output)
  }

  # Main output handler in the wrapper function
  output = prepare_output_bgmCompare(
    out = out,
    observations = observations,
    num_categories = num_categories,
    is_ordinal_variable = ordinal_variable,
    num_groups = num_groups,
    group = sorted_group,
    iter = iter,
    warmup = warmup,
    main_effect_indices = main_effect_indices,
    pairwise_effect_indices = pairwise_effect_indices,
    data_columnnames = if(is.null(colnames(x))) paste0("Variable ", seq_len(ncol(x))) else colnames(x),
    difference_selection = difference_selection,
    difference_prior = difference_prior,
    difference_selection_alpha = beta_bernoulli_alpha,
    difference_selection_beta = beta_bernoulli_beta,
    inclusion_probability = model$inclusion_probability_difference,
    pairwise_scale = pairwise_scale,
    difference_scale = difference_scale,
    update_method = update_method,
    target_accept = target_accept,
    nuts_max_depth = nuts_max_depth,
    hmc_num_leapfrogs = hmc_num_leapfrogs,
    learn_mass_matrix = learn_mass_matrix,
    num_chains = chains, projection = projection
  )

  if(update_method == "nuts") {
    nuts_diag = summarize_nuts_diagnostics(out, nuts_max_depth = nuts_max_depth)
    output$nuts_diag = nuts_diag
  }

  return(output)
}
