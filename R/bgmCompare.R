#' Bayesian Variable Selection or Bayesian Estimation for Differences in Markov Random Fields
#'
#' @description
#' The \code{bgmCompare} function estimates the pseudoposterior distribution of
#' the parameters of a Markov Random Field (MRF) model for mixed binary and ordinal
#' variables, as well as differences in pairwise interactions and category thresholds
#' across groups. Groups are assumed to be \code{G} independent samples.
#'
#' @details
#' This function models group differences in Markov Random Fields using a Bayesian
#' framework. It supports binary and ordinal variables, and includes options for
#' Bayesian variable selection on the differences in both pairwise interactions
#' and threshold parameters. Key components are described in the sections below.
#'
#' @section Pairwise Interactions:
#' Pairwise interactions between variables \code{i} and \code{j} are modeled as:
#' \deqn{\boldsymbol{\theta}_{ij} = \phi_{ij} + \boldsymbol{\delta}_{ij},}{
#' \boldsymbol{\theta}_{ij} = \phi_{ij} + \boldsymbol{\delta}_{ij},}
#' where:
#'
#' \itemize{
#'  \item \eqn{\boldsymbol{\theta}_{ij}}{\boldsymbol{\theta}_{ij}} is the vector of pairwise interaction parameters of length \code{G}.
#'  \item \eqn{\phi_{ij}}{\phi_{ij}} is the overall pairwise interaction (nuisance parameter).
#'  \item \eqn{\boldsymbol{\delta}_{ij}}{\boldsymbol{\delta}_{ij}} represents group-specific differences constrained to sum to zero for identification.
#' }
#'
#' @section Ordinal Variables:
#' The function supports two types of ordinal variables:
#'
#' \strong{Regular ordinal variables}:
#' Introduce a threshold parameter for each category except the lowest, modeled as:
#' \deqn{\boldsymbol{\mu}_{ic} = \tau_{ic} + \boldsymbol{\epsilon}_{ic},}{\boldsymbol{\mu}_{ic} = \tau_{ic} + \boldsymbol{\epsilon}_{ic},}
#' where:
#' \itemize{
#'  \item \eqn{\tau_{ic}}{\tau_{ic}} denotes an overall effect (nuisance parameter).
#'  \item \eqn{\boldsymbol{\epsilon}_{ic}}{\boldsymbol{\epsilon}_{ic}} represents group-specific differences constrained to sum to zero.
#' }
#'
#' \strong{Blume-Capel ordinal variables}:
#' Assume a specific reference category and score responses based on distance to it:
#' \deqn{\boldsymbol{\mu}_{ic} = (\tau_{i1} + \boldsymbol{\epsilon}_{i1}) \cdot c + (\tau_{i2} + \boldsymbol{\epsilon}_{i2}) \cdot (c - r)^2,}{
#' \boldsymbol{\mu}_{ic} = (\tau_{i1} + \boldsymbol{\epsilon}_{i1}) \cdot c + (\tau_{i2} + \boldsymbol{\epsilon}_{i2}) \cdot (c - r)^2,}
#' where:
#' \itemize{
#'  \item \code{r} is the reference category.
#'  \item \eqn{\tau_{i1}}{\tau_{i1}} and \eqn{\tau_{i2}}{\tau_{i2}} are nuisance parameters.
#'  \item \eqn{\boldsymbol{\epsilon}_{i1}}{\boldsymbol{\epsilon}_{i1}} and \eqn{\boldsymbol{\epsilon}_{i2}}{\boldsymbol{\epsilon}_{i2}} represent group-specific differences.
#' }
#'
#' @section Variable Selection:
#' Bayesian variable selection enables testing of parameter differences or equivalence across groups. Independent spike-and-slab priors are applied to difference parameters:
#' \itemize{
#'  \item \strong{Bernoulli Model}: Assigns a fixed probability to parameter inclusion.
#'  \item \strong{Beta-Bernoulli Model}: Incorporates a beta prior to model inclusion probabilities.
#' }
#'
#' @section Gibbs Sampling:
#'
#' Parameters are estimated using a Metropolis-within-Gibbs sampling scheme.
#' When \code{difference_selection = TRUE}, the algorithm runs \code{2 * burnin} warmup iterations:
#' \itemize{
#'   \item First half without difference selection.
#'   \item Second half with edge selection enabled.
#' }
#' This warmup strategy improves stability of adaptive Metropolis-Hastings proposals and starting values.
#'
#' @section Saving Options:
#' Users can store sampled states for parameters (\code{main_effects}, \code{pairwise_effects}, \code{indicator}) during Gibbs sampling. Enabling these options (\code{save_main}, \code{save_pairwise}, \code{save_indicator}) increases output size and memory usage, so use them judiciously.
#'
#' @param x Data frame or matrix with binary and ordinal responses. Regular ordinal variables should be coded as integers starting from 0. Missing categories are collapsed for regular ordinal variables but retained for Blume-Capel variables.
#' @param y A data frame or matrix similar to \code{x}, used for two-group designs. \code{x} contains Group 1 data, and \code{y} contains Group 2 data. Ignored for multi-group designs.
#' @param g Group membership vector for rows in \code{x}. Required for multi-group designs; ignored if \code{y} is provided.
#' @param difference_selection Logical. If \code{TRUE}, the function models the inclusion or exclusion of parameter differences. Default: \code{TRUE}.
#' @param save_main,save_pairwise,save_indicator Logical. Enable saving sampled states for \code{main_effects}, \code{pairwise_effects}, and \code{indicator}, respectively. Default: \code{FALSE}.
#' @param main_difference_model Character. Specifies how to handle threshold differences when categories are unmatched. Options: \code{"Collapse"}, \code{"Free"}. The option "Collapse" collapses categories unobserved in one or more groups. The option "Free" option estimates thresholds separately for each group and does not model their difference. Default: \code{"Free"}.
#' @param variable_type Character or vector. Specifies the type of variables in \code{x} (\code{"ordinal"} or \code{"blume-capel"}). Default: \code{"ordinal"}.
#' @param reference_category Integer or vector. Reference category for Blume-Capel variables. Required if there is at least one Blume-Capel variable.
#' @param pairwise_difference_scale Double. Scale parameter for the Cauchy prior on pairwise differences. Default: \code{1}.
#' @param main_difference_scale Double. Scale parameter for the Cauchy prior on threshold differences. Default: \code{1}.
#' @param pairwise_difference_prior,main_difference_prior Character. Specifies the inclusion probability model (\code{"Bernoulli"} or \code{"Beta-Bernoulli"}). Default: \code{"Bernoulli"}.
#' @param pairwise_difference_probability A numeric value or a \eqn{p \times p} matrix specifying the prior inclusion probability of a pairwise difference in the Bernoulli model. A single value applies the same probability to all pairs, while a matrix allows for edge-specific probabilities. Default: 0.5 for equal prior probability for inclusion and exclusion.
#' @param main_difference_probability A numeric value or a length-\eqn{p} vector specifying the prior inclusion probability of a threshold difference in the Bernoulli model. A single value applies the same probability to all variables, while a vector allows for variable-specific probabilities. Default: 0.5 to indicate no prior preference.
#' @param iter,burnin Integer. Number of Gibbs iterations (\code{iter}) and burn-in iterations (\code{burnin}). Defaults: \code{iter = 1e4}, \code{burnin = 1e3}.
#' @param na_action Character. Specifies handling of missing data. \code{"listwise"} deletes rows with missing values; \code{"impute"} imputes values during Gibbs sampling. Default: \code{"listwise"}.
#' @param display_progress Logical. Show progress bar during computation. Default: \code{TRUE}.
#' @param threshold_alpha,threshold_beta Double. Shape parameters for the beta-prime prior on nuisance threshold parameters.
#' @param interaction_scale Double. Scale of the Cauchy prior for nuisance pairwise interactions. Default: \code{2.5}.
#' @param main_beta_bernoulli_alpha,main_beta_bernoulli_beta Double. Shape parameters for the Beta-Bernoulli prior on threshold differences.
#' @param pairwise_beta_bernoulli_alpha,pairwise_beta_bernoulli_beta Double. Shape parameters for the Beta-Bernoulli prior on pairwise differences.
#' @param save Logical. If true, sampled states for all parameters are returned. Deprecated.
#' @param save_main,save_pairwise,save_indicator Logical. Enable saving sampled states for `main_effects`, `pairwise_effects`, and `indicator`, respectively. Default: `FALSE`.
#'
#' @return A list containing the posterior means and, optionally, sampled states based on the \code{save_*} options. The returned components include:
#' \itemize{
#'  \item \code{posterior_mean_main}, \code{posterior_mean_pairwise}, and \code{posterior_mean_indicator} for posterior means.
#'  \item If saving options are enabled, the list also includes:
#'    \itemize{
#'      \item \code{raw_samples_main} – sampled states of main effects.
#'      \item \code{raw_samples_pairwise} – sampled states of pairwise effects.
#'      \item \code{raw_samples_indicator} – sampled states of inclusion indicators.
#'    }
#' }
#' In addition to the results of the analysis, the output lists some of the
#' arguments of its call. This is useful for post-processing the results.
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
    interaction_scale = 2.5,
    threshold_alpha = 0.5,
    threshold_beta = 0.5,
    iter = 1e3,
    burnin = 1e3,
    na_action = c("listwise", "impute"),
    display_progress =  c("per-chain", "total", "none"),
    update_method = c("nuts", "adaptive-metropolis", "hamiltonian-mc"),
    target_accept,
    hmc_num_leapfrogs = 100,
    nuts_max_depth = 10,
    learn_mass_matrix = FALSE,
    chains = 4,
    cores = parallel::detectCores(),
    seed = NULL
) {
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
      target_accept = 0.60
    }
  }

  # Check and preprocess data
  x = data_check(x, "x")
  if (hasArg(y)) {
    y = data_check(y, "y")
    if (ncol(x) != ncol(y)) stop("x and y must have the same number of columns.")
  }

  if(!hasArg(y) & !hasArg(group_indicator))
    stop(paste0("For multi-group designs, the bgmCompare function requires input for\n",
                "either y (group 2 data) or group_indicator (group indicator)."))

  # Validate group indicators
  if (!hasArg(y) && hasArg(group_indicator)) {
    group_indicator = as.vector(group_indicator)
    if (anyNA(group_indicator)) stop("group_indicator cannot contain missing values.")
    if (length(group_indicator) != nrow(x)) stop("Length of group_indicator must match number of rows in x.")
  }

  # Model and preprocessing
  if(!hasArg(y))
    y = NULL
  if(!hasArg(group_indicator))
    group_indicator = NULL

  model = check_compare_model(
    x = x, y = y, g = group_indicator, difference_selection = difference_selection,
    variable_type = variable_type, baseline_category = baseline_category,
    difference_scale = difference_scale, difference_prior = difference_prior,
    difference_probability = difference_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    interaction_scale = interaction_scale, threshold_alpha = threshold_alpha,
    threshold_beta = threshold_beta
  )

  x = model$x
  group = model$group
  ordinal_variable = model$variable_bool
  baseline_category = model$baseline_category

  # Check Gibbs input
  check_positive_integer(iter, "iter")
  check_non_negative_integer(burnin, "burnin")

  # Check na_action
  na_action_input = na_action
  na_action = try(match.arg(na_action), silent = TRUE)
  if (inherits(na_action, "try-error")) {
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

  # Compute sufficient statistics for pairwise interactions
  pairwise_stats = compute_pairwise_stats(
    x, group
  )


  # Index vector used to sample interactions in a random order -----------------
  Index = matrix(0, nrow = num_interactions, ncol = 3)
  counter = 0
  for(variable1 in 1:(num_variables - 1)) {
    for(variable2 in (variable1 + 1):num_variables) {
      counter =  counter + 1
      Index[counter, ] = c(counter, variable1 - 1, variable2 - 1)
    }
  }

  # Gibbs sampling
  # Prepare indices for main and pairwise effects
  main_effect_indices = matrix(NA, nrow = num_variables, ncol = 2)
  for (variable in seq_len(num_variables)) {
    if (variable > 1) {
      main_effect_indices[variable, 1] = 1 + main_effect_indices[variable - 1, 2]
    } else {
      main_effect_indices[variable, 1] = 0  # C++ starts at zero
    }
    if (ordinal_variable[variable]) {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + max(num_categories[variable, ]) - 1
    } else {
      main_effect_indices[variable, 2] = main_effect_indices[variable, 1] + 1
    }
  }

  pairwise_effect_indices = matrix(NA, nrow = num_variables, ncol = num_variables)
  tel = 0
  for (v1 in seq_len(num_variables - 1)) {
    for (v2 in seq((v1 + 1), num_variables)) {
      pairwise_effect_indices[v1, v2] = tel
      pairwise_effect_indices[v2, v1] = tel
      tel = tel + 1  # C++ starts at zero
    }
  }

  # Compute group-level data
  num_groups = length(unique(group))
  group_indices = matrix(NA, nrow = num_groups, ncol = 2)

  # Align observations with sorted group
  observations = x
  sorted_group = sort(group)
  for (g in unique(group)) {
    observations[which(sorted_group == g), ] = x[which(group == g), ]
    group_indices[g, 1] = min(which(sorted_group == g)) - 1  # C++ starts at zero
    group_indices[g, 2] = max(which(sorted_group == g)) - 1  # C++ starts at zero
  }

  # Compute projection matrix for group differences
  one = matrix(1, nrow = num_groups, ncol = num_groups)
  V = diag(num_groups) - one / num_groups
  projection = eigen(V)$vectors[, -num_groups]
  if (num_groups == 2) {
    projection = matrix(projection, ncol = 1) / sqrt(2)
  }

  #Setting the seed
  if (missing(seed) || is.null(seed)) {
    # Draw a random seed if none provided
    seed = sample.int(.Machine$integer.max, 1)
  }

  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
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
    num_categories = num_categories[, 1],
    main_alpha = threshold_alpha,
    main_beta = threshold_beta,
    pairwise_scale = interaction_scale,
    difference_scale = difference_scale,
    difference_selection_alpha = beta_bernoulli_alpha,
    difference_selection_beta = beta_bernoulli_beta,
    difference_prior = model$difference_prior, iter = iter, burnin = burnin,
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

  # Main output handler in the wrapper function
  # output = prepare_output_bgmCompare2(
  #   out = out, ...
  # )

  return(out)
}