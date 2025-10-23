#' Bayesian Estimation or Edge Selection for Markov Random Fields
#'
#' @description
#' The \code{bgm} function estimates the pseudoposterior distribution of
#' category thresholds (main effects) and pairwise interaction parameters of a
#' Markov Random Field (MRF) model for binary and/or ordinal variables.
#' Optionally, it performs Bayesian edge selection using spike-and-slab
#' priors to infer the network structure.
#'
#' @details
#' This function models the joint distribution of binary and ordinal variables
#' using a Markov Random Field, with support for edge selection through Bayesian
#' variable selection. The statistical foundation of the model is described in
#' \insertCite{MarsmanVandenBerghHaslbeck_2024;textual}{bgms}, where the ordinal
#' MRF model and its Bayesian estimation procedure were first introduced. While
#' the implementation in \pkg{bgms} has since been extended and updated (e.g.,
#' alternative priors, parallel chains, HMC/NUTS warmup), it builds on that
#' original framework.
#'
#' Key components of the model are described in the sections below.
#'
#' @seealso \code{vignette("intro", package = "bgms")} for a worked example.
#'
#' @section Ordinal Variables:
#' The function supports two types of ordinal variables:
#'
#' \strong{Regular ordinal variables}:
#' Assigns a category threshold parameter to each response category except the
#' lowest. The model imposes no additional constraints on the distribution of
#' category responses.
#'
#' \strong{Blume-Capel ordinal variables}:
#' Assume a baseline category (e.g., a “neutral” response) and score responses
#' by distance from this baseline. Category thresholds are modeled as:
#'
#' \deqn{\mu_{c} = \alpha \cdot c + \beta \cdot (c - b)^2}
#'
#' where:
#' \itemize{
#'   \item \eqn{\mu_{c}}: category threshold for category \eqn{c}
#'   \item \eqn{\alpha}: linear trend across categories
#'   \item \eqn{\beta}: preference toward or away from the baseline
#'    \itemize{
#'      \item If \eqn{\beta < 0}, the model favors responses near the baseline
#'      category;
#'      \item if \eqn{\beta > 0}, it favors responses farther away (i.e.,
#'      extremes).
#'    }
#'   \item \eqn{b}: baseline category
#' }
#'
#'
#' @section Edge Selection:
#' When \code{edge_selection = TRUE}, the function performs Bayesian variable
#' selection on the pairwise interactions (edges) in the MRF using
#' spike-and-slab priors.
#'
#' Supported priors for edge inclusion:
#' \itemize{
#'   \item \strong{Bernoulli}: Fixed inclusion probability across edges.
#'   \item \strong{Beta-Bernoulli}: Inclusion probability is assigned a Beta
#'   prior distribution.
#'   \item \strong{Stochastic-Block}: Cluster-based edge priors with Beta,
#'   Dirichlet, and Poisson hyperpriors.
#' }
#'
#' All priors operate via binary indicator variables controlling the inclusion
#' or exclusion of each edge in the MRF.
#'
#' @section Prior Distributions:
#'
#' \itemize{
#'   \item \strong{Pairwise effects}: Modeled with a Cauchy (slab) prior.
#'   \item \strong{Main effects}: Modeled using a beta-prime
#'   distribution.
#'   \item \strong{Edge indicators}: Use either a Bernoulli, Beta-Bernoulli, or
#'   Stochastic-Block prior (as above).
#' }
#'
#' @section Sampling Algorithms and Warmup:
#'
#' Parameters are updated within a Gibbs framework, but the conditional
#' updates can be carried out using different algorithms:
#' \itemize{
#'   \item \strong{Adaptive Metropolis–Hastings}: Componentwise random–walk
#'     updates for main effects and pairwise effects. Proposal standard
#'     deviations are adapted during burn–in via Robbins–Monro updates
#'     toward a target acceptance rate.
#'
#'   \item \strong{Hamiltonian Monte Carlo (HMC)}: Joint updates of all
#'     parameters using fixed–length leapfrog trajectories. Step size is
#'     tuned during warmup via dual–averaging; the diagonal mass matrix can
#'     also be adapted if \code{learn_mass_matrix = TRUE}.
#'
#'   \item \strong{No–U–Turn Sampler (NUTS)}: An adaptive extension of HMC
#'     that dynamically chooses trajectory lengths. Warmup uses a staged
#'     adaptation schedule (fast–slow–fast) to stabilize step size and, if
#'     enabled, the mass matrix.
#' }
#'
#' When \code{edge_selection = TRUE}, updates of edge–inclusion indicators
#' are carried out with Metropolis–Hastings steps. These are switched on
#' after the core warmup phase, ensuring that graph updates occur only once
#' the samplers’ tuning parameters (step size, mass matrix, proposal SDs)
#' have stabilized.
#'
#' After warmup, adaptation is disabled. Step size and mass matrix are
#' fixed at their learned values, and proposal SDs remain constant.
#'
#' @section Warmup and Adaptation:
#'
#' The warmup procedure in \code{bgm} is based on the multi–stage adaptation
#' schedule used in Stan \insertCite{stan-manual}{bgms}. Warmup iterations are
#' split into several phases:
#'
#' \itemize{
#'   \item \strong{Stage 1 (fast adaptation)}: A short initial interval
#'     where only step size (for HMC/NUTS) is adapted, allowing the chain
#'     to move quickly toward the typical set.
#'
#'   \item \strong{Stage 2 (slow windows)}: A sequence of expanding,
#'     memoryless windows where both step size and, if
#'     \code{learn_mass_matrix = TRUE}, the diagonal mass matrix are
#'     adapted. Each window ends with a reset of the dual–averaging scheme
#'     for improved stability.
#'
#'   \item \strong{Stage 3a (final fast interval)}: A short interval at the
#'     end of the core warmup where the step size is adapted one final time.
#'
#'   \item \strong{Stage 3b (proposal–SD tuning)}: Only active when
#'     \code{edge_selection = TRUE} under HMC/NUTS. In this phase,
#'     Robbins–Monro adaptation of proposal standard deviations is
#'     performed for the Metropolis steps used in edge–selection moves.
#'
#'   \item \strong{Stage 3c (graph selection warmup)}: Also only relevant
#'     when \code{edge_selection = TRUE}. At the start of this phase, a
#'     random graph structure is initialized, and Metropolis–Hastings
#'     updates for edge inclusion indicators are switched on.
#' }
#'
#' When \code{edge_selection = FALSE}, the total number of warmup iterations
#' equals the user–specified \code{burnin}. When \code{edge_selection = TRUE}
#' and \code{update_method} is \code{"nuts"} or \code{"hamiltonian-mc"},
#' the schedule automatically appends additional Stage-3b and Stage-3c
#' intervals, so the total warmup is strictly greater than the requested
#' \code{burnin}.
#'
#' After all warmup phases, the sampler transitions to the sampling phase
#' with adaptation disabled. Step size and mass matrix (for HMC/NUTS) are
#' fixed at their learned values, and proposal SDs remain constant.
#'
#' This staged design improves stability of proposals and ensures that both
#' local parameters (step size) and global parameters (mass matrix, proposal
#' SDs) are tuned before collecting posterior samples.
#'
#' For adaptive Metropolis–Hastings runs, step size and mass matrix
#' adaptation are not relevant. Proposal SDs are tuned continuously during
#' burn–in using Robbins–Monro updates, without staged fast/slow intervals.
#'
#' @section Missing Data:
#'
#' If \code{na_action = "listwise"}, observations with missing values are
#' removed.
#' If \code{na_action = "impute"}, missing values are imputed during Gibbs
#' sampling.
#'
#' @param x A data frame or matrix with \code{n} rows and \code{p} columns
#'   containing binary and ordinal responses. Variables are automatically
#'   recoded to non-negative integers (\code{0, 1, ..., m}). For regular
#'   ordinal variables, unobserved categories are collapsed; for
#'   Blume–Capel variables, all categories are retained.
#'
#' @param variable_type Character or character vector. Specifies the type of
#'   each variable in \code{x}. Allowed values: \code{"ordinal"} or
#'   \code{"blume-capel"}. Binary variables are automatically treated as
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
#'   interaction parameters. Default: \code{2.5}.
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
#'   priors. Must be positive. Defaults: \code{beta_bernoulli_alpha = 1} and
#'   \code{beta_bernoulli_beta = 1}.
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
#' @param display_progress Logical. Whether to show a progress bar during
#'   sampling. Default: \code{TRUE}.
#'
#' @param update_method Character. Specifies how the MCMC sampler updates
#'   the model parameters:
#'   \describe{
#'     \item{"adaptive-metropolis"}{Componentwise adaptive Metropolis–Hastings
#'       with Robbins–Monro proposal adaptation.}
#'     \item{"hamiltonian-mc"}{Hamiltonian Monte Carlo with fixed path length
#'       (number of leapfrog steps set by \code{hmc_num_leapfrogs}).}
#'     \item{"nuts"}{The No-U-Turn Sampler, an adaptive form of HMC with
#'       dynamically chosen trajectory lengths.}
#'   }
#'   Default: \code{"nuts"}.
#'
#' @param target_accept Numeric between 0 and 1. Target acceptance rate for
#'   the sampler. Defaults are set automatically if not supplied:
#'   \code{0.44} for adaptive Metropolis, \code{0.65} for HMC,
#'   and \code{0.60} for NUTS.
#'
#' @param hmc_num_leapfrogs Integer. Number of leapfrog steps for Hamiltonian
#'   Monte Carlo. Must be positive. Default: \code{100}.
#'
#' @param nuts_max_depth Integer. Maximum tree depth in NUTS. Must be positive.
#'   Default: \code{10}.
#'
#' @param learn_mass_matrix Logical. If \code{TRUE}, adapt a diagonal mass
#'   matrix during warmup (HMC/NUTS only). If \code{FALSE}, use the identity
#'   matrix. Default: \code{FALSE}.
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
#'   Deprecated arguments as of **bgms 0.1.6.0**.
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
#'     (mean, sd, MCSE, ESS, Rhat) for category threshold parameters.
#'   \item \code{posterior_summary_pairwise}: Data frame with posterior
#'     summaries for pairwise interaction parameters.
#'   \item \code{posterior_summary_indicator}: Data frame with posterior
#'     summaries for edge inclusion indicators (if \code{edge_selection = TRUE}).
#'
#'   \item \code{posterior_mean_main}: Matrix of posterior mean thresholds
#'     (rows = variables, cols = categories or parameters).
#'   \item \code{posterior_mean_pairwise}: Symmetric matrix of posterior mean
#'     pairwise interaction strengths.
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
#' fit = bgm(x = Wenchuan[, 1:5])
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
    pairwise_scale = 2.5,
    main_alpha = 0.5,
    main_beta = 0.5,
    edge_selection = TRUE,
    edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
    inclusion_probability = 0.5,
    beta_bernoulli_alpha = 1,
    beta_bernoulli_beta = 1,
    dirichlet_alpha = 1,
    lambda = 1,
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
    interaction_scale,
    burnin,
    save,
    threshold_alpha,
    threshold_beta
) {
  if (hasArg(interaction_scale)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(interaction_scale =)", "bgm(pairwise_scale =)")
    if (!hasArg(pairwise_scale)) {
      pairwise_scale = interaction_scale
    }
  }

  if (hasArg(burnin)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(burnin =)", "bgm(warmup =)")
    if (!hasArg(warmup)) {
      warmup = burnin
    }
  }

  if (hasArg(save)) {
    lifecycle::deprecate_warn("0.1.6.0", "bgm(save =)")
  }

  if (hasArg(threshold_alpha) || hasArg(threshold_beta)) {
    lifecycle::deprecate_warn("0.1.6.0",
                              "bgm(threshold_alpha =, threshold_beta =)",
                              "bgm(main_alpha =, main_beta =)"
    )
    if (!hasArg(main_alpha)) main_alpha = threshold_alpha
    if (!hasArg(main_beta)) main_beta = threshold_beta
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

  #Check data input ------------------------------------------------------------
  if(!inherits(x, what = "matrix") && !inherits(x, what = "data.frame"))
    stop("The input x needs to be a matrix or dataframe.")
  if(inherits(x, what = "data.frame"))
    x = data.matrix(x)
  if(ncol(x) < 2)
    stop("The matrix x should have more than one variable (columns).")
  if(nrow(x) < 2)
    stop("The matrix x should have more than one observation (rows).")

  #Check model input -----------------------------------------------------------
  model = check_model(x = x,
                      variable_type = variable_type,
                      baseline_category = baseline_category,
                      pairwise_scale = pairwise_scale,
                      main_alpha = main_alpha,
                      main_beta = main_beta,
                      edge_selection = edge_selection,
                      edge_prior = edge_prior,
                      inclusion_probability = inclusion_probability,
                      beta_bernoulli_alpha = beta_bernoulli_alpha,
                      beta_bernoulli_beta = beta_bernoulli_beta,
                      dirichlet_alpha = dirichlet_alpha,
                      lambda = lambda)

  # ----------------------------------------------------------------------------
  # The vector variable_type is now coded as boolean.
  # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
  # ----------------------------------------------------------------------------
  variable_bool = model$variable_bool
  # ----------------------------------------------------------------------------

  baseline_category = model$baseline_category
  edge_selection = model$edge_selection
  edge_prior = model$edge_prior
  inclusion_probability = model$inclusion_probability

  #Check Gibbs input -----------------------------------------------------------
  check_positive_integer(iter, "iter")
  check_non_negative_integer(warmup, "warmup")
  if(warmup < 1e3)
    warning("The warmup parameter is set to a low value. This may lead to unreliable results. Reset to a minimum of 1000 iterations.")
  warmup = max(warmup, 1e3) # Set minimum warmup to 1000 iterations

  check_positive_integer(hmc_num_leapfrogs, "hmc_num_leapfrogs")
  hmc_num_leapfrogs = max(hmc_num_leapfrogs, 1) # Set minimum hmc_num_leapfrogs to 1

  check_positive_integer(nuts_max_depth, "nuts_max_depth")
  nuts_max_depth = max(nuts_max_depth, 1) # Set minimum nuts_max_depth to 1

  #Check na_action -------------------------------------------------------------
  na_action_input = na_action
  na_action = try(match.arg(na_action), silent = TRUE)
  if(inherits(na_action, what = "try-error"))
    stop(paste0("The na_action argument should equal listwise or impute, not ",
                na_action_input,
                "."))

  #Check display_progress ------------------------------------------------------
  progress_type = progress_type_from_display_progress(display_progress)

  #Format the data input -------------------------------------------------------
  data = reformat_data(x = x,
                       na_action = na_action,
                       variable_bool = variable_bool,
                       baseline_category = baseline_category)
  x = data$x
  num_categories = data$num_categories
  missing_index = data$missing_index
  na_impute = data$na_impute
  baseline_category = data$baseline_category

  num_variables = ncol(x)
  num_interactions = num_variables * (num_variables - 1) / 2
  num_thresholds = sum(num_categories)

  # Starting value of model matrix ---------------------------------------------
  indicator = matrix(1,
                 nrow = num_variables,
                 ncol = num_variables)


  #Starting values of interactions and thresholds (posterior mode) -------------
  interactions = matrix(0, nrow = num_variables, ncol = num_variables)
  thresholds = matrix(0, nrow = num_variables, ncol = max(num_categories))

  #Precompute the number of observations per category for each variable --------
  counts_per_category = matrix(0,
                     nrow = max(num_categories) + 1,
                     ncol = num_variables)
  for(variable in 1:num_variables) {
    for(category in 0:num_categories[variable]) {
      counts_per_category[category + 1, variable] = sum(x[, variable] == category)
    }
  }

  #Precompute the sufficient statistics for the two Blume-Capel parameters -----
  blume_capel_stats = matrix(0, nrow = 2, ncol = num_variables)
  if(any(!variable_bool)) {
    # Ordinal (variable_bool == TRUE) or Blume-Capel (variable_bool == FALSE)
    bc_vars = which(!variable_bool)
    for(i in bc_vars) {
      blume_capel_stats[1, i] = sum(x[, i] - baseline_category[i])
      blume_capel_stats[2, i] = sum((x[, i] - baseline_category[i]) ^ 2)
      x[, i] = x[, i] - baseline_category[i]
    }
  }
  pairwise_stats = t(x) %*% x

  # Index matrix used in the c++ functions  ------------------------------------
  interaction_index_matrix = matrix(0,
                 nrow = num_variables * (num_variables - 1) / 2,
                 ncol = 3)
  cntr = 0
  for(variable1 in 1:(num_variables - 1)) {
    for(variable2 in (variable1 + 1):num_variables) {
      cntr =  cntr + 1
      interaction_index_matrix[cntr, 1] = cntr - 1
      interaction_index_matrix[cntr, 2] = variable1 - 1
      interaction_index_matrix[cntr, 3] = variable2 - 1
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

  #Setting the seed
  if (missing(seed) || is.null(seed)) {
    # Draw a random seed if none provided
    seed = sample.int(.Machine$integer.max, 1)
  }

  if (!is.numeric(seed) || length(seed) != 1 || is.na(seed) || seed < 0) {
    stop("Argument 'seed' must be a single non-negative integer.")
  }

  seed <- as.integer(seed)

  out = run_bgm_parallel(
    observations = x, num_categories = num_categories,
    pairwise_scale = pairwise_scale, edge_prior = edge_prior,
    inclusion_probability = inclusion_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    dirichlet_alpha = dirichlet_alpha, lambda = lambda,
    interaction_index_matrix = interaction_index_matrix, iter = iter,
    warmup = warmup, counts_per_category = counts_per_category,
    blume_capel_stats = blume_capel_stats,
    main_alpha = main_alpha, main_beta = main_beta,
    na_impute = na_impute, missing_index = missing_index,
    is_ordinal_variable = variable_bool,
    baseline_category = baseline_category, edge_selection = edge_selection,
    update_method = update_method,
    pairwise_effect_indices = pairwise_effect_indices,
    target_accept = target_accept, pairwise_stats = pairwise_stats,
    hmc_num_leapfrogs = hmc_num_leapfrogs, nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix, num_chains = chains,
    nThreads = cores, seed = seed, progress_type = progress_type
  )

  userInterrupt = any(vapply(out, FUN = `[[`, FUN.VALUE = logical(1L), "userInterrupt"))
  if (userInterrupt) {
    warning("Stopped sampling after user interrupt, results are likely uninterpretable.")
    # Try to prepare output, but catch any errors
    output <- tryCatch(
      prepare_output_bgm(
        out = out, x = x, num_categories = num_categories, iter = iter,
        data_columnnames = if (is.null(colnames(x))) paste0("Variable ", seq_len(ncol(x))) else colnames(x),
        is_ordinal_variable = variable_bool,
        warmup = warmup, pairwise_scale = pairwise_scale,
        main_alpha = main_alpha, main_beta = main_beta,
        na_action = na_action, na_impute = na_impute,
        edge_selection = edge_selection, edge_prior = edge_prior, inclusion_probability = inclusion_probability,
        beta_bernoulli_alpha = beta_bernoulli_alpha, beta_bernoulli_beta = beta_bernoulli_beta,
        dirichlet_alpha = dirichlet_alpha, lambda = lambda,
        variable_type = variable_type,
        update_method = update_method,
        target_accept = target_accept,
        hmc_num_leapfrogs = hmc_num_leapfrogs,
        nuts_max_depth = nuts_max_depth,
        learn_mass_matrix = learn_mass_matrix,
        num_chains = chains
      ),
      error = function(e) {
        list(partial = out, error = conditionMessage(e))
      },
      warning = function(w) {
        # still salvage what we can
        list(partial = out, warning = conditionMessage(w))
      }
    )
    return(output)
  }

  # Main output handler in the wrapper function
  output = prepare_output_bgm (
    out = out, x = x, num_categories = num_categories, iter = iter,
    data_columnnames = if (is.null(colnames(x))) paste0("Variable ", seq_len(ncol(x))) else colnames(x),
    is_ordinal_variable = variable_bool,
    warmup = warmup, pairwise_scale = pairwise_scale,
    main_alpha = main_alpha, main_beta = main_beta,
    na_action = na_action, na_impute = na_impute,
    edge_selection = edge_selection, edge_prior = edge_prior, inclusion_probability = inclusion_probability,
    beta_bernoulli_alpha = beta_bernoulli_alpha,
    beta_bernoulli_beta = beta_bernoulli_beta,
    dirichlet_alpha = dirichlet_alpha, lambda = lambda,
    variable_type = variable_type,
    update_method = update_method,
    target_accept = target_accept,
    hmc_num_leapfrogs = hmc_num_leapfrogs,
    nuts_max_depth = nuts_max_depth,
    learn_mass_matrix = learn_mass_matrix,
    num_chains = chains
  )

  if (update_method == "nuts") {
    nuts_diag = summarize_nuts_diagnostics(out, nuts_max_depth = nuts_max_depth)
    output$nuts_diag = nuts_diag
  }

  # -------------------------------------------------------------------
  # TODO: REMOVE after easybgm >= 0.2.2 is on CRAN
  # Compatibility shim for easybgm <= 0.2.1
  # -------------------------------------------------------------------
  if ("easybgm" %in% loadedNamespaces()) {
    ebgm_version <- utils::packageVersion("easybgm")
    if (ebgm_version <= "0.2.1") {
      warning("bgms is running in compatibility mode for easybgm (<= 0.2.1). ",
              "This will be removed once easybgm >= 0.2.2 is on CRAN.")

      # Add legacy variables to output
      output$arguments$save <- TRUE
      if (edge_selection) {
        output$indicator <- extract_indicators(output)
      }
      output$interactions <- extract_pairwise_interactions(output)
      output$thresholds   <- extract_category_thresholds(output)
    }
  }

  return(output)
}