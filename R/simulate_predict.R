# ==============================================================================
#   MRF Simulation and Prediction Functions
#
#   This file contains:
#     - simulate_mrf(): Standalone MRF data simulation
#     - mrfSampler(): Deprecated wrapper for simulate_mrf()
#     - simulate.bgms(): S3 method for simulating from fitted models
#     - predict.bgms(): S3 method for conditional probability prediction
# ==============================================================================


#' Simulate Observations from a Markov Random Field
#'
#' @description
#' `simulate_mrf()` generates observations from a Markov Random Field using
#' user-specified parameters. For ordinal and Blume-Capel variables, observations
#' are generated via Gibbs sampling. For continuous variables (Gaussian graphical
#' model), observations are drawn directly from the multivariate normal
#' distribution implied by the precision matrix.
#'
#' @details
#' \strong{Ordinal / Blume-Capel variables:}
#' The Gibbs sampler is initiated with random values from the response options,
#' after which it proceeds by simulating states for each variable from its full
#' conditional distribution given the other variable states.
#'
#' \strong{Continuous variables (GGM):}
#' Observations are drawn from \eqn{N(\mu, \Omega^{-1})}{N(mu, Omega^{-1})}
#' where \eqn{\Omega}{Omega} is the precision matrix specified via
#' `pairwise` and \eqn{\mu}{mu} is the means vector specified via `main`.
#' No Gibbs sampling is needed; `iter` is ignored.
#'
#' There are two modeling options for the category thresholds. The default
#' option assumes that the category thresholds are free, except that the first
#' threshold is set to zero for identification. The user then only needs to
#' specify the thresholds for the remaining response categories. This option is
#' useful for any type of ordinal variable and gives the user the most freedom
#' in specifying their model.
#'
#' The Blume-Capel option is specifically designed for ordinal variables that
#' have a special type of baseline_category category, such as the neutral
#' category in a Likert scale. The Blume-Capel model specifies the following
#' quadratic model for the threshold parameters:
#' \deqn{\mu_{\text{c}} = \alpha \times (\text{c} - \text{r}) + \beta \times (\text{c} - \text{r})^2,}{{\mu_{\text{c}} = \alpha \times (\text{c} - \text{r}) + \beta \times (\text{c} - \text{r})^2,}}
#' where \eqn{\mu_{\text{c}}}{\mu_{\text{c}}} is the threshold for category c
#' (which now includes zero), \eqn{\alpha}{\alpha} offers a linear trend
#' across categories (increasing threshold values if
#' \eqn{\alpha > 0}{\alpha > 0} and decreasing threshold values if
#' \eqn{\alpha <0}{\alpha <0}), if \eqn{\beta < 0}{\beta < 0}, it offers an
#' increasing penalty for responding in a category further away from the
#' baseline_category category r, while \eqn{\beta > 0}{\beta > 0} suggests a
#' preference for responding in the baseline_category category.
#'
#' @param num_states The number of observations to be generated.
#'
#' @param num_variables The number of variables in the MRF.
#'
#' @param num_categories Either a positive integer or a vector of positive
#' integers of length \code{num_variables}. The number of response categories on top
#' of the base category: \code{num_categories = 1} generates binary states.
#' Only used for ordinal and Blume-Capel variables; ignored when
#' \code{variable_type = "continuous"}.
#'
#' @param pairwise A symmetric \code{num_variables} by \code{num_variables} matrix.
#' For ordinal and Blume-Capel variables, this contains the pairwise interaction
#' parameters; only the off-diagonal elements are used. For continuous variables,
#' this is the precision matrix \eqn{\Omega}{Omega} (including diagonal) and
#' must be positive definite.
#'
#' @param main For ordinal and Blume-Capel variables: a
#' \code{num_variables} by \code{max(num_categories)} matrix of category
#' thresholds. The elements in row \code{i} indicate the thresholds of
#' variable \code{i}. If \code{num_categories} is a vector, only the first
#' \code{num_categories[i]} elements are used in row \code{i}. If the Blume-Capel
#' model is used for the category thresholds for variable \code{i}, then row
#' \code{i} requires two values (details below); the first is
#' \eqn{\alpha}{\alpha}, the linear contribution of the Blume-Capel model and
#' the second is \eqn{\beta}{\beta}, the quadratic contribution.
#' For continuous variables: a numeric vector of length \code{num_variables}
#' containing the means \eqn{\mu}{mu} for each variable. Defaults to zeros
#' if not supplied or if all values are zero.
#'
#' @param variable_type What kind of variables are simulated? Can be a single
#' character string specifying the variable type of all \code{p} variables at
#' once or a vector of character strings of length \code{p} specifying the type
#' for each variable separately. Currently, bgm supports ``ordinal'',
#' ``blume-capel'', and ``continuous''. Binary variables are automatically
#' treated as ``ordinal''. Ordinal and Blume-Capel variables can be mixed
#' freely, but continuous variables cannot be mixed with ordinal or Blume-Capel
#' variables. When \code{variable_type = "continuous"}, the function simulates
#' from a Gaussian graphical model.
#' Defaults to \code{variable_type = "ordinal"}.
#'
#' @param baseline_category An integer vector of length \code{num_variables} specifying the
#' baseline_category category that is used for the Blume-Capel model (details below).
#' Can be any integer value between \code{0} and \code{num_categories} (or
#' \code{num_categories[i]}).
#'
#' @param iter The number of iterations used by the Gibbs sampler
#' (ordinal/Blume-Capel variables only). The function provides the last state
#' of the Gibbs sampler as output. Ignored for continuous variables.
#' By default set to \code{1e3}.
#'
#' @param seed Optional integer seed for reproducibility. If \code{NULL},
#' a seed is generated from R's random number generator (so \code{set.seed()}
#' can be used before calling this function).
#'
#' @return A \code{num_states} by \code{num_variables} matrix of simulated
#' observations. For ordinal/Blume-Capel variables, entries are non-negative
#' integers. For continuous variables, entries are real-valued.
#'
#' @examples
#' # Generate responses from a network of five binary and ordinal variables.
#' num_variables = 5
#' num_categories = sample(1:5, size = num_variables, replace = TRUE)
#'
#' Pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
#' Pairwise[2, 1] = Pairwise[4, 1] = Pairwise[3, 2] =
#'   Pairwise[5, 2] = Pairwise[5, 4] = .25
#' Pairwise = Pairwise + t(Pairwise)
#' Main = matrix(0, nrow = num_variables, ncol = max(num_categories))
#'
#' x = simulate_mrf(
#'   num_states = 1e3,
#'   num_variables = num_variables,
#'   num_categories = num_categories,
#'   pairwise = Pairwise,
#'   main = Main
#' )
#'
#' # Generate responses from a network of 2 ordinal and 3 Blume-Capel variables.
#' num_variables = 5
#' num_categories = 4
#'
#' Pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
#' Pairwise[2, 1] = Pairwise[4, 1] = Pairwise[3, 2] =
#'   Pairwise[5, 2] = Pairwise[5, 4] = .25
#' Pairwise = Pairwise + t(Pairwise)
#'
#' Main = matrix(NA, num_variables, num_categories)
#' Main[, 1] = -1
#' Main[, 2] = -1
#' Main[3, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#' Main[5, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
#'
#' x = simulate_mrf(
#'   num_states = 1e3,
#'   num_variables = num_variables,
#'   num_categories = num_categories,
#'   pairwise = Pairwise,
#'   main = Main,
#'   variable_type = c("b", "b", "o", "b", "o"),
#'   baseline_category = 2
#' )
#'
#' # Generate responses from a Gaussian graphical model (GGM) with 4 variables.
#' num_variables = 4
#'
#' # Precision matrix (symmetric, positive definite)
#' Omega = diag(c(1, 1.2, 0.8, 1.5))
#' Omega[2, 1] = Omega[1, 2] = 0.3
#' Omega[3, 1] = Omega[1, 3] = 0.3
#' Omega[4, 2] = Omega[2, 4] = -0.2
#'
#' x = simulate_mrf(
#'   num_states = 500,
#'   num_variables = num_variables,
#'   pairwise = Omega,
#'   variable_type = "continuous"
#' )
#'
#' @seealso \code{\link{simulate.bgms}} for simulating from a fitted model.
#' @family prediction
#'
#' @export
simulate_mrf = function(num_states,
                        num_variables,
                        num_categories,
                        pairwise,
                        main,
                        variable_type = "ordinal",
                        baseline_category,
                        iter = 1e3,
                        seed = NULL) {
  # Check num_states, num_variables ---------------------------------------------
  check_positive_integer(num_states, "num_states")
  check_positive_integer(num_variables, "num_variables")

  # Check variable specification -----------------------------------------------
  vt = validate_variable_types(
    variable_type = variable_type,
    num_variables = num_variables,
    allow_continuous = TRUE,
    caller = "simulate_mrf"
  )
  variable_type = vt$variable_type
  is_continuous = vt$is_continuous

  # Blume-Capel binary guard (simulation-specific: uses num_categories)
  if(any(variable_type == "blume-capel")) {
    bc_binary = variable_type == "blume-capel" & num_categories < 2
    if(any(bc_binary)) {
      stop(paste0(
        "The Blume-Capel model only works for ordinal variables with more than two \n",
        "response options. But variables ",
        paste(which(bc_binary), collapse = ", "),
        " are binary variables."
      ))
    }
  }

  # ===========================================================================
  #   Continuous (GGM) path — direct multivariate normal sampling
  # ===========================================================================
  if(is_continuous) {
    # Check pairwise (full precision matrix, including diagonal)
    if(!inherits(pairwise, what = "matrix")) {
      pairwise = as.matrix(pairwise)
    }
    # NAs indicate excluded edges (zero precision)
    pairwise[is.na(pairwise)] = 0
    if(!isSymmetric(pairwise)) {
      stop("The matrix 'pairwise' needs to be symmetric.")
    }
    if(nrow(pairwise) != num_variables) {
      stop("The matrix 'pairwise' needs to have 'num_variables' rows and columns.")
    }
    if(any(diag(pairwise) <= 0)) {
      stop("The diagonal of the precision matrix 'pairwise' must be positive.")
    }

    precision = pairwise

    # Handle means (from 'main', default to zero)
    if(missing(main)) {
      means = rep(0, num_variables)
    } else {
      means = as.numeric(main)
      if(length(means) != num_variables) {
        stop(paste0(
          "'main' must have ", num_variables,
          " elements (one mean per variable), but has ",
          length(means), "."
        ))
      }
      if(any(!is.finite(means))) {
        stop("All elements of 'main' must be finite.")
      }
    }

    # Handle seed
    seed = check_seed(seed)

    x = sample_ggm_direct(
      num_states = num_states,
      precision = precision,
      means = means,
      seed = seed
    )

    return(x)
  }

  # ===========================================================================
  #   Ordinal / Blume-Capel path — Gibbs sampling
  # ===========================================================================
  check_positive_integer(iter, "iter")

  # Check num_categories --------------------------------------------------------
  if(length(num_categories) == 1) {
    if(num_categories <= 0 ||
      abs(num_categories - round(num_categories)) > .Machine$double.eps) {
      stop("``num_categories'' needs be a (vector of) positive integer(s).")
    }
    num_categories = rep(num_categories, num_variables)
  } else {
    for(variable in 1:num_variables) {
      if(num_categories[variable] <= 0 ||
        abs(num_categories[variable] - round(num_categories[variable])) >
          .Machine$double.eps) {
        stop(paste("For variable", variable, "``num_categories'' was not a positive integer."))
      }
    }
  }

  # Check the baseline_category for Blume-Capel variables ---------------------
  if(any(variable_type == "blume-capel")) {
    if(length(baseline_category) == 1) {
      baseline_category = rep(baseline_category, num_variables)
    }
    if(any(baseline_category < 0) || any(abs(baseline_category - round(baseline_category)) > .Machine$double.eps)) {
      stop(paste0(
        "For variables ",
        which(baseline_category < 0),
        " ``baseline_category'' was either negative or not integer."
      ))
    }
    if(any(baseline_category - num_categories > 0)) {
      stop(paste0(
        "For variables ",
        which(baseline_category - num_categories > 0),
        " the ``baseline_category'' category was larger than the maximum category value."
      ))
    }
  }

  # Check pairwise ---------------------------------------------------------
  if(!inherits(pairwise, what = "matrix")) {
    pairwise = as.matrix(pairwise)
  }
  if(!isSymmetric(pairwise)) {
    stop("The matrix ``pairwise'' needs to be symmetric.")
  }
  if(nrow(pairwise) != num_variables) {
    stop("The matrix ``pairwise'' needs to have ``num_variables'' rows and columns.")
  }

  # Check the threshold values -------------------------------------------------
  if(!inherits(main, what = "matrix")) {
    if(max(num_categories) == 1) {
      if(length(main) == num_variables) {
        main = matrix(main, ncol = 1)
      } else {
        stop(paste0(
          "The matrix ``main'' has ",
          length(main),
          " elements, but requires",
          num_variables,
          "."
        ))
      }
    } else {
      stop("``main'' needs to be a matrix.")
    }
  }

  if(nrow(main) != num_variables) {
    stop("The matrix ``main'' needs to be have ``num_variables'' rows.")
  }

  for(variable in 1:num_variables) {
    if(variable_type[variable] != "blume-capel") {
      if(anyNA(main[variable, 1:num_categories[variable]])) {
        tmp = which(is.na(main[variable, 1:num_categories[variable]]))

        string = paste(tmp, sep = ",")

        stop(paste0(
          "The matrix ``main'' contains NA(s) for variable ",
          variable,
          " in category \n",
          "(categories) ",
          paste(which(is.na(main[variable, 1:num_categories[variable]])), collapse = ", "),
          ", where a numeric value is needed."
        ))
      }
      if(ncol(main) > num_categories[variable]) {
        if(!anyNA(main[variable, (num_categories[variable] + 1):ncol(main)])) {
          warning(paste0(
            "The matrix ``main'' contains numeric values for variable ",
            variable,
            " for category \n",
            "(categories, i.e., columns) exceding the maximum of ",
            num_categories[variable],
            ". These values will \n",
            "be ignored."
          ))
        }
      }
    } else {
      if(anyNA(main[variable, 1:2])) {
        stop(paste0(
          "The Blume-Capel model is chosen for the category thresholds of variable ",
          variable,
          ". \n",
          "This model has two parameters that need to be placed in columns 1 and 2, row \n",
          variable,
          ", of the ``main'' input matrix. Currently, there are NA(s) in these \n",
          "entries, where a numeric value is needed."
        ))
      }
      if(ncol(main) > 2) {
        if(!anyNA(main[variable, 3:ncol(main)])) {
          warning(paste0(
            "The Blume-Capel model is chosen for the category thresholds of variable ",
            variable,
            ". \n",
            "This model has two parameters that need to be placed in columns 1 and 2, row \n",
            variable,
            ", of the ``main'' input matrix. However, there are numeric values \n",
            "in higher categories. These values will be ignored."
          ))
        }
      }
    }
  }

  for(variable in 1:num_variables) {
    if(variable_type[variable] != "blume-capel") {
      for(category in 1:num_categories[variable]) {
        if(!is.finite(main[variable, category])) {
          stop(paste(
            "The threshold parameter for variable", variable, "and category",
            category, "is NA or not finite."
          ))
        }
      }
    } else {
      if(!is.finite(main[variable, 1])) {
        stop(paste0(
          "The alpha parameter for the Blume-Capel model for variable ",
          variable,
          " is NA \n",
          " or not finite."
        ))
      }
      if(!is.finite(main[variable, 2])) {
        stop(paste0(
          "The beta parameter for the Blume-Capel model for variable",
          variable,
          "is NA \n",
          " or not finite."
        ))
      }
    }
  }

  # Handle seed ----------------------------------------------------------------
  seed = check_seed(seed)

  # The Gibbs sampler ----------------------------------------------------------
  if(!any(variable_type == "blume-capel")) {
    x = sample_omrf_gibbs(
      num_states = num_states,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      iter = iter,
      seed = seed
    )
  } else {
    x = sample_bcomrf_gibbs(
      num_states = num_states,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      variable_type_r = variable_type,
      baseline_category = baseline_category,
      iter = iter,
      seed = seed
    )
  }

  return(x)
}


# ==============================================================================
#   mrfSampler() - Deprecated Wrapper for simulate_mrf()
# ==============================================================================

#' Sample observations from the ordinal MRF
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `mrfSampler()` was renamed to [simulate_mrf()] as of bgms 0.1.6.3 to
#' follow the package's naming conventions.
#'
#' @inheritParams simulate_mrf
#'
#' @return A matrix of simulated observations (see [simulate_mrf()]).
#'
#' @seealso [simulate_mrf()] for the current function.
#'
#' @keywords internal
#' @export
mrfSampler = function(num_states,
                      num_variables,
                      num_categories,
                      pairwise,
                      main,
                      variable_type = "ordinal",
                      baseline_category,
                      iter = 1e3,
                      seed = NULL) {
  lifecycle::deprecate_warn("0.1.6.3", "mrfSampler()", "simulate_mrf()")

  simulate_mrf(
    num_states = num_states,
    num_variables = num_variables,
    num_categories = num_categories,
    pairwise = pairwise,
    main = main,
    variable_type = variable_type,
    baseline_category = baseline_category,
    iter = iter,
    seed = seed
  )
}


# ==============================================================================
#   simulate.bgms() - S3 Method for Simulating from Fitted Models
# ==============================================================================

#' Simulate Data from a Fitted bgms Model
#'
#' @description
#' Generates new observations from the Markov Random Field model using the
#' estimated parameters from a fitted \code{bgms} object.
#'
#' @param object An object of class \code{bgms}.
#' @param nsim Number of observations to simulate. Default: \code{500}.
#' @param seed Optional random seed for reproducibility.
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters (faster,
#'       single simulation).}
#'     \item{\code{"posterior-sample"}}{Sample from posterior draws, producing
#'       one dataset per draw (accounts for parameter uncertainty). This method
#'       uses parallel processing when \code{cores > 1}.}
#'   }
#' @param ndraws Number of posterior draws to use when
#'   \code{method = "posterior-sample"}. If \code{NULL}, uses all available draws.
#' @param iter Number of Gibbs iterations for equilibration before collecting
#'   samples. Default: \code{1000}.
#' @param cores Number of CPU cores for parallel execution when
#'   \code{method = "posterior-sample"}.
#'   Default: \code{parallel::detectCores()}.
#' @param display_progress Character string specifying the type of progress bar.
#'   Options: \code{"per-chain"}, \code{"total"}, \code{"none"}.
#'   Default: \code{"per-chain"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' If \code{method = "posterior-mean"}: A matrix with \code{nsim} rows and
#' \code{p} columns containing simulated observations.
#'
#' If \code{method = "posterior-sample"}: A list of matrices, one per posterior
#' draw, each with \code{nsim} rows and \code{p} columns.
#'
#' @details
#' This function uses the estimated interaction and threshold parameters to
#' generate new data via Gibbs sampling. When \code{method = "posterior-sample"},
#' parameter uncertainty is propagated to the simulated data by using different
#' posterior draws. Parallel processing is available for this method via the
#' \code{cores} argument.
#'
#' @seealso \code{\link{predict.bgms}} for computing conditional probabilities,
#'   \code{\link{simulate_mrf}} for simulation with user-specified parameters.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a model
#' fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#'
#' # Simulate 100 new observations using posterior means
#' new_data = simulate(fit, nsim = 100)
#'
#' # Simulate with parameter uncertainty (10 datasets)
#' new_data_list = simulate(fit, nsim = 100, method = "posterior-sample", ndraws = 10)
#'
#' # Use parallel processing for faster simulation
#' new_data_list = simulate(fit,
#'   nsim = 100, method = "posterior-sample",
#'   ndraws = 100, cores = 2
#' )
#' }
#'
#' @importFrom stats simulate
#' @export
simulate.bgms = function(object,
                         nsim = 500,
                         seed = NULL,
                         method = c("posterior-mean", "posterior-sample"),
                         ndraws = NULL,
                         iter = 1000,
                         cores = parallel::detectCores(),
                         display_progress = c("per-chain", "total", "none"),
                         ...) {
  method = match.arg(method)
  progress_type = progress_type_from_display_progress(display_progress)

  # Validate cores
  check_positive_integer(cores, "cores")
  cores = as.integer(cores)

  # Setting the seed
  seed = check_seed(seed)

  # Extract model information

  arguments = extract_arguments(object)
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  variable_type = arguments$variable_type
  data_columnnames = arguments$data_columnnames

  # Handle variable_type

  if(length(variable_type) == 1) {
    variable_type = rep(variable_type, num_variables)
  }

  # Get baseline_category (for Blume-Capel variables)
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  # ============================================================================
  #   GGM (continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_continuous)) {
    return(simulate_bgms_ggm(
      object = object,
      nsim = nsim,
      seed = seed,
      method = method,
      ndraws = ndraws,
      num_variables = num_variables,
      data_columnnames = data_columnnames,
      cores = cores,
      progress_type = progress_type
    ))
  }

  # ============================================================================
  #   Mixed MRF (discrete + continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_mixed)) {
    return(simulate_bgms_mixed(
      object = object,
      nsim = nsim,
      seed = seed,
      method = method,
      ndraws = ndraws,
      arguments = arguments,
      iter = iter,
      cores = cores,
      progress_type = progress_type
    ))
  }

  # ============================================================================
  #   OMRF (ordinal / Blume-Capel) path
  # ============================================================================

  if(method == "posterior-mean") {
    # Use posterior mean parameters
    pairwise = object$posterior_mean_pairwise
    main = object$posterior_mean_main

    # Set R's RNG for simulate_mrf
    if(!is.null(seed)) set.seed(seed)

    # Call simulate_mrf
    result = simulate_mrf(
      num_states = nsim,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      variable_type = variable_type,
      baseline_category = baseline_category,
      iter = iter
    )

    colnames(result) = data_columnnames
    return(result)
  } else {
    # Use posterior samples with parallel processing
    pairwise_samples = do.call(rbind, object$raw_samples$pairwise)
    main_samples = do.call(rbind, object$raw_samples$main)

    total_draws = nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws = total_draws
    }
    ndraws = min(ndraws, total_draws)

    # Sample which draws to use
    if(!is.null(seed)) set.seed(seed)
    draw_indices = sample.int(total_draws, ndraws)

    # Call parallel C++ function
    results = run_simulation_parallel(
      pairwise_samples = pairwise_samples,
      main_samples = main_samples,
      draw_indices = as.integer(draw_indices),
      num_states = as.integer(nsim),
      num_variables = as.integer(num_variables),
      num_categories = as.integer(num_categories),
      variable_type_r = variable_type,
      baseline_category = as.integer(baseline_category),
      iter = as.integer(iter),
      nThreads = cores,
      seed = seed,
      progress_type = progress_type
    )

    # Add column names
    for(i in seq_along(results)) {
      colnames(results[[i]]) = data_columnnames
    }

    return(results)
  }
}


# ==============================================================================
#   simulate.bgmCompare() - S3 Method for Simulating from Group-Comparison Models
# ==============================================================================

#' Simulate Data from a Fitted bgmCompare Model
#'
#' @description
#' Generates new observations from the Markov Random Field model for a
#' specified group using the estimated parameters from a fitted
#' \code{bgmCompare} object.
#'
#' @param object An object of class \code{bgmCompare}.
#' @param nsim Number of observations to simulate. Default: \code{500}.
#' @param seed Optional random seed for reproducibility.
#' @param group Integer specifying which group to simulate from (1 to
#'   number of groups). Required argument.
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters (faster,
#'       single simulation).}
#'   }
#' @param iter Number of Gibbs iterations for equilibration before collecting
#'   samples. Default: \code{1000}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix with \code{nsim} rows and \code{p} columns containing
#'   simulated observations for the specified group.
#'
#' @details
#' Group-specific parameters are obtained by applying the projection matrix
#' to convert baseline parameters and differences into group-level estimates:
#' \code{group_param = baseline + projection[group, ] \%*\% differences}.
#'
#' The function then uses these group-specific interaction and threshold
#' parameters to generate new data via Gibbs sampling.
#'
#' @seealso \code{\link{simulate.bgms}} for simulating from single-group models,
#'   \code{\link{predict.bgmCompare}} for computing conditional probabilities.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a comparison model
#' x = Boredom[Boredom$language == "fr", 2:6]
#' y = Boredom[Boredom$language != "fr", 2:6]
#' fit = bgmCompare(x, y, chains = 2)
#'
#' # Simulate 100 observations from group 1
#' new_data_g1 = simulate(fit, nsim = 100, group = 1)
#'
#' # Simulate 100 observations from group 2
#' new_data_g2 = simulate(fit, nsim = 100, group = 2)
#' }
#'
#' @export
simulate.bgmCompare = function(object,
                               nsim = 500,
                               seed = NULL,
                               group,
                               method = c("posterior-mean"),
                               iter = 1000,
                               ...) {
  method = match.arg(method)

  # Validate group argument
  if(missing(group)) {
    stop("Argument 'group' is required. Specify which group to simulate from (1 to num_groups).")
  }

  arguments = extract_arguments(object)
  num_groups = arguments$num_groups

  if(!is.numeric(group) || length(group) != 1 || is.na(group) ||
    group < 1 || group > num_groups) {
    stop(sprintf("Argument 'group' must be an integer between 1 and %d.", num_groups))
  }
  group = as.integer(group)

  # Setting the seed
  seed = check_seed(seed)

  # Extract model information
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal_variable
  data_columnnames = arguments$data_columnnames
  projection = arguments$projection # [num_groups x (num_groups-1)]

  # Determine variable_type from is_ordinal
  variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  # Get baseline_category (for Blume-Capel variables)
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  if(method == "posterior-mean") {
    # Extract group-specific parameters using projection
    group_params = extract_group_params(object)

    main_group = group_params$main_effects_groups[, group]
    pairwise_group = group_params$pairwise_effects_groups[, group]

    # Reconstruct threshold matrix
    max_cats = max(num_categories)
    main = matrix(NA_real_, nrow = num_variables, ncol = max_cats)

    pos = 1
    for(v in seq_len(num_variables)) {
      if(is_ordinal[v]) {
        k = num_categories[v]
        main[v, 1:k] = main_group[pos:(pos + k - 1)]
        pos = pos + k
      } else {
        # Blume-Capel: 2 parameters
        main[v, 1:2] = main_group[pos:(pos + 1)]
        pos = pos + 2
      }
    }

    # Reconstruct interaction matrix
    pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
    pairwise[lower.tri(pairwise)] = pairwise_group
    pairwise = pairwise + t(pairwise)

    # Set R's RNG for simulate_mrf
    set.seed(seed)

    # Call simulate_mrf
    result = simulate_mrf(
      num_states = nsim,
      num_variables = num_variables,
      num_categories = num_categories,
      pairwise = pairwise,
      main = main,
      variable_type = variable_type,
      baseline_category = baseline_category,
      iter = iter
    )

    colnames(result) = data_columnnames
    return(result)
  }
}


# ==============================================================================
#   predict.bgms() - S3 Method for Conditional Probability Prediction
# ==============================================================================

#' Predict Conditional Probabilities from a Fitted bgms Model
#'
#' @description
#' Computes conditional probability distributions for one or more variables
#' given the observed values of other variables in the data.
#'
#' @param object An object of class \code{bgms}.
#' @param newdata A matrix or data frame with \code{n} rows and \code{p} columns
#'   containing the observed data. Must have the same variables (columns) as
#'   the original data used to fit the model.
#' @param variables Which variables to predict. Can be:
#'   \itemize{
#'     \item A character vector of variable names
#'     \item An integer vector of column indices
#'     \item \code{NULL} (default) to predict all variables
#'   }
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{\code{"probabilities"}}{Return the full conditional probability
#'       distribution for each variable and observation.}
#'     \item{\code{"response"}}{Return the predicted category (mode of the
#'       conditional distribution).}
#'   }
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters.}
#'     \item{\code{"posterior-sample"}}{Average predictions over posterior draws.}
#'   }
#' @param ndraws Number of posterior draws to use when
#'   \code{method = "posterior-sample"}. If \code{NULL}, uses all available draws.
#' @param seed Optional random seed for reproducibility when
#'   \code{method = "posterior-sample"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' \strong{Ordinal models:}
#'
#' For \code{type = "probabilities"}: A named list with one element per
#' predicted variable. Each element is a matrix with \code{n} rows and
#' \code{num_categories + 1} columns containing \eqn{P(X_j = c | X_{-j})}{P(X_j = c | X_-j)} for each
#' observation and category.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing predicted categories.
#'
#' When \code{method = "posterior-sample"}, probabilities are averaged over
#' posterior draws, and an attribute \code{"sd"} is included containing the
#' standard deviation across draws.
#'
#' \strong{GGM (continuous) models:}
#'
#' For \code{type = "probabilities"}: A named list with one element per
#' predicted variable. Each element is a matrix with \code{n} rows and
#' 2 columns (\code{"mean"} and \code{"sd"}) containing the conditional
#' Gaussian parameters \eqn{E(X_j | X_{-j})}{E(X_j | X_{-j})} and
#' \eqn{\text{SD}(X_j | X_{-j})}{SD(X_j | X_{-j})}.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing conditional means.
#'
#' When \code{method = "posterior-sample"}, conditional parameters are
#' averaged over posterior draws, and an attribute \code{"sd"} is included.
#'
#' @details
#' For each observation, the function computes the conditional distribution
#' of the target variable(s) given the observed values of all other variables.
#' This is the same conditional distribution used internally by the Gibbs
#' sampler.
#'
#' For GGM (continuous) models, the conditional distribution of
#' \eqn{X_j | X_{-j}}{X_j | X_{-j}} is Gaussian with mean
#' \eqn{-\omega_{jj}^{-1} \sum_{k \neq j} \omega_{jk} x_k}{-omega_jj^{-1} sum_{k != j} omega_jk x_k}
#' and variance \eqn{\omega_{jj}^{-1}}{omega_jj^{-1}}, where \eqn{\Omega}{Omega}
#' is the precision matrix.
#'
#' @seealso \code{\link{simulate.bgms}} for generating new data from the model.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a model
#' fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#'
#' # Compute conditional probabilities for all variables
#' probs = predict(fit, newdata = Wenchuan[1:10, 1:5])
#'
#' # Predict the first variable only
#' probs_v1 = predict(fit, newdata = Wenchuan[1:10, 1:5], variables = 1)
#'
#' # Get predicted categories
#' pred_class = predict(fit, newdata = Wenchuan[1:10, 1:5], type = "response")
#' }
#'
#' @importFrom stats predict
#' @export
predict.bgms = function(object,
                        newdata,
                        variables = NULL,
                        type = c("probabilities", "response"),
                        method = c("posterior-mean", "posterior-sample"),
                        ndraws = NULL,
                        seed = NULL,
                        ...) {
  type = match.arg(type)
  method = match.arg(method)

  # Setting the seed (for R's RNG used by sample.int for draw selection)
  if(!is.null(seed)) {
    seed = check_seed(seed)
    set.seed(seed)
  }

  # Validate newdata
  if(missing(newdata)) {
    stop("Argument 'newdata' is required. Provide the data for which to compute predictions.")
  }

  if(!inherits(newdata, "matrix") && !inherits(newdata, "data.frame")) {
    stop("'newdata' must be a matrix or data frame.")
  }

  if(inherits(newdata, "data.frame")) {
    newdata = data.matrix(newdata)
  }

  # Extract model information
  arguments = extract_arguments(object)
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  variable_type = arguments$variable_type
  data_columnnames = arguments$data_columnnames

  # Validate dimensions

  if(ncol(newdata) != num_variables) {
    stop(paste0(
      "'newdata' must have ", num_variables, " columns (same as fitted model), ",
      "but has ", ncol(newdata), "."
    ))
  }

  # Handle variable_type
  if(length(variable_type) == 1) {
    variable_type = rep(variable_type, num_variables)
  }

  # Get baseline_category
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  # Convert variable_type to is_ordinal logical vector
  is_ordinal = variable_type != "blume-capel"

  # Determine which variables to predict
  if(is.null(variables)) {
    predict_vars = seq_len(num_variables)
  } else if(is.character(variables)) {
    predict_vars = match(variables, data_columnnames)
    if(anyNA(predict_vars)) {
      stop("Variable names not found: ", paste(variables[is.na(predict_vars)], collapse = ", "))
    }
  } else {
    predict_vars = as.integer(variables)
    if(any(predict_vars < 1 | predict_vars > num_variables)) {
      stop("Variable indices must be between 1 and ", num_variables)
    }
  }

  # ============================================================================
  #   GGM (continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_continuous)) {
    return(predict_bgms_ggm(
      object = object,
      newdata = newdata,
      predict_vars = predict_vars,
      data_columnnames = data_columnnames,
      num_variables = num_variables,
      type = type,
      method = method,
      ndraws = ndraws
    ))
  }

  # ============================================================================
  #   Mixed MRF (discrete + continuous) path
  # ============================================================================
  if(isTRUE(arguments$is_mixed)) {
    return(predict_bgms_mixed(
      object = object,
      newdata = newdata,
      predict_vars = predict_vars,
      arguments = arguments,
      type = type,
      method = method,
      ndraws = ndraws
    ))
  }

  # ============================================================================
  #   OMRF (ordinal) path
  # ============================================================================

  # Recode data to 0-based integers (matching what bgm() does)
  newdata_recoded = recode_data_for_prediction(newdata, num_categories, is_ordinal)

  if(method == "posterior-mean") {
    # Use posterior mean parameters
    pairwise = object$posterior_mean_pairwise
    main = object$posterior_mean_main

    probs = compute_conditional_probs(
      observations = newdata_recoded,
      predict_vars = predict_vars - 1L, # C++ uses 0-based indexing
      pairwise = pairwise,
      main = main,
      num_categories = num_categories,
      variable_type = variable_type,
      baseline_category = baseline_category
    )

    # Add names
    names(probs) = data_columnnames[predict_vars]
    for(v in seq_along(probs)) {
      var_idx = predict_vars[v]
      n_cats = num_categories[var_idx] + 1
      colnames(probs[[v]]) = paste0("cat_", 0:(n_cats - 1))
    }
  } else {
    # Use posterior samples
    pairwise_samples = do.call(rbind, object$raw_samples$pairwise)
    main_samples = do.call(rbind, object$raw_samples$main)

    total_draws = nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws = total_draws
    }
    ndraws = min(ndraws, total_draws)

    draw_indices = sample.int(total_draws, ndraws)

    # Collect probabilities from each draw
    all_probs = vector("list", ndraws)

    for(i in seq_len(ndraws)) {
      idx = draw_indices[i]

      # Reconstruct interaction matrix
      pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
      pairwise[lower.tri(pairwise)] = pairwise_samples[idx, ]
      pairwise = pairwise + t(pairwise)

      # Reconstruct threshold matrix
      main = reconstruct_main(
        main_samples[idx, ],
        num_variables,
        num_categories,
        variable_type
      )

      all_probs[[i]] = compute_conditional_probs(
        observations = newdata_recoded,
        predict_vars = predict_vars - 1L,
        pairwise = pairwise,
        main = main,
        num_categories = num_categories,
        variable_type = variable_type,
        baseline_category = baseline_category
      )
    }

    # Average over draws
    probs = vector("list", length(predict_vars))
    probs_sd = vector("list", length(predict_vars))
    names(probs) = data_columnnames[predict_vars]
    names(probs_sd) = data_columnnames[predict_vars]

    for(v in seq_along(predict_vars)) {
      # Stack probabilities from all draws: n x categories x ndraws
      var_probs = lapply(all_probs, `[[`, v)
      prob_array = array(unlist(var_probs),
        dim = c(nrow(newdata), ncol(var_probs[[1]]), ndraws)
      )

      probs[[v]] = apply(prob_array, c(1, 2), mean)
      probs_sd[[v]] = apply(prob_array, c(1, 2), sd)

      var_idx = predict_vars[v]
      n_cats = num_categories[var_idx] + 1
      colnames(probs[[v]]) = paste0("cat_", 0:(n_cats - 1))
      colnames(probs_sd[[v]]) = paste0("cat_", 0:(n_cats - 1))
    }

    attr(probs, "sd") = probs_sd
  }

  if(type == "response") {
    # Return predicted categories (mode)
    pred_matrix = sapply(probs, function(p) {
      apply(p, 1, which.max) - 1L # Convert to 0-based category
    })
    if(is.vector(pred_matrix)) {
      pred_matrix = matrix(pred_matrix, ncol = 1)
    }
    colnames(pred_matrix) = data_columnnames[predict_vars]
    return(pred_matrix)
  }

  return(probs)
}


# ==============================================================================
#   predict.bgmCompare() - S3 Method for Group-Comparison Models
# ==============================================================================

#' Predict Conditional Probabilities from a Fitted bgmCompare Model
#'
#' @description
#' Computes conditional probability distributions for one or more variables
#' given the observed values of other variables in the data, using
#' group-specific parameters from a \code{bgmCompare} model.
#'
#' @param object An object of class \code{bgmCompare}.
#' @param newdata A matrix or data frame with \code{n} rows and \code{p} columns
#'   containing the observed data. Must have the same variables (columns) as
#'   the original data used to fit the model.
#' @param group Integer specifying which group's parameters to use for
#'   prediction (1 to number of groups). Required argument.
#' @param variables Which variables to predict. Can be:
#'   \itemize{
#'     \item A character vector of variable names
#'     \item An integer vector of column indices
#'     \item \code{NULL} (default) to predict all variables
#'   }
#' @param type Character string specifying the type of prediction:
#'   \describe{
#'     \item{\code{"probabilities"}}{Return the full conditional probability
#'       distribution for each variable and observation.}
#'     \item{\code{"response"}}{Return the predicted category (mode of the
#'       conditional distribution).}
#'   }
#' @param method Character string specifying which parameter estimates to use:
#'   \describe{
#'     \item{\code{"posterior-mean"}}{Use posterior mean parameters.}
#'   }
#' @param ... Additional arguments (currently ignored).
#'
#' @return
#' For \code{type = "probabilities"}: A named list with one element per
#' predicted variable. Each element is a matrix with \code{n} rows and
#' \code{num_categories + 1} columns containing \eqn{P(X_j = c | X_{-j})}{P(X_j = c | X_-j)}
#' for each observation and category.
#'
#' For \code{type = "response"}: A matrix with \code{n} rows and
#' \code{length(variables)} columns containing predicted categories.
#'
#' @details
#' Group-specific parameters are obtained by applying the projection matrix
#' to convert baseline parameters and differences into group-level estimates.
#' The function then computes the conditional distribution of target variables
#' given the observed values of all other variables.
#'
#' @seealso \code{\link{predict.bgms}} for predicting from single-group models,
#'   \code{\link{simulate.bgmCompare}} for simulating from group-comparison models.
#' @family prediction
#'
#' @examples
#' \donttest{
#' # Fit a comparison model
#' x = Boredom[Boredom$language == "fr", 2:6]
#' y = Boredom[Boredom$language != "fr", 2:6]
#' fit = bgmCompare(x, y, chains = 2)
#'
#' # Predict conditional probabilities using group 1 parameters
#' probs_g1 = predict(fit, newdata = x[1:10, ], group = 1)
#'
#' # Predict responses using group 2 parameters
#' pred_g2 = predict(fit, newdata = y[1:10, ], group = 2, type = "response")
#' }
#'
#' @export
predict.bgmCompare = function(object,
                              newdata,
                              group,
                              variables = NULL,
                              type = c("probabilities", "response"),
                              method = c("posterior-mean"),
                              ...) {
  type = match.arg(type)
  method = match.arg(method)

  # Validate group argument
  if(missing(group)) {
    stop("Argument 'group' is required. Specify which group's parameters to use (1 to num_groups).")
  }

  arguments = extract_arguments(object)
  num_groups = arguments$num_groups

  if(!is.numeric(group) || length(group) != 1 || is.na(group) ||
    group < 1 || group > num_groups) {
    stop(sprintf("Argument 'group' must be an integer between 1 and %d.", num_groups))
  }
  group = as.integer(group)

  # Validate newdata
  if(missing(newdata)) {
    stop("Argument 'newdata' is required. Provide the data for which to compute predictions.")
  }

  if(!inherits(newdata, "matrix") && !inherits(newdata, "data.frame")) {
    stop("'newdata' must be a matrix or data frame.")
  }

  if(inherits(newdata, "data.frame")) {
    newdata = data.matrix(newdata)
  }

  # Extract model information
  num_variables = arguments$num_variables
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal_variable
  data_columnnames = arguments$data_columnnames
  projection = arguments$projection

  # Validate dimensions
  if(ncol(newdata) != num_variables) {
    stop(paste0(
      "'newdata' must have ", num_variables, " columns (same as fitted model), ",
      "but has ", ncol(newdata), "."
    ))
  }

  # Determine variable_type from is_ordinal
  variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  # Get baseline_category (for Blume-Capel variables)
  baseline_category = arguments$baseline_category
  if(is.null(baseline_category)) {
    baseline_category = rep(0L, num_variables)
  }

  # Determine which variables to predict
  if(is.null(variables)) {
    predict_vars = seq_len(num_variables)
  } else if(is.character(variables)) {
    predict_vars = match(variables, data_columnnames)
    if(anyNA(predict_vars)) {
      stop("Variable names not found: ", paste(variables[is.na(predict_vars)], collapse = ", "))
    }
  } else {
    predict_vars = as.integer(variables)
    if(any(predict_vars < 1 | predict_vars > num_variables)) {
      stop("Variable indices must be between 1 and ", num_variables)
    }
  }

  # Recode data to 0-based integers
  newdata_recoded = recode_data_for_prediction(newdata, num_categories, is_ordinal)

  if(method == "posterior-mean") {
    # Extract group-specific parameters using projection
    group_params = extract_group_params(object)

    main_group = group_params$main_effects_groups[, group]
    pairwise_group = group_params$pairwise_effects_groups[, group]

    # Reconstruct threshold matrix
    max_cats = max(num_categories)
    main = matrix(NA_real_, nrow = num_variables, ncol = max_cats)

    pos = 1
    for(v in seq_len(num_variables)) {
      if(is_ordinal[v]) {
        k = num_categories[v]
        main[v, 1:k] = main_group[pos:(pos + k - 1)]
        pos = pos + k
      } else {
        # Blume-Capel: 2 parameters
        main[v, 1:2] = main_group[pos:(pos + 1)]
        pos = pos + 2
      }
    }

    # Reconstruct interaction matrix
    pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
    pairwise[lower.tri(pairwise)] = pairwise_group
    pairwise = pairwise + t(pairwise)

    probs = compute_conditional_probs(
      observations = newdata_recoded,
      predict_vars = predict_vars - 1L, # C++ uses 0-based indexing
      pairwise = pairwise,
      main = main,
      num_categories = num_categories,
      variable_type = variable_type,
      baseline_category = baseline_category
    )

    # Add names
    names(probs) = data_columnnames[predict_vars]
    for(v in seq_along(probs)) {
      var_idx = predict_vars[v]
      n_cats = num_categories[var_idx] + 1
      colnames(probs[[v]]) = paste0("cat_", 0:(n_cats - 1))
    }
  }

  if(type == "response") {
    # Return predicted categories (mode)
    pred_matrix = sapply(probs, function(p) {
      apply(p, 1, which.max) - 1L # Convert to 0-based category
    })
    if(is.vector(pred_matrix)) {
      pred_matrix = matrix(pred_matrix, ncol = 1)
    }
    colnames(pred_matrix) = data_columnnames[predict_vars]
    return(pred_matrix)
  }

  return(probs)
}


# ==============================================================================
#   Helper Functions
# ==============================================================================

# Helper function to reconstruct threshold matrix from flat vector
reconstruct_main = function(main_vec, num_variables, num_categories, variable_type) {
  if(length(variable_type) == 1) {
    variable_type = rep(variable_type, num_variables)
  }

  max_cats = max(num_categories)
  main = matrix(NA, nrow = num_variables, ncol = max_cats)

  pos = 1
  for(v in seq_len(num_variables)) {
    if(variable_type[v] != "blume-capel") {
      k = num_categories[v]
      main[v, 1:k] = main_vec[pos:(pos + k - 1)]
      pos = pos + k
    } else {
      main[v, 1:2] = main_vec[pos:(pos + 1)]
      pos = pos + 2
    }
  }

  return(main)
}


# Helper function to recode data for prediction
recode_data_for_prediction = function(x, num_categories, is_ordinal) {
  x = as.matrix(x)
  num_variables = ncol(x)

  for(v in seq_len(num_variables)) {
    if(is_ordinal[v]) {
      # For ordinal variables, ensure values are in 0:num_categories[v]
      x[, v] = as.integer(x[, v])
      if(min(x[, v], na.rm = TRUE) > 0) {
        # Shift to 0-based if necessary
        x[, v] = x[, v] - min(x[, v], na.rm = TRUE)
      }
    }
  }

  return(x)
}


# ==============================================================================
#   GGM Prediction Helpers
# ==============================================================================

# Reconstruct the full precision matrix from posterior mean components.
#
# @param posterior_mean_pairwise p x p symmetric matrix with off-diagonal
#   precision elements (diagonal is zero).
# @param posterior_mean_main p x 1 matrix with diagonal precision elements
#   (column named "precision_diag").
#
# @return p x p precision matrix (Omega).
reconstruct_precision = function(posterior_mean_pairwise, posterior_mean_main) {
  omega = posterior_mean_pairwise
  # Excluded edges (NA) have zero precision
  omega[is.na(omega)] = 0
  diag(omega) = as.numeric(posterior_mean_main)
  return(omega)
}


# Reconstruct precision matrix from a single posterior draw.
#
# @param pairwise_vec Vector of p*(p-1)/2 off-diagonal precision elements
#   (lower-triangle order).
# @param main_vec Vector of p diagonal precision elements.
# @param p Number of variables.
#
# @return p x p precision matrix (Omega).
reconstruct_precision_from_draw = function(pairwise_vec, main_vec, p) {
  omega = matrix(0, nrow = p, ncol = p)
  omega[lower.tri(omega)] = pairwise_vec
  omega = omega + t(omega)
  diag(omega) = main_vec
  return(omega)
}


# GGM prediction implementation (called from predict.bgms).
#
# @param object Fitted bgms object (GGM).
# @param newdata n x p numeric matrix of observed continuous data.
# @param predict_vars Integer vector of 1-based variable indices to predict.
# @param data_columnnames Character vector of variable names.
# @param num_variables Number of variables p.
# @param type "probabilities" or "response".
# @param method "posterior-mean" or "posterior-sample".
# @param ndraws Number of posterior draws (NULL = all).
#
# @return See predict.bgms() documentation for GGM return format.
predict_bgms_ggm = function(object, newdata, predict_vars, data_columnnames,
                            num_variables,
                            type, method, ndraws) {
  # Center newdata by its own column means
  newdata_means = colMeans(newdata)
  newdata_centered = sweep(newdata, 2, newdata_means)

  if(method == "posterior-mean") {
    # Reconstruct precision matrix from posterior means
    omega = reconstruct_precision(
      object$posterior_mean_pairwise,
      object$posterior_mean_main
    )

    result = compute_conditional_ggm(
      observations = newdata_centered,
      predict_vars = predict_vars - 1L,
      precision = omega
    )

    # Add names and shift conditional means back to original scale
    names(result) = data_columnnames[predict_vars]
    for(v in seq_along(result)) {
      colnames(result[[v]]) = c("mean", "sd")
      result[[v]][, "mean"] = result[[v]][, "mean"] + newdata_means[predict_vars[v]]
    }
  } else {
    # Use posterior samples
    pairwise_samples = do.call(rbind, object$raw_samples$pairwise)
    main_samples = do.call(rbind, object$raw_samples$main)

    total_draws = nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws = total_draws
    }
    ndraws = min(ndraws, total_draws)

    draw_indices = sample.int(total_draws, ndraws)

    # Collect predictions from each draw
    all_preds = vector("list", ndraws)

    for(i in seq_len(ndraws)) {
      idx = draw_indices[i]

      omega = reconstruct_precision_from_draw(
        pairwise_vec = pairwise_samples[idx, ],
        main_vec = main_samples[idx, ],
        p = num_variables
      )

      preds = compute_conditional_ggm(
        observations = newdata_centered,
        predict_vars = predict_vars - 1L,
        precision = omega
      )

      # Shift conditional means back to original scale
      for(v in seq_along(predict_vars)) {
        preds[[v]][, 1] = preds[[v]][, 1] + newdata_means[predict_vars[v]]
      }

      all_preds[[i]] = preds
    }

    # Average over draws
    result = vector("list", length(predict_vars))
    result_sd = vector("list", length(predict_vars))
    names(result) = data_columnnames[predict_vars]
    names(result_sd) = data_columnnames[predict_vars]

    for(v in seq_along(predict_vars)) {
      # Stack predictions: n x 2 x ndraws
      var_preds = lapply(all_preds, `[[`, v)
      pred_array = array(unlist(var_preds),
        dim = c(nrow(newdata), 2, ndraws)
      )

      result[[v]] = apply(pred_array, c(1, 2), mean)
      result_sd[[v]] = apply(pred_array, c(1, 2), sd)

      colnames(result[[v]]) = c("mean", "sd")
      colnames(result_sd[[v]]) = c("mean", "sd")
    }

    attr(result, "sd") = result_sd
  }

  if(type == "response") {
    # Return conditional means
    pred_matrix = sapply(result, function(m) m[, "mean"])
    if(is.vector(pred_matrix)) {
      pred_matrix = matrix(pred_matrix, ncol = 1)
    }
    colnames(pred_matrix) = data_columnnames[predict_vars]
    return(pred_matrix)
  }

  return(result)
}


# ==============================================================================
#   GGM Simulation Helpers
# ==============================================================================

# GGM simulation implementation (called from simulate.bgms).
#
# @param object Fitted bgms object (GGM).
# @param nsim Number of observations to simulate.
# @param seed Random seed.
# @param method "posterior-mean" or "posterior-sample".
# @param ndraws Number of posterior draws (NULL = all).
# @param num_variables Number of variables p.
# @param data_columnnames Character vector of variable names.
# @param cores Number of parallel threads.
# @param progress_type Integer progress type (0/1/2).
#
# @return See simulate.bgms() documentation.
simulate_bgms_ggm = function(object, nsim, seed, method, ndraws,
                             num_variables, data_columnnames,
                             cores, progress_type) {
  if(method == "posterior-mean") {
    # Reconstruct precision matrix: inject diagonal from posterior_mean_main
    precision = reconstruct_precision(
      object$posterior_mean_pairwise,
      object$posterior_mean_main
    )

    # Call simulate_mrf with variable_type = "continuous"
    result = simulate_mrf(
      num_states = nsim,
      num_variables = num_variables,
      pairwise = precision,
      main = rep(0, num_variables),
      variable_type = "continuous",
      seed = seed
    )

    colnames(result) = data_columnnames
    return(result)
  } else {
    # Use posterior samples with parallel processing
    pairwise_samples = do.call(rbind, object$raw_samples$pairwise)
    main_samples = do.call(rbind, object$raw_samples$main)

    total_draws = nrow(pairwise_samples)
    if(is.null(ndraws)) {
      ndraws = total_draws
    }
    ndraws = min(ndraws, total_draws)

    # Sample which draws to use
    if(!is.null(seed)) set.seed(seed)
    draw_indices = sample.int(total_draws, ndraws)

    # Call parallel C++ function for GGM
    results = run_ggm_simulation_parallel(
      pairwise_samples = pairwise_samples,
      main_samples = main_samples,
      draw_indices = as.integer(draw_indices),
      num_states = as.integer(nsim),
      num_variables = as.integer(num_variables),
      means = rep(0, num_variables),
      nThreads = cores,
      seed = seed,
      progress_type = progress_type
    )

    # Add column names
    for(i in seq_along(results)) {
      colnames(results[[i]]) = data_columnnames
    }

    return(results)
  }
}


# ==============================================================================
#   Mixed MRF Simulation Helper
# ==============================================================================

# ------------------------------------------------------------------
# simulate_bgms_mixed
# ------------------------------------------------------------------
# Simulation implementation for mixed MRF models (called from simulate.bgms).
#
# @param object     Fitted bgms object (mixed MRF).
# @param nsim       Number of observations to simulate.
# @param seed       Random seed.
# @param method     "posterior-mean" or "posterior-sample".
# @param ndraws     Number of posterior draws (for posterior-sample).
# @param arguments  Output of extract_arguments().
# @param iter       Gibbs burn-in iterations.
# @param cores      Number of threads.
# @param progress_type  Progress bar type.
#
# Returns: matrix (posterior-mean) or list of matrices (posterior-sample).
# ------------------------------------------------------------------
simulate_bgms_mixed = function(object, nsim, seed, method, ndraws,
                               arguments, iter, cores, progress_type) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  data_columnnames = arguments$data_columnnames
  disc_idx = arguments$discrete_indices
  cont_idx = arguments$continuous_indices
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal
  baseline_category_disc = arguments$baseline_category

  disc_variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  bc = integer(p)
  for(s in seq_len(p)) {
    if(is_ordinal[s]) {
      bc[s] = 0L
    } else {
      bc[s] = as.integer(baseline_category_disc[s])
    }
  }

  if(method == "posterior-mean") {
    params = reconstruct_mixed_params_from_means(object, arguments)

    seed = check_seed(seed)

    result = sample_mixed_mrf_gibbs(
      num_states = as.integer(nsim),
      Kxx_r = params$Kxx,
      Kxy_r = params$Kxy,
      Kyy_r = params$Kyy,
      mux_r = params$mux,
      muy_r = params$muy,
      num_categories_r = as.integer(num_categories),
      variable_type_r = disc_variable_type,
      baseline_category_r = as.integer(bc),
      iter = as.integer(iter),
      seed = seed
    )

    out = combine_mixed_result(result, disc_idx, cont_idx, data_columnnames)
    return(out)

  } else {
    sample_info = split_mixed_raw_samples(object, arguments)

    total_draws = sample_info$total_draws
    if(is.null(ndraws)) ndraws = total_draws
    ndraws = min(ndraws, total_draws)

    if(!is.null(seed)) set.seed(seed)
    draw_indices = sample.int(total_draws, ndraws)

    results = run_mixed_simulation_parallel(
      mux_samples = sample_info$mux_samples,
      kxx_samples = sample_info$kxx_samples,
      muy_samples = sample_info$muy_samples,
      kyy_samples = sample_info$kyy_samples,
      kxy_samples = sample_info$kxy_samples,
      draw_indices = as.integer(draw_indices),
      num_states = as.integer(nsim),
      p = as.integer(p),
      q = as.integer(q),
      num_categories = as.integer(num_categories),
      variable_type_r = disc_variable_type,
      baseline_category = as.integer(bc),
      iter = as.integer(iter),
      nThreads = cores,
      seed = check_seed(seed),
      progress_type = progress_type
    )

    for(i in seq_along(results)) {
      results[[i]] = combine_mixed_result(
        results[[i]], disc_idx, cont_idx, data_columnnames
      )
    }

    return(results)
  }
}


# ==============================================================================
#   Mixed MRF Prediction Helper
# ==============================================================================

# ------------------------------------------------------------------
# predict_bgms_mixed
# ------------------------------------------------------------------
# Prediction implementation for mixed MRF models (called from predict.bgms).
#
# @param object       Fitted bgms object (mixed MRF).
# @param newdata      n x (p+q) matrix of observed data.
# @param predict_vars 1-based indices into the combined variable list.
# @param arguments    Output of extract_arguments().
# @param type         "probabilities" or "response".
# @param method       "posterior-mean" or "posterior-sample".
# @param ndraws       Number of posterior draws (for posterior-sample).
#
# Returns: Named list of prediction matrices.
# ------------------------------------------------------------------
predict_bgms_mixed = function(object, newdata, predict_vars, arguments,
                              type, method, ndraws) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  data_columnnames = arguments$data_columnnames
  disc_idx = arguments$discrete_indices
  cont_idx = arguments$continuous_indices
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal
  baseline_category_disc = arguments$baseline_category

  disc_variable_type = ifelse(is_ordinal, "ordinal", "blume-capel")

  bc = integer(p)
  for(s in seq_len(p)) {
    if(is_ordinal[s]) {
      bc[s] = 0L
    } else {
      bc[s] = as.integer(baseline_category_disc[s])
    }
  }

  # Split newdata into discrete and continuous parts
  x_data = as.matrix(newdata[, disc_idx, drop = FALSE])
  storage.mode(x_data) = "integer"
  y_data = as.matrix(newdata[, cont_idx, drop = FALSE])
  storage.mode(y_data) = "double"

  # Map user predict_vars (1-based, original order) to internal 0-based indices
  # Internal layout: [discrete_0..p-1, continuous_p..p+q-1]
  internal_predict_vars = integer(length(predict_vars))
  for(k in seq_along(predict_vars)) {
    orig_col = predict_vars[k]
    disc_match = match(orig_col, disc_idx)
    cont_match = match(orig_col, cont_idx)
    if(!is.na(disc_match)) {
      internal_predict_vars[k] = disc_match - 1L
    } else if(!is.na(cont_match)) {
      internal_predict_vars[k] = p + cont_match - 1L
    } else {
      stop("Variable index ", orig_col,
           " not found in discrete or continuous indices.")
    }
  }

  compute_one_draw = function(Kxx, Kxy, Kyy, mux, muy) {
    compute_conditional_mixed(
      x_observations = x_data,
      y_observations = y_data,
      predict_vars = as.integer(internal_predict_vars),
      Kxx = Kxx,
      Kxy = Kxy,
      Kyy = Kyy,
      mux = mux,
      muy = muy,
      num_categories = as.integer(num_categories),
      variable_type = disc_variable_type,
      baseline_category = as.integer(bc)
    )
  }

  if(method == "posterior-mean") {
    params = reconstruct_mixed_params_from_means(object, arguments)
    raw_result = compute_one_draw(
      params$Kxx, params$Kxy, params$Kyy, params$mux, params$muy
    )

    probs = format_mixed_predictions(
      raw_result, predict_vars, internal_predict_vars,
      p, num_categories, data_columnnames
    )
  } else {
    sample_info = split_mixed_raw_samples(object, arguments)
    total_draws = sample_info$total_draws
    if(is.null(ndraws)) ndraws = total_draws
    ndraws = min(ndraws, total_draws)

    draw_indices = sample.int(total_draws, ndraws)

    all_results = vector("list", ndraws)
    for(i in seq_len(ndraws)) {
      params_i = reconstruct_mixed_params_from_row(
        sample_info, draw_indices[i], p, q, num_categories, is_ordinal
      )
      all_results[[i]] = compute_one_draw(
        params_i$Kxx, params_i$Kxy, params_i$Kyy, params_i$mux, params_i$muy
      )
    }

    # Average predictions across draws
    num_pv = length(predict_vars)
    probs = vector("list", num_pv)
    probs_sd = vector("list", num_pv)
    names(probs) = data_columnnames[predict_vars]
    names(probs_sd) = data_columnnames[predict_vars]

    for(k in seq_len(num_pv)) {
      # Stack all draws into an array: n x ncol x ndraws
      var_preds = lapply(all_results, `[[`, k)
      pred_array = array(unlist(var_preds),
        dim = c(nrow(var_preds[[1]]), ncol(var_preds[[1]]), ndraws)
      )

      probs[[k]] = apply(pred_array, c(1, 2), mean)
      probs_sd[[k]] = apply(pred_array, c(1, 2), sd)
    }

    probs = format_mixed_predictions(
      probs, predict_vars, internal_predict_vars,
      p, num_categories, data_columnnames
    )
    names(probs_sd) = names(probs)
    attr(probs, "sd") = probs_sd
  }

  if(type == "response") {
    return(format_mixed_response(probs, predict_vars, internal_predict_vars,
                                 p, data_columnnames))
  }

  return(probs)
}


# ==============================================================================
#   Mixed MRF Internal Helpers
# ==============================================================================

# ------------------------------------------------------------------
# reconstruct_mixed_params_from_means
# ------------------------------------------------------------------
# Reconstruct Kxx, Kxy, Kyy, mux, muy from posterior mean summaries.
#
# @param object     Fitted bgms object (mixed MRF).
# @param arguments  Output of extract_arguments().
#
# Returns: list with Kxx, Kxy, Kyy, mux, muy.
# ------------------------------------------------------------------
reconstruct_mixed_params_from_means = function(object, arguments) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  disc_idx = arguments$discrete_indices
  cont_idx = arguments$continuous_indices
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal

  pmat = object$posterior_mean_pairwise

  Kxx = matrix(0, p, p)
  for(i in seq_len(p)) {
    for(j in seq_len(p)) {
      if(i != j) Kxx[i, j] = pmat[disc_idx[i], disc_idx[j]]
    }
  }

  Kxy = matrix(0, p, q)
  for(i in seq_len(p)) {
    for(j in seq_len(q)) {
      Kxy[i, j] = pmat[disc_idx[i], cont_idx[j]]
    }
  }

  Kyy = matrix(0, q, q)
  for(i in seq_len(q)) {
    for(j in seq_len(q)) {
      if(i != j) Kyy[i, j] = pmat[cont_idx[i], cont_idx[j]]
    }
  }
  for(j in seq_len(q)) {
    Kyy[j, j] = object$posterior_mean_main$continuous[j, "precision"]
  }

  mux = object$posterior_mean_main$discrete
  mux[is.na(mux)] = 0

  muy = object$posterior_mean_main$continuous[, "mean"]

  list(Kxx = Kxx, Kxy = Kxy, Kyy = Kyy, mux = mux, muy = muy)
}


# ------------------------------------------------------------------
# split_mixed_raw_samples
# ------------------------------------------------------------------
# Split raw main and pairwise sample matrices into separate component
# matrices for the C++ parallel simulation worker.
#
# @param object     Fitted bgms object (mixed MRF).
# @param arguments  Output of extract_arguments().
#
# Returns: list with mux_samples, kxx_samples, muy_samples,
#   kyy_samples, kxy_samples, total_draws.
# ------------------------------------------------------------------
split_mixed_raw_samples = function(object, arguments) {
  p = arguments$num_discrete
  q = arguments$num_continuous
  num_categories = arguments$num_categories
  is_ordinal = arguments$is_ordinal

  main_all = do.call(rbind, object$raw_samples$main)
  pairwise_all = do.call(rbind, object$raw_samples$pairwise)
  total_draws = nrow(main_all)

  # Main layout: [mux_flat | muy | kyy_diag]
  num_mux = sum(ifelse(is_ordinal, num_categories, 2L))
  mux_cols = seq_len(num_mux)
  muy_cols = num_mux + seq_len(q)
  kyy_diag_cols = num_mux + q + seq_len(q)

  mux_samples = main_all[, mux_cols, drop = FALSE]
  muy_samples = main_all[, muy_cols, drop = FALSE]
  kyy_diag_values = main_all[, kyy_diag_cols, drop = FALSE]

  # Pairwise layout: [kxx_ut | kyy_offdiag | kxy]
  nxx = as.integer(p * (p - 1) / 2)
  nyy_offdiag = as.integer(q * (q - 1) / 2)
  nxy = as.integer(p * q)

  kyy_off_end = nxx + nyy_offdiag
  kxy_end = kyy_off_end + nxy

  kxx_samples = if(nxx > 0) {
    pairwise_all[, seq_len(nxx), drop = FALSE]
  } else {
    matrix(0, nrow = total_draws, ncol = 0)
  }

  kyy_offdiag_values = if(nyy_offdiag > 0) {
    pairwise_all[, (nxx + 1):kyy_off_end, drop = FALSE]
  } else {
    matrix(0, nrow = total_draws, ncol = 0)
  }

  kxy_samples = if(nxy > 0) {
    pairwise_all[, (kyy_off_end + 1):kxy_end, drop = FALSE]
  } else {
    matrix(0, nrow = total_draws, ncol = 0)
  }

  # Combine Kyy diagonal and off-diagonal into upper-triangle format
  # C++ expects column-major upper-triangle including diagonal
  nyy_total = as.integer(q * (q + 1) / 2)
  kyy_samples = matrix(0, nrow = total_draws, ncol = nyy_total)
  diag_pos = 0L
  offdiag_pos = 0L
  out_pos = 0L
  for(col in seq_len(q)) {
    for(row in col:q) {
      out_pos = out_pos + 1L
      if(row == col) {
        diag_pos = diag_pos + 1L
        kyy_samples[, out_pos] = kyy_diag_values[, diag_pos]
      } else {
        offdiag_pos = offdiag_pos + 1L
        kyy_samples[, out_pos] = kyy_offdiag_values[, offdiag_pos]
      }
    }
  }

  list(
    mux_samples = mux_samples,
    kxx_samples = kxx_samples,
    muy_samples = muy_samples,
    kyy_samples = kyy_samples,
    kxy_samples = kxy_samples,
    total_draws = total_draws
  )
}


# ------------------------------------------------------------------
# reconstruct_mixed_params_from_row
# ------------------------------------------------------------------
# Reconstruct Kxx, Kxy, Kyy, mux, muy from a single row of split
# sample matrices (used by predict posterior-sample).
#
# @param sample_info  Output of split_mixed_raw_samples().
# @param row_idx      1-based row index.
# @param p            Number of discrete variables.
# @param q            Number of continuous variables.
# @param num_categories  Categories per discrete variable.
# @param is_ordinal   Logical vector.
#
# Returns: list with Kxx, Kxy, Kyy, mux, muy.
# ------------------------------------------------------------------
reconstruct_mixed_params_from_row = function(sample_info, row_idx,
                                             p, q, num_categories,
                                             is_ordinal) {
  mux_vec = sample_info$mux_samples[row_idx, ]
  num_params_disc = ifelse(is_ordinal, num_categories, 2L)
  max_cats = max(num_params_disc)
  mux = matrix(0, nrow = p, ncol = max_cats)
  pos = 1L
  for(s in seq_len(p)) {
    nc = num_params_disc[s]
    mux[s, seq_len(nc)] = mux_vec[pos:(pos + nc - 1L)]
    pos = pos + nc
  }

  Kxx = matrix(0, p, p)
  if(p > 1) {
    kxx_vec = sample_info$kxx_samples[row_idx, ]
    idx = 0L
    for(col in seq_len(p - 1)) {
      for(row in (col + 1):p) {
        idx = idx + 1L
        Kxx[row, col] = kxx_vec[idx]
        Kxx[col, row] = kxx_vec[idx]
      }
    }
  }

  muy = sample_info$muy_samples[row_idx, ]

  kyy_vec = sample_info$kyy_samples[row_idx, ]
  Kyy = matrix(0, q, q)
  idx = 0L
  for(col in seq_len(q)) {
    for(row in col:q) {
      idx = idx + 1L
      Kyy[row, col] = kyy_vec[idx]
      if(row != col) Kyy[col, row] = kyy_vec[idx]
    }
  }

  Kxy = matrix(0, p, q)
  if(p > 0 && q > 0) {
    kxy_vec = sample_info$kxy_samples[row_idx, ]
    idx = 0L
    for(s in seq_len(p)) {
      for(j in seq_len(q)) {
        idx = idx + 1L
        Kxy[s, j] = kxy_vec[idx]
      }
    }
  }

  list(Kxx = Kxx, Kxy = Kxy, Kyy = Kyy, mux = mux, muy = muy)
}


# ------------------------------------------------------------------
# combine_mixed_result
# ------------------------------------------------------------------
# Combine $x (n x p integer) and $y (n x q double) into a single
# n x (p+q) matrix in the original column order.
#
# @param result    List with $x and $y matrices.
# @param disc_idx  Original column indices for discrete variables.
# @param cont_idx  Original column indices for continuous variables.
# @param colnames  Original data column names.
#
# Returns: n x (p+q) numeric matrix.
# ------------------------------------------------------------------
combine_mixed_result = function(result, disc_idx, cont_idx, colnames) {
  n = nrow(result$x)
  num_vars = length(disc_idx) + length(cont_idx)
  out = matrix(NA_real_, nrow = n, ncol = num_vars)
  out[, disc_idx] = result$x
  out[, cont_idx] = result$y
  colnames(out) = colnames
  out
}


# ------------------------------------------------------------------
# format_mixed_predictions
# ------------------------------------------------------------------
# Add names and column labels to C++ prediction output.
#
# @param raw_result         List from C++ compute_conditional_mixed.
# @param predict_vars       1-based user-facing variable indices.
# @param internal_predict_vars  0-based internal indices.
# @param p                  Number of discrete variables.
# @param num_categories     Categories per discrete variable.
# @param data_columnnames   Original data column names.
#
# Returns: Named list of prediction matrices.
# ------------------------------------------------------------------
format_mixed_predictions = function(raw_result, predict_vars,
                                    internal_predict_vars, p,
                                    num_categories, data_columnnames) {
  probs = raw_result
  names(probs) = data_columnnames[predict_vars]

  for(k in seq_along(predict_vars)) {
    int_idx = internal_predict_vars[k]
    if(int_idx < p) {
      s = int_idx + 1L
      n_cats = num_categories[s] + 1
      colnames(probs[[k]]) = paste0("cat_", 0:(n_cats - 1))
    } else {
      colnames(probs[[k]]) = c("mean", "sd")
    }
  }

  probs
}


# ------------------------------------------------------------------
# format_mixed_response
# ------------------------------------------------------------------
# Convert probability predictions to point predictions for mixed models.
#
# @param probs              Named list of prediction matrices.
# @param predict_vars       1-based user-facing variable indices.
# @param internal_predict_vars  0-based internal indices.
# @param p                  Number of discrete variables.
# @param data_columnnames   Original data column names.
#
# Returns: n x length(predict_vars) matrix of predicted values.
# ------------------------------------------------------------------
format_mixed_response = function(probs, predict_vars,
                                 internal_predict_vars, p,
                                 data_columnnames) {
  n = nrow(probs[[1]])
  out = matrix(NA_real_, nrow = n, ncol = length(predict_vars))
  colnames(out) = data_columnnames[predict_vars]

  for(k in seq_along(predict_vars)) {
    int_idx = internal_predict_vars[k]
    if(int_idx < p) {
      out[, k] = apply(probs[[k]], 1, which.max) - 1L
    } else {
      out[, k] = probs[[k]][, 1]
    }
  }

  out
}
