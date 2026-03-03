# ==============================================================================
# Extractor Functions - S3 Generics
# ==============================================================================
# These are internal S3 generics. See the method documentation for details:
# - extract_arguments.bgms, extract_arguments.bgmCompare
# - extract_indicators.bgms, extract_indicators.bgmCompare, etc.
# ==============================================================================

# S3 generic for extract_arguments
extract_arguments = function(bgms_object) {
  UseMethod("extract_arguments")
}

#' Extract Model Arguments from bgms Objects
#'
#' Retrieves the arguments used when fitting a model with [bgm()].
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A named list containing all arguments passed to [bgm()], including
#'   data dimensions, prior settings, and MCMC configuration.
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_arguments.bgms = function(bgms_object) {
  if(is.null(bgms_object$arguments)) {
    stop("Fit object predates bgms version 0.1.3. Upgrade the model output.")
  }
  return(bgms_object$arguments)
}

#' Extract Model Arguments from bgmCompare Objects
#'
#' Retrieves the arguments used when fitting a model with [bgmCompare()].
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A named list containing all arguments passed to [bgmCompare()],
#'   including data dimensions, prior settings, and MCMC configuration.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_arguments.bgmCompare = function(bgms_object) {
  if(is.null(bgms_object$arguments)) {
    stop("Fit object predates bgms version 0.1.3. Upgrade the model output.")
  }
  return(bgms_object$arguments)
}

# S3 generic for extract_indicators
extract_indicators = function(bgms_object) {
  UseMethod("extract_indicators")
}

#' Extract Edge Indicator Samples from bgms Objects
#'
#' Retrieves posterior samples of edge inclusion indicators from a
#' model fitted with [bgm()] using edge selection.
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A matrix with one row per post-warmup iteration and one column
#'   per edge, containing binary (0/1) indicator samples.
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_indicators.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$edge_selection)) {
    stop("To access edge indicators, the model must be run with edge_selection = TRUE.")
  }

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$raw_samples$indicator)) {
    indicators_list = bgms_object$raw_samples$indicator
    indicator_samples = do.call(rbind, indicators_list)
    param_names = bgms_object$raw_samples$parameter_names$indicator
    stopifnot("parameter_names$indicator missing in fit object" = !is.null(param_names))
    colnames(indicator_samples) = param_names
    return(indicator_samples)
  }

  # Deprecated format (0.1.4–0.1.5): $indicator stored at top level
  if(!is.null(bgms_object$indicator)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$indicator)
  }

  # Defunct format (pre-0.1.4): $gamma field
  lifecycle::deprecate_stop(
    "0.1.4",
    I("The '$gamma' field is defunct; please refit with bgms >= 0.1.6.0")
  )
}

#' Extract Difference Indicator Samples from bgmCompare Objects
#'
#' Retrieves posterior samples of difference indicators from a model
#' fitted with [bgmCompare()] using difference selection.
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A matrix with one row per post-warmup iteration, containing
#'   binary indicator samples for main-effect and pairwise differences.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_indicators.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To access difference indicators, the model must be run with difference_selection = TRUE.")
  }

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$raw_samples$indicator)) {
    indicator_samples = do.call(rbind, bgms_object$raw_samples$indicator)
    param_names = bgms_object$raw_samples$parameter_names$indicators
    if(!is.null(param_names)) {
      colnames(indicator_samples) = param_names
    }
    return(indicator_samples)
  }

  # Deprecated format (0.1.4–0.1.5): $pairwise_difference_indicator at top level
  if(!is.null(bgms_object$pairwise_difference_indicator)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$pairwise_difference_indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$pairwise_difference_indicator)
  }

  stop("No indicator samples found in fit object.")
}

# S3 generic for extract_posterior_inclusion_probabilities
extract_posterior_inclusion_probabilities = function(bgms_object) {
  UseMethod("extract_posterior_inclusion_probabilities")
}

#' Extract Posterior Inclusion Probabilities from bgms Objects
#'
#' Computes the posterior probability of edge inclusion from a model
#' fitted with [bgm()] using edge selection.
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A symmetric p x p matrix of posterior inclusion
#'   probabilities, with variable names as row and column names.
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_posterior_inclusion_probabilities.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$edge_selection)) {
    stop("To estimate posterior inclusion probabilities, run bgm() with edge_selection = TRUE.")
  }

  # Handle legacy field name (no_variables → num_variables in 0.1.6.0)
  num_vars = arguments$num_variables %||% arguments$no_variables
  data_columnnames = arguments$data_columnnames

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$raw_samples$indicator)) {
    indicator_samples = extract_indicators(bgms_object)
    edge_means = colMeans(indicator_samples)
  } else if(!is.null(bgms_object$indicator)) {
    # Deprecated format (0.1.4–0.1.5): $indicator at top level
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    edge_means = colMeans(bgms_object$indicator)
  } else {
    # Defunct format (pre-0.1.4)
    lifecycle::deprecate_stop(
      "0.1.4.2",
      I("The '$gamma' field is defunct; please refit with bgms >= 0.1.6.0")
    )
  }

  pip_matrix = matrix(0, num_vars, num_vars)
  pip_matrix[lower.tri(pip_matrix)] = edge_means
  pip_matrix = pip_matrix + t(pip_matrix)

  colnames(pip_matrix) = data_columnnames
  rownames(pip_matrix) = data_columnnames

  return(pip_matrix)
}


# S3 generic for extract_sbm
extract_sbm = function(bgms_object) {
  UseMethod("extract_sbm")
}

#' Extract Stochastic Block Model Summaries from bgms Objects
#'
#' Retrieves posterior summaries from a model fitted with [bgm()] using
#' the Stochastic Block prior on edge inclusion.
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A list with `posterior_num_blocks`, `posterior_mean_allocations`,
#'   `posterior_mode_allocations`, and `posterior_mean_coclustering_matrix`.
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_sbm.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$edge_selection)) {
    stop("To extract SBM summaries, run bgm() with edge_selection = TRUE.")
  }
  if(!identical(arguments$edge_prior, "Stochastic-Block")) {
    stop(paste0(
      "edge_prior must be 'Stochastic-Block' (got '",
      as.character(arguments$edge_prior), "')."
    ))
  }

  return(list(
    posterior_num_blocks               = bgms_object$posterior_num_blocks,
    posterior_mean_allocations         = bgms_object$posterior_mean_allocations,
    posterior_mode_allocations         = bgms_object$posterior_mode_allocations,
    posterior_mean_coclustering_matrix = bgms_object$posterior_mean_coclustering_matrix
  ))
}


#' Extract Posterior Inclusion Probabilities from bgmCompare Objects
#'
#' Computes the posterior probability of difference inclusion from a
#' model fitted with [bgmCompare()] using difference selection.
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A symmetric p x p matrix where diagonal entries
#'   are main-effect inclusion probabilities and off-diagonal entries are
#'   pairwise difference inclusion probabilities.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_posterior_inclusion_probabilities.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To estimate posterior inclusion probabilities, run bgmCompare() with difference_selection = TRUE.")
  }

  var_names = arguments$data_columnnames
  # Handle legacy field name (no_variables → num_variables in 0.1.6.0)
  num_variables = as.integer(arguments$num_variables %||% arguments$no_variables)

  # ---- helper: combine chains into [iter, chain, param]
  to_array3d = function(xlist) {
    stopifnot(length(xlist) >= 1)
    mats = lapply(xlist, as.matrix)
    niter = nrow(mats[[1]])
    nparam = ncol(mats[[1]])
    arr = array(NA_real_, dim = c(niter, length(mats), nparam))
    for(c in seq_along(mats)) arr[, c, ] = mats[[c]]
    arr
  }

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$raw_samples$indicator)) {
    array3d_ind = to_array3d(bgms_object$raw_samples$indicator)
    mean_ind = apply(array3d_ind, 3, mean)

    # reconstruct VxV matrix using the sampler’s interleaved order:
    # (1,1),(1,2),...,(1,V),(2,2),...,(2,V),...,(V,V)
    V = num_variables
    stopifnot(length(mean_ind) == V * (V + 1L) / 2L)

    ind_mat = matrix(0,
      nrow = V, ncol = V,
      dimnames = list(var_names, var_names)
    )
    pos = 1L
    for(i in seq_len(V)) {
      # diagonal (main indicator)
      ind_mat[i, i] = mean_ind[pos]
      pos = pos + 1L
      if(i < V) {
        for(j in (i + 1L):V) {
          val = mean_ind[pos]
          pos = pos + 1L
          ind_mat[i, j] = val
          ind_mat[j, i] = val
        }
      }
    }

    rownames(ind_mat) = arguments$data_columnnames
    colnames(ind_mat) = arguments$data_columnnames
    return(ind_mat)
  }

  # Deprecated format (0.1.4–0.1.5): $pairwise_difference_indicator at top level
  if(!is.null(bgms_object$pairwise_difference_indicator)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$pairwise_difference_indicator' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    edge_means = colMeans(bgms_object$pairwise_difference_indicator)
    V = num_variables
    ind_mat = matrix(0, nrow = V, ncol = V)
    ind_mat[lower.tri(ind_mat)] = edge_means
    ind_mat = ind_mat + t(ind_mat)
    dimnames(ind_mat) = list(var_names, var_names)
    return(ind_mat)
  }

  stop("No indicator samples found in fit object.")
}

# S3 generic for extract_indicator_priors
extract_indicator_priors = function(bgms_object) {
  UseMethod("extract_indicator_priors")
}

#' Extract Indicator Prior Structure from bgms Objects
#'
#' Retrieves the prior specification used for edge indicators in a
#' model fitted with [bgm()].
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A named list describing the prior structure, including the prior type
#'   and any hyperparameters.
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_indicator_priors.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  if(!isTRUE(arguments$edge_selection)) stop("No edge selection performed.")

  switch(arguments$edge_prior,
    "Bernoulli" = list(type = "Bernoulli", prior_inclusion_probability = arguments$inclusion_probability),
    "Beta-Bernoulli" = list(type = "Beta-Bernoulli", alpha = arguments$beta_bernoulli_alpha, beta = arguments$beta_bernoulli_beta),
    "Stochastic-Block" = list(
      type = "Stochastic-Block",
      beta_bernoulli_alpha = arguments$beta_bernoulli_alpha,
      beta_bernoulli_beta = arguments$beta_bernoulli_beta,
      dirichlet_alpha = arguments$dirichlet_alpha
    )
  )
}


#' Extract Indicator Prior Structure from bgmCompare Objects
#'
#' Retrieves the prior specification used for difference indicators in a
#' model fitted with [bgmCompare()].
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A named list describing the prior structure, including the prior type
#'   and any hyperparameters.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_indicator_priors.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("The model ran without selection, so there are no indicator priors specified.")
  }

  return(arguments$difference_prior)
}


# S3 generic for extract_pairwise_interactions
extract_pairwise_interactions = function(bgms_object) {
  UseMethod("extract_pairwise_interactions")
}

#' Extract Pairwise Interaction Samples from bgms Objects
#'
#' Retrieves posterior samples of pairwise interaction parameters from
#' a model fitted with [bgm()].
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A matrix with one row per post-warmup iteration and one column per
#'   edge, containing posterior samples of interaction strengths.
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_pairwise_interactions.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  # Handle legacy field name (no_variables → num_variables in 0.1.6.0)
  num_vars = arguments$num_variables %||% arguments$no_variables
  var_names = arguments$data_columnnames

  # Current format (0.1.6.0+): raw samples
  if(!is.null(bgms_object$raw_samples)) {
    mats = bgms_object$raw_samples$pairwise
    mat = do.call(rbind, mats)

    edge_names = character()
    for(i in 1:(num_vars - 1)) {
      for(j in (i + 1):num_vars) {
        edge_names = c(edge_names, paste0(var_names[i], "-", var_names[j]))
      }
    }

    dimnames(mat) = list(paste0("iter", seq_len(nrow(mat))), edge_names)
    return(mat)
  }

  # Deprecated format (0.1.4–0.1.5): $interactions
  if(!is.null(bgms_object$interactions)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$interactions' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    edge_means = colMeans(bgms_object$interactions)
    mat = matrix(0, nrow = num_vars, ncol = num_vars)
    mat[lower.tri(mat)] = edge_means
    mat = mat + t(mat)
    dimnames(mat) = list(var_names, var_names)
    return(mat)
  }

  # Defunct format (pre-0.1.4)
  lifecycle::deprecate_stop(
    "0.1.4.2",
    I("The '$pairwise_effects' field is defunct; please refit with bgms >= 0.1.6.0")
  )
}


#' Extract Pairwise Interaction Samples from bgmCompare Objects
#'
#' Retrieves posterior samples of baseline pairwise interaction
#' parameters from a model fitted with [bgmCompare()].
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A matrix with one row per post-warmup iteration and one column per
#'   edge, containing posterior samples of baseline interaction strengths.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_pairwise_interactions.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$raw_samples$pairwise)) {
    pairwise_samples = do.call(rbind, bgms_object$raw_samples$pairwise)

    num_vars = bgms_object$arguments$num_variables
    num_pairs = num_vars * (num_vars - 1) / 2

    pairwise_samples = pairwise_samples[, 1:num_pairs]
    colnames(pairwise_samples) = bgms_object$raw_samples$parameter_names$pairwise_baseline

    return(pairwise_samples)
  }

  # Deprecated format (0.1.4–0.1.5): $interactions at top level
  if(!is.null(bgms_object$interactions)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$interactions' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$interactions)
  }

  stop("No pairwise interaction samples found in fit object.")
}

# S3 generic for extract_category_thresholds
extract_category_thresholds = function(bgms_object) {
  UseMethod("extract_category_thresholds")
}

#' Extract Category Threshold Estimates from bgms Objects
#'
#' Retrieves posterior mean category threshold parameters from a model
#' fitted with [bgm()].
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A matrix with one row per variable and one column per category
#'   threshold, containing posterior means.
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_category_thresholds.bgms = function(bgms_object) {
  arguments = extract_arguments(bgms_object)
  var_names = arguments$data_columnnames

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$posterior_summary_main)) {
    vec = bgms_object$posterior_summary_main[, "mean"]
    # Handle legacy field name (no_variables → num_variables in 0.1.6.0)
    num_vars = arguments$num_variables %||% arguments$no_variables
    variable_type = arguments$variable_type
    if(length(variable_type) == 1) {
      variable_type = rep(variable_type, num_vars)
    }
    num_cats = arguments$num_categories
    max_cats = max(num_cats)
    mat = matrix(NA_real_, nrow = num_vars, ncol = max_cats)
    rownames(mat) = var_names
    pos = 1
    for(v in seq_len(num_vars)) {
      if(variable_type[v] == "ordinal") {
        k = num_cats[v]
        mat[v, 1:k] = vec[pos:(pos + k - 1)]
        pos = pos + k
      } else {
        mat[v, 1:2] = vec[pos:(pos + 1)]
        pos = pos + 2
      }
    }
    return(mat)
  }

  # Deprecated format (0.1.4–0.1.5): $thresholds
  if(!is.null(bgms_object$thresholds)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$thresholds' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    means = colMeans(bgms_object$thresholds)
    # For binary variables in 0.1.4.x, there's 1 threshold per variable
    mat = matrix(means, nrow = length(means), ncol = 1)
    rownames(mat) = var_names
    return(mat)
  }

  # Defunct format (pre-0.1.4)
  lifecycle::deprecate_stop(
    "0.1.4.2",
    I("The '$main_effects' field is defunct; please refit with bgms >= 0.1.6.0")
  )
}

#' Extract Category Threshold Samples from bgmCompare Objects
#'
#' Retrieves posterior samples of baseline main-effect (threshold)
#' parameters from a model fitted with [bgmCompare()].
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A matrix with one row per post-warmup iteration, containing
#'   posterior samples of baseline threshold parameters.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_category_thresholds.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$raw_samples$main)) {
    main_samples = do.call(rbind, bgms_object$raw_samples$main)

    num_vars = bgms_object$arguments$num_variables
    num_main = length(bgms_object$raw_samples$parameter_names$main_baseline)

    main_samples = main_samples[, 1:num_main]
    colnames(main_samples) = bgms_object$raw_samples$parameter_names$main_baseline

    return(main_samples)
  }

  # Deprecated format (0.1.4–0.1.5): $thresholds or $thresholds_gr1/$thresholds_gr2 at top level
  if(!is.null(bgms_object$thresholds)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$thresholds' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(bgms_object$thresholds)
  }

  # Alternative deprecated format (0.1.4.1+): $thresholds_gr1, $thresholds_gr2
  if(!is.null(bgms_object$thresholds_gr1)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$thresholds_gr*' fields are deprecated; please refit with bgms >= 0.1.6.0")
    )
    # Combine the two groups' thresholds
    return(cbind(bgms_object$thresholds_gr1, bgms_object$thresholds_gr2))
  }

  stop("No category threshold samples found in fit object.")
}

# S3 generic for extract_group_params
extract_group_params = function(bgms_object) {
  UseMethod("extract_group_params")
}

#' Extract Group-Specific Parameters from bgmCompare Objects
#'
#' Computes group-specific parameter estimates by combining baseline
#' parameters and group differences from a model fitted with [bgmCompare()].
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A list with `main_effects_groups` (main effects per group) and
#'   `pairwise_effects_groups` (pairwise effects per group).
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_group_params.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  # Current format (0.1.6.0+)
  if(!is.null(bgms_object$raw_samples$main)) {
    return(.extract_group_params_current(bgms_object, arguments))
  }

  # Deprecated format (0.1.4–0.1.5): separate fields for baseline and differences
  if(!is.null(bgms_object$interactions) && !is.null(bgms_object$pairwise_difference)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The legacy bgmCompare format is deprecated; please refit with bgms >= 0.1.6.0")
    )
    return(.extract_group_params_legacy(bgms_object, arguments))
  }

  stop("No group parameter samples found in fit object.")
}

# Helper for current format (0.1.6+)
.extract_group_params_current = function(bgms_object, arguments) {
  var_names = arguments$data_columnnames
  num_categories = as.integer(arguments$num_categories)
  is_ordinal = as.logical(arguments$is_ordinal_variable)
  num_groups = as.integer(arguments$num_groups)
  num_variables = as.integer(arguments$num_variables)
  projection = arguments$projection # [num_groups x (num_groups-1)]

  # ---- helper: combine chains into [iter, chain, param]
  to_array3d = function(xlist) {
    stopifnot(length(xlist) >= 1)
    mats = lapply(xlist, as.matrix)
    niter = nrow(mats[[1]])
    nparam = ncol(mats[[1]])
    arr = array(NA_real_, dim = c(niter, length(mats), nparam))
    for(c in seq_along(mats)) arr[, c, ] = mats[[c]]
    arr
  }

  # ============================================================
  # ---- main effects ----
  array3d_main = to_array3d(bgms_object$raw_samples$main)
  mean_main = apply(array3d_main, 3, mean)

  stopifnot(length(mean_main) %% num_groups == 0L)
  num_main = as.integer(length(mean_main) / num_groups)

  main_mat = matrix(mean_main, nrow = num_main, ncol = num_groups, byrow = FALSE)

  # row names in sampler row order
  rownames(main_mat) = unlist(lapply(seq_len(num_variables), function(v) {
    if(is_ordinal[v]) {
      paste0(var_names[v], "(c", seq_len(num_categories[v]), ")")
    } else {
      c(
        paste0(var_names[v], "(linear)"),
        paste0(var_names[v], "(quadratic)")
      )
    }
  }))
  colnames(main_mat) = c("baseline", paste0("diff", seq_len(num_groups - 1L)))

  # group-specific main effects: baseline + P %*% diffs
  main_effects_groups = matrix(NA_real_, nrow = num_main, ncol = num_groups)
  for(r in seq_len(num_main)) {
    baseline = main_mat[r, 1]
    diffs = main_mat[r, -1, drop = TRUE]
    main_effects_groups[r, ] = baseline + as.vector(projection %*% diffs)
  }
  rownames(main_effects_groups) = rownames(main_mat)
  colnames(main_effects_groups) = paste0("group", seq_len(num_groups))

  # ============================================================
  # ---- pairwise effects ----
  array3d_pair = to_array3d(bgms_object$raw_samples$pairwise)
  mean_pair = apply(array3d_pair, 3, mean)

  stopifnot(length(mean_pair) %% num_groups == 0L)
  num_pair = as.integer(length(mean_pair) / num_groups)

  pairwise_mat = matrix(mean_pair, nrow = num_pair, ncol = num_groups, byrow = FALSE)

  # row names in sampler row order (upper-tri i<j)
  pair_names = character()
  if(num_variables >= 2L) {
    for(i in 1L:(num_variables - 1L)) {
      for(j in (i + 1L):num_variables) {
        pair_names = c(pair_names, paste0(var_names[i], "-", var_names[j]))
      }
    }
  }
  rownames(pairwise_mat) = pair_names
  colnames(pairwise_mat) = c("baseline", paste0("diff", seq_len(num_groups - 1L)))

  # group-specific pairwise effects
  pairwise_effects_groups = matrix(NA_real_, nrow = num_pair, ncol = num_groups)
  for(r in seq_len(num_pair)) {
    baseline = pairwise_mat[r, 1]
    diffs = pairwise_mat[r, -1, drop = TRUE]
    pairwise_effects_groups[r, ] = baseline + as.vector(projection %*% diffs)
  }
  rownames(pairwise_effects_groups) = rownames(pairwise_mat)
  colnames(pairwise_effects_groups) = paste0("group", seq_len(num_groups))

  return(list(
    main_effects_groups = main_effects_groups,
    pairwise_effects_groups = pairwise_effects_groups
  ))
}

# Helper for legacy format (0.1.4–0.1.5)
# v0.1.4.x only supported 2 groups with parameterization:
#   group1 = baseline + diff, group2 = baseline - diff
.extract_group_params_legacy = function(bgms_object, arguments) {
  var_names = arguments$data_columnnames
  # Handle legacy field name (no_variables → num_variables in 0.1.6.0)
  num_variables = as.integer(arguments$num_variables %||% arguments$no_variables)

  # v0.1.4 format: baseline interactions and differences are separate
  # $interactions: [iter x n_pairs] baseline pairwise effects
  # $pairwise_difference: [iter x n_pairs] pairwise differences
  # $thresholds or $thresholds_gr1/$thresholds_gr2: main effects
  # $main_difference: [iter x n_vars] main differences

  # Compute posterior means
  mean_interactions = colMeans(bgms_object$interactions)
  mean_pairwise_diff = colMeans(bgms_object$pairwise_difference)

  # Get thresholds (handles both v0.1.4 and v0.1.4.1+ formats)
  if(!is.null(bgms_object$thresholds)) {
    mean_thresholds = colMeans(bgms_object$thresholds)
  } else if(!is.null(bgms_object$thresholds_gr1)) {
    # v0.1.4.1+ stored group-specific thresholds directly
    mean_thresholds_gr1 = colMeans(bgms_object$thresholds_gr1)
    mean_thresholds_gr2 = colMeans(bgms_object$thresholds_gr2)
    # Return directly since we have group-specific values
    main_effects_groups = cbind(mean_thresholds_gr1, mean_thresholds_gr2)
    colnames(main_effects_groups) = c("group1", "group2")
    rownames(main_effects_groups) = var_names

    pairwise_effects_groups = cbind(
      mean_interactions + mean_pairwise_diff,
      mean_interactions - mean_pairwise_diff
    )
    colnames(pairwise_effects_groups) = c("group1", "group2")

    # Row names for pairs
    pair_names = character()
    if(num_variables >= 2L) {
      for(i in 1L:(num_variables - 1L)) {
        for(j in (i + 1L):num_variables) {
          pair_names = c(pair_names, paste0(var_names[i], "-", var_names[j]))
        }
      }
    }
    rownames(pairwise_effects_groups) = pair_names

    return(list(
      main_effects_groups = main_effects_groups,
      pairwise_effects_groups = pairwise_effects_groups
    ))
  } else {
    stop("No threshold samples found in legacy fit object.")
  }

  mean_main_diff = colMeans(bgms_object$main_difference)

  # v0.1.4 parameterization: group1 = baseline + diff, group2 = baseline - diff
  main_effects_groups = cbind(
    mean_thresholds + mean_main_diff,
    mean_thresholds - mean_main_diff
  )
  colnames(main_effects_groups) = c("group1", "group2")
  rownames(main_effects_groups) = var_names

  pairwise_effects_groups = cbind(
    mean_interactions + mean_pairwise_diff,
    mean_interactions - mean_pairwise_diff
  )
  colnames(pairwise_effects_groups) = c("group1", "group2")

  # Row names for pairs
  pair_names = character()
  if(num_variables >= 2L) {
    for(i in 1L:(num_variables - 1L)) {
      for(j in (i + 1L):num_variables) {
        pair_names = c(pair_names, paste0(var_names[i], "-", var_names[j]))
      }
    }
  }
  rownames(pairwise_effects_groups) = pair_names

  return(list(
    main_effects_groups = main_effects_groups,
    pairwise_effects_groups = pairwise_effects_groups
  ))
}

#' Deprecated: Use extract_indicators instead
#' @param bgms_object A bgms or bgmCompare object.
#' @keywords internal
#' @export
extract_edge_indicators = function(bgms_object) {
  lifecycle::deprecate_warn("0.1.4.2", "extract_edge_indicators()", "extract_indicators()")
  extract_indicators(bgms_object)
}

#' Deprecated: Use extract_category_thresholds instead
#' @param bgms_object A bgms or bgmCompare object.
#' @keywords internal
#' @export
extract_pairwise_thresholds = function(bgms_object) {
  lifecycle::deprecate_warn("0.1.4.2", "extract_pairwise_thresholds()", "extract_category_thresholds()")
  extract_category_thresholds(bgms_object)
}


# ------------------------------------------------------------------------------
# extract_rhat() - R-hat Convergence Diagnostics
# ------------------------------------------------------------------------------

# S3 generic for extract_rhat
extract_rhat = function(bgms_object) {
  UseMethod("extract_rhat")
}

#' Extract R-hat Diagnostics from bgms Objects
#'
#' Retrieves R-hat convergence diagnostics for all parameters from a
#' model fitted with [bgm()].
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A named list with R-hat values for each parameter type present in
#'   the model (e.g., `main`, `pairwise`, `indicator`).
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_rhat.bgms = function(bgms_object) {
  result = list()

  # Main effect Rhat
  if(!is.null(bgms_object$posterior_summary_main)) {
    result$main = bgms_object$posterior_summary_main$Rhat
    names(result$main) = rownames(bgms_object$posterior_summary_main)
  }

  # Pairwise interaction Rhat
  if(!is.null(bgms_object$posterior_summary_pairwise)) {
    result$pairwise = bgms_object$posterior_summary_pairwise$Rhat
    names(result$pairwise) = rownames(bgms_object$posterior_summary_pairwise)
  }

  # Indicator Rhat (if edge selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = bgms_object$posterior_summary_indicator$Rhat
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}

#' Extract R-hat Diagnostics from bgmCompare Objects
#'
#' Retrieves R-hat convergence diagnostics for all parameters from a
#' model fitted with [bgmCompare()].
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A named list with R-hat values for each parameter type present in
#'   the model.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_rhat.bgmCompare = function(bgms_object) {
  result = list()

  # Main baseline Rhat
  if(!is.null(bgms_object$posterior_summary_main_baseline)) {
    result$main_baseline = bgms_object$posterior_summary_main_baseline$Rhat
    names(result$main_baseline) = rownames(bgms_object$posterior_summary_main_baseline)
  }

  # Main differences Rhat
  if(!is.null(bgms_object$posterior_summary_main_differences)) {
    result$main_differences = bgms_object$posterior_summary_main_differences$Rhat
    names(result$main_differences) = rownames(bgms_object$posterior_summary_main_differences)
  }

  # Pairwise baseline Rhat
  if(!is.null(bgms_object$posterior_summary_pairwise_baseline)) {
    result$pairwise_baseline = bgms_object$posterior_summary_pairwise_baseline$Rhat
    names(result$pairwise_baseline) = rownames(bgms_object$posterior_summary_pairwise_baseline)
  }

  # Pairwise differences Rhat
  if(!is.null(bgms_object$posterior_summary_pairwise_differences)) {
    result$pairwise_differences = bgms_object$posterior_summary_pairwise_differences$Rhat
    names(result$pairwise_differences) = rownames(bgms_object$posterior_summary_pairwise_differences)
  }

  # Indicator Rhat (if difference selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = bgms_object$posterior_summary_indicator$Rhat
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}


# ------------------------------------------------------------------------------
# extract_ess() - Effective Sample Size
# ------------------------------------------------------------------------------

# S3 generic for extract_ess
extract_ess = function(bgms_object) {
  UseMethod("extract_ess")
}

#' Extract Effective Sample Size from bgms Objects
#'
#' Retrieves effective sample size estimates for all parameters from a
#' model fitted with [bgm()].
#'
#' @param bgms_object An object of class `bgms` returned by [bgm()].
#'
#' @return A named list with ESS values for each parameter type present in
#'   the model (e.g., `main`, `pairwise`, `indicator`).
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family extractors
#'
#' @exportS3Method
extract_ess.bgms = function(bgms_object) {
  result = list()

  # Main effect ESS
  if(!is.null(bgms_object$posterior_summary_main)) {
    result$main = bgms_object$posterior_summary_main$n_eff
    names(result$main) = rownames(bgms_object$posterior_summary_main)
  }

  # Pairwise interaction ESS
  if(!is.null(bgms_object$posterior_summary_pairwise)) {
    result$pairwise = bgms_object$posterior_summary_pairwise$n_eff
    names(result$pairwise) = rownames(bgms_object$posterior_summary_pairwise)
  }

  # Indicator ESS (if edge selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = bgms_object$posterior_summary_indicator$n_eff
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}

#' Extract Effective Sample Size from bgmCompare Objects
#'
#' Retrieves effective sample size estimates for all parameters from a
#' model fitted with [bgmCompare()].
#'
#' @param bgms_object An object of class `bgmCompare` returned by [bgmCompare()].
#'
#' @return A named list with ESS values for each parameter type present in
#'   the model.
#'
#' @seealso [bgmCompare()], [summary.bgmCompare()], [coef.bgmCompare()]
#' @family extractors
#'
#' @exportS3Method
extract_ess.bgmCompare = function(bgms_object) {
  result = list()

  # Main baseline ESS
  if(!is.null(bgms_object$posterior_summary_main_baseline)) {
    result$main_baseline = bgms_object$posterior_summary_main_baseline$n_eff
    names(result$main_baseline) = rownames(bgms_object$posterior_summary_main_baseline)
  }

  # Main differences ESS
  if(!is.null(bgms_object$posterior_summary_main_differences)) {
    result$main_differences = bgms_object$posterior_summary_main_differences$n_eff
    names(result$main_differences) = rownames(bgms_object$posterior_summary_main_differences)
  }

  # Pairwise baseline ESS
  if(!is.null(bgms_object$posterior_summary_pairwise_baseline)) {
    result$pairwise_baseline = bgms_object$posterior_summary_pairwise_baseline$n_eff
    names(result$pairwise_baseline) = rownames(bgms_object$posterior_summary_pairwise_baseline)
  }

  # Pairwise differences ESS
  if(!is.null(bgms_object$posterior_summary_pairwise_differences)) {
    result$pairwise_differences = bgms_object$posterior_summary_pairwise_differences$n_eff
    names(result$pairwise_differences) = rownames(bgms_object$posterior_summary_pairwise_differences)
  }

  # Indicator ESS (if difference selection was used)
  if(!is.null(bgms_object$posterior_summary_indicator)) {
    result$indicator = bgms_object$posterior_summary_indicator$n_eff
    names(result$indicator) = rownames(bgms_object$posterior_summary_indicator)
  }

  if(length(result) == 0) {
    stop("No posterior summary information found in this object.")
  }

  return(result)
}
