##' Extractor Functions for bgms Objects
#'
#' These functions extract various components from objects returned by the `bgm()` function,
#' such as edge indicators, posterior inclusion probabilities, and parameter summaries.
#'
#' @param bgms_object An object of class `bgms` or `bgmCompare`.
#'
#' @section Functions:
#' - `extract_arguments()` – Extract model arguments
#' - `extract_indicators()` – Get sampled edge indicators
#' - `extract_posterior_inclusion_probabilities()` – Posterior edge inclusion probabilities
#' - `extract_pairwise_interactions()` – Posterior mean of pairwise interactions
#' - `extract_category_thresholds()` – Posterior mean of category thresholds
#' - `extract_indicator_priors()` – Prior structure used for edge indicators
#' - `extract_sbm()` – Extract stochastic block model parameters (if applicable)
#' - `extract_rhat()` – Extract R-hat convergence diagnostics
#' - `extract_ess()` – Extract effective sample size estimates
#'
#' @name extractor_functions
#' @title Extractor Functions for bgms Objects
#' @keywords internal
NULL

#' @name extractor_functions
#' @export
extract_arguments = function(bgms_object) {
  UseMethod("extract_arguments")
}

#' @rdname extractor_functions
#' @export
extract_arguments.bgms = function(bgms_object) {
  if(!inherits(bgms_object, "bgms")) stop("Object must be of class 'bgms'.")
  if(is.null(bgms_object$arguments)) {
    stop("Fit object predates bgms version 0.1.3. Upgrade the model output.")
  }
  return(bgms_object$arguments)
}

#' @rdname extractor_functions
#' @export
extract_arguments.bgmCompare = function(bgms_object) {
  if(is.null(bgms_object$arguments)) {
    stop(paste0(
      "Extractor functions have been defined for bgms versions 0.1.3 and up but not \n",
      "for older versions. The current fit object predates version 0.1.3."
    ))
  } else {
    return(bgms_object$arguments)
  }
}

#' @rdname extractor_functions
#' @export
extract_indicators = function(bgms_object) {
  UseMethod("extract_indicators")
}

#' @rdname extractor_functions
#' @details
#' Internally, indicator samples were stored in `$gamma` (pre-0.1.4, now defunct)
#' and `$indicator` (0.1.4–0.1.5, deprecated). As of **bgms 0.1.6.0**, they are
#' stored in `$raw_samples$indicator`.
#' @export
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
    # In 0.1.4.x, $indicator is a matrix (iter x edges), not a list
    if(is.matrix(bgms_object$indicator)) {
      return(bgms_object$indicator)
    }
    # Fallback for list format (if any intermediate versions used this)
    indicator_samples = do.call(rbind, bgms_object$indicator)
    return(indicator_samples)
  }

  # Defunct format (pre-0.1.4): $gamma field
  if(!is.null(bgms_object$gamma)) {
    lifecycle::deprecate_stop(
      "0.1.4",
      I("The '$gamma' field is defunct; please refit with bgms >= 0.1.6.0")
    )
  }

  stop("No indicator samples found in this object.")
}

#' @rdname extractor_functions
#' @export
extract_indicators.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To access difference indicators, the model must be run with difference_selection = TRUE.")
  }

  indicators_list = bgms_object$raw_samples$indicator
  if(is.null(indicators_list)) {
    stop("No indicator samples found in this object.")
  }

  indicator_samples = do.call(rbind, indicators_list)
  param_names = bgms_object$raw_samples$parameter_names$indicators
  if(!is.null(param_names)) {
    colnames(indicator_samples) = param_names
  }
  return(indicator_samples)
}

#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities = function(bgms_object) {
  UseMethod("extract_posterior_inclusion_probabilities")
}

#' @rdname extractor_functions
#' @details
#' Posterior inclusion probabilities are computed from edge indicators.
#'
#' Internally, indicator samples were stored in `$gamma` (pre-0.1.4, now defunct)
#' and `$indicator` (0.1.4–0.1.5, deprecated). As of **bgms 0.1.6.0**, they are
#' stored in `$raw_samples$indicator`.
#' @export
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
    # Check data structure: matrix (raw samples) vs vector (pre-computed means)
    if(is.matrix(bgms_object$indicator) && nrow(bgms_object$indicator) > 1) {
      edge_means = colMeans(bgms_object$indicator)
    } else {
      edge_means = as.vector(bgms_object$indicator)
    }
  } else if(!is.null(bgms_object$gamma)) {
    # Defunct format (pre-0.1.4)
    lifecycle::deprecate_stop(
      "0.1.4.2",
      I("The '$gamma' field is defunct; please refit with bgms >= 0.1.6.0")
    )
  } else {
    stop("No indicator data found to compute posterior inclusion probabilities.")
  }

  pip_matrix = matrix(0, num_vars, num_vars)
  pip_matrix[lower.tri(pip_matrix)] = edge_means
  pip_matrix = pip_matrix + t(pip_matrix)

  colnames(pip_matrix) = data_columnnames
  rownames(pip_matrix) = data_columnnames

  return(pip_matrix)
}


#' @rdname extractor_functions
#' @export
extract_sbm = function(bgms_object) {
  UseMethod("extract_sbm")
}

#' @rdname extractor_functions
#' @export
extract_sbm.bgms = function(bgms_object) {
  if(!inherits(bgms_object, "bgms")) stop("Object must be of class 'bgms'.")

  # Checks
  ver = try(utils::packageVersion("bgms"), silent = TRUE)
  if(inherits(ver, "try-error") || is.na(ver)) {
    stop("Could not determine 'bgms' package version.")
  }
  if(utils::compareVersion(as.character(ver), "0.1.6.0") < 0) {
    stop(paste0(
      "Extractor functions for the SBM prior are defined for bgms version 0.1.6.0. ",
      "The current installed version is ", as.character(ver), "."
    ))
  }

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

  posterior_num_blocks = try(bgms_object$posterior_num_blocks, silent = TRUE)
  posterior_mean_allocations = try(bgms_object$posterior_mean_allocations, silent = TRUE)
  posterior_mode_allocations = try(bgms_object$posterior_mode_allocations, silent = TRUE)
  posterior_mean_coclustering_matrix = try(bgms_object$posterior_mean_coclustering_matrix, silent = TRUE)

  if(inherits(posterior_num_blocks, "try-error")) posterior_num_blocks = NULL
  if(inherits(posterior_mean_allocations, "try-error")) posterior_mean_allocations = NULL
  if(inherits(posterior_mode_allocations, "try-error")) posterior_mode_allocations = NULL
  if(inherits(posterior_mean_coclustering_matrix, "try-error")) posterior_mean_coclustering_matrix = NULL

  if(is.null(posterior_num_blocks) ||
    is.null(posterior_mean_allocations) ||
    is.null(posterior_mode_allocations) ||
    is.null(posterior_mean_coclustering_matrix)) {
    stop(paste0(
      "SBM summaries not found in this object. Missing one or more of: ",
      "posterior_num_blocks, posterior_mean_allocations, ",
      "posterior_mode_allocations, posterior_mean_coclustering_matrix."
    ))
  }


  return(list(
    posterior_num_blocks               = posterior_num_blocks,
    posterior_mean_allocations         = posterior_mean_allocations,
    posterior_mode_allocations         = posterior_mode_allocations,
    posterior_mean_coclustering_matrix = posterior_mean_coclustering_matrix
  ))
}


#' @rdname extractor_functions
#' @export
extract_posterior_inclusion_probabilities.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("To estimate posterior inclusion probabilities, run bgmCompare() with difference_selection = TRUE.")
  }

  var_names = arguments$data_columnnames
  num_categories = as.integer(arguments$num_categories)
  is_ordinal = as.logical(arguments$is_ordinal_variable)
  num_groups = as.integer(arguments$num_groups)
  num_variables = as.integer(arguments$num_variables)
  projection = arguments$projection # [num_groups x (num_groups-1)]

  # ---- helper: combine chains into [iter, chain, param], robust to vectors/1-col
  to_array3d = function(xlist) {
    if(is.null(xlist)) {
      return(NULL)
    }
    stopifnot(length(xlist) >= 1)
    mats = lapply(xlist, function(x) {
      m = as.matrix(x)
      if(is.null(dim(m))) m = matrix(m, ncol = 1L)
      m
    })
    niter = nrow(mats[[1]])
    nparam = ncol(mats[[1]])
    arr = array(NA_real_, dim = c(niter, length(mats), nparam))
    for(c in seq_along(mats)) arr[, c, ] = mats[[c]]
    arr
  }

  array3d_ind = to_array3d(bgms_object$raw_samples$indicator)
  if(!is.null(array3d_ind)) {
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
    indicators = ind_mat
    rownames(indicators) = arguments$data_columnnames
    colnames(indicators) = arguments$data_columnnames
  } else {
    indicators = NULL
  }

  return(indicators)
}

#' @rdname extractor_functions
#' @export
extract_indicator_priors = function(bgms_object) {
  UseMethod("extract_indicator_priors")
}

#' @rdname extractor_functions
#' @export
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


#' @rdname extractor_functions
#' @export
extract_indicator_priors.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(!isTRUE(arguments$difference_selection)) {
    stop("The model ran without selection, so there are no indicator priors specified.")
  }

  return(arguments$difference_prior)
}


#' @rdname extractor_functions
#' @export
extract_pairwise_interactions = function(bgms_object) {
  UseMethod("extract_pairwise_interactions")
}

#' @rdname extractor_functions
#' @details
#' Pairwise interactions were previously stored in `$pairwise_effects` (pre-0.1.4, now
#' defunct) and `$posterior_mean_pairwise` (0.1.4–0.1.5, deprecated). As of **bgms
#' 0.1.6.0**, they are stored in `$raw_samples$pairwise` (raw samples) and
#' `$posterior_summary_pairwise` (summaries).
#' @export
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

    dimnames(mat) = list(paste0("iter", 1:nrow(mat)), edge_names)
    return(mat)
  }

  # Current format fallback (0.1.6.0+): summary only
  if(!is.null(bgms_object$posterior_summary_pairwise)) {
    vec = bgms_object$posterior_summary_pairwise[, "mean"]
    mat = matrix(0, nrow = num_vars, ncol = num_vars)
    mat[lower.tri(mat)] = vec
    mat = mat + t(mat)
    dimnames(mat) = list(var_names, var_names)
    return(mat)
  }

  # Deprecated format (0.1.4–0.1.5): $interactions or $posterior_mean_pairwise
  if(!is.null(bgms_object$interactions)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$interactions' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    # In 0.1.4.x, $interactions is a matrix (iter x edges) of raw samples
    if(is.matrix(bgms_object$interactions) && nrow(bgms_object$interactions) > 1) {
      edge_means = colMeans(bgms_object$interactions)
    } else {
      edge_means = as.vector(bgms_object$interactions)
    }
    mat = matrix(0, nrow = num_vars, ncol = num_vars)
    mat[lower.tri(mat)] = edge_means
    mat = mat + t(mat)
    dimnames(mat) = list(var_names, var_names)
    return(mat)
  }
  if(!is.null(bgms_object$posterior_mean_pairwise)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$posterior_mean_pairwise' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    mat = bgms_object$posterior_mean_pairwise
    dimnames(mat) = list(var_names, var_names)
    return(mat)
  }

  # Defunct format (pre-0.1.4)
  if(!is.null(bgms_object$pairwise_effects)) {
    lifecycle::deprecate_stop(
      "0.1.4.2",
      I("The '$pairwise_effects' field is defunct; please refit with bgms >= 0.1.6.0")
    )
  }

  stop("No pairwise interaction effects found in the object.")
}


#' @rdname extractor_functions
#' @export
extract_pairwise_interactions.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(is.null(bgms_object$raw_samples$pairwise)) {
    stop("No raw samples found for the pairwise effects in the object.")
  }

  pairwise_list = bgms_object$raw_samples$pairwise
  pairwise_samples = do.call(rbind, pairwise_list)

  num_vars = bgms_object$arguments$num_variables
  num_pairs = num_vars * (num_vars - 1) / 2

  pairwise_samples = pairwise_samples[, 1:num_pairs]
  colnames(pairwise_samples) = bgms_object$raw_samples$parameter_names$pairwise_baseline

  return(pairwise_samples)
}

#' @rdname extractor_functions
#' @export
extract_category_thresholds = function(bgms_object) {
  UseMethod("extract_category_thresholds")
}

#' @rdname extractor_functions
#' @details
#' Category thresholds were previously stored in `$main_effects` (pre-0.1.4, now defunct)
#' and `$posterior_mean_main` (0.1.4–0.1.5, deprecated). As of **bgms 0.1.6.0**, they
#' are stored in `$posterior_summary_main`.
#' @export
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

  # Deprecated format (0.1.4–0.1.5): $thresholds or $posterior_mean_main
  if(!is.null(bgms_object$thresholds)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$thresholds' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    # In 0.1.4.x, $thresholds is a matrix (iter x vars) of raw samples
    if(is.matrix(bgms_object$thresholds) && nrow(bgms_object$thresholds) > 1) {
      means = colMeans(bgms_object$thresholds)
    } else {
      means = as.vector(bgms_object$thresholds)
    }
    # For binary variables in 0.1.4.x, there's 1 threshold per variable
    mat = matrix(means, nrow = length(means), ncol = 1)
    rownames(mat) = var_names
    return(mat)
  }
  if(!is.null(bgms_object$posterior_mean_main)) {
    lifecycle::deprecate_warn(
      "0.1.6.0",
      I("The '$posterior_mean_main' field is deprecated; please refit with bgms >= 0.1.6.0")
    )
    mat = bgms_object$posterior_mean_main
    rownames(mat) = var_names
    return(mat)
  }

  # Defunct format (pre-0.1.4)
  if(!is.null(bgms_object$main_effects)) {
    lifecycle::deprecate_stop(
      "0.1.4.2",
      I("The '$main_effects' field is defunct; please refit with bgms >= 0.1.6.0")
    )
  }

  stop("No threshold or main effects found in the object.")
}

#' @rdname extractor_functions
#' @export
extract_category_thresholds.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  if(is.null(bgms_object$raw_samples$main)) {
    stop("No raw samples found for the main effects in the object.")
  }

  main_list = bgms_object$raw_samples$main
  main_samples = do.call(rbind, main_list)

  num_vars = bgms_object$arguments$num_variables
  num_main = length(bgms_object$raw_samples$parameter_names$main_baseline)

  main_samples = main_samples[, 1:num_main]
  colnames(main_samples) = bgms_object$raw_samples$parameter_names$main_baseline

  return(main_samples)
}

#' @rdname extractor_functions
#' @export
extract_group_params = function(bgms_object) {
  UseMethod("extract_group_params")
}

#' @rdname extractor_functions
#' @export
extract_group_params.bgmCompare = function(bgms_object) {
  arguments = extract_arguments(bgms_object)

  var_names = arguments$data_columnnames
  num_categories = as.integer(arguments$num_categories)
  is_ordinal = as.logical(arguments$is_ordinal_variable)
  num_groups = as.integer(arguments$num_groups)
  num_variables = as.integer(arguments$num_variables)
  projection = arguments$projection # [num_groups x (num_groups-1)]

  # ---- helper: combine chains into [iter, chain, param], robust to vectors/1-col
  to_array3d = function(xlist) {
    if(is.null(xlist)) {
      return(NULL)
    }
    stopifnot(length(xlist) >= 1)
    mats = lapply(xlist, function(x) {
      m = as.matrix(x)
      if(is.null(dim(m))) m = matrix(m, ncol = 1L)
      m
    })
    niter = nrow(mats[[1]])
    nparam = ncol(mats[[1]])
    arr = array(NA_real_, dim = c(niter, length(mats), nparam))
    for(c in seq_along(mats)) arr[, c, ] = mats[[c]]
    arr
  }

  # ============================================================
  # ---- main effects ----
  array3d_main = to_array3d(bgms_object$raw_samples$main)
  stopifnot(!is.null(array3d_main))
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
  stopifnot(!is.null(array3d_pair))
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

#' @rdname extractor_functions
#' @export
extract_edge_indicators = function(bgms_object) {
  lifecycle::deprecate_warn("0.1.4.2", "extract_edge_indicators()", "extract_indicators()")
  extract_indicators(bgms_object)
}

#' @rdname extractor_functions
#' @export
extract_pairwise_thresholds = function(bgms_object) {
  lifecycle::deprecate_warn("0.1.4.2", "extract_pairwise_thresholds()", "extract_category_thresholds()")
  extract_category_thresholds(bgms_object)
}


# ------------------------------------------------------------------------------
# extract_rhat() - R-hat Convergence Diagnostics
# ------------------------------------------------------------------------------

#' @rdname extractor_functions
#' @export
extract_rhat = function(bgms_object) {

  UseMethod("extract_rhat")
}

#' @rdname extractor_functions
#' @export
extract_rhat.bgms = function(bgms_object) {
  if(!inherits(bgms_object, "bgms")) stop("Object must be of class 'bgms'.")

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

#' @rdname extractor_functions
#' @export
extract_rhat.bgmCompare = function(bgms_object) {
  if(!inherits(bgms_object, "bgmCompare")) stop("Object must be of class 'bgmCompare'.")

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

#' @rdname extractor_functions
#' @export
extract_ess = function(bgms_object) {
  UseMethod("extract_ess")
}

#' @rdname extractor_functions
#' @export
extract_ess.bgms = function(bgms_object) {
  if(!inherits(bgms_object, "bgms")) stop("Object must be of class 'bgms'.")

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

#' @rdname extractor_functions
#' @export
extract_ess.bgmCompare = function(bgms_object) {
  if(!inherits(bgms_object, "bgmCompare")) stop("Object must be of class 'bgmCompare'.")

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

