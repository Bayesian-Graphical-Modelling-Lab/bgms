# R/bgmCompare-methods.R

#' @name print.bgmCompare
#' @title Print method for `bgmCompare` objects
#' @description Minimal console output for `bgmCompare` fit objects.
#' @param x An object of class `bgmCompare`.
#' @param ... Ignored.
#' @export
print.bgmCompare = function(x, ...) {
  arguments = extract_arguments(x)

  # Model type
  if(isTRUE(arguments$difference_selection)) {
    prior_msg = switch(as.character(arguments$difference_prior),
      "Bernoulli" = "Bayesian Difference Selection (Bernoulli prior on inclusion)",
      "Beta-Bernoulli" = "Bayesian Difference Selection (Beta-Bernoulli prior on inclusion)",
      "Bayesian Difference Selection"
    )
    cat(prior_msg, "\n")
  } else {
    cat("Bayesian Estimation (multi-group)\n")
  }

  # Dataset info
  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  if(!is.null(arguments$num_groups)) {
    cat(paste0(" Number of groups: ", arguments$num_groups, "\n"))
  }
  if(!is.null(arguments$num_cases)) {
    # In our prepare_output_bgmCompare() we stored total cases in num_cases.
    if(isTRUE(arguments$na_impute)) {
      cat(paste0(" Number of cases: ", arguments$num_cases, " (missings imputed)\n"))
    } else {
      cat(paste0(" Number of cases: ", arguments$num_cases, "\n"))
    }
  }

  # Iterations and chains
  if(!is.null(arguments$num_chains)) {
    total_iter = arguments$iter * arguments$num_chains
    cat(paste0(" Number of post-burnin MCMC iterations: ", total_iter, "\n"))
    cat(paste0(" Number of MCMC chains: ", arguments$num_chains, "\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, "\n"))
  }

  cat("Use the `summary()` function for posterior summaries and diagnostics.\n")
  cat("See the `easybgm` package for additional summaries and plotting.\n")
}


#' @name summary.bgmCompare
#' @title Summary method for `bgmCompare` objects
#'
#' @description Returns posterior summaries and diagnostics for a fitted `bgmCompare` model.
#'
#' @param object An object of class `bgmCompare`.
#' @param ... Currently ignored.
#'
#' @return An object of class `summary.bgmCompare` with posterior summaries.
#' @export
summary.bgmCompare = function(object, ...) {
  arguments = extract_arguments(object)

  if(!is.null(object$posterior_summary_main_baseline) &&
    !is.null(object$posterior_summary_pairwise_baseline)) {
    out = list(
      main      = object$posterior_summary_main_baseline,
      pairwise  = object$posterior_summary_pairwise_baseline
    )

    if(!is.null(object$posterior_summary_indicator)) {
      out$indicator = object$posterior_summary_indicator
    }

    if(!is.null(object$posterior_summary_main_differences)) {
      out$main_diff = object$posterior_summary_main_differences
    }

    if(!is.null(object$posterior_summary_pairwise_differences)) {
      out$pairwise_diff = object$posterior_summary_pairwise_differences
    }

    out$arguments = arguments
    class(out) = "summary.bgmCompare"
    return(out)
  }

  message(
    "No summary statistics available for this model object.\n",
    "Try fitting the model again using the latest bgms version,\n",
    "or use the `easybgm` package for diagnostic summaries and plotting."
  )
  invisible(NULL)
}


#' @export
print.summary.bgmCompare = function(x, digits = 3, ...) {
  cat("Posterior summaries from Bayesian grouped MRF estimation (bgmCompare):\n\n")

  print_df = function(df, digits) {
    df2 = df
    if(ncol(df2) > 1) {
      df2[, -1] = lapply(df2[, -1, drop = FALSE], round, digits = digits)
    }
    print(head(df2, 6))
  }

  if(!is.null(x$main)) {
    cat("Category thresholds:\n")
    print_df(x$main, digits)
    if(nrow(x$main) > 6) {
      cat("... (use `summary(fit)$main` to see full output)\n")
    }
    cat("\n")
  }

  if(!is.null(x$pairwise)) {
    cat("Pairwise interactions:\n")
    print_df(x$pairwise, digits)
    if(nrow(x$pairwise) > 6) {
      cat("... (use `summary(fit)$pairwise` to see full output)\n")
    }
    cat("\n")
  }

  if(!is.null(x$indicator)) {
    cat("Inclusion probabilities:\n")
    ind <- head(x$indicator, 6)

    # round only numeric columns
    ind[] <- lapply(ind, function(col) {
      if(is.numeric(col)) {
        round(col, digits)
      } else {
        col
      }
    })

    # replace NA with empty string for printing
    ind[] <- lapply(ind, function(col) {
      ifelse(is.na(col), "", col)
    })

    print(ind, row.names = FALSE)
    if(nrow(x$indicator) > 6) {
      cat("... (use `summary(fit)$indicator` to see full output)\n")
    }
    cat("Note: NA values are suppressed in the print table. They occur when an indicator\n")
    cat("was constant (all 0 or all 1) across all iterations, so sd/mcse/n_eff/Rhat\n")
    cat("are undefined; `summary(fit)$indicator` still contains the NA values.\n\n")
  }

  if(!is.null(x$main_diff)) {
    cat("Group differences (main effects):\n")

    maind <- head(x$main_diff, 6)

    # Only round numeric columns
    is_num <- vapply(maind, is.numeric, logical(1))
    maind[is_num] <- lapply(
      maind[is_num],
      function(col) ifelse(is.na(col), "", round(col, digits))
    )

    print(maind, row.names = FALSE)

    if(nrow(x$main_diff) > 6) {
      cat("... (use `summary(fit)$main_diff` to see full output)\n")
    }

    if(!is.null(x$indicator)) {
      cat("Note: NA values are suppressed in the print table. They occur here when an\n")
      cat("indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;\n")
      cat("`summary(fit)$main_diff` still contains the NA values.\n")
    }
    cat("\n")
  }

  if(!is.null(x$pairwise_diff)) {
    cat("Group differences (pairwise effects):\n")

    pairwised <- head(x$pairwise_diff, 6)

    # Only round numeric columns
    is_num <- vapply(pairwised, is.numeric, logical(1))
    pairwised[is_num] <- lapply(
      pairwised[is_num],
      function(col) ifelse(is.na(col), "", round(col, digits))
    )

    print(pairwised, row.names = FALSE)

    if(nrow(x$pairwise_diff) > 6) {
      cat("... (use `summary(fit)$pairwise_diff` to see full output)\n")
    }

    if(!is.null(x$indicator)) {
      cat("Note: NA values are suppressed in the print table. They occur here when an\n")
      cat("indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;\n")
      cat("`summary(fit)$pairwise_diff` still contains the NA values.\n")
    }
    cat("\n")
  }

  cat("Use `summary(fit)$<component>` to access full results.\n")
  cat("See the `easybgm` package for other summary and plotting tools.\n")
}


#' @title Extract Coefficients from a bgmCompare Object
#' @name coef.bgmCompare
#' @description Returns posterior means for raw parameters (baseline + differences)
#' and group-specific effects from a \code{bgmCompare} fit, as well as inclusion indicators.
#'
#' @param object An object of class \code{bgmCompare}.
#' @param ... Ignored.
#'
#' @return A list with components:
#' \describe{
#'   \item{main_effects_raw}{Posterior means of the raw main-effect parameters
#'   (variables x [baseline + differences]).}
#'   \item{pairwise_effects_raw}{Posterior means of the raw pairwise-effect parameters
#'   (pairs x [baseline + differences]).}
#'   \item{main_effects_groups}{Posterior means of group-specific main effects
#'   (variables x groups), computed as baseline plus projected differences.}
#'   \item{pairwise_effects_groups}{Posterior means of group-specific pairwise effects
#'   (pairs x groups), computed as baseline plus projected differences.}
#'   \item{indicators}{Posterior mean inclusion probabilities as a symmetric matrix,
#'   with diagonals corresponding to main effects and off-diagonals to pairwise effects.}
#' }
#' @export
coef.bgmCompare <- function(object, ...) {
  args <- extract_arguments(object)

  var_names <- args$data_columnnames
  num_categories <- as.integer(args$num_categories)
  is_ordinal <- as.logical(args$is_ordinal_variable)
  num_groups <- as.integer(args$num_groups)
  num_variables <- as.integer(args$num_variables)
  projection <- args$projection # [num_groups x (num_groups-1)]

  # ---- helper: combine chains into [iter, chain, param], robust to vectors/1-col
  to_array3d <- function(xlist) {
    if(is.null(xlist)) {
      return(NULL)
    }
    stopifnot(length(xlist) >= 1)
    mats <- lapply(xlist, function(x) {
      m <- as.matrix(x)
      if(is.null(dim(m))) m <- matrix(m, ncol = 1L)
      m
    })
    niter <- nrow(mats[[1]])
    nparam <- ncol(mats[[1]])
    arr <- array(NA_real_, dim = c(niter, length(mats), nparam))
    for(c in seq_along(mats)) arr[, c, ] <- mats[[c]]
    arr
  }

  # ============================================================
  # ---- main effects ----
  array3d_main <- to_array3d(object$raw_samples$main)
  stopifnot(!is.null(array3d_main))
  mean_main <- apply(array3d_main, 3, mean)

  stopifnot(length(mean_main) %% num_groups == 0L)
  num_main <- as.integer(length(mean_main) / num_groups)

  main_mat <- matrix(mean_main, nrow = num_main, ncol = num_groups, byrow = FALSE)

  # row names in sampler row order
  rownames(main_mat) <- unlist(lapply(seq_len(num_variables), function(v) {
    if(is_ordinal[v]) {
      paste0(var_names[v], "(c", seq_len(num_categories[v]), ")")
    } else {
      c(
        paste0(var_names[v], "(linear)"),
        paste0(var_names[v], "(quadratic)")
      )
    }
  }))
  colnames(main_mat) <- c("baseline", paste0("diff", seq_len(num_groups - 1L)))

  # group-specific main effects: baseline + P %*% diffs
  main_effects_groups <- matrix(NA_real_, nrow = num_main, ncol = num_groups)
  for(r in seq_len(num_main)) {
    baseline <- main_mat[r, 1]
    diffs <- main_mat[r, -1, drop = TRUE]
    main_effects_groups[r, ] <- baseline + as.vector(projection %*% diffs)
  }
  rownames(main_effects_groups) <- rownames(main_mat)
  colnames(main_effects_groups) <- paste0("group", seq_len(num_groups))

  # ============================================================
  # ---- pairwise effects ----
  array3d_pair <- to_array3d(object$raw_samples$pairwise)
  stopifnot(!is.null(array3d_pair))
  mean_pair <- apply(array3d_pair, 3, mean)

  stopifnot(length(mean_pair) %% num_groups == 0L)
  num_pair <- as.integer(length(mean_pair) / num_groups)

  pairwise_mat <- matrix(mean_pair, nrow = num_pair, ncol = num_groups, byrow = FALSE)

  # row names in sampler row order (upper-tri i<j)
  pair_names <- character()
  if(num_variables >= 2L) {
    for(i in 1L:(num_variables - 1L)) {
      for(j in (i + 1L):num_variables) {
        pair_names <- c(pair_names, paste0(var_names[i], "-", var_names[j]))
      }
    }
  }
  rownames(pairwise_mat) <- pair_names
  colnames(pairwise_mat) <- c("baseline", paste0("diff", seq_len(num_groups - 1L)))

  # group-specific pairwise effects
  pairwise_effects_groups <- matrix(NA_real_, nrow = num_pair, ncol = num_groups)
  for(r in seq_len(num_pair)) {
    baseline <- pairwise_mat[r, 1]
    diffs <- pairwise_mat[r, -1, drop = TRUE]
    pairwise_effects_groups[r, ] <- baseline + as.vector(projection %*% diffs)
  }
  rownames(pairwise_effects_groups) <- rownames(pairwise_mat)
  colnames(pairwise_effects_groups) <- paste0("group", seq_len(num_groups))

  # ============================================================
  # ---- indicators (present only if selection was used) ----
  indicators <- NULL
  array3d_ind <- to_array3d(object$raw_samples$indicator)
  if(!is.null(array3d_ind)) {
    mean_ind <- apply(array3d_ind, 3, mean)

    # reconstruct VxV matrix using the sampler’s interleaved order:
    # (1,1),(1,2),...,(1,V),(2,2),...,(2,V),...,(V,V)
    V <- num_variables
    stopifnot(length(mean_ind) == V * (V + 1L) / 2L)

    ind_mat <- matrix(0,
      nrow = V, ncol = V,
      dimnames = list(var_names, var_names)
    )
    pos <- 1L
    for(i in seq_len(V)) {
      # diagonal (main indicator)
      ind_mat[i, i] <- mean_ind[pos]
      pos <- pos + 1L
      if(i < V) {
        for(j in (i + 1L):V) {
          val <- mean_ind[pos]
          pos <- pos + 1L
          ind_mat[i, j] <- val
          ind_mat[j, i] <- val
        }
      }
    }
    indicators <- ind_mat
  }

  # ============================================================
  # ---- return both raw + group-specific ----
  list(
    main_effects_raw        = main_mat,
    pairwise_effects_raw    = pairwise_mat,
    main_effects_groups     = main_effects_groups,
    pairwise_effects_groups = pairwise_effects_groups,
    indicators              = indicators
  )
}
