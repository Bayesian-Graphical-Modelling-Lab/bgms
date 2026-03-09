#' @name print.bgms
#' @title Print method for `bgms` objects
#'
#' @description Minimal console output for `bgms` fit objects.
#'
#' @param x An object of class `bgms`.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' print(fit)
#' }
#'
#' @seealso [bgm()], [summary.bgms()], [coef.bgms()]
#' @family posterior-methods
#'
#' @export
print.bgms = function(x, ...) {
  arguments = extract_arguments(x)

  # Estimation method
  if(isTRUE(arguments$edge_selection)) {
    prior_msg = switch(arguments$edge_prior,
      "Bernoulli" = "Bayesian Edge Selection using a Bernoulli prior on edge inclusion",
      "Beta-Bernoulli" = "Bayesian Edge Selection using a Beta-Bernoulli prior on edge inclusion",
      "Stochastic-Block" = "Bayesian Edge Selection using a Stochastic Block prior on edge inclusion",
      "Bayesian Edge Selection"
    )
    cat(prior_msg, "\n")
  } else {
    cat("Bayesian Estimation\n")
  }

  # Model type
  mt = arguments$model_type
  if(!is.null(mt)) {
    mt_label = switch(mt,
      ggm = "GGM (Gaussian Graphical Model)",
      omrf = "OMRF (Ordinal Markov Random Field)",
      mixed_mrf = sprintf(
        "Mixed MRF (%d discrete, %d continuous)",
        arguments$num_discrete, arguments$num_continuous
      ),
      mt
    )
    cat(paste0(" Model: ", mt_label, "\n"))
  }

  # Dataset info
  cat(paste0(" Number of variables: ", arguments$num_variables, "\n"))
  if(isTRUE(arguments$standardize)) {
    cat(" Prior standardization: enabled\n")
  }
  if(isTRUE(arguments$na_impute)) {
    cat(paste0(" Number of cases: ", arguments$num_cases, " (missings imputed)\n"))
  } else {
    cat(paste0(" Number of cases: ", arguments$num_cases, "\n"))
  }

  # Iterations and chains
  if(!is.null(arguments$num_chains)) {
    total_iter = arguments$iter * arguments$num_chains
    cat(paste0(" Number of post-burnin MCMC iterations: ", total_iter, "\n"))
    cat(paste0(" Number of MCMC chains: ", arguments$num_chains, "\n"))
  } else {
    cat(paste0(" Number of post-burnin MCMC iterations: ", arguments$iter, "\n"))
  }

  cat("Use the `summary()` function for posterior summaries and chain diagnostics.\n")
  cat("See the `easybgm` package for summary and plotting tools.\n")
  invisible(x)
}


#' @name summary.bgms
#' @title Summary method for `bgms` objects
#'
#' @description Returns posterior summaries and diagnostics for a fitted `bgms` model.
#'
#' @param object An object of class `bgms`.
#' @param ... Currently ignored.
#'
#' @return An object of class `summary.bgms` with posterior summaries.
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' summary(fit)
#' }
#'
#' @seealso [bgm()], [print.bgms()], [coef.bgms()]
#' @family posterior-methods
#'
#' @export
summary.bgms = function(object, ...) {
  arguments = extract_arguments(object)

  if(!is.null(object$posterior_summary_main) && !is.null(object$posterior_summary_pairwise)) {
    mt = arguments$model_type
    main_label = switch(mt,
      ggm       = "Precision diagonal (quadratic effects):",
      omrf      = "Category thresholds:",
      mixed_mrf = "Main effects (discrete thresholds and continuous means/precisions):",
      "Main effects:"
    )
    out = list(
      main = object$posterior_summary_main,
      pairwise = object$posterior_summary_pairwise
    )
    attr(out, "main_label") = main_label

    if(!is.null(object$posterior_summary_indicator)) {
      out$indicator = object$posterior_summary_indicator
    }

    if(!is.null(object$posterior_summary_pairwise_allocations)) {
      out$allocations = object$posterior_summary_pairwise_allocations
      out$mean_allocations = object$posterior_mean_allocations
      out$mode_allocations = object$posterior_mode_allocations
      out$num_blocks = object$posterior_num_blocks
    }

    class(out) = "summary.bgms"
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
print.summary.bgms = function(x, digits = 3, ...) {
  cat("Posterior summaries from Bayesian estimation:\n\n")

  if(!is.null(x$main)) {
    main_label = attr(x, "main_label") %||% "Main effects:"
    cat(main_label, "\n")
    print(round(head(x$main, 6), digits = digits))
    if(nrow(x$main) > 6) cat("... (use `summary(fit)$main` to see full output)\n")
    cat("\n")
  }

  if(!is.null(x$pairwise)) {
    cat("Pairwise interactions:\n")
    pair = head(x$pairwise, 6)
    pair[] = lapply(pair, function(col) ifelse(is.na(col), "", round(col, digits)))
    print(pair)
    # print(round(head(x$pairwise, 6), digits = digits))
    if(nrow(x$pairwise) > 6) cat("... (use `summary(fit)$pairwise` to see full output)\n")
    if(!is.null(x$indicator)) {
      cat("Note: NA values are suppressed in the print table. They occur here when an \n")
      cat("indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;\n")
      cat("`summary(fit)$pairwise` still contains the NA values.\n")
    }
    cat("\n")
  }

  if(!is.null(x$indicator)) {
    cat("Inclusion probabilities:\n")
    ind = head(x$indicator, 6)
    ind[] = lapply(ind, function(col) ifelse(is.na(col), "", round(col, digits)))
    print(ind)
    if(nrow(x$indicator) > 6) cat("... (use `summary(fit)$indicator` to see full output)\n")
    cat("Note: NA values are suppressed in the print table. They occur when an indicator\n")
    cat("was constant (all 0 or all 1) across all iterations, so sd/mcse/n_eff/Rhat\n")
    cat("are undefined; `summary(fit)$indicator` still contains the NA values.\n\n")
  }

  if(!is.null(x$allocations)) {
    cat("Pairwise node co-clustering proportion:\n")
    print(round(head(x$allocations, 6), digits = digits))
    if(nrow(x$allocations) > 6) cat("... (use `summary(fit)$allocations` to see full output)\n")
    cat("\n")
  }

  if(!is.null(x$mean_allocations)) {
    cat("Mean posterior node allocation vector:\n")
    print(round(head(x$mean_allocations, 6), digits = digits))
    cat("Mode posterior node allocation vector:\n")
    print(round(head(x$mode_allocations, 6), digits = digits))
    cat("\n")
  }

  if(!is.null(x$num_blocks)) {
    cat("Number of blocks and their posterior probability :\n")
    print(round(head(x$num_blocks, 6), digits = digits))
    if(nrow(x$num_blocks) > 6) cat("... (use `summary(fit)$num_blocks` to see full output)\n")
    cat("\n")
  }

  cat("Use `summary(fit)$<component>` to access full results.\n")
  cat("See the `easybgm` package for other summary and plotting tools.\n")
}


#' @title Extract Coefficients from a bgms Object
#' @name coef.bgms
#' @description Returns the posterior mean main effects, pairwise effects, and edge inclusion indicators from a \code{bgms} model fit.
#'
#' @param object An object of class \code{bgms}.
#' @param ... Ignored.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{main}{Posterior mean of the main-effect parameters. For GGM models
#'     this contains the precision matrix diagonal (quadratic effects, not main
#'     effects). For OMRF models this is a numeric matrix (p x max_categories)
#'     of category thresholds. For mixed MRF models this is a list with
#'     `$discrete` (p x max_categories matrix) and `$continuous` (q x 2 matrix
#'     with columns "mean" and "precision").}
#'   \item{pairwise}{Posterior mean of the pairwise interaction matrix.}
#'   \item{indicator}{Posterior mean of the edge inclusion indicators (if available).}
#' }
#'
#' @examples
#' \donttest{
#' fit = bgm(x = Wenchuan[, 1:3])
#' coef(fit)
#' }
#'
#' @seealso [bgm()], [print.bgms()], [summary.bgms()]
#' @family posterior-methods
#'
#' @export
coef.bgms = function(object, ...) {
  out = list(
    main = object$posterior_mean_main,
    pairwise = object$posterior_mean_pairwise
  )
  if(!is.null(object$posterior_mean_indicator)) {
    out$indicator = object$posterior_mean_indicator
  }

  if(!is.null(object$posterior_mean_allocations)) {
    out$mean_allocations = object$posterior_mean_allocations
    out$mode_allocations = object$posterior_mode_allocations
    out$num_blocks = object$posterior_num_blocks
  }

  return(out)
}


.warning_state = new.env(parent = emptyenv())
.warning_state$issued = FALSE

warning_once = function(msg) {
  if(!.warning_state$issued) {
    warning(msg, call. = FALSE)
    .warning_state$issued = TRUE
  }
}
