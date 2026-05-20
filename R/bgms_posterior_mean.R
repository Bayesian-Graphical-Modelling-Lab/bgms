#' Sign-corrected posterior means for a bgms fit
#'
#' Computes posterior means using Lyne (2015)'s sign-corrected ergodic
#' estimator. For a functional \eqn{f(K)} of the precision matrix,
#'
#' \deqn{E[f(K) \mid Y] \;\approx\; \frac{\sum_i \mathrm{sign}_i \cdot f(K_i)}{\sum_i \mathrm{sign}_i}}
#'
#' where \eqn{\mathrm{sign}_i} is the sign of the Russian-Roulette
#' estimator of \eqn{1/Z(\Gamma)} at iteration \eqn{i}. The correction
#' is required for the ergodic averages to converge to the true
#' posterior expectation when running the hierarchical-spec GGM
#' (\code{graph_prior_spec = "hierarchical"}) under settings that
#' permit signed V values.
#'
#' Under operational tunings (low \eqn{\delta}, reasonably large
#' \eqn{\kappa}) sign is identically \eqn{+1} and the correction
#' collapses to the plain posterior mean. Sign correction matters
#' primarily at high \eqn{\delta} or aggressive \eqn{\kappa}. The
#' diagnostic vector lives at \code{fit$v_ratio_diagnostics$sign}.
#'
#' For joint-spec fits, or for any fit that has no
#' \code{v_ratio_diagnostics} field, this function returns the
#' \code{posterior_mean_*} fields unchanged.
#'
#' @param fit  A bgm() fit object.
#' @return A list with components:
#'   \describe{
#'     \item{\code{main}}{Sign-corrected posterior mean of main effects.
#'       \code{NULL} for GGM (which has no main effects).}
#'     \item{\code{pairwise}}{Sign-corrected posterior mean of pairwise
#'       associations, as a symmetric \eqn{p \times p} matrix. For GGM
#'       this is on the association scale \eqn{-K_{ij}/2}.}
#'     \item{\code{indicator}}{Sign-corrected posterior mean of edge
#'       indicators \eqn{\gamma_{ij}}, as a symmetric \eqn{p \times p}
#'       matrix. Present only when the fit used edge selection.}
#'     \item{\code{residual_variance}}{Sign-corrected posterior mean
#'       of the residual variance \eqn{1/K_{ii}}. Present only for GGM
#'       fits.}
#'   }
#' @export
bgms_posterior_mean = function(fit) {
  diag = fit$v_ratio_diagnostics

  # No sign data → identity fall-through to existing posterior means.
  # (Joint-spec fits, or hierarchical fits before F2 wiring landed.)
  if(is.null(diag)) {
    return(list(
      main              = fit$posterior_mean_main,
      pairwise          = fit$posterior_mean_pairwise,
      indicator         = fit$posterior_mean_indicator,
      residual_variance = fit$posterior_mean_residual_variance
    ))
  }

  raw = fit$raw_samples
  signs = unlist(diag$sign)
  sum_signs = sum(signs)
  if(sum_signs == 0) {
    stop(
      "All sign(V) entries sum to zero; the sign-corrected posterior ",
      "mean is undefined. This usually indicates that the V/RR estimator ",
      "is operating outside its convergence regime - try refitting with ",
      "a larger kappa (passed via z_ratio_tuning).",
      call. = FALSE
    )
  }

  # Pool per-chain samples and compute sign-weighted column means in one shot.
  # Chain order in `signs` must match chain order in `raw$*` (both unlisted
  # / rbind'd over the same lapply order).
  weighted_col_means = function(samples_per_chain) {
    if(is.null(samples_per_chain)) return(NULL)
    pooled = do.call(rbind, samples_per_chain)
    colSums(pooled * signs) / sum_signs
  }

  # Whether this is a GGM fit (the only model class that produces sign data
  # currently). Heuristic: GGM has residual_variance, OMRF/mixed do not.
  is_continuous = !is.null(fit$posterior_mean_residual_variance)
  num_variables = nrow(fit$posterior_mean_pairwise)

  out = list()

  # Pairwise → symmetric matrix in association scale (GGM: precision * -0.5).
  pair_means = weighted_col_means(raw$pairwise)
  pairwise_mat = matrix(0, num_variables, num_variables,
                        dimnames = dimnames(fit$posterior_mean_pairwise))
  pairwise_mat[lower.tri(pairwise_mat)] = pair_means
  pairwise_mat = pairwise_mat + t(pairwise_mat)
  if(is_continuous) pairwise_mat = -0.5 * pairwise_mat
  out$pairwise = pairwise_mat

  # Main: GGM has no main effects; OMRF/mixed don't (yet) support
  # hierarchical V-correction. Carry the existing field through for shape
  # parity with the field of the same name on the fit object. Use
  # `out["main"] = list(NULL)` so the slot name survives even when the
  # value is NULL (plain `out$main = NULL` would drop the slot).
  out["main"] = list(fit$posterior_mean_main)

  # Indicator (when edge selection was active).
  if(!is.null(raw$indicator)) {
    ind_means = weighted_col_means(raw$indicator)
    indicator_mat = matrix(0, num_variables, num_variables,
                           dimnames = dimnames(fit$posterior_mean_indicator))
    indicator_mat[lower.tri(indicator_mat)] = ind_means
    indicator_mat = indicator_mat + t(indicator_mat)
    out$indicator = indicator_mat
  } else {
    out["indicator"] = list(NULL)
  }

  # Residual variance for GGM: sign-correct on the per-sample functional
  # 1/K_ii (matches the per-sample-inversion form used by the plain field
  # to avoid Jensen bias).
  if(is_continuous && !is.null(raw$main)) {
    inv_main_per_chain = lapply(raw$main, function(m) 1 / m)
    rv = weighted_col_means(inv_main_per_chain)
    names(rv) = names(fit$posterior_mean_residual_variance)
    out$residual_variance = rv
  } else {
    out["residual_variance"] = list(NULL)
  }

  out
}
