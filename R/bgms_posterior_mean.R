#' Posterior means for a bgms fit
#'
#' Returns the posterior-mean fields of a \code{bgm()} fit in a flat
#' named list. Provided as a stable accessor that does not depend on S7
#' property names.
#'
#' For GGM and mixed MRF fits, \code{pairwise} is the partial-correlation
#' style association \eqn{-K_{ij}/2} (i.e. the value already stored in
#' \code{fit$posterior_mean_pairwise}, not the raw precision entry).
#'
#' @param fit  A bgm() fit object.
#' @return A list with components:
#'   \describe{
#'     \item{\code{main}}{Posterior mean of main effects. \code{NULL}
#'       for GGM (which has no main effects).}
#'     \item{\code{pairwise}}{Posterior mean of pairwise associations,
#'       symmetric \eqn{p \times p}.}
#'     \item{\code{indicator}}{Posterior mean of edge indicators
#'       \eqn{\gamma_{ij}}. Present only when edge selection was active.}
#'     \item{\code{residual_variance}}{Posterior mean of the residual
#'       variance \eqn{1/K_{ii}}. Present only for GGM fits.}
#'   }
#' @export
bgms_posterior_mean = function(fit) {
  list(
    main              = fit$posterior_mean_main,
    pairwise          = fit$posterior_mean_pairwise,
    indicator         = fit$posterior_mean_indicator,
    residual_variance = fit$posterior_mean_residual_variance
  )
}
