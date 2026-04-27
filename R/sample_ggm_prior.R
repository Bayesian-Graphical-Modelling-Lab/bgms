#' Sample from the GGM Prior
#'
#' Draws precision matrices \eqn{K} from the prior of a Gaussian graphical
#' model using the same constrained NUTS sampler that drives \code{\link{bgm}}
#' for continuous data. The likelihood is omitted (\eqn{n = 0},
#' \eqn{S = 0}), so the chain targets the prior alone.
#'
#' Off-diagonals are placed on the association scale
#' \eqn{K_{yy,ij} = -K_{ij}/2} and assigned the supplied
#' \code{interaction_prior}. A \code{normal_prior(scale = s)} therefore
#' constrains \eqn{K_{yy}} with standard deviation \eqn{s}, equivalent to a
#' \eqn{\textrm{Normal}(0, 2s)} prior on \eqn{K_{ij}} itself. The
#' diagonals \eqn{K_{ii}} are drawn from the supplied
#' \code{precision_scale_prior}. When \code{edge_indicators} is supplied,
#' off-diagonals at excluded positions are constrained to zero throughout
#' the chain.
#'
#' @param p Integer. Dimension of the precision matrix (\eqn{p \ge 2}).
#' @param n_samples Integer. Number of post-warmup draws to keep.
#' @param n_warmup Integer. NUTS warmup iterations. Default \code{1000}.
#' @param interaction_prior A \code{bgms_parameter_prior} for the
#'   off-diagonal entries. Use \code{\link{cauchy_prior}()} or
#'   \code{\link{normal_prior}()}; \code{\link{beta_prime_prior}()} is not
#'   supported here. Default: \code{cauchy_prior(scale = 2.5)}.
#' @param precision_scale_prior A \code{bgms_scale_prior} for the diagonal
#'   entries of \eqn{K}. Use \code{\link{gamma_prior}()} or
#'   \code{\link{exponential_prior}()}. Default: \code{gamma_prior(1, 1)}.
#' @param step_size Positive numeric. Initial NUTS step size used to seed
#'   dual-averaging adaptation. Default \code{0.1}.
#' @param max_depth Integer. Maximum NUTS tree depth. Default \code{10}.
#' @param seed Integer. RNG seed for the chain. Default \code{1L}.
#' @param verbose Logical. If \code{TRUE} (default), print a progress bar.
#' @param edge_indicators Optional integer \eqn{p \times p} matrix with
#'   \code{1} = edge included, \code{0} = excluded. Must be symmetric with
#'   \code{1}s on the diagonal. Default: full graph (all edges included).
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{K_offdiag}}{Numeric matrix of size
#'       \code{n_samples} x \code{p * (p - 1) / 2} containing the upper-triangle
#'       off-diagonal entries of \eqn{K} for each draw, in column-major order
#'       \eqn{(K_{12}, K_{13}, K_{23}, K_{14}, \ldots)}. Excluded edges are
#'       returned as \code{0}.}
#'     \item{\code{K_diag}}{Numeric matrix of size
#'       \code{n_samples} x \code{p} containing the diagonal entries
#'       \eqn{K_{11}, \ldots, K_{pp}}.}
#'     \item{\code{offdiag_names}}{Character vector of length
#'       \code{p * (p - 1) / 2} naming the columns of \code{K_offdiag}
#'       (e.g. \code{"K_1_2"}).}
#'     \item{\code{diag_names}}{Character vector of length \code{p} naming
#'       the columns of \code{K_diag}.}
#'     \item{\code{step_size}}{The (initial) NUTS step size used.}
#'     \item{\code{edge_indicators}}{The integer edge-indicator matrix used
#'       (full graph if not supplied).}
#'   }
#'
#' @seealso \code{\link{cauchy_prior}}, \code{\link{normal_prior}},
#'   \code{\link{gamma_prior}}, \code{\link{exponential_prior}},
#'   \code{\link{bgm}}
#'
#' @examples
#' \donttest{
#' # Default Cauchy(0, 2.5) off-diagonal, Gamma(1, 1) diagonal, p = 4.
#' draws = sample_ggm_prior(
#'   p = 4, n_samples = 200, n_warmup = 200,
#'   verbose = FALSE
#' )
#' dim(draws$K_offdiag) # 200 x 6
#' colnames(draws$K_offdiag) = draws$offdiag_names
#' head(draws$K_offdiag)
#'
#' # Sparser graph: drop the (1, 4) edge.
#' E = matrix(1L, 4, 4)
#' E[1, 4] = E[4, 1] = 0L
#' draws = sample_ggm_prior(
#'   p = 4, n_samples = 200, n_warmup = 200,
#'   edge_indicators = E, verbose = FALSE
#' )
#' colnames(draws$K_offdiag) = draws$offdiag_names
#' all(draws$K_offdiag[, "K_1_4"] == 0) # TRUE
#' }
#' @export
sample_ggm_prior = function(
  p,
  n_samples,
  n_warmup = 1000L,
  interaction_prior = cauchy_prior(scale = 2.5),
  precision_scale_prior = gamma_prior(shape = 1, rate = 1),
  step_size = 0.1,
  max_depth = 10L,
  seed = 1L,
  verbose = TRUE,
  edge_indicators = NULL
) {
  validate_positive_integer(p, "p", min_value = 2L)
  validate_positive_integer(n_samples, "n_samples", min_value = 1L)
  validate_positive_integer(n_warmup, "n_warmup", min_value = 0L)
  validate_positive_integer(max_depth, "max_depth", min_value = 1L)
  validate_finite_scalar(step_size, "step_size", positive = TRUE)
  validate_positive_integer(seed, "seed", min_value = 0L)
  if(!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.")
  }

  ip = unpack_interaction_prior(interaction_prior)
  if(identical(ip$interaction_prior_type, "beta-prime")) {
    stop(
      "beta_prime_prior() is not supported for 'interaction_prior' in ",
      "sample_ggm_prior(). Use cauchy_prior() or normal_prior()."
    )
  }
  sp = unpack_scale_prior(precision_scale_prior)

  edge_indicators = validate_ggm_prior_edge_indicators(edge_indicators, p)

  sample_ggm_prior_cpp(
    p                        = as.integer(p),
    n_samples                = as.integer(n_samples),
    n_warmup                 = as.integer(n_warmup),
    pairwise_scale           = ip$pairwise_scale,
    interaction_prior_type   = ip$interaction_prior_type,
    scale_prior_type         = sp$scale_prior_type,
    gamma_shape              = sp$scale_shape,
    gamma_rate               = sp$scale_rate,
    step_size                = step_size,
    max_depth                = as.integer(max_depth),
    seed                     = as.integer(seed),
    verbose                  = verbose,
    edge_indicators_nullable = edge_indicators
  )
}


# Internal helpers -------------------------------------------------------------

validate_positive_integer = function(x, name, min_value = 1L) {
  if(!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite integer.", name))
  }
  if(x != as.integer(x)) {
    stop(sprintf("'%s' must be an integer (got %s).", name, format(x)))
  }
  if(x < min_value) {
    stop(sprintf("'%s' must be >= %d.", name, as.integer(min_value)))
  }
  invisible(as.integer(x))
}

validate_finite_scalar = function(x, name, positive = FALSE) {
  if(!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("'%s' must be a single finite numeric.", name))
  }
  if(positive && x <= 0) {
    stop(sprintf("'%s' must be positive.", name))
  }
  invisible(x)
}

validate_ggm_prior_edge_indicators = function(edge_indicators, p) {
  if(is.null(edge_indicators)) {
    return(NULL)
  }
  if(!is.matrix(edge_indicators) ||
    nrow(edge_indicators) != p || ncol(edge_indicators) != p) {
    stop("'edge_indicators' must be a p x p matrix.")
  }
  if(any(is.na(edge_indicators))) {
    stop("'edge_indicators' must not contain NA values.")
  }
  vals = as.integer(edge_indicators)
  if(any(!vals %in% c(0L, 1L))) {
    stop("'edge_indicators' must contain only 0 or 1.")
  }
  E = matrix(vals, nrow = p, ncol = p)
  if(!isTRUE(all.equal(E, t(E)))) {
    stop("'edge_indicators' must be symmetric.")
  }
  if(!all(diag(E) == 1L)) {
    stop("'edge_indicators' must have 1s on the diagonal.")
  }
  E
}
