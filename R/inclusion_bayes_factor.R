# ==============================================================================
# inclusion_bayes_factor() — extractor for per-edge inclusion Bayes factors
# ==============================================================================
#
# For a GGM fit with edge_selection = TRUE, returns posterior PIP, prior PIP,
# and inclusion BF for each potential edge. The prior PIP is estimated
# empirically by running a short prior-only chain that targets the joint
# prior p(K, Gamma) on the same V_ij as the data fit (the C++ prior-only
# flag zeros n_ and suf_stat_ after V_ij has been cached).
#
# This is the joint-specification BF: ratio of posterior odds to prior
# odds, both conditional on the same V_ij(Y). Works for any continuous
# interaction prior (cauchy, normal, graphical_g, ...) — the prior-only
# chain machinery is prior-agnostic.
#
# Results are cached on fit$cache so subsequent calls are free unless
# `recompute = TRUE`.
# ==============================================================================


# ------------------------------------------------------------------
# Accessors
# ------------------------------------------------------------------
get_fit_spec = function(fit) {
  if(inherits(fit, "S7_object")) {
    S7::prop(fit, ".bgm_spec")
  } else {
    .subset2(fit, ".bgm_spec")
  }
}


# ------------------------------------------------------------------
# run_prior_only_chain (internal)
# ------------------------------------------------------------------
# Runs a chain on the joint prior p(K, Gamma) using the fit's existing
# spec. Likelihood is muted via C++ `prior_only = TRUE`, which zeroes
# n_ and suf_stat_ after V_ij has been cached (GG-prior path) or after
# the model is constructed (non-GG-prior path — V_ij isn't built so
# zeroing is harmless). Returns the per-edge inclusion proportions in
# the same upper-triangle row-major order as posterior_mean_indicator.
# ------------------------------------------------------------------
run_prior_only_chain = function(spec, iter = 2000L, warmup = 500L,
                                seed_offset = 0L) {
  stopifnot(inherits(spec, "bgm_spec"))
  if(!identical(spec$model_type, "ggm")) {
    stop("inclusion_bayes_factor() currently supports GGM fits only ",
         "(continuous data). Got model_type = '", spec$model_type, "'.")
  }
  prior_spec = spec
  prior_spec$prior$prior_only = TRUE
  prior_spec$sampler$iter   = as.integer(iter)
  prior_spec$sampler$warmup = as.integer(warmup)
  prior_spec$sampler$chains = 1L
  prior_spec$sampler$cores  = 1L
  prior_spec$sampler$seed   = as.integer(spec$sampler$seed + seed_offset)
  prior_spec$sampler$progress_type = 0L
  prior_spec$sampler$display_progress = "none"
  # Force adaptive-Metropolis for the prior-only chain. NUTS uses the
  # gradient_engine on the precision block and adds no value for PIP
  # estimation with the likelihood muted.
  prior_spec$sampler$update_method = "adaptive-metropolis"
  prior_spec$sampler$target_accept = 0.44

  raw = run_sampler(prior_spec)
  if(length(raw) == 0L) {
    stop("inclusion_bayes_factor(): empty prior-only chain output.")
  }
  ind = raw[[1L]]$indicator_samples
  if(is.null(ind)) {
    stop("inclusion_bayes_factor(): prior chain did not return ",
         "indicator samples — was edge_selection = TRUE on the fit?")
  }

  # ind is (p(p+1)/2) x n_iter with main-then-pair layout matching the
  # data chain. Reduce to off-diagonals (i < j) only, in the same order
  # build_output_bgm uses for posterior_mean_indicator.
  num_variables = spec$data$num_variables
  offdiag_idx = integer(num_variables * (num_variables - 1L) / 2L)
  pos = 0L; oi = 0L
  for(i in seq_len(num_variables)) {
    for(j in i:num_variables) {
      pos = pos + 1L
      if(i != j) {
        oi = oi + 1L
        offdiag_idx[oi] = pos
      }
    }
  }
  rowMeans(ind[offdiag_idx, , drop = FALSE])
}


# ------------------------------------------------------------------
# Edge-name helper
# ------------------------------------------------------------------
# Returns the p*(p-1)/2 character vector of edge names in the same
# row-major upper-triangle order used by build_output_bgm.
gg_edge_names = function(spec) {
  vars = spec$data$data_columnnames
  p = length(vars)
  out = character(p * (p - 1L) / 2L)
  k = 0L
  for(i in seq_len(p - 1L)) {
    for(j in (i + 1L):p) {
      k = k + 1L
      out[k] = paste0(vars[i], "-", vars[j])
    }
  }
  out
}


# ------------------------------------------------------------------
# inclusion_bayes_factor() — main extractor
# ------------------------------------------------------------------

#' Per-edge inclusion Bayes factors for a bgm() fit
#'
#' For a GGM fit run with `edge_selection = TRUE`, returns the
#' posterior inclusion probability, the empirical prior inclusion
#' probability, and the inclusion Bayes factor for each potential edge.
#' The prior PIP is estimated from a short prior-only chain that
#' samples from the joint prior on `(K, Gamma)` using the same `V_ij`
#' as the data fit. The resulting BF is the *joint-specification* BF
#' (ratio of posterior odds to prior odds, both conditional on
#' `V_ij(Y)`); it is consistent across interaction-prior families
#' because both PIPs come from the same chain machinery.
#'
#' Results are cached on `fit$cache` so repeated calls are free unless
#' `recompute = TRUE`.
#'
#' @param fit A `bgms` fit from `bgm()` with `edge_selection = TRUE`
#'   on continuous data (GGM).
#' @param iter Integer. Number of post-warmup iterations for the
#'   prior-only chain. Default `2000`.
#' @param warmup Integer. Warmup length for the prior-only chain.
#'   Default `500`.
#' @param recompute Logical. If `TRUE`, recompute the prior PIPs even
#'   if a cached result is present. Default `FALSE`.
#'
#' @return A data frame with one row per potential off-diagonal edge
#'   and columns
#'   \describe{
#'     \item{`edge`}{Character. Edge name (e.g. `"V1-V2"`).}
#'     \item{`posterior_inclusion`}{Numeric. Posterior PIP from the
#'       fitted chain.}
#'     \item{`prior_inclusion`}{Numeric. Empirical prior PIP from the
#'       prior-only chain.}
#'     \item{`bf`}{Numeric. `(post / (1-post)) / (prior / (1-prior))`.
#'       `Inf` when `posterior -> 1` with `prior < 1`.}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- bgm(X, variable_type = "continuous",
#'            interaction_prior = graphical_g_prior(),
#'            edge_selection = TRUE)
#' inclusion_bayes_factor(fit)
#' }
#'
#' @export
inclusion_bayes_factor = function(fit, iter = 2000L, warmup = 500L,
                                  recompute = FALSE) {
  spec = get_fit_spec(fit)
  if(is.null(spec)) {
    stop("inclusion_bayes_factor(): fit has no embedded .bgm_spec. ",
         "Re-fit with the current bgms version.")
  }
  if(!identical(spec$model_type, "ggm")) {
    stop("inclusion_bayes_factor() currently supports GGM fits only ",
         "(continuous data). Got model_type = '", spec$model_type, "'.")
  }
  if(!isTRUE(spec$prior$edge_selection)) {
    stop("inclusion_bayes_factor() requires a fit run with ",
         "edge_selection = TRUE.")
  }

  cache = get_fit_cache(fit)

  # Posterior PIPs come from the existing fit object (no recomputation).
  post_mat = if(inherits(fit, "S7_object")) {
    S7::prop(fit, "posterior_mean_indicator")
  } else {
    .subset2(fit, "posterior_mean_indicator")
  }
  post_pip = post_mat[upper.tri(post_mat)]

  prior_pip = if(!recompute && !is.null(cache) &&
                 !is.null(cache$inclusion_bf_prior_pips) &&
                 length(cache$inclusion_bf_prior_pips) == length(post_pip)) {
    cache$inclusion_bf_prior_pips
  } else {
    p = run_prior_only_chain(spec, iter = iter, warmup = warmup)
    if(!is.null(cache)) cache$inclusion_bf_prior_pips = p
    p
  }

  post_odds  = post_pip  / (1 - post_pip)
  prior_odds = prior_pip / (1 - prior_pip)
  bf = post_odds / prior_odds

  data.frame(
    edge                = gg_edge_names(spec),
    posterior_inclusion = post_pip,
    prior_inclusion     = prior_pip,
    bf                  = bf,
    stringsAsFactors    = FALSE
  )
}
