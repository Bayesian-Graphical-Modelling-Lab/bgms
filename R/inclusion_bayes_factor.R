# ==============================================================================
# Inclusion Bayes-factor support for Graphical G-prior fits
# ==============================================================================
#
# Surfaces per-edge inclusion Bayes factors on the bgm() return for fits
# using interaction_prior = graphical_g_prior(). The marginal prior P(γ_ij)
# under the joint specification carries an implicit Z̃(Γ) weighting that
# cancels in the ratio P(γ_ij | data) / P(γ_ij) provided BOTH probabilities
# come from the same chain machinery. We therefore estimate the prior PIP
# empirically by running a second, prior-only chain on the same V_ij as
# the data fit (achieved by zeroing n_ and suf_stat_ after V_ij has been
# cached at enable_gg_prior() time).
#
# The prior-only chain is run lazily on the first summary()/coef() call
# and cached on the fit's cache environment, so subsequent calls are free.
# ==============================================================================


# ------------------------------------------------------------------
# get_fit_spec
# ------------------------------------------------------------------
get_fit_spec = function(fit) {
  if(inherits(fit, "S7_object")) {
    S7::prop(fit, ".bgm_spec")
  } else {
    .subset2(fit, ".bgm_spec")
  }
}


# ------------------------------------------------------------------
# fit_uses_graphical_g_prior
# ------------------------------------------------------------------
fit_uses_graphical_g_prior = function(fit) {
  spec = get_fit_spec(fit)
  if(is.null(spec)) return(FALSE)
  identical(spec$prior$interaction_prior_type, "graphical_g")
}


# ------------------------------------------------------------------
# run_prior_only_chain
# ------------------------------------------------------------------
# Runs an MCMC chain on the prior alone (likelihood muted via the C++
# `prior_only` flag) using the fit's existing spec. Returns the per-edge
# inclusion proportions averaged over the chain.
#
# @param spec        A validated bgm_spec (currently GGM-only).
# @param iter        Number of post-warmup iterations (default 2000).
# @param warmup      Warmup length (default 500).
# @param seed_offset Added to spec$sampler$seed for reproducibility.
#
# Returns: a numeric vector of length p*(p-1)/2 with per-edge prior PIPs,
#   in upper-triangle row-major order (same layout as the fit's
#   posterior_mean_indicator off-diagonal entries).
# ------------------------------------------------------------------
run_prior_only_chain = function(spec, iter = 2000L, warmup = 500L,
                                seed_offset = 0L) {
  stopifnot(inherits(spec, "bgm_spec"))
  if(!identical(spec$model_type, "ggm")) {
    stop("run_prior_only_chain: only GGM specs are supported.")
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
  # gradient_engine on the precision block, and exercising NUTS with a
  # muted likelihood under edge selection adds no value for PIP
  # estimation (the slab + diag + edge prior alone drive the chain).
  prior_spec$sampler$update_method = "adaptive-metropolis"
  prior_spec$sampler$target_accept = 0.44

  raw = run_sampler(prior_spec)
  if(length(raw) == 0L) {
    stop("run_prior_only_chain: empty raw output.")
  }
  ind = raw[[1L]]$indicator_samples
  if(is.null(ind)) {
    stop("run_prior_only_chain: prior chain did not return indicator samples.")
  }

  # ind is (p(p+1)/2) x n_iter with main-then-pair layout matching the data
  # chain. Reduce to off-diagonals (i < j) only, in the same column order
  # that build_output_bgm uses for posterior_mean_indicator.
  num_variables = spec$data$num_variables
  diag_idx = integer(num_variables)
  offdiag_idx = integer(num_variables * (num_variables - 1L) / 2L)
  pos = 0L; di = 0L; oi = 0L
  for(i in seq_len(num_variables)) {
    for(j in i:num_variables) {
      pos = pos + 1L
      if(i == j) {
        di = di + 1L
        diag_idx[di] = pos
      } else {
        oi = oi + 1L
        offdiag_idx[oi] = pos
      }
    }
  }
  rowMeans(ind[offdiag_idx, , drop = FALSE])
}


# ------------------------------------------------------------------
# ensure_inclusion_bayes_factor
# ------------------------------------------------------------------
# Computes prior PIPs (and inclusion BFs) for a fit and caches them on
# the fit's cache environment. Idempotent: only runs the prior-only
# chain on first call.
#
# Active only when the fit uses graphical_g_prior. For other interaction
# priors the marginal P(γ_ij) = inclusion_probability is analytic and
# does not require a prior-only chain — those cases are handled directly
# in summary.bgms without going through here.
# ------------------------------------------------------------------
ensure_inclusion_bayes_factor = function(fit) {
  if(!fit_uses_graphical_g_prior(fit)) {
    return(invisible(NULL))
  }
  cache = get_fit_cache(fit)
  if(is.null(cache)) return(invisible(NULL))
  if(!is.null(cache$gg_prior_pips)) return(invisible(NULL))

  spec = get_fit_spec(fit)
  if(is.null(spec)) return(invisible(NULL))
  if(!isTRUE(spec$prior$edge_selection)) return(invisible(NULL))

  prior_pips = run_prior_only_chain(spec)
  cache$gg_prior_pips = prior_pips
  invisible(NULL)
}


# ------------------------------------------------------------------
# inclusion_bayes_factor
# ------------------------------------------------------------------
# BF_10(γ_ij) = (post / (1 - post)) / (prior / (1 - prior)).
# Returns Inf when post -> 1 with prior < 1, NaN when both endpoints
# saturate.
# ------------------------------------------------------------------
inclusion_bayes_factor = function(post_pip, prior_pip) {
  post_odds  = post_pip  / (1 - post_pip)
  prior_odds = prior_pip / (1 - prior_pip)
  post_odds / prior_odds
}
