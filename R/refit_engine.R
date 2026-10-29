# ==============================================================================
# Refit engine for prior_sensitivity_check(): refits of a fit at other slab
# scales, the per-fit statistics the check classifies on, and the per-fit
# convergence diagnostics it reports.
# ==============================================================================
#
# A refit reuses the original fit's spec with the slab scale changed and a fresh
# seed. An ordinal (omrf) refit starts each chain from the original fit's final
# per-chain state, and a NUTS refit also from its adapted step size and
# diagonal metric; the GGM and mixed samplers take no initial state and start
# cold.
# ==============================================================================


# ------------------------------------------------------------------
# extract_warm_state
# ------------------------------------------------------------------
# Per-chain final state of a fit, in the exact vectorization the C++ sampler
# consumes (set_storage_vectorized_parameters + set_vectorized_indicator_
# parameters). For omrf the storage vector is c(main effects, all pairwise
# upper-triangle row-major); the indicator vector is the pairwise upper triangle.
#
# @param fit  A fitted bgms object (omrf) with edge selection.
#
# The returned object is sampler-agnostic: `parameters` is universal (every
# sampler consumes it), while the tuning fields are per method and optional --
# `step_sizes` and `inv_mass` are the NUTS metric (NULL/NA for a non-NUTS
# origin). A refit sampler that does not recognize a field simply ignores it.
# Warm starts use a dense all-edges-active graph, so no indicator state is
# carried (warm-starting a sparse graph would desynchronize the active-parameter
# dimension; see refit_at_scale).
#
# Returns: list(parameters = per-chain numeric vectors,
#               step_sizes  = per-chain NUTS step sizes (NA if none),
#               inv_mass    = per-chain NUTS diagonal metric (NULL if none)).
# ------------------------------------------------------------------
extract_warm_state = function(fit) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  if(!identical(spec$model_type, "omrf")) {
    stop(
      "Warm-started refits are currently implemented for ordinal (omrf) ",
      "models only; got model type '", spec$model_type, "'."
    )
  }
  nchains = length(raw$pairwise)
  parameters = vector("list", nchains)
  for(c in seq_len(nchains)) {
    ni = nrow(raw$pairwise[[c]])
    # c(main, all pairwise) is the storage vector, row for the last iteration.
    parameters[[c]] = as.numeric(c(raw$main[[c]][ni, ], raw$pairwise[[c]][ni, ]))
  }
  # Per-chain final NUTS step size (NA for non-NUTS fits) and diagonal metric
  # (NULL for non-NUTS fits), for warm-starting a NUTS refit's tuning state. A
  # fit saved before these properties existed lacks them entirely; treat that as
  # absent tuning state (the refit falls back to a cold metric/step size).
  get_prop = function(name) {
    tryCatch(S7::prop(fit, name), error = function(e) {
      tryCatch(fit[[name]], error = function(e2) NULL)
    })
  }
  step_sizes = get_prop("refit_step_sizes")
  inv_mass = get_prop("refit_inv_mass")
  list(parameters = parameters, step_sizes = step_sizes, inv_mass = inv_mass)
}


# ------------------------------------------------------------------
# refit_at_scale
# ------------------------------------------------------------------
# One refit of `fit` at absolute slab scale `scale`, reusing the original spec
# with that scale, the given run length, a warm start where one is supplied,
# and a fresh seed. Returns a fitted bgms object.
#
# @param fit          Original fitted bgms object.
# @param scale        Absolute interaction slab scale for this refit.
# @param warm_state   Output of extract_warm_state(fit), or NULL for a cold refit.
# @param warmup       Warmup iterations for the refit.
# @param iter         Sampling iterations for the refit.
# @param seed         Seed for the refit.
# @param cores        Threads for the refit (chains run in parallel internally).
# @param sampler      Optional update-method override ("nuts",
#                     "adaptive-metropolis", "gibbs"); NULL keeps the fit's.
# @param show_progress  Show the sampler's native progress bar for this refit.
# @param scale_field  Name of the prior field the scale is written to:
#                     "pairwise_scale" for bgm(), "difference_scale" for
#                     bgmCompare().
# @param diagonal_rate  Raw rate for the precision diagonal at this scale, or
#                     NULL for a model without a precision diagonal; see
#                     diagonal_eta().
# ------------------------------------------------------------------
refit_at_scale = function(fit, scale, warm_state, warmup, iter, seed,
                          cores = 1L, sampler = NULL, show_progress = FALSE,
                          scale_field = "pairwise_scale",
                          diagonal_rate = NULL) {
  spec = get_fit_spec(fit)

  spec$prior[[scale_field]] = scale
  # The spec stores the diagonal rate already resolved to the raw frame, so a
  # new slab scale leaves it alone unless the caller asks otherwise.
  if(!is.null(diagonal_rate)) {
    spec$prior$scale_rate = diagonal_rate
  }

  if(!is.null(sampler)) {
    spec$sampler$update_method = sampler
  }
  spec$sampler$warmup = as.integer(warmup)
  spec$sampler$iter = as.integer(iter)
  spec$sampler$seed = as.integer(seed)
  # The refit bypasses validate_sampler(), so it applies the core guard itself.
  spec$sampler$cores = normalize_parallel_cores(cores)
  # Render the sampler's own native display (a bar per chain, with the
  # warmup/sampling stage and ETA) for this refit when asked; otherwise stay
  # silent. The manager redraws its block in place, so a refit's chains do not
  # pile up as they run.
  spec$sampler$display_progress = if(show_progress) "per-chain" else "none"
  spec$sampler$progress_type = if(show_progress) 2L else 0L
  spec$sampler$progress_callback = NULL

  # Warm-start the continuous parameters only, on the default dense
  # (all-edges-active) start. Warm-starting a sparse indicator configuration
  # would desynchronize the sampler's active-parameter dimension from the full
  # storage dimension; a dense start matches the cold-start dimension exactly,
  # and the warmup re-settles the graph (which shifts with the scale anyway).
  #
  # A NUTS refit additionally starts from the origin fit's per-chain tuning
  # state (step size + diagonal metric). Non-NUTS refits ignore these fields.
  spec$initial_state = if(is.null(warm_state)) {
    NULL
  } else {
    state = list(parameters = warm_state$parameters)
    if(identical(spec$sampler$update_method, "nuts")) {
      if(!is.null(warm_state$step_sizes) && all(is.finite(warm_state$step_sizes))) {
        state$step_sizes = warm_state$step_sizes
      }
      if(!is.null(warm_state$inv_mass)) {
        state$inv_mass = warm_state$inv_mass
      }
    }
    state
  }

  raw = run_sampler(spec)
  build_output(spec, raw)
}


# ------------------------------------------------------------------
# diagonal_eta
# ------------------------------------------------------------------
# The standardized rate eta of the prior on the precision diagonal, held fixed
# across the sweep. At slab scale s the sampler runs that prior with raw rate
# eta / s.
#
# @param spec  The original fit's spec.
#
# Returns: eta, or NA for a model without a precision diagonal. A fit whose
# diagonal prior was given as a raw rate stops with an error.
# ------------------------------------------------------------------
diagonal_eta = function(spec) {
  prior = spec$prior
  has_diagonal = !is.null(prior$scale_rate) && !is.na(prior$scale_rate)
  if(!has_diagonal) {
    return(NA_real_)
  }
  if(is.null(prior$scale_eta) || is.na(prior$scale_eta)) {
    stop(
      "prior_sensitivity_check(): this fit gives the prior on the precision ",
      "diagonal as a raw rate. The check holds the standardized rate eta fixed ",
      "while the slab scale moves and has no sweep that holds a raw rate ",
      "fixed. Give the diagonal prior as eta (exponential_prior(eta = ) or ",
      "gamma_prior(eta = )), or refit at a few slab scales directly; see ",
      "'Fits with a raw-rate diagonal prior' in ?prior_sensitivity_check.",
      call. = FALSE
    )
  }
  prior$scale_eta
}


# ------------------------------------------------------------------
# resolve_refit_sampler
# ------------------------------------------------------------------
# Resolve the requested refit sampler against the original fit's update method.
# "same-as-fit" inherits; otherwise the named method is used. Returns the method
# and whether to recommend NUTS for speed (an AM/Gibbs original inherited as-is).
# ------------------------------------------------------------------
resolve_refit_sampler = function(fit, refit_sampler) {
  fit_method = get_fit_spec(fit)$sampler$update_method
  method = if(identical(refit_sampler, "same-as-fit")) fit_method else refit_sampler
  valid = c("nuts", "adaptive-metropolis", "gibbs")
  if(!method %in% valid) {
    stop(
      "'refit_sampler' must be \"same-as-fit\" or one of ",
      paste(sQuote(valid), collapse = ", "), "; got ", sQuote(method), "."
    )
  }
  # Recommend NUTS only when the user inherited a slower same-as-fit method.
  recommend_nuts = identical(refit_sampler, "same-as-fit") && method != "nuts"
  list(method = method, recommend_nuts = recommend_nuts, fit_method = fit_method)
}


# ------------------------------------------------------------------
# data_preferred_scale
# ------------------------------------------------------------------
# The scale the data prefer, read off the ORIGINAL fit's draws with no refit:
# s_hat is the root-mean-square of the included interactions, with an approximate
# 1/sqrt(2m) log-scale standard error (m = mean number of included edges). For
# the GGM the slab sits on Kyy = -0.5 * Omega, matching the reweighting path.
#
# Returns: list(s_hat, lo, hi, log_sd, m). NA-valued if no edge is ever included.
# ------------------------------------------------------------------
data_preferred_scale = function(fit) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  slab_factor = if(identical(spec$model_type, "ggm")) -0.5 else 1
  theta = slab_factor * do.call(rbind, raw$pairwise)
  gamma = do.call(rbind, raw$indicator)
  incl = gamma == 1
  n_incl = sum(incl)
  if(n_incl == 0L) {
    return(list(
      s_hat = NA_real_, lo = NA_real_, hi = NA_real_,
      log_sd = NA_real_, m = 0
    ))
  }
  s_hat = sqrt(mean(theta[incl]^2))
  m = mean(rowSums(incl))
  log_sd = 1 / sqrt(2 * max(m, 1))
  list(
    s_hat = s_hat, lo = s_hat * exp(-1.96 * log_sd),
    hi = s_hat * exp(1.96 * log_sd), log_sd = log_sd, m = m
  )
}


# ------------------------------------------------------------------
# verdict_from_lbf
# ------------------------------------------------------------------
# Map a natural log inclusion Bayes factor to a verdict at threshold
# lthr = log(t).
# ------------------------------------------------------------------
verdict_from_lbf = function(lbf, lthr) {
  out = rep("undecided", length(lbf))
  out[lbf >= lthr] = "presence"
  out[lbf <= -lthr] = "absence"
  out[is.na(lbf)] = NA_character_
  out
}


# ------------------------------------------------------------------
# refit_edge_stats
# ------------------------------------------------------------------
# Per-edge quantities from one refit: the Rao-Blackwellized inclusion
# probability (canonical), its natural log inclusion Bayes factor and verdict,
# the within-chain MCSE of the log BF (RB machinery, PR #182), the between-chain
# half-band (2 * SEM of the grand mean from chain spread, propagated to the log
# BF), per-chain verdict unanimity, and the zero-flip mask. Edge order is the
# native row-major upper triangle.
# ------------------------------------------------------------------
refit_edge_stats = function(fit, evidence_threshold) {
  raw = get_raw_samples(fit)
  # bgm() stores `indicator`, bgmCompare() stores `indicators`; read both
  # exactly.
  pnm = raw$parameter_names
  enm = pnm[["indicator"]] %||% pnm[["indicators"]] %||% pnm[["pairwise"]]
  M = length(raw$rb_inclusion)
  pc = sapply(raw$rb_inclusion, colMeans) # n_edges x M
  if(is.null(dim(pc))) pc = matrix(pc, ncol = M)
  pbar = rowMeans(pc)

  prior_q = if(inherits(fit, "bgmCompare")) {
    # The difference prior is exchangeable across difference indicators, so one
    # inclusion probability covers the main-effect and the pairwise families
    # alike. Recycled to the indicator count below.
    rep(difference_prior_inclusion(fit), length(enm))
  } else {
    prm = extract_prior_inclusion_probabilities(fit)
    prm[upper.tri(prm)][order_upper_tri_rowmajor(nrow(prm))]
  }
  prior_odds = prior_q / (1 - prior_q)

  ind = fit@posterior_summary_indicator
  mcse_pip = ind[, "mcse"]
  lthr = log(evidence_threshold)

  lbf = function(p) log((p / (1 - p)) / prior_odds)
  dfac = 1 / pmax(pbar * (1 - pbar), 1e-6) # d log BF / d p
  lbf_bar = lbf(pbar)
  mcse_lbf = mcse_pip * dfac
  se_between = apply(pc, 1, stats::sd) / sqrt(M)
  band_half = 2 * se_between * dfac

  verdict = verdict_from_lbf(lbf_bar, lthr)
  vc = apply(pc, 2, function(p) verdict_from_lbf(lbf(p), lthr)) # n_edges x M
  if(is.null(dim(vc))) vc = matrix(vc, ncol = M)
  unanimous = apply(vc, 1, function(v) length(unique(v[!is.na(v)])) <= 1L)
  zeroflip = apply(pc, 1, function(p) all(p <= 0) || all(p >= 1))

  list(
    edge = enm, pip = pbar, prior_odds = prior_odds, lbf = lbf_bar,
    verdict = verdict, mcse_lbf = mcse_lbf, band_half = band_half,
    unanimous = unanimous, zeroflip = zeroflip
  )
}


# ------------------------------------------------------------------
# finite_reduce
# ------------------------------------------------------------------
# Reduce the finite entries of a diagnostic vector, with a stated value for
# the case where there are none. min() and max() answer that case with a
# warning and an infinity, which reads as a value rather than as an absent
# one.
#
# @param values  Numeric vector, possibly empty or all NA.
# @param fn      The reduction.
# @param empty   Value to return when nothing is finite.
#
# Returns: the reduction, or `empty`.
# ------------------------------------------------------------------
finite_reduce = function(values, fn, empty) {
  values = values[is.finite(values)]
  if(length(values) == 0L) empty else fn(values)
}


# ------------------------------------------------------------------
# refit_diagnostics_row
# ------------------------------------------------------------------
# The convergence diagnostics of one fit, as the check reports them in $grid:
# the largest split R-hat and the smallest effective sample size over the
# continuous parameters (main effects and pairwise interactions; for a
# bgmCompare fit the baseline and difference blocks of both), the smallest
# effective sample size over the inclusion indicators, and for NUTS the
# smallest E-BFMI and the largest first/second-half energy variance ratio over
# the chains. NA where a quantity is not computed for the fit.
# ------------------------------------------------------------------
refit_diagnostics_row = function(fit) {
  rhats = numeric(0)
  esss = numeric(0)
  blocks = if(inherits(fit, "bgmCompare")) {
    c(
      "posterior_summary_main_baseline", "posterior_summary_main_differences",
      "posterior_summary_pairwise_baseline",
      "posterior_summary_pairwise_differences"
    )
  } else {
    c("posterior_summary_main", "posterior_summary_pairwise")
  }
  for(name in blocks) {
    block = tryCatch(S7::prop(fit, name), error = function(e) NULL)
    if(is.null(block)) next
    rhats = c(rhats, block[, "Rhat"])
    esss = c(esss, block[, "n_eff"])
  }
  ind = fit@posterior_summary_indicator
  wc = tryCatch(fit@nuts_diag$warmup_check, error = function(e) NULL)
  wi = tryCatch(fit@nuts_diag$summary$warmup_incomplete, error = function(e) NULL)
  list(
    rhat_max = finite_reduce(rhats, max, NA_real_),
    ess_min = finite_reduce(esss, min, NA_real_),
    inclusion_ess_min = finite_reduce(ind[, "n_eff"], min, NA_real_),
    ebfmi_min = if(is.null(wc)) {
      NA_real_
    } else {
      finite_reduce(c(wc$ebfmi_first_half, wc$ebfmi_second_half), min, NA_real_)
    },
    energy_variance_ratio_max = if(is.null(wc)) {
      NA_real_
    } else {
      finite_reduce(wc$var_ratio, max, NA_real_)
    },
    warmup_incomplete = if(is.null(wi)) NA else isTRUE(wi)
  )
}
