# ==============================================================================
# prior_sensitivity_check: scale-robustness of edge inclusion verdicts
# ==============================================================================
#
# Refit the model at a few multipliers of the fit's chosen slab scale, each
# refit at the original fit's run length, and compare every edge's verdict
# across the fits. The refits and the per-fit statistics live in
# refit_engine.R; the classification is a function of those statistics alone
# (sensitivity_classification).
# ==============================================================================


# ------------------------------------------------------------------
# verdict_from_bf
# ------------------------------------------------------------------
# Map an inclusion Bayes factor to an evidence verdict.
#
# @param bf         Numeric inclusion Bayes factor(s).
# @param threshold  Upper evidence threshold; lower is 1 / threshold.
#
# Returns: character vector of "presence" / "absence" / "undecided".
# ------------------------------------------------------------------
verdict_from_bf = function(bf, threshold) {
  out = rep("undecided", length(bf))
  out[bf >= threshold] = "presence"
  out[bf <= 1 / threshold] = "absence"
  out[is.na(bf)] = NA_character_
  out
}


#' @title Prior Sensitivity of Inclusion Verdicts
#'
#' @description
#' Refits the model at a few multipliers of the slab scale that gates the
#' indicators, each refit at the original fit's run length, and compares every
#' indicator's inclusion verdict across the fits.
#'
#' On a \code{\link{bgm}()} fit with edge selection the unit is the edge
#' indicator and the swept prior is the interaction slab. On a
#' \code{\link{bgmCompare}()} fit with difference selection the unit is the
#' difference indicator and the swept prior is the difference slab
#' (\code{difference_scale}), which covers the pairwise and the main-effect
#' difference families alike, since \code{bgmCompare()} gives them one scale.
#'
#' @details
#' \strong{The fits.} The model is refit at every multiplier in
#' \code{scale_multipliers} of the chosen scale \eqn{s_0}, with the original
#' fit's \code{iter} and \code{warmup} unless \code{iter} and \code{warmup}
#' are given. The original fit is the \code{1x} fit and is not refit. An
#' ordinal (omrf) refit starts each chain from the original fit's final
#' per-chain state, and a NUTS refit also from its adapted step size and mass
#' matrix; continuous (GGM) and mixed-MRF refits start cold. Every quantity
#' the check reports for a fit is that fit's own Rao-Blackwellized inclusion
#' probability, its natural log inclusion Bayes factor and verdict, and the
#' Monte Carlo standard error of the log Bayes factor from the fit's own
#' summary; the \code{1x} column is the original fit's reported analysis.
#'
#' \strong{Classification.} An edge's reference verdict is the original fit's.
#' An edge is \code{stable} when every refit's verdict equals it. Otherwise
#' the change in log Bayes factor from the original fit to each refit whose
#' verdict differs is compared with twice the combined Monte Carlo standard
#' error of the two fits, \eqn{2\sqrt{\mathrm{MCSE}_s^2 + \mathrm{MCSE}_1^2}}:
#' the edge is \code{beyond-monte-carlo-error} when the change exceeds it at
#' some such refit and \code{within-monte-carlo-error} when it does not at
#' any. \code{$edges$moved_at} is the multiplier nearest 1 at which the
#' verdict differs. \code{stability_lower} and \code{stability_upper} are the
#' smallest and largest multiplier in the contiguous run of fits around
#' \code{1x} that hold the reference verdict.
#'
#' \strong{Sufficiency at the chosen scale.} An edge whose per-chain verdicts
#' disagree in the original fit, or whose log Bayes factor lies within
#' \code{max(between-chain band, 2 * MCSE)} of a verdict threshold there, is
#' marked \code{insufficient}: the fit cannot certify the verdict the check
#' starts from, and the printed report lists the edge as not certifiable
#' rather than in a mover category.
#'
#' \strong{Convergence.} Each fit's largest split-\eqn{\hat R} over the
#' continuous parameters, smallest effective sample size over the continuous
#' parameters and over the indicators, and for NUTS the smallest E-BFMI and
#' the largest energy variance ratio, are in \code{$grid}. The classification
#' uses every fit; the printed report names a fit whose largest
#' split-\eqn{\hat R} exceeds 1.01.
#'
#' \strong{Models with a precision diagonal.} A continuous or mixed fit has a
#' Gamma or exponential prior on the precision diagonal with a standardized
#' rate \eqn{\eta}; at slab scale \eqn{s} the sampler runs that prior with raw
#' rate \eqn{\eta / s}. The refits hold \eqn{\eta} fixed and move the raw rate
#' with the slab scale, and the printed report gives the \eqn{\eta} held
#' fixed. A fit whose diagonal prior was given as a raw \code{rate} stops with
#' an error; see the section below.
#'
#' @section Fits with a raw-rate diagonal prior:
#' A diagonal prior given as \code{exponential_prior(rate = )} or
#' \code{gamma_prior(rate = )} fixes the raw rate independently of the slab
#' scale. The check has no sweep that holds a raw rate fixed. The table of
#' inclusion probabilities at a few slab scales with the raw rate held at the
#' fitted value is a loop over \code{\link{bgm}()}:
#' \preformatted{
#' scales = c(0.5, 1, 2)
#' fits = lapply(scales, function(s) {
#'   bgm(x,
#'       interaction_prior = normal_prior(scale = s),
#'       precision_scale_prior = exponential_prior(rate = 2),
#'       seed = 1)
#' })
#' pip = sapply(fits, function(f) {
#'   m = extract_posterior_inclusion_probabilities(f)
#'   m[lower.tri(m)]
#' })
#' colnames(pip) = paste0(scales, "x")
#' }
#' \code{\link{verdicts}()} on each fit gives the verdicts at every scale.
#'
#' @param bgms_object A fitted \code{bgms} object from \code{\link{bgm}()} run
#'   with \code{edge_selection = TRUE}, or a \code{bgmCompare} object from
#'   \code{\link{bgmCompare}()} run with \code{difference_selection = TRUE}.
#' @param scale_multipliers Numeric vector of positive multipliers of the
#'   chosen scale at which the model is refit. Default \code{c(0.5, 2)}. The
#'   multiplier \code{1} is the original fit, always included and never refit;
#'   each other multiplier is one refit.
#' @param evidence_threshold Positive numeric. Inclusion Bayes factor threshold
#'   for a presence verdict; \code{1 / evidence_threshold} is the absence
#'   threshold. Default: \code{10}.
#' @param refit_sampler One of \code{"same-as-fit"} (default; inherit the
#'   original fit's update method) or an explicit \code{"nuts"},
#'   \code{"adaptive-metropolis"}, or \code{"gibbs"}. Every comparison the
#'   check makes is between a refit and the original fit, so the refit sampler
#'   should be the original fit's unless the original fit is to be rerun with
#'   the other sampler too.
#' @param iter,warmup Integer sampling and warmup iterations per refit, or
#'   \code{NULL} (default) to use the original fit's.
#' @param include_preferred_scale Logical. Add a refit at the data-preferred
#'   scale \eqn{\hat s}. Default: \code{FALSE}.
#' @param cores Integer thread count for each refit's chains. Default: the
#'   original fit's core count.
#' @param seed Integer base seed for the refits. Default: \code{1}.
#' @param keep_fits Logical. Retain the full refit objects in the result (for
#'   power users); the default keeps only per-fit summaries. Default:
#'   \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, print each refit's raw sampler
#'   notes (energy, tree-depth, dropped chains) live as it runs. By default
#'   these are captured, not printed, and the text stays available in
#'   \code{$refit_diagnostics}. Default: \code{FALSE}.
#'
#' @return An object of class \code{"bgms_prior_sensitivity"}: a list with
#'   \describe{
#'   \item{\code{edges}}{One row per indicator: \code{edge}, the prior
#'     inclusion probability, the original fit's \code{chosen_scale_pip},
#'     \code{chosen_scale_log_bf}, \code{chosen_scale_mcse} (of the log Bayes
#'     factor) and \code{chosen_scale_verdict}; \code{stability_lower} and
#'     \code{stability_upper}; \code{mover} (\code{stable},
#'     \code{within-monte-carlo-error}, \code{beyond-monte-carlo-error}, or
#'     \code{NA} for an indicator without a verdict) and \code{moved_at};
#'     \code{insufficient} with its two subcauses \code{insufficient_noisy}
#'     (log Bayes factor within Monte Carlo error of a threshold) and
#'     \code{insufficient_disagree} (chains disagree on the verdict);
#'     \code{saturated} (the indicator never flipped in the original fit);
#'     then one \code{verdict_x<m>}, \code{log_bf_x<m>} and \code{mcse_x<m>}
#'     column per multiplier \code{m}, from that fit's own statistics.}
#'   \item{\code{grid}}{One row per fit: \code{multiplier}, \code{scale},
#'     \code{original_fit}, \code{iter}, \code{warmup}, \code{seed}, the
#'     diagnostics \code{rhat_max}, \code{ess_min}, \code{inclusion_ess_min},
#'     \code{ebfmi_min}, \code{energy_variance_ratio_max} and
#'     \code{warmup_incomplete}, \code{sampler_notes}, and \code{seconds}.}
#'   \item{\code{verdict_by_scale}}{Fit-by-indicator matrix of verdicts.}
#'   \item{\code{preferred_scale}}{The scale the estimated interactions imply:
#'     \code{s_hat} is the root-mean-square of the included interactions over
#'     the original fit's draws, with \code{lo} and \code{hi} an approximate
#'     95\% interval from a \eqn{1/\sqrt{2m}} log-scale standard error, \eqn{m}
#'     the mean number of included edges. \code{NA} where the swept prior is
#'     not the one on the pairwise interactions.}
#'   \item{\code{refit_diagnostics}}{The captured raw sampler notes per fit.}
#'   \item{\code{scale_multipliers}, \code{chosen_scale}, \code{unit},
#'     \code{diagonal_eta}, \code{evidence_threshold}, \code{refit_sampler},
#'     \code{warm}, \code{model_type}, \code{runtime_seconds},
#'     \code{fits}}{The multipliers (1 included) and the chosen scale, the
#'     unit the check reported on, the standardized diagonal rate held fixed
#'     (\code{NA} for a model without a precision diagonal), the settings, and
#'     the refit objects when \code{keep_fits = TRUE}.}
#'   }
#'
#' @seealso \code{\link{bgm}()}, \code{\link{verdicts}()},
#'   \code{\link{extract_posterior_inclusion_probabilities}()}
#' @family diagnostics
#' @export
#' @examples
#' \donttest{
#' fit = bgm(Wenchuan[, 1:6], chains = 2)
#' ps = prior_sensitivity_check(fit)
#' ps
#' }
prior_sensitivity_check = function(bgms_object,
                                   scale_multipliers = c(0.5, 2),
                                   evidence_threshold = 10,
                                   refit_sampler = "same-as-fit",
                                   iter = NULL,
                                   warmup = NULL,
                                   include_preferred_scale = FALSE,
                                   cores = NULL,
                                   seed = 1L,
                                   keep_fits = FALSE,
                                   verbose = FALSE) {
  UseMethod("prior_sensitivity_check")
}

#' @inheritParams prior_sensitivity_check
#' @exportS3Method
#' @noRd
prior_sensitivity_check.bgms = function(bgms_object,
                                        scale_multipliers = c(0.5, 2),
                                        evidence_threshold = 10,
                                        refit_sampler = "same-as-fit",
                                        iter = NULL,
                                        warmup = NULL,
                                        include_preferred_scale = FALSE,
                                        cores = NULL,
                                        seed = 1L,
                                        keep_fits = FALSE,
                                        verbose = FALSE) {
  spec = get_fit_spec(bgms_object)
  if(is.null(spec) || !isTRUE(spec$prior$edge_selection)) {
    stop(
      "prior_sensitivity_check() needs edge selection. Refit with ",
      "edge_selection = TRUE."
    )
  }
  prior_sensitivity_engine(bgms_object,
    unit = single_network_unit(spec),
    scale_multipliers = scale_multipliers,
    evidence_threshold = evidence_threshold,
    refit_sampler = refit_sampler, iter = iter, warmup = warmup,
    include_preferred_scale = include_preferred_scale, cores = cores,
    seed = seed, keep_fits = keep_fits, verbose = verbose
  )
}

#' @inheritParams prior_sensitivity_check
#' @exportS3Method
#' @noRd
prior_sensitivity_check.bgmCompare = function(bgms_object,
                                              scale_multipliers = c(0.5, 2),
                                              evidence_threshold = 10,
                                              refit_sampler = "same-as-fit",
                                              iter = NULL,
                                              warmup = NULL,
                                              include_preferred_scale = FALSE,
                                              cores = NULL,
                                              seed = 1L,
                                              keep_fits = FALSE,
                                              verbose = FALSE) {
  spec = get_fit_spec(bgms_object)
  if(is.null(spec) || !isTRUE(spec$prior$difference_selection)) {
    stop(
      "prior_sensitivity_check() needs difference selection. Refit with ",
      "bgmCompare(difference_selection = TRUE); without it every difference ",
      "is in the model and there is no inclusion Bayes factor to compare."
    )
  }
  if(is.na(difference_prior_inclusion(bgms_object))) {
    stop(
      "A stochastic-block difference prior has no single marginal inclusion ",
      "probability, so the check would report posterior odds rather than ",
      "Bayes factors. Refit with bernoulli_prior() or beta_bernoulli_prior() ",
      "to compare the difference verdicts across scales."
    )
  }
  prior_sensitivity_engine(bgms_object,
    unit = difference_unit(spec),
    scale_multipliers = scale_multipliers,
    evidence_threshold = evidence_threshold,
    refit_sampler = refit_sampler, iter = iter, warmup = warmup,
    include_preferred_scale = include_preferred_scale, cores = cores,
    seed = seed, keep_fits = keep_fits, verbose = verbose
  )
}


# ------------------------------------------------------------------
# single_network_unit / difference_unit
# ------------------------------------------------------------------
# What the check sweeps and what it reports on. bgm() compares edge indicators
# across the interaction slab scale; bgmCompare() compares difference
# indicators across the difference slab scale. Everything else -- the refits,
# the classification, the reporting -- is shared.
#
# Returns: list(
#   scale_field  Prior field the refits rescale.
#   chosen_scale The fit's own value of it.
#   noun         Singular name of one reported unit.
#   nouns        Its plural.
#   headline     The question the printed report opens with.
#   preferred    Whether a data-preferred scale can be estimated.
# )
# ------------------------------------------------------------------
single_network_unit = function(spec) {
  list(
    scale_field = "pairwise_scale",
    chosen_scale = spec$prior$pairwise_scale,
    noun = "edge", nouns = "edges",
    headline = "are the edge verdicts robust to the slab scale?",
    preferred = TRUE
  )
}

difference_unit = function(spec) {
  list(
    scale_field = "difference_scale",
    chosen_scale = spec$prior$difference_scale,
    noun = "difference", nouns = "differences",
    headline = "are the difference verdicts robust to the difference scale?",
    preferred = FALSE
  )
}


# ------------------------------------------------------------------
# prior_sensitivity_engine
# ------------------------------------------------------------------
# The shared machinery behind both methods: runs the refits, takes each fit's
# statistics and diagnostics, and hands the statistics to
# sensitivity_classification(). `unit` names the prior field to sweep and the
# vocabulary of the report; see single_network_unit and difference_unit.
# ------------------------------------------------------------------
prior_sensitivity_engine = function(bgms_object,
                                    unit,
                                    scale_multipliers,
                                    evidence_threshold,
                                    refit_sampler,
                                    iter,
                                    warmup,
                                    include_preferred_scale,
                                    cores,
                                    seed,
                                    keep_fits,
                                    verbose) {
  if(!is.numeric(evidence_threshold) || length(evidence_threshold) != 1L ||
    evidence_threshold <= 1) {
    stop("'evidence_threshold' must be a single number greater than 1.")
  }
  if(!is.numeric(scale_multipliers) || any(!is.finite(scale_multipliers)) ||
    any(scale_multipliers <= 0)) {
    stop("'scale_multipliers' must be positive multipliers of the chosen scale.")
  }

  spec = get_fit_spec(bgms_object)
  chosen_scale = unit$chosen_scale
  # The multipliers apply to this scale, so a fit whose swept prior has no
  # scale (a beta-prime slab carries none) has nothing to sweep.
  if(!is.numeric(chosen_scale) || length(chosen_scale) != 1L ||
    !is.finite(chosen_scale) || chosen_scale <= 0) {
    stop(
      "The fit carries no usable ", unit$scale_field, " (got ",
      format(chosen_scale), "), so there is no scale for the check to sweep. ",
      "This is what a scale-free slab such as beta_prime_prior() leaves ",
      "behind; refit with cauchy_prior() or normal_prior() to compare scales."
    )
  }
  lthr = log(evidence_threshold)
  if(is.null(cores)) cores = spec$sampler$chains

  # Standardized rate of the precision-diagonal prior, held fixed across the
  # sweep; NA for a model without a precision diagonal.
  eta = diagonal_eta(spec)

  rs = resolve_refit_sampler(bgms_object, refit_sampler)

  # An ordinal refit starts from the original fit's final state; the GGM and
  # mixed samplers take no initial state.
  is_omrf = identical(spec$model_type, "omrf")
  warm_state = if(is_omrf) extract_warm_state(bgms_object) else NULL

  # The scale the estimated interactions imply, read off the original fit.
  # Only defined where the swept prior is the one on the pairwise interactions;
  # the difference slab has no comparable plug-in estimate.
  preferred = if(isTRUE(unit$preferred)) {
    data_preferred_scale(bgms_object)
  } else {
    list(s_hat = NA_real_, lo = NA_real_, hi = NA_real_)
  }

  # --- The fits: 1x is the original fit, every other multiplier a refit -------
  multipliers = sort(unique(c(scale_multipliers, 1)))
  if(include_preferred_scale && is.finite(preferred$s_hat)) {
    multipliers = sort(unique(c(multipliers, preferred$s_hat / chosen_scale)))
  }
  n_fit = length(multipliers)
  s0_idx = which(multipliers == 1)
  refit_idx = setdiff(seq_len(n_fit), s0_idx)
  if(length(refit_idx) == 0L) {
    stop("'scale_multipliers' must include at least one multiplier other than 1.")
  }
  # Every refit has the original fit's run length unless the caller sets one.
  run_iter = as.integer(if(is.null(iter)) spec$sampler$iter else iter)
  run_warmup = as.integer(if(is.null(warmup)) spec$sampler$warmup else warmup)
  refit_cores = min(as.integer(cores), spec$sampler$chains)

  refits = vector("list", length(refit_idx))
  walls = numeric(length(refit_idx))
  # Per-refit sampler chatter (energy/tree-depth notes, dropped-chain warnings)
  # is captured, not printed: a user watching the check cannot act on per-refit
  # alarms mid-run, and $grid carries the diagnostics. The captured text stays
  # inspectable in $refit_diagnostics, and verbose = TRUE prints it live.
  refit_notes = vector("list", length(refit_idx))
  run_one_refit = function(k, show_progress = FALSE) {
    s = multipliers[refit_idx[k]] * chosen_scale
    refit_at_scale(
      bgms_object,
      scale = s,
      warm_state = warm_state, warmup = run_warmup, iter = run_iter,
      seed = seed + k, cores = refit_cores, sampler = rs$method,
      show_progress = show_progress,
      scale_field = unit$scale_field,
      diagonal_rate = if(is.na(eta)) NULL else eta / s
    )
  }
  # Each refit renders the sampler's own per-chain display, with the whole
  # refit set announced up front. A refit's own cat() diagnostics (NUTS
  # settle-in notes and the like) are captured, not printed, unless the user
  # asked to see them; the progress bar is not gated by that option.
  if(!verbose) {
    old_bgms_verbose = options(bgms.verbose = FALSE)
    on.exit(options(old_bgms_verbose), add = TRUE)
  }
  show_bar = interactive() && !verbose
  if(show_bar) {
    refit_scales = sprintf("%.2gx", multipliers[refit_idx])
    n_chains = length(get_raw_samples(bgms_object)$pairwise)
    message(sprintf(
      paste0(
        "Refitting at %d scale%s (%s the chosen scale).\n",
        "Each refit runs %d chain%s for %d warmup and %d sampling iterations."
      ),
      length(refit_scales), if(length(refit_scales) == 1L) "" else "s",
      paste(refit_scales, collapse = ", "),
      n_chains, if(n_chains == 1L) "" else "s", run_warmup, run_iter
    ))
  }
  for(k in seq_along(refit_idx)) {
    t0 = Sys.time()
    warns_k = character(0)
    fit_k = NULL
    absorb_warning = function(w) {
      warns_k <<- c(warns_k, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
    if(show_bar) {
      # The native display writes to the console, so it must not be captured;
      # announce which refit is starting, then let its per-chain bars animate.
      message(sprintf("\n[%d/%d] %.2gx scale", k, length(refit_idx), multipliers[refit_idx[k]]))
      fit_k = withCallingHandlers(
        run_one_refit(k, show_progress = TRUE),
        warning = absorb_warning
      )
      notes = warns_k
    } else if(verbose) {
      cat_k = utils::capture.output({
        fit_k = withCallingHandlers(run_one_refit(k), warning = absorb_warning)
      })
      notes = c(cat_k[nzchar(cat_k)], warns_k)
    } else {
      fit_k = withCallingHandlers(run_one_refit(k), warning = absorb_warning)
      notes = warns_k
    }
    refits[[k]] = fit_k
    refit_notes[[k]] = notes
    walls[k] = as.numeric(Sys.time() - t0, units = "secs")
    if(verbose && length(notes) > 0) {
      cat(sprintf("[%.2gx refit] sampler notes:\n", multipliers[refit_idx[k]]))
      cat(paste0("  ", notes), sep = "\n")
    }
  }

  # --- Per-fit statistics and diagnostics -------------------------------------
  fits = vector("list", n_fit)
  fits[[s0_idx]] = bgms_object
  for(k in seq_along(refit_idx)) fits[[refit_idx[k]]] = refits[[k]]

  fit_stats = lapply(fits, refit_edge_stats, evidence_threshold = evidence_threshold)
  diag_rows = lapply(fits, refit_diagnostics_row)

  cls = sensitivity_classification(fit_stats, multipliers, s0_idx, lthr)

  # --- Fit table: one row per fit ---------------------------------------------
  secs = numeric(n_fit)
  secs[refit_idx] = walls
  seeds = rep(NA_integer_, n_fit)
  seeds[refit_idx] = seed + seq_along(refit_idx)
  seeds[s0_idx] = spec$sampler$seed
  grid_notes = vector("list", n_fit)
  grid_notes[[s0_idx]] = character(0)
  for(k in seq_along(refit_idx)) grid_notes[[refit_idx[k]]] = refit_notes[[k]]
  names(grid_notes) = sprintf("%gx", multipliers)
  diag_col = function(name, type) vapply(diag_rows, `[[`, type, name)
  original = seq_len(n_fit) == s0_idx
  grid = data.frame(
    multiplier = multipliers,
    scale = multipliers * chosen_scale,
    original_fit = original,
    iter = ifelse(original, spec$sampler$iter, run_iter),
    warmup = ifelse(original, spec$sampler$warmup, run_warmup),
    seed = seeds,
    rhat_max = diag_col("rhat_max", numeric(1)),
    ess_min = diag_col("ess_min", numeric(1)),
    inclusion_ess_min = diag_col("inclusion_ess_min", numeric(1)),
    ebfmi_min = diag_col("ebfmi_min", numeric(1)),
    energy_variance_ratio_max = diag_col("energy_variance_ratio_max", numeric(1)),
    warmup_incomplete = diag_col("warmup_incomplete", logical(1)),
    sampler_notes = lengths(grid_notes) > 0,
    seconds = secs,
    row.names = NULL
  )

  structure(
    list(
      edges = cls$edges,
      grid = grid,
      scale_multipliers = multipliers,
      chosen_scale = chosen_scale,
      verdict_by_scale = cls$verdict_by_scale,
      preferred_scale = preferred,
      unit = unit,
      diagonal_eta = eta,
      evidence_threshold = evidence_threshold,
      refit_sampler = rs$method,
      warm = !is.null(warm_state),
      model_type = spec$model_type,
      runtime_seconds = sum(walls),
      refit_diagnostics = grid_notes,
      fits = if(keep_fits) refits else NULL
    ),
    class = "bgms_prior_sensitivity"
  )
}


# ------------------------------------------------------------------
# sensitivity_classification
# ------------------------------------------------------------------
# The classification, from the per-fit statistics alone.
#
# @param fit_stats    List over fits of refit_edge_stats() output; the element
#                     at s0_idx is the original fit's.
# @param multipliers  The fits' multipliers of the chosen scale, sorted, 1
#                     included.
# @param s0_idx       Position of 1 in `multipliers`.
# @param lthr         log(evidence_threshold).
#
# Returns: list(edges = the per-indicator table, verdict_by_scale = fit x
#   indicator verdict matrix).
# ------------------------------------------------------------------
sensitivity_classification = function(fit_stats, multipliers, s0_idx, lthr) {
  s0 = fit_stats[[s0_idx]]
  edge_names = s0$edge
  n_edges = length(edge_names)
  n_fit = length(multipliers)
  prior_odds = s0$prior_odds

  stat_mat = function(field, type) {
    m = vapply(fit_stats, function(st) as.vector(st[[field]]), vector(type, n_edges))
    # vapply drops to a vector when there is one indicator; keep the
    # fit-by-indicator shape.
    matrix(m, nrow = n_edges, ncol = n_fit)
  }
  lbf_mat = t(stat_mat("lbf", "numeric")) # n_fit x n_edges
  mcse_mat = t(stat_mat("mcse_lbf", "numeric"))
  verdict_mat = t(stat_mat("verdict", "character"))
  # A saturated indicator has no finite MCSE; its error enters as zero, and
  # its log BF is infinite, so any change from it reads as beyond error.
  mcse0 = mcse_mat
  mcse0[!is.finite(mcse0)] = 0

  # --- Sufficiency at the chosen scale, from the original fit -----------------
  band = s0$band_half
  band[is.na(band)] = 0
  m0 = s0$mcse_lbf
  m0[is.na(m0)] = 0
  hw = pmax(band, 2 * m0)
  straddle = (!s0$zeroflip) &
    (((s0$lbf - hw < lthr) & (s0$lbf + hw > lthr)) |
      ((s0$lbf - hw < -lthr) & (s0$lbf + hw > -lthr)))
  disagree = (!s0$zeroflip) & !s0$unanimous
  insufficient = (!s0$zeroflip) & (disagree | straddle)
  insufficient[is.na(insufficient)] = FALSE
  insufficient_noisy = straddle & insufficient
  insufficient_noisy[is.na(insufficient_noisy)] = FALSE
  insufficient_disagree = disagree & insufficient
  insufficient_disagree[is.na(insufficient_disagree)] = FALSE

  # --- Mover category, where the verdict first differs, stability bracket -----
  ref_v = s0$verdict
  mover = rep(NA_character_, n_edges)
  moved_at = rep(NA_real_, n_edges)
  stability_lower = rep(NA_real_, n_edges)
  stability_upper = rep(NA_real_, n_edges)
  # Fits in order of distance from the chosen scale, for moved_at.
  by_distance = order(abs(log(multipliers)))
  for(e in seq_len(n_edges)) {
    if(is.na(ref_v[e])) next
    v = verdict_mat[, e]
    differs = !is.na(v) & v != ref_v[e]
    if(!any(differs)) {
      mover[e] = "stable"
    } else {
      change = abs(lbf_mat[, e] - lbf_mat[s0_idx, e])
      error = 2 * sqrt(mcse0[, e]^2 + mcse0[s0_idx, e]^2)
      beyond = differs & (is.na(change) | change > error)
      mover[e] = if(any(beyond)) "beyond-monte-carlo-error" else "within-monte-carlo-error"
      moved_at[e] = multipliers[by_distance[which(differs[by_distance])[1]]]
    }
    lo = s0_idx
    while(lo > 1L && !is.na(v[lo - 1L]) && v[lo - 1L] == ref_v[e]) lo = lo - 1L
    hi = s0_idx
    while(hi < n_fit && !is.na(v[hi + 1L]) && v[hi + 1L] == ref_v[e]) hi = hi + 1L
    stability_lower[e] = multipliers[lo]
    stability_upper[e] = multipliers[hi]
  }

  edges = data.frame(
    edge = edge_names,
    prior_inclusion_probability = prior_odds / (1 + prior_odds),
    chosen_scale_pip = s0$pip,
    chosen_scale_log_bf = s0$lbf,
    chosen_scale_mcse = s0$mcse_lbf,
    chosen_scale_verdict = ref_v,
    stability_lower = stability_lower,
    stability_upper = stability_upper,
    mover = mover,
    moved_at = moved_at,
    insufficient = insufficient,
    insufficient_noisy = insufficient_noisy,
    insufficient_disagree = insufficient_disagree,
    saturated = s0$zeroflip,
    stringsAsFactors = FALSE, row.names = NULL
  )
  vcols = as.data.frame(t(verdict_mat), stringsAsFactors = FALSE)
  names(vcols) = paste0("verdict_x", multipliers)
  lcols = as.data.frame(t(lbf_mat))
  names(lcols) = paste0("log_bf_x", multipliers)
  mcols = as.data.frame(t(mcse_mat))
  names(mcols) = paste0("mcse_x", multipliers)
  edges = cbind(edges, vcols, lcols, mcols)

  list(edges = edges, verdict_by_scale = verdict_mat)
}


# ------------------------------------------------------------------
# order_upper_tri_rowmajor
# ------------------------------------------------------------------
# Permutation that reorders the column-major upper-triangle extraction
# matrix[upper.tri(matrix)] into the row-major edge order used by the
# raw pairwise / indicator columns and edge names.
#
# @param p  Number of variables.
#
# Returns: integer permutation of length p(p-1)/2.
# ------------------------------------------------------------------
order_upper_tri_rowmajor = function(p) {
  # Column-major upper triangle visits (i, j) with j outer, i < j inner.
  cm = which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  # Row-major target order sorts by row then column.
  order(cm[, "row"], cm[, "col"])
}


# ==============================================================================
# Methods for bgms_prior_sensitivity
# ==============================================================================

#' @title Print a Prior Sensitivity Check
#'
#' @description
#' Prints the verdict counts at each scale, the number of edges whose verdict
#' held at every scale, changed within Monte Carlo error, changed beyond it,
#' or could not be certified at the chosen scale, the edges that changed
#' beyond Monte Carlo error by name, the fits whose largest split R-hat
#' exceeds 1.01, and how the chosen scale compares with the scale the
#' estimated interactions imply.
#'
#' @param x A \code{bgms_prior_sensitivity} object.
#' @param max_rows Integer. Maximum edges to name in the table of edges that
#'   changed beyond Monte Carlo error; the rest are counted and left to
#'   \code{$edges}. Default: \code{10}.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{prior_sensitivity_check}()}
#' @family diagnostics
#' @export
print.bgms_prior_sensitivity = function(x, max_rows = 10L, ...) {
  edges = x$edges
  n_edges = nrow(edges)
  multipliers = x$scale_multipliers
  mlab = sprintf("%.2gx", multipliers)

  w = min(getOption("width", 80L), 80L)
  wrap = function(...) writeLines(strwrap(paste0(...), width = w))
  list_of = function(v) {
    if(length(v) <= 1L) {
      return(paste(v))
    }
    paste0(paste(v[-length(v)], collapse = ", "), " and ", v[length(v)])
  }

  unit = x$unit %||% list(
    noun = "edge", nouns = "edges",
    headline = "are the edge verdicts robust to the slab scale?"
  )

  cat(sprintf("Prior sensitivity check: %s\n", unit$headline))
  wrap(sprintf(
    "Verdicts at an inclusion Bayes factor of %g from fits at %s the chosen scale; the 1x fit is the original fit. %d %s.",
    x$evidence_threshold, list_of(mlab), n_edges, unit$nouns
  ))
  eta = x$diagonal_eta %||% NA_real_
  if(is.finite(eta)) {
    wrap(sprintf(
      "The refits move the precision diagonal's prior with the slab scale, holding its standardized rate eta at %.3g.",
      eta
    ))
  }
  cat("\n")

  # An indicator the sampler never updated has no verdict to be robust or
  # sensitive; bgmCompare() leaves the main-effect difference indicators there
  # unless main_difference_selection = TRUE.
  no_verdict = is.na(edges$chosen_scale_verdict)
  if(any(no_verdict)) {
    wrap(sprintf(
      "%d of the %d %s were never updated by the sampler and carry no verdict at any scale; they are excluded from the counts below.",
      sum(no_verdict), n_edges, unit$nouns
    ))
    cat("\n")
    edges = edges[!no_verdict, , drop = FALSE]
    x$verdict_by_scale = x$verdict_by_scale[, !no_verdict, drop = FALSE]
    n_edges = nrow(edges)
  }

  # Verdict counts at each scale.
  cat("Verdict counts by scale:\n")
  vt = apply(x$verdict_by_scale, 1, function(v) {
    c(
      presence = sum(v == "presence", na.rm = TRUE),
      undecided = sum(v == "undecided", na.rm = TRUE),
      absence = sum(v == "absence", na.rm = TRUE)
    )
  })
  vt = matrix(vt, nrow = 3L, dimnames = list(c("presence", "undecided", "absence"), mlab))
  out = utils::capture.output(print(vt))
  cat(paste0("  ", out, collapse = "\n"), "\n", sep = "")
  cat("\n")

  # One category per edge: not certifiable at the chosen scale first, then the
  # mover category.
  uncert = edges$insufficient %in% TRUE
  held = !uncert & edges$mover %in% "stable"
  within = !uncert & edges$mover %in% "within-monte-carlo-error"
  beyond = !uncert & edges$mover %in% "beyond-monte-carlo-error"
  counts = c(sum(held), sum(within), sum(beyond), sum(uncert))
  labels = c(
    "verdict held at every scale",
    "changed, within Monte Carlo error",
    "changed, beyond Monte Carlo error",
    "not certifiable at the chosen scale"
  )
  cat(sprintf("  %-38s %4d\n", labels, counts), sep = "")
  cat("\n")

  # The edges whose verdict changed beyond Monte Carlo error, with their
  # verdict at each scale and the scale nearest 1x at which it differs.
  if(sum(beyond) > 0) {
    cat(sprintf(
      "%d %s changed verdict beyond Monte Carlo error:\n",
      sum(beyond), if(sum(beyond) == 1L) unit$noun else unit$nouns
    ))
    idx = which(beyond)
    show = utils::head(idx, max_rows)
    vmat = t(x$verdict_by_scale[, show, drop = FALSE])
    tab = cbind(edges$edge[show], vmat, sprintf("%.2gx", edges$moved_at[show]))
    colnames(tab) = c(unit$noun, mlab, "differs from")
    rownames(tab) = rep("", nrow(tab))
    print(tab, quote = FALSE, print.gap = 2)
    if(length(idx) > length(show)) {
      cat(sprintf("  and %d more; see $edges.\n", length(idx) - length(show)))
    }
    cat("\n")
  }

  # The edges the original fit cannot certify, and why.
  if(sum(uncert) > 0) {
    nm = edges$edge[uncert]
    shown = utils::head(nm, max_rows)
    tail_txt = if(length(nm) > length(shown)) {
      sprintf(", and %d more (see $edges)", length(nm) - length(shown))
    } else {
      ""
    }
    one = length(nm) == 1L
    noisy = any(edges$insufficient_noisy[uncert])
    disagree = any(edges$insufficient_disagree[uncert])
    cause = if(noisy && disagree) {
      "the Bayes factor lies within Monte Carlo error of a threshold or the chains disagree on the verdict"
    } else if(disagree) {
      "the chains disagree on the verdict"
    } else {
      "the Bayes factor lies within Monte Carlo error of a threshold"
    }
    wrap(sprintf(
      "%d %s %s not certifiable at the chosen scale: %s%s. In the original fit %s, so a rerun with another seed could change %s verdict without any change of prior.",
      length(nm), if(one) unit$noun else unit$nouns, if(one) "is" else "are",
      paste(shown, collapse = ", "), tail_txt, cause,
      if(one) "its" else "the"
    ))
    cat("\n")
  }

  # Fits whose largest split R-hat exceeds 1.01.
  g = x$grid
  if(!is.null(g$rhat_max)) {
    high = which(is.finite(g$rhat_max) & g$rhat_max > 1.01)
    if(length(high) > 0) {
      lab = ifelse(g$original_fit[high], "1x original", sprintf("%.2gx", g$multiplier[high]))
      wrap(sprintf(
        "Largest split R-hat above 1.01: %s.",
        list_of(sprintf("%s fit %s", lab, format(round(g$rhat_max[high], 3))))
      ))
      cat("\n")
    }
  }

  # How the chosen scale compares with the scale the estimated interactions
  # imply.
  ps = x$preferred_scale
  if(!is.null(ps) && is.finite(ps$s_hat)) {
    ratio = x$chosen_scale / ps$s_hat
    direction = if(ratio > 2) {
      " A slab wider than the weights it describes moves near-threshold verdicts toward absence."
    } else if(ratio < 0.5) {
      " A slab narrower than the weights it describes moves near-threshold verdicts toward presence."
    } else {
      ""
    }
    wrap(sprintf(
      "Note: the slab scale is %.3g, and the estimated edge weights imply a scale of %.3g (approximate 95%% interval %.3g to %.3g); the slab is %.2g times that.%s",
      x$chosen_scale, ps$s_hat, ps$lo, ps$hi, ratio, direction
    ))
    cat("\n")
  }

  secs = if(x$runtime_seconds < 10) {
    sprintf("%.1f s", x$runtime_seconds)
  } else {
    sprintf("%.0f s", x$runtime_seconds)
  }
  n_refit = sum(!g$original_fit)
  start_txt = if(isTRUE(x$warm)) "started from the original fit's final state" else "cold starts"
  it = g$iter[!g$original_fit][1]
  wu = g$warmup[!g$original_fit][1]
  cat(sprintf(
    "Refits:  %d %s refit%s (%d warmup, %d sampling iterations each), %s, %s total.\n",
    n_refit, x$refit_sampler, if(n_refit == 1L) "" else "s", wu, it, start_txt, secs
  ))
  cat("See ?prior_sensitivity_check for the construction.\n")
  invisible(x)
}
