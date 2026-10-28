# ==============================================================================
# Anchored sensitivity curve: importance reweighting between fixed-scale fits.
# ==============================================================================
#
# A fit at fixed slab scale s_a is reweighted to a nearby scale s with per-draw
# slab-density ratios over the included edges (and the precision-diagonal prior
# ratio for GGM and mixed fits). The likelihood cancels in the ratio; the
# posterior normalizing constants cancel through self-normalization. The
# self-normalized average of the indicator draws is the inclusion probability
# at s; its Monte Carlo error comes from batch means within chains; the
# importance ESS and the maximum normalized weight describe the weight
# concentration. The curve pools the anchors that reach a grid point by inverse
# variance and takes each anchor fit's own Rao-Blackwellized value at the
# anchor's grid position. Design and validation: dev/plans/active/
# 2026-10-22_sensitivity-engine_DESIGN.md; dev/audit/checks/2026-10-reweighting-reach/.
# ==============================================================================


# ------------------------------------------------------------------
# anchor_draws
# ------------------------------------------------------------------
# Per-chain draws of a fit in the slab frame, for reweighting.
#
# @param fit  A fitted bgms object.
#
# Returns: list(theta = per-chain iter x parameters matrices in the frame the
#   slab prior applies to (GGM pairwise values carry the -0.5 factor),
#   gamma = per-chain iter x parameters indicator matrices aligned with theta,
#   indicator = per-chain iter x indicators matrices, the unit the inclusion
#   probability is reported on (identical to gamma for bgm()),
#   family = "normal" or "cauchy").
# ------------------------------------------------------------------
anchor_draws = function(fit) {
  spec = get_fit_spec(fit)
  raw = get_raw_samples(fit)
  if(identical(spec$model_type, "compare")) {
    return(compare_anchor_draws(fit, spec, raw))
  }
  # GGM pairwise draws are K_ij; the slab prior is evaluated at -K_ij / 2
  # (ggm_model.cpp:313). The mixed sampler stores every pairwise block on the
  # slab's own scale: the discrete and cross blocks as is
  # (mixed_mrf_gradient.cpp:503, 521), the continuous block as -Omega / 2
  # (mixed_mrf_model.cpp:534, mixed_mrf_gradient.cpp:574).
  slab_factor = if(identical(spec$model_type, "ggm")) -0.5 else 1
  theta = lapply(raw$pairwise, function(m) slab_factor * m)
  gamma = raw$indicator
  family = spec$prior$interaction_prior_type %||% "normal"
  list(
    theta = theta, gamma = gamma, indicator = gamma,
    family = tolower(family),
    diagonal = anchor_diagonal(spec, raw)
  )
}


# ------------------------------------------------------------------
# anchor_diagonal
# ------------------------------------------------------------------
# The Gamma prior on the precision diagonal, in the frame the weight ratio
# needs it, or NULL for a model that has no such diagonal.
#
# The refits set the diagonal rate to eta / s, and the weight ratio carries
# the diagonal prior as well as the slab. The sampler evaluates that prior at
# K_jj / 2
# (ggm_model.cpp:320, mixed_mrf_gradient.cpp:559). The two model types store
# the diagonal on different scales:
#
#   ggm        raw main draws are K_jj (build_precision_from_draw());
#              K_jj / 2 = 0.5 * draw.
#   mixed_mrf  the "(precision diag)" columns are the diagonal of
#              pairwise_effects_continuous_, stored as -Omega / 2 for the whole
#              continuous block (mixed_mrf_model.cpp:534); K_jj / 2 = -draw.
#
# @param spec  The fit's spec.
# @param raw   Its raw samples.
#
# Returns: list(sum = per-chain vectors of sum_j K_jj / 2, n = number of
#   diagonal elements, shape = the Gamma shape alpha), or NULL.
# ------------------------------------------------------------------
anchor_diagonal = function(spec, raw) {
  is_ggm = identical(spec$model_type, "ggm")
  cols = if(is_ggm) {
    seq_len(ncol(raw$main[[1]]))
  } else if(identical(spec$model_type, "mixed_mrf")) {
    which(endsWith(raw$parameter_names[["main"]], "(precision diag)"))
  } else {
    integer(0)
  }
  if(length(cols) == 0L) {
    return(NULL)
  }
  to_half_precision = if(is_ggm) 0.5 else -1
  list(
    sum = lapply(raw$main, function(m) {
      to_half_precision * rowSums(m[, cols, drop = FALSE])
    }),
    n = length(cols),
    shape = spec$prior$scale_shape %||% 1
  )
}


# ------------------------------------------------------------------
# compare_anchor_draws
# ------------------------------------------------------------------
# Anchor draws for a bgmCompare() fit, where the swept prior is the difference
# slab and one indicator can gate several parameters.
#
# A pairwise difference indicator gates that pair's difference in every group
# contrast; a main-effect difference indicator gates the whole block of that
# variable's threshold differences, again in every contrast. The prior is
# evaluated once per gated parameter, so the importance weight sums over
# parameters while the inclusion probability is reported per indicator.
#
# @param fit   A fitted bgmCompare object.
# @param spec  Its spec.
# @param raw   Its raw samples.
#
# Returns: as anchor_draws().
# ------------------------------------------------------------------
compare_anchor_draws = function(fit, spec, raw) {
  arguments = extract_arguments(fit)
  names_all = raw$parameter_names
  num_variables = as.integer(
    arguments[["num_variables"]] %||% arguments[["no_variables"]]
  )

  num_pairs = length(names_all$pairwise_baseline)
  num_contrasts = length(names_all$pairwise_diff) / num_pairs
  num_main_baseline = length(names_all$main_baseline)

  # Main parameters per variable: one per threshold for an ordinal variable,
  # a linear and a quadratic term for a Blume-Capel one. The field is named
  # per model type -- is_ordinal_variable on a compare fit, is_ordinal on the
  # single-network ones -- and is read exactly, not by partial match.
  is_ordinal = arguments[["is_ordinal_variable"]] %||% arguments[["is_ordinal"]]
  block = ifelse(is_ordinal, as.integer(arguments[["num_categories"]]), 2L)

  # Indicators run over the upper triangle with the diagonal: (v, v) is that
  # variable's main-effect difference, (i, j) the pair's.
  idx = compare_indicator_index(num_variables)
  main_indicator = which(idx[, 1] == idx[, 2])
  pair_indicator = which(idx[, 1] != idx[, 2])

  # Which indicator gates each difference parameter, in the column order the
  # draws carry: pairwise differences contrast by contrast, then main ones.
  pairwise_owner = rep(pair_indicator, times = num_contrasts)
  main_owner = rep(rep(main_indicator, times = block), times = num_contrasts)
  owner = c(pairwise_owner, main_owner)

  pairwise_diff_cols = num_pairs + seq_len(num_pairs * num_contrasts)
  main_diff_cols = num_main_baseline + seq_len(num_main_baseline * num_contrasts)

  theta = vector("list", length(raw$pairwise))
  gamma = vector("list", length(raw$pairwise))
  for(c in seq_along(raw$pairwise)) {
    theta[[c]] = cbind(
      raw$pairwise[[c]][, pairwise_diff_cols, drop = FALSE],
      raw$main[[c]][, main_diff_cols, drop = FALSE]
    )
    gamma[[c]] = raw$indicator[[c]][, owner, drop = FALSE]
  }

  family = spec$prior$difference_prior_type %||% "cauchy"
  list(
    theta = theta, gamma = gamma, indicator = raw$indicator,
    family = tolower(family)
  )
}


# ------------------------------------------------------------------
# anchor_log_weights
# ------------------------------------------------------------------
# Per-draw log importance weights from anchor scale s_a to each target
# scale, for one chain. Only included edges contribute; the weight is
# the slab-density ratio summed over them.
#
# The refits also move the Gamma prior on the precision diagonal, rate =
# eta / s; the target differs from the anchor there too and the weight carries
# that ratio as well. For Gamma(alpha,
# eta / s) on each of the n diagonal values x_j = K_jj / 2, the per-draw term is
#
#   n * alpha * log((eta / s) / (eta / s_a)) - (eta / s - eta / s_a) * S_t
#     = -n * alpha * log(s / s_a) - eta * (1 / s - 1 / s_a) * S_t,
#
# with S_t = sum_j K_jj / 2 at draw t. The leading term is the same for every
# draw and cancels in the self-normalized average; the S_t term reweights.
#
# @param theta    iter x edges matrix, slab frame.
# @param gamma    iter x edges indicator matrix.
# @param family   "normal" or "cauchy".
# @param s_a      Anchor slab scale.
# @param s_grid   Target slab scales.
# @param diagonal The anchor's anchor_diagonal() record, or NULL; `sum` is
#                 this chain's per-draw S_t.
# @param eta      Standardized diagonal rate held fixed by the sweep, or NULL
#                 when the sweep leaves the diagonal prior alone.
#
# Returns: iter x length(s_grid) matrix of log weights.
# ------------------------------------------------------------------
anchor_log_weights = function(theta, gamma, family, s_a, s_grid,
                              diagonal = NULL, eta = NULL) {
  # The reweighting identity is written per slab family; a scale-free family
  # (beta-prime) has no s to sweep and any other one has a different density
  # ratio.
  if(!family %in% c("normal", "cauchy")) {
    stop(
      "The anchored sensitivity curve supports a normal or Cauchy slab on ",
      "the swept parameters; this fit uses a '", family, "' slab, for which ",
      "there is no slab scale to reweight across."
    )
  }
  incl = gamma == 1
  m_t = rowSums(incl)
  if(identical(family, "normal")) {
    a_t = rowSums((theta * incl)^2)
    lw = vapply(s_grid, function(s) {
      -m_t * log(s / s_a) - 0.5 * a_t * (1 / s^2 - 1 / s_a^2)
    }, numeric(length(m_t)))
  } else {
    base = log1p((theta / s_a)^2)
    lw = vapply(s_grid, function(s) {
      -m_t * log(s / s_a) -
        rowSums(incl * (log1p((theta / s)^2) - base))
    }, numeric(length(m_t)))
  }
  if(is.null(dim(lw))) lw = matrix(lw, nrow = 1L)
  # Only when the sweep moved the diagonal prior.
  if(!is.null(eta) && !is.null(diagonal)) {
    s_t = diagonal$sum
    dlw = vapply(s_grid, function(s) {
      -diagonal$n * diagonal$shape * log(s / s_a) -
        eta * (1 / s - 1 / s_a) * s_t
    }, numeric(length(s_t)))
    if(is.null(dim(dlw))) dlw = matrix(dlw, nrow = 1L)
    lw = lw + dlw
  }
  lw
}


# ------------------------------------------------------------------
# anchor_reweight
# ------------------------------------------------------------------
# Reweighted inclusion probabilities at each target scale from one anchor
# fit's draws, with their Monte Carlo error and weight diagnostics; pooled
# over chains, and per chain for the unanimity check.
#
# The Monte Carlo error is estimated by batch means on the estimator's
# influence process. With p_hat the self-normalized estimate and w_bar the mean
# weight, each block j of b consecutive draws within a chain gives
#
#   H_j = (1 / b) * sum_{i in j} w_i (g_i - p_hat) / w_bar,
#
# and var(p_hat) = var(H_j) / n_blocks over the blocks of all chains. The
# blocks use the global p_hat and w_bar, since a block with a larger total
# weight contributes more to the full estimate. The variance is bounded below
# by the smoothed-count variance p~ (1 - p~) / ESS with
# p~ = (ESS p_hat + 1/2) / (ESS + 1), which is a regularization for edges whose
# indicator never flips (every H_j is zero there), not an estimate.
#
# @param draws    Output of anchor_draws() for the anchor fit.
# @param s_a      Anchor slab scale.
# @param s_grid   Target slab scales (length P).
# @param eta      Standardized diagonal rate the sweep holds fixed, or NULL
#                 for a model without a precision diagonal.
# @param block    Block length in draws for the batch means; shortened so
#                 that every chain has at least two blocks.
#
# Returns: list(pip = P x E pooled estimate, var = P x E Monte Carlo
#   variance, ess = length-P importance ESS, max_weight = length-P largest
#   normalized weight, n = number of draws over all chains,
#   chain_pip = list of P x E per chain, chain_ess = C x P).
# ------------------------------------------------------------------
anchor_reweight = function(draws, s_a, s_grid, eta = NULL, block = 300L) {
  n_chain = length(draws$theta)
  diag = draws$diagonal
  lw = lapply(seq_len(n_chain), function(c) {
    anchor_log_weights(
      draws$theta[[c]], draws$gamma[[c]], draws$family, s_a, s_grid,
      diagonal = if(is.null(diag)) NULL else list(sum = diag$sum[[c]], n = diag$n, shape = diag$shape),
      eta = eta
    )
  })
  # The weight sums over gated parameters; the inclusion probability is
  # reported on the indicators, which are the same thing for bgm() and a
  # coarser unit for bgmCompare().
  indicator = draws$indicator %||% draws$gamma
  n_grid = length(s_grid)
  n_edge = ncol(indicator[[1]])
  n_draw = nrow(indicator[[1]])
  n_total = n_draw * n_chain
  b = max(1L, min(as.integer(block), n_draw %/% 2L))
  n_block = n_draw %/% b
  rows = seq_len(n_block * b)
  blk = rep(seq_len(n_block), each = b)

  chain_pip = vector("list", n_chain)
  chain_ess = matrix(NA_real_, n_chain, n_grid)
  pip = matrix(NA_real_, n_grid, n_edge)
  var = matrix(NA_real_, n_grid, n_edge)
  ess = numeric(n_grid)
  max_weight = numeric(n_grid)
  gamma_all = do.call(rbind, indicator)
  for(p in seq_len(n_grid)) {
    lw_all = unlist(lapply(lw, function(m) m[, p]))
    lmax = max(lw_all)
    w_all = exp(lw_all - lmax)
    w_sum = sum(w_all)
    ess[p] = w_sum^2 / sum(w_all^2)
    max_weight[p] = max(w_all) / w_sum
    p_hat = colSums(w_all * gamma_all) / w_sum
    pip[p, ] = p_hat
    w_bar = w_sum / n_total
    H = NULL
    for(c in seq_len(n_chain)) {
      w_c = exp(lw[[c]][rows, p] - lmax)
      g_c = indicator[[c]][rows, , drop = FALSE]
      infl = w_c * (g_c - rep(p_hat, each = length(rows))) / w_bar
      H = rbind(H, rowsum(infl, blk, reorder = FALSE) / b)
    }
    v_bm = if(nrow(H) >= 2L) apply(H, 2, stats::var) / nrow(H) else rep(NA_real_, n_edge)
    p_s = (ess[p] * p_hat + 0.5) / (ess[p] + 1)
    v_floor = p_s * (1 - p_s) / ess[p]
    var[p, ] = pmax(v_bm, v_floor, na.rm = TRUE)
  }
  for(c in seq_len(n_chain)) {
    cp = matrix(NA_real_, n_grid, n_edge)
    for(p in seq_len(n_grid)) {
      w = exp(lw[[c]][, p] - max(lw[[c]][, p]))
      chain_ess[c, p] = sum(w)^2 / sum(w^2)
      cp[p, ] = colSums(w * indicator[[c]]) / sum(w)
    }
    chain_pip[[c]] = cp
  }
  list(
    pip = pip, var = var, ess = ess, max_weight = max_weight, n = n_total,
    chain_pip = chain_pip, chain_ess = chain_ess
  )
}


# ------------------------------------------------------------------
# assemble_curve
# ------------------------------------------------------------------
# The display curve from the anchors' reweighted estimates.
#
# An anchor contributes to a grid point when it passed the convergence gate
# and its importance ESS there is at least max(ess_floor, ess_fraction * n),
# n its number of draws. Contributing anchors are pooled by inverse Monte
# Carlo variance (anchor_reweight()$var), per edge; the pooled variance is the
# reciprocal of the summed weights. A point no anchor reaches is NA.
#
# A point is disputed for an edge when the contributing anchors differ both
# by Cochran's Q above the q_level quantile of chi-squared on k - 1 degrees of
# freedom and by more than `tolerance` in log Bayes factor between the two
# most different anchors. A disputed point takes the nearest contributing
# anchor's estimate (nearest in log scale) and is recorded in `disagree`.
#
# At an anchor's own grid position the curve takes that fit's
# Rao-Blackwellized estimate and its reported Monte Carlo error
# (`anchor_exact`), nothing pooled.
#
# @param reweights     List over anchors of anchor_reweight() output.
# @param usable        Logical per anchor (convergence gate).
# @param ess_floor     Minimum importance ESS in draws.
# @param anchor_index  Grid position of each anchor.
# @param ess_fraction  Minimum importance ESS as a fraction of the draws.
# @param anchor_exact  List over anchors of list(pip, mcse_pip), each length
#                      E, or NULL to leave anchor positions pooled.
# @param prior_odds    Per-edge prior inclusion odds, for the log Bayes factor
#                      in the dispute check; NULL skips that condition.
# @param tolerance     Log Bayes factor difference for a dispute.
# @param q_level       Cochran's Q quantile for a dispute.
# @param anchor_mult   Anchor scale multipliers (for the nearest anchor).
# @param grid_mult     Grid scale multipliers (for the nearest anchor).
#
# Returns: list(pip = P x E, var = P x E, chain_pip = list C of P x E,
#   anchor_used = highest-ESS contributing anchor per point (NA if masked),
#   ess = length-P best ESS, disagree = P x E logical, n_contrib = length-P,
#   floor = per-anchor ESS floor applied).
# ------------------------------------------------------------------
assemble_curve = function(reweights, usable, ess_floor, anchor_index = NULL,
                          ess_fraction = 0, anchor_exact = NULL,
                          prior_odds = NULL, tolerance = 0.5 * log(10),
                          q_level = 0.99, anchor_mult = NULL, grid_mult = NULL) {
  n_grid = length(reweights[[1]]$ess)
  n_edge = ncol(reweights[[1]]$pip)
  n_chain = length(reweights[[1]]$chain_pip)
  n_anchor = length(reweights)
  ess_mat = vapply(reweights, `[[`, numeric(n_grid), "ess")
  if(is.null(dim(ess_mat))) ess_mat = matrix(ess_mat, nrow = 1L)
  n_draws = vapply(reweights, function(r) as.numeric(r$n %||% NA_real_), numeric(1))
  n_draws[is.na(n_draws)] = 0
  floor_a = pmax(ess_floor, ess_fraction * n_draws)
  pass = ess_mat >= matrix(floor_a, n_grid, n_anchor, byrow = TRUE)
  pass[, !usable] = FALSE

  # Distance between grid points and anchors, for the nearest anchor.
  dist = if(!is.null(anchor_mult) && !is.null(grid_mult)) {
    abs(outer(log(grid_mult), log(anchor_mult), "-"))
  } else {
    abs(outer(seq_len(n_grid), anchor_index %||% seq_len(n_anchor), "-"))
  }
  clamp = function(p) pmin(pmax(p, 1e-6), 1 - 1e-6)

  pip = matrix(NA_real_, n_grid, n_edge)
  var = matrix(NA_real_, n_grid, n_edge)
  disagree = matrix(FALSE, n_grid, n_edge)
  chain_pip = replicate(n_chain, matrix(NA_real_, n_grid, n_edge), simplify = FALSE)
  anchor_used = rep(NA_integer_, n_grid)
  best_ess = rep(NA_real_, n_grid)
  n_contrib = integer(n_grid)

  for(p in seq_len(n_grid)) {
    contrib = which(pass[p, ])
    if(length(contrib) == 0L) next
    n_contrib[p] = length(contrib)
    best_ess[p] = max(ess_mat[p, contrib])
    anchor_used[p] = contrib[which.max(ess_mat[p, contrib])]
    P_a = matrix(vapply(contrib, function(a) reweights[[a]]$pip[p, ], numeric(n_edge)), nrow = n_edge)
    V_a = matrix(vapply(contrib, function(a) reweights[[a]]$var[p, ], numeric(n_edge)), nrow = n_edge)
    W = 1 / V_a
    wsum = rowSums(W)
    pooled = rowSums(W * P_a) / wsum
    pooled_var = 1 / wsum
    disp = rep(FALSE, n_edge)
    if(length(contrib) > 1L) {
      Q = rowSums(W * (P_a - pooled)^2)
      stat = Q > stats::qchisq(q_level, df = length(contrib) - 1L)
      practical = if(is.null(prior_odds)) {
        TRUE
      } else {
        lbf = log((clamp(P_a) / (1 - clamp(P_a))) / prior_odds)
        apply(lbf, 1, function(v) diff(range(v))) > tolerance
      }
      disp = stat & practical
    }
    nearest = contrib[which.min(dist[p, contrib])]
    pip[p, ] = ifelse(disp, reweights[[nearest]]$pip[p, ], pooled)
    var[p, ] = ifelse(disp, reweights[[nearest]]$var[p, ], pooled_var)
    disagree[p, ] = disp
    for(cc in seq_len(n_chain)) {
      cnum = numeric(n_edge)
      for(k in seq_along(contrib)) {
        cnum = cnum + W[, k] * reweights[[contrib[k]]]$chain_pip[[cc]][p, ]
      }
      chain_pip[[cc]][p, ] = ifelse(disp, reweights[[nearest]]$chain_pip[[cc]][p, ], cnum / wsum)
    }
  }

  # An anchor's own grid position carries that fit's own estimate.
  for(a in seq_along(anchor_index)) {
    idx = anchor_index[a]
    if(!usable[a] || !pass[idx, a]) next
    anchor_used[idx] = a
    if(!is.null(anchor_exact) && !is.null(anchor_exact[[a]])) {
      ex = anchor_exact[[a]]
      pip[idx, ] = ex$pip
      v = ex$mcse_pip^2
      v[!is.finite(v)] = reweights[[a]]$var[idx, !is.finite(v)]
      var[idx, ] = v
      disagree[idx, ] = FALSE
      for(cc in seq_len(n_chain)) chain_pip[[cc]][idx, ] = reweights[[a]]$chain_pip[[cc]][idx, ]
    }
  }
  list(
    pip = pip, var = var, chain_pip = chain_pip, anchor_used = anchor_used,
    ess = best_ess, disagree = disagree, n_contrib = n_contrib, floor = floor_a
  )
}
