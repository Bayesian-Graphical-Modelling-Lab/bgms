# --------------------------------------------------------------------------- #
# SBC for GGM (no edge selection)
#
# Simulation-based calibration checks that posterior ranks are uniform
# when data is generated from the prior predictive distribution. This
# validates the NUTS and Adaptive-Metropolis samplers for the GGM with
#
# p = 3, R = 200 replications, L = 999 posterior draws per replication.
# Prior: Cauchy(0, 2.5) off-diagonal, Gamma(1, 1) diagonal. PD
# constraint enforced by rejection sampling.
#
# Edge selection is tested separately via parameter recovery (PR.2).
# Omitting edge selection here avoids the tied-rank issue inherent
# in spike-and-slab SBC (discrete indicators + point-mass at zero
# make standard uniformity tests unreliable).
#
# Uniformity tested by KS test at alpha = 0.01 per parameter.
# Global chi-squared test as a fallback.
#
# Gated behind BGMS_RUN_SLOW_TESTS.
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

skip_unless_slow_sbc = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run SBC tests"
  )
}


# ---- Prior sampler -----------------------------------------------------------

# Draw a symmetric PD precision matrix from the GGM prior (no edge
# selection):
#   K_ij ~ Cauchy(0, scale)  (off-diagonal, i < j, symmetrised)
#   K_ii ~ Gamma(1, 1)
# Rejection-samples until the result is positive definite.
draw_prior_K = function(p, scale = 2.5, max_tries = 10000) {
  for(attempt in seq_len(max_tries)) {
    K = matrix(0, p, p)

    # Off-diagonal (upper triangle)
    for(i in seq_len(p - 1)) {
      for(j in (i + 1):p) {
        K[i, j] = rcauchy(1, 0, scale)
        K[j, i] = K[i, j]
      }
    }

    # Diagonal
    for(i in seq_len(p)) {
      K[i, i] = rgamma(1, shape = 1, rate = 1)
    }

    # Check positive definiteness
    ev = eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if(all(ev > 1e-8)) {
      return(K)
    }
  }
  stop("draw_prior_K: failed to draw PD matrix in ", max_tries, " attempts")
}


# ---- Rank computation --------------------------------------------------------

# Compute SBC rank: number of posterior draws less than the true value.
# Returns a named vector of ranks (one per parameter).
compute_sbc_ranks = function(K_true, p, fit, thin_idx = NULL) {
  # Off-diagonal precision entries (raw precision scale)
  pw_samples = do.call(rbind, fit$raw_samples$pairwise)
  if(!is.null(thin_idx)) pw_samples = pw_samples[thin_idx, , drop = FALSE]

  # Diagonal precision entries
  main_samples = do.call(rbind, fit$raw_samples$main)
  if(!is.null(thin_idx)) main_samples = main_samples[thin_idx, , drop = FALSE]

  ranks = numeric(0)
  names_out = character(0)

  # Off-diagonal K entries
  col_idx = 0
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      col_idx = col_idx + 1
      true_k = K_true[i, j]
      ranks = c(ranks, sum(pw_samples[, col_idx] < true_k))
      names_out = c(names_out, paste0("K_", i, j))
    }
  }

  # Diagonal K entries
  for(i in seq_len(p)) {
    true_k = K_true[i, i]
    ranks = c(ranks, sum(main_samples[, i] < true_k))
    names_out = c(names_out, paste0("K_", i, i))
  }

  names(ranks) = names_out
  ranks
}


# ---- SBC test ----------------------------------------------------------------

test_that("SBC: GGM NUTS produces uniform ranks (p=3, no edge selection)", {
  skip_unless_slow_sbc()

  p = 3
  n = 100
  R = 200
  L = 999
  scale = 2.5

  set.seed(2026)

  # Pre-draw all prior K matrices to separate randomness
  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K(p, scale)
  }

  # Storage: 3 off-diagonal + 3 diagonal = 6 parameters
  n_off = p * (p - 1) / 2
  n_params = n_off + p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    K_true = prior_draws[[r]]
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L, warmup = 1000, chains = 1,
      edge_selection = FALSE, update_method = "nuts",
      pairwise_scale = scale,
      display_progress = "none", seed = 2026L + r
    )

    ranks[r, ] = compute_sbc_ranks(K_true, p, fit)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  # Per-parameter KS test: ranks / (L + 1) ~ Uniform(0, 1)
  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  # At alpha=0.01 with 6 parameters, allow at most 1 false positive
  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS: %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  # Global chi-squared on aggregated ranks (20 bins)
  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared p=%.4f", chisq_p)
  )
})


# ---- SBC test: Adaptive Metropolis ------------------------------------------

test_that("SBC: GGM MH produces uniform ranks (p=3, no edge selection)", {
  skip_unless_slow_sbc()

  p = 3
  n = 100
  R = 200
  L = 999
  thin = 5
  L_raw = L * thin
  scale = 2.5

  set.seed(2027)

  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K(p, scale)
  }

  n_off = p * (p - 1) / 2
  n_params = n_off + p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    K_true = prior_draws[[r]]
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L_raw, warmup = 5000, chains = 1,
      edge_selection = FALSE, update_method = "adaptive-metropolis",
      pairwise_scale = scale,
      display_progress = "none", seed = 2027L + r
    )

    # Thin to reduce autocorrelation in MH chains
    thin_idx = seq(1, L_raw, by = thin)
    ranks[r, ] = compute_sbc_ranks(K_true, p, fit, thin_idx = thin_idx)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS: %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared (MH) p=%.4f", chisq_p)
  )
})


# ---- Prior sampler with edge selection ---------------------------------------

# Draw a precision matrix from the spike-and-slab GGM prior:
#   gamma_ij ~ Bernoulli(0.5)
#   K_ij | gamma_ij=1 ~ Cauchy(0, scale); K_ij | gamma_ij=0 = 0
#   K_ii ~ Gamma(1, 1)
# Rejection-samples until K is positive definite.
draw_prior_K_es = function(p, scale = 2.5, inclusion_prob = 0.5,
                           max_tries = 50000) {
  for(attempt in seq_len(max_tries)) {
    K = matrix(0, p, p)
    gamma = matrix(0L, p, p)

    for(i in seq_len(p - 1)) {
      for(j in (i + 1):p) {
        if(runif(1) < inclusion_prob) {
          gamma[i, j] = 1L
          gamma[j, i] = 1L
          K[i, j] = rcauchy(1, 0, scale)
          K[j, i] = K[i, j]
        }
      }
    }

    for(i in seq_len(p)) {
      K[i, i] = rgamma(1, shape = 1, rate = 1)
    }

    ev = eigen(K, symmetric = TRUE, only.values = TRUE)$values
    if(all(ev > 1e-8)) {
      return(list(K = K, gamma = gamma))
    }
  }
  stop("draw_prior_K_es: failed to draw PD matrix in ", max_tries, " attempts")
}


# ---- Diagonal rank computation for edge-selection SBC ------------------------

# Compute SBC ranks for diagonal K_ii only (avoids tied-rank issues
# from the spike-and-slab on off-diagonals).
compute_sbc_ranks_diag = function(K_true, p, fit, thin_idx = NULL) {
  main_samples = do.call(rbind, fit$raw_samples$main)
  if(!is.null(thin_idx)) main_samples = main_samples[thin_idx, , drop = FALSE]

  ranks = numeric(p)
  for(i in seq_len(p)) {
    ranks[i] = sum(main_samples[, i] < K_true[i, i])
  }
  names(ranks) = paste0("K_", seq_len(p), seq_len(p))
  ranks
}


# ---- SBC test: Edge selection (MH, diagonal elements) ------------------------

test_that("SBC: GGM MH produces uniform diagonal ranks (p=3, edge selection)", {
  skip_unless_slow_sbc()

  p = 3
  n = 100
  R = 200
  L = 999
  thin = 5
  L_raw = L * thin
  scale = 2.5

  set.seed(2028)

  prior_draws = vector("list", R)
  for(r in seq_len(R)) {
    prior_draws[[r]] = draw_prior_K_es(p, scale)
  }

  n_params = p
  ranks = matrix(NA_real_, nrow = R, ncol = n_params)

  for(r in seq_len(R)) {
    draw = prior_draws[[r]]
    K_true = draw$K
    Sigma = solve(K_true)
    X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    dat = as.data.frame(X)
    colnames(dat) = paste0("V", seq_len(p))

    fit = bgm(dat,
      variable_type = "continuous",
      iter = L_raw, warmup = 5000, chains = 1,
      edge_selection = TRUE, update_method = "adaptive-metropolis",
      pairwise_scale = scale,
      display_progress = "none", seed = 2028L + r
    )

    # Thin to reduce autocorrelation
    thin_idx = seq(1, L_raw, by = thin)
    ranks[r, ] = compute_sbc_ranks_diag(K_true, p, fit, thin_idx = thin_idx)
    if(r == 1) colnames(ranks) = names(ranks[1, ])
  }

  n_fail_ks = 0
  for(j in seq_len(ncol(ranks))) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1
  }

  max_fail = max(1, ceiling(n_params * 0.01 * 2))
  expect_true(n_fail_ks <= max_fail,
    info = sprintf(
      "SBC KS (edge selection, diag): %d/%d parameters failed (limit %d)",
      n_fail_ks, n_params, max_fail
    )
  )

  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
    breaks = seq(0, 1, length.out = 21),
    include.lowest = TRUE
  )
  counts = tabulate(bins, nbins = 20)
  chisq_p = chisq.test(counts)$p.value
  expect_true(chisq_p > 0.001,
    info = sprintf("SBC global chi-squared p=%.4f (edge selection, diag)", chisq_p)
  )
})
