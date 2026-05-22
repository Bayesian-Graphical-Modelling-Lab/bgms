# --------------------------------------------------------------------------- #
# Conditional SBC for the Graphical G-prior on GGM
#
# Standard SBC fails for the GG-prior because V_ij = n / (4·S̄_ii·S̄_jj)
# depends on the data: the prior used to draw K_true and the prior used
# in the bgm fit on the simulated data have different V_ij values, so
# ranks are not uniform by construction.
#
# The conditional SBC fixes V_ij at a reference V_ref (computed once
# from some X_ref) and routes BOTH the generator and the fitter
# through that same matrix via the V_ij_external knob added to
# graphical_g_prior(). Generator and fitter now share the same prior,
# so under a correct sampler the SBC ranks ARE uniform and the
# inclusion-PIP consistency E[P(γ_ij | y)] = P(γ_ij | prior) holds.
#
# Two tests:
#   1. Rank uniformity on diagonals of K (continuous, no point mass).
#      Off-diagonals are skipped because spike-and-slab tied ranks
#      break standard uniformity tests.
#   2. Inclusion-PIP consistency: the across-replication mean
#      posterior PIP equals the prior PIP (estimated from a much
#      longer independent prior-only chain so the comparison reference
#      is tight). Joint MC SE accounts for uncertainty in BOTH the
#      across-rep mean and the prior-chain estimate.
#
# Gated behind BGMS_RUN_SLOW_TESTS.
# --------------------------------------------------------------------------- #


skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run SBC tests"
  )
}


# Reconstruct a p x p symmetric K from upper-triangle (i <= j, row-major)
# raw samples produced by bgm() under variable_type = "continuous".
reconstruct_K = function(upper_row, p) {
  K = matrix(0, p, p)
  idx = 1L
  for(i in seq_len(p)) {
    for(j in i:p) {
      K[i, j] = upper_row[idx]
      K[j, i] = upper_row[idx]
      idx = idx + 1L
    }
  }
  K
}


# Build the upper-triangle (i < j) row-major off-diagonal index used to
# pull edges out of bgm()'s indicator_samples matrix.
offdiag_index = function(p) {
  out = integer(p * (p - 1L) / 2L)
  pos = 0L; oi = 0L
  for(i in seq_len(p)) for(j in i:p) {
    pos = pos + 1L
    if(i != j) { oi = oi + 1L; out[oi] = pos }
  }
  out
}


# Compute V_ij = n / (4 * S_ii * S_jj) from a data matrix.
compute_V_ij = function(X) {
  n = nrow(X)
  p = ncol(X)
  S_diag = colSums(X^2)
  V = matrix(0, p, p)
  for(i in seq_len(p)) for(j in seq_len(p)) {
    if(i != j) V[i, j] = n / (4 * S_diag[i] * S_diag[j])
  }
  V
}


test_that("conditional SBC: GG-prior diagonal ranks are uniform (p=3)", {
  skip_unless_slow()

  # Diagonal SBC: with V_ij fixed at V_ref on both sides, ranks of
  # K_true diagonal entries within posterior samples should be uniform.

  p = 3L
  n = 100L
  R = 200L
  L = 999L
  thin = 5L

  set.seed(2026)
  X_ref = MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  V_ref = compute_V_ij(X_ref)

  prior_obj = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                 a0 = 1, b0 = 1, g_init = 1,
                                 V_ij_external = V_ref)

  setup_fit = bgm(x = X_ref, variable_type = "continuous",
                  interaction_prior = prior_obj,
                  edge_selection = TRUE, delta = NULL,
                  iter = 1L, warmup = 1L, chains = 1L, cores = 1L,
                  update_method = "adaptive-metropolis",
                  seed = 2026L, display_progress = "none")

  # Generator chain (thinned).
  spec = setup_fit$.bgm_spec
  spec$prior$prior_only = TRUE
  spec$sampler$iter = as.integer(R * thin)
  spec$sampler$warmup = 2000L
  spec$sampler$chains = 1L
  spec$sampler$update_method = "adaptive-metropolis"
  spec$sampler$target_accept = 0.44
  spec$sampler$progress_type = 0L
  gen_raw = run_sampler(spec)
  K_true_samples = gen_raw[[1L]]$samples

  thin_idx = seq(thin, R * thin, by = thin)
  K_true_list = lapply(thin_idx, function(it) {
    reconstruct_K(K_true_samples[, it], p)
  })

  ranks = matrix(NA_integer_, nrow = R, ncol = p)
  for(r in seq_len(R)) {
    K_true = K_true_list[[r]]
    eig = tryCatch(eigen(K_true, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NULL)
    if(is.null(eig) || min(eig) <= 1e-8) next
    Sigma = solve(K_true)
    y_star = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    fit_r = bgm(x = y_star, variable_type = "continuous",
                interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                                       a0 = 1, b0 = 1, g_init = 1,
                                                       V_ij_external = V_ref),
                edge_selection = TRUE, delta = NULL,
                iter = L, warmup = 1000L, chains = 1L, cores = 1L,
                update_method = "adaptive-metropolis",
                seed = 2026L + r, display_progress = "none")
    main_samples = do.call(rbind, fit_r$raw_samples$main)
    for(i in seq_len(p)) ranks[r, i] = sum(main_samples[, i] < K_true[i, i])
  }

  ranks = ranks[stats::complete.cases(ranks), , drop = FALSE]
  expect_true(nrow(ranks) > R * 0.9,
              info = sprintf("Dropped too many reps: %d / %d kept.",
                             nrow(ranks), R))

  # Per-parameter KS test (alpha = 0.01) — allow up to 1 false positive.
  n_fail = sum(vapply(seq_len(p), function(j) {
    u = ranks[, j] / (L + 1)
    suppressWarnings(stats::ks.test(u, "punif")$p.value) <= 0.01
  }, logical(1L)))
  expect_true(n_fail <= 1L,
              info = sprintf("KS failures: %d / %d at alpha=0.01.", n_fail, p))

  # Pooled chi-squared on aggregated ranks.
  all_u = as.vector(ranks) / (L + 1)
  counts = tabulate(cut(all_u, seq(0, 1, length.out = 21),
                        include.lowest = TRUE), nbins = 20L)
  chisq_p = suppressWarnings(stats::chisq.test(counts)$p.value)
  expect_true(chisq_p > 0.001,
              info = sprintf("Pooled chi-squared p = %.4f", chisq_p))
})


test_that("conditional SBC: GG-prior inclusion PIPs are calibrated (p=10)", {
  skip_unless_slow()

  # Inclusion-PIP consistency: E[P(γ_ij | y)] across SBC replications
  # equals the prior P(γ_ij | V_ref) when generator and fitter share
  # V_ref. Earlier versions of this test used a 1000-iter prior chain
  # as the reference and reported huge "biases" — those were prior-
  # chain MC noise wearing posterior-chain clothes. We now use a 50k-
  # iter prior chain and a joint MC SE that includes both sources.

  p = 10L
  n = 200L
  R = 200L
  L = 999L
  thin = 5L
  n_prior_iter = 50000L

  set.seed(2026)
  X_ref = MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  V_ref = compute_V_ij(X_ref)

  prior_obj = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                 a0 = 1, b0 = 1, g_init = 1,
                                 V_ij_external = V_ref)
  setup_fit = bgm(x = X_ref, variable_type = "continuous",
                  interaction_prior = prior_obj,
                  edge_selection = TRUE, delta = NULL,
                  iter = 1L, warmup = 1L, chains = 1L, cores = 1L,
                  update_method = "adaptive-metropolis",
                  seed = 2026L, display_progress = "none")

  # Long independent prior-only chain to estimate the prior PIPs tightly.
  spec_prior = setup_fit$.bgm_spec
  spec_prior$prior$prior_only = TRUE
  spec_prior$sampler$iter = n_prior_iter
  spec_prior$sampler$warmup = 5000L
  spec_prior$sampler$chains = 1L
  spec_prior$sampler$update_method = "adaptive-metropolis"
  spec_prior$sampler$target_accept = 0.44
  spec_prior$sampler$progress_type = 0L
  prior_raw = run_sampler(spec_prior)
  ind_samples = prior_raw[[1L]]$indicator_samples
  off_idx = offdiag_index(p)
  prior_pip = rowMeans(ind_samples[off_idx, , drop = FALSE])

  # Generator chain (separate seed) producing thinned K_true draws.
  spec_gen = spec_prior
  spec_gen$sampler$iter = as.integer(R * thin)
  spec_gen$sampler$warmup = 3000L
  spec_gen$sampler$seed = 99L
  gen_raw = run_sampler(spec_gen)
  thin_idx = seq(thin, R * thin, by = thin)
  K_true_list = lapply(thin_idx, function(it) {
    reconstruct_K(gen_raw[[1L]]$samples[, it], p)
  })

  # Fitter loop: record posterior PIPs per replication.
  post_pips = matrix(NA_real_, R, p * (p - 1L) / 2L)
  for(r in seq_len(R)) {
    K_true = K_true_list[[r]]
    eig = tryCatch(eigen(K_true, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NULL)
    if(is.null(eig) || min(eig) <= 1e-8) next
    Sigma = solve(K_true)
    y_star = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    fit_r = bgm(x = y_star, variable_type = "continuous",
                interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                                       a0 = 1, b0 = 1, g_init = 1,
                                                       V_ij_external = V_ref),
                edge_selection = TRUE, delta = NULL,
                iter = L, warmup = 1500L, chains = 1L, cores = 1L,
                update_method = "adaptive-metropolis",
                seed = 2026L + r, display_progress = "none")
    m = fit_r$posterior_mean_indicator
    post_pips[r, ] = m[upper.tri(m)]
  }

  post_pips = post_pips[stats::complete.cases(post_pips), , drop = FALSE]
  expect_true(nrow(post_pips) > R * 0.9,
              info = sprintf("Dropped too many reps: %d / %d kept.",
                             nrow(post_pips), R))

  mean_post = colMeans(post_pips)
  diff = mean_post - prior_pip
  # Joint MC SE: across-rep variance of mean_post + binomial-ish variance
  # of prior_pip from a chain of effective length ~ n_prior_iter / acf_lag.
  # acf_lag is a conservative guess (50 for indicator chains).
  n_eff_prior = n_prior_iter / 50
  mc_se_post = apply(post_pips, 2, sd) / sqrt(nrow(post_pips))
  mc_se_prior = sqrt(prior_pip * (1 - prior_pip) / n_eff_prior)
  mc_se = sqrt(mc_se_post^2 + mc_se_prior^2)
  zscore = diff / mc_se

  n_edges = length(zscore)
  # Under exact calibration: # |z| > 2 ~ 0.05 * n_edges, # |z| > 3 ~ 0.003 * n_edges.
  # Allow generous tail counts (each test is alpha=0.05 family-wise).
  expect_true(sum(abs(zscore) > 2) <= 0.20 * n_edges,
              info = sprintf("Too many edges with |z| > 2: %d / %d",
                             sum(abs(zscore) > 2), n_edges))
  expect_true(sum(abs(zscore) > 3) <= 0.05 * n_edges,
              info = sprintf("Too many edges with |z| > 3: %d / %d",
                             sum(abs(zscore) > 3), n_edges))
  # Median absolute z should be near 0.67 (the median of a standard half-normal).
  expect_lt(median(abs(zscore)), 1.5)
})
