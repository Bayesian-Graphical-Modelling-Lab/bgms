# --------------------------------------------------------------------------- #
# Row-block conjugate Gibbs within-step (update_method = "Gibbs").
#
# One block per feature increment landed on feat/row-block-gibbs-within, so a
# regression in any single commit is isolated to its block:
#
#   a9068d9  Normal slab, alpha = 1, delta = 0  -> exact conjugate draw anchor
#   04e8122  within_step_kind dispatch wiring   -> update_method = "Gibbs" routes
#   041e77d  delta != 0 (xi shape shift)        -> AM agreement at delta = 0.5
#   64f952e  Gamma(alpha != 1) independent-MH   -> AM agreement at alpha = 2
#   c193dec  Cauchy slab via omega mixture      -> AM agreement under cauchy_prior
#   d3ca987  merge into update_method           -> covered by the dispatch block
#
# The conjugate-draw block is the strongest, lowest-variance anchor: it calls
# the C++ test entry ggm_test_row_block_gibbs directly and checks the draw
# against its closed-form Gamma marginal. The remaining blocks compare Gibbs
# vs adaptive-Metropolis posterior means on a *fixed complete graph*
# (edge_selection = FALSE). Both samplers target the identical K | G posterior,
# so their means must agree within MCMC error -- this catches a wrong target
# in the Gibbs kernel that a "runs without error" smoke test would miss.
# --------------------------------------------------------------------------- #


# ---- Exact conjugate anchor (commit a9068d9) -------------------------------- #

test_that("row-block Gibbs draws K_ii from its closed-form Gamma at prior-only, q=0", {
  # Empty graph -> every row has zero active neighbours, so the (beta, xi)
  # block collapses to gamma = xi alone. At n = 0 (prior-only) the marginal is
  # xi ~ Gamma(shape = n/2 + delta + 1 = 1, rate = S_ii/2 + gamma_rate/2 =
  # gamma_rate/2). With gamma_rate = 2 that is Gamma(1, 1), mean 1.
  p = 4L
  edges = matrix(0L, p, p)
  S = matrix(0, p, p)          # n = 0: sufficient statistics are irrelevant
  beta0 = 2.0
  n_sweeps = 5000L

  out = ggm_test_row_block_gibbs(
    suf_stat = S, n = 0L, edge_indicators = edges,
    pairwise_scale = 1.0, gamma_shape = 1.0, gamma_rate = beta0,
    n_sweeps = n_sweeps, seed = 42L
  )

  burn = 201:n_sweeps          # forget the identity init
  K_ii = sapply(seq_len(p), function(i) out$K_samples[i, i, burn])

  # Moments: Gamma(shape = 1, rate = beta0 / 2) -> mean 2/beta0, var 4/beta0^2.
  theory_mean = 2.0 / beta0
  theory_var  = 4.0 / beta0^2
  tol_mean = 4 * sqrt(theory_var / length(burn))   # 4 SE, generous for CI
  expect_true(
    all(abs(colMeans(K_ii) - theory_mean) < tol_mean),
    info = sprintf("means: %s (theory %.4f, tol %.4f)",
                   paste(sprintf("%.4f", colMeans(K_ii)), collapse = " "),
                   theory_mean, tol_mean)
  )

  # Distributional: KS against the analytic Gamma, Bonferroni over the p rows.
  ks_p = vapply(seq_len(p), function(i) {
    suppressWarnings(
      stats::ks.test(K_ii[, i], "pgamma", shape = 1, rate = beta0 / 2)$p.value
    )
  }, numeric(1L))
  expect_gt(min(ks_p), 0.01 / p)
})


# ---- Gibbs-vs-AM agreement helpers ------------------------------------------ #

# Posterior means of the precision matrix (diagonal + off-diagonal) from a
# fixed-graph bgm() fit. With edge_selection = FALSE the graph stays complete,
# so indicators are all 1 and pairwise samples are the raw K_ij.
ggm_K_means = function(Y, update_method, interaction_prior, alpha, delta,
                       iter, warmup, seed = 1L) {
  fit = bgm(
    Y, variable_type = "continuous",
    interaction_prior = interaction_prior,
    precision_scale_prior = gamma_prior(shape = alpha, rate = 1),
    delta = delta,
    prior_factorization = "joint",
    edge_selection = FALSE,
    iter = iter, warmup = warmup,
    update_method = update_method,
    chains = 1L, cores = 1L, seed = seed,
    display_progress = "none", verbose = FALSE
  )
  raw = S7::prop(fit, "raw_samples")
  list(
    diag = colMeans(raw$main[[1L]]),
    off  = colMeans(raw$pairwise[[1L]])
  )
}

# Shared test fixture: a small dense GGM, centered as bgm() expects.
ggm_agreement_data = function(seed = 7L, p = 6L, n = 200L) {
  set.seed(seed)
  K_true = diag(p) + 0.3 * (abs(row(diag(p)) - col(diag(p))) == 1)
  Y = MASS::mvrnorm(n, mu = rep(0, p), Sigma = solve(K_true))
  Y = scale(Y, center = TRUE, scale = FALSE)
  colnames(Y) = paste0("V", seq_len(p))
  Y
}

# Assert two K-mean lists agree within MCMC error.
expect_K_agreement = function(a, b, tol_abs, tol_rel) {
  d_diag = abs(a$diag - b$diag)
  d_off  = abs(a$off - b$off)
  max_abs = max(d_diag, d_off)
  rel = mean(c(d_diag, d_off)) / mean(abs(c(b$diag, b$off)))
  expect_lt(max_abs, tol_abs)
  expect_lt(rel, tol_rel)
}


# ---- Dispatch wiring + Normal alpha=1 delta=0 (commits 04e8122 / a9068d9) --- #

test_that("update_method = 'Gibbs' routes and agrees with AM (Normal, alpha=1, delta=0)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  args = list(Y = Y, interaction_prior = normal_prior(scale = 1),
              alpha = 1, delta = 0, iter = 2500L, warmup = 600L)
  am   = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "Gibbs"))
  expect_K_agreement(gibb, am, tol_abs = 0.06, tol_rel = 0.025)
})


# ---- delta != 0 (commit 041e77d) -------------------------------------------- #

test_that("Gibbs agrees with AM at delta = 0.5 (xi shape shift)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  args = list(Y = Y, interaction_prior = normal_prior(scale = 1),
              alpha = 1, delta = 0.5, iter = 2500L, warmup = 600L)
  am   = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "Gibbs"))
  expect_K_agreement(gibb, am, tol_abs = 0.06, tol_rel = 0.025)
})


# ---- Gamma(alpha != 1) (commit 64f952e) ------------------------------------- #

test_that("Gibbs agrees with AM at alpha = 2 (independent-MH on gamma)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  args = list(Y = Y, interaction_prior = normal_prior(scale = 1),
              alpha = 2, delta = 0, iter = 2500L, warmup = 600L)
  am   = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "Gibbs"))
  # alpha != 1 adds an MH correction on the diagonal; allow a touch more slack.
  expect_K_agreement(gibb, am, tol_abs = 0.07, tol_rel = 0.035)
})


# ---- Cauchy slab via omega scale-mixture (commit c193dec) ------------------- #

test_that("Gibbs agrees with AM under a Cauchy slab (omega scale-mixture)", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  Y = ggm_agreement_data()
  # Cauchy mixing is slower than Normal: longer chain, looser tolerance.
  args = list(Y = Y, interaction_prior = cauchy_prior(scale = 1),
              alpha = 1, delta = 0, iter = 4000L, warmup = 1000L)
  am   = do.call(ggm_K_means, c(args, update_method = "adaptive-metropolis"))
  gibb = do.call(ggm_K_means, c(args, update_method = "Gibbs"))
  expect_K_agreement(gibb, am, tol_abs = 0.10, tol_rel = 0.05)
})
