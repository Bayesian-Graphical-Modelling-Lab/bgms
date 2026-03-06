# =============================================================================
# test-mixed-mrf-sampling.R — Phase B.5 sampler correctness tests
# =============================================================================
# Validates that do_one_metropolis_step() in MixedMRFModel produces correct
# posterior samples by comparing bgms against the mixedGM reference sampler.
#
# Because no proposal adaptation is implemented yet (Phase F), acceptance
# rates are low with the default proposal SD = 1.0. Tests use correlation
# thresholds rather than absolute agreement.

skip_if_not_installed("mixedGM")

# ---------------------------------------------------------------------------
# Helper: extract parameter matrices from bgms sample matrix
# ---------------------------------------------------------------------------
# The vectorization order in bgms is:
#   1. mux: per-variable, C_s thresholds each (for ordinal)
#   2. Kxx: upper-triangular, row-major
#   3. muy: q means
#   4. Kyy: upper-triangle including diagonal
#   5. Kxy: p*q entries, row-major
extract_bgms_estimates = function(samples, p, q, num_cats_bgms) {
  S = colMeans(samples)
  idx = 1L

  # 1. mux (per-variable: C_s thresholds each)
  max_cat = max(num_cats_bgms)
  mux = matrix(0, p, max_cat)
  for(s in 1:p) {
    for(c in seq_len(num_cats_bgms[s])) {
      mux[s, c] = S[idx]; idx = idx + 1L
    }
  }

  # 2. Kxx upper-tri
  Kxx = matrix(0, p, p)
  for(i in 1:(p - 1)) for(j in (i + 1):p) {
    Kxx[i, j] = Kxx[j, i] = S[idx]; idx = idx + 1L
  }

  # 3. muy
  muy = S[idx:(idx + q - 1)]; idx = idx + q

  # 4. Kyy upper-tri + diag
  Kyy = matrix(0, q, q)
  for(i in 1:q) for(j in i:q) {
    Kyy[i, j] = Kyy[j, i] = S[idx]; idx = idx + 1L
  }

  # 5. Kxy row-major
  Kxy = matrix(S[idx:(idx + p * q - 1)], nrow = p, ncol = q, byrow = TRUE)

  list(mux = mux, Kxx = Kxx, muy = muy, Kyy = Kyy, Kxy = Kxy)
}


# ===========================================================================
# Test 1: ordinal simulation recovery (3-category, p=3, q=2, n=1000)
# ===========================================================================
test_that("conditional PL sampler recovers parameters (ordinal simulation)", {
  skip_on_cran()

  p = 3L; q = 2L; n = 1000L

  Kxx_true = matrix(c(0, 0.4, -0.2,
                       0.4, 0, 0.3,
                       -0.2, 0.3, 0), 3, 3)
  Kxy_true = matrix(c(0.2, -0.1, 0.15,
                       0.1, 0.3, -0.2), nrow = p, ncol = q)
  Kyy_true = matrix(c(1.3, 0.2, 0.2, 1.1), 2, 2)
  mux_true = matrix(c(0.5, -0.5, 0.3,
                       -0.3, -1.0, -0.5), nrow = p, ncol = 2)
  muy_true = c(0.5, -0.3)

  set.seed(456)
  sim = mixedGM::mixed_gibbs_generate(
    n = n, Kxx = Kxx_true, Kxy = Kxy_true, Kyy = Kyy_true,
    mux = mux_true, muy = muy_true,
    num_categories = rep(3L, p), n_burnin = 1000
  )

  # --- run mixedGM (with Robbins-Monro adaptation) ---
  set.seed(42)
  mgm = mixedGM::mixed_sampler(
    x = sim$x, y = sim$y, num_categories = rep(3L, p),
    n_warmup = 4000L, n_samples = 8000L,
    edge_selection = FALSE, verbose = FALSE,
    pseudolikelihood = "conditional"
  )
  mgm_est = mgm$final_parameters

  # --- run bgms (no adaptation) ---
  num_cats_bgms = rep(2L, p)  # max category index
  pq = p + q
  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  bgms_out = bgms:::test_mixed_mrf_sampler(
    sim$x, sim$y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, FALSE, "conditional",
    4000L, 8000L, 42L
  )
  bgms_est = extract_bgms_estimates(bgms_out$samples, p, q, num_cats_bgms)

  # --- correlations bgms vs mixedGM (r >= 0.95) ---
  kxx_ut = upper.tri(bgms_est$Kxx)
  expect_gt(
    cor(bgms_est$Kxx[kxx_ut], mgm_est$Kxx[kxx_ut]), 0.95,
    label = "Kxx bgms-mgm correlation"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kxy), as.vector(mgm_est$Kxy)), 0.95,
    label = "Kxy bgms-mgm correlation"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kyy), as.vector(mgm_est$Kyy)), 0.95,
    label = "Kyy bgms-mgm correlation"
  )
  expect_gt(
    cor(as.vector(bgms_est$mux), as.vector(mgm_est$mux)), 0.95,
    label = "mux bgms-mgm correlation"
  )

  # --- signs agree with truth ---
  expect_true(
    all(sign(bgms_est$Kxx[kxx_ut]) == sign(Kxx_true[kxx_ut])),
    label = "Kxx signs match truth"
  )
})


# ===========================================================================
# Test 2: real data cross-validation (Wenchuan, ordinal + continuous)
# ===========================================================================
test_that("conditional PL sampler agrees with mixedGM on Wenchuan data", {
  skip_on_cran()

  data(Wenchuan, package = "bgms")
  W = na.omit(Wenchuan)

  p = 4L; q = 2L
  x = as.matrix(W[, 1:p]) - 1L   # convert 1-5 to 0-4
  y = as.matrix(W[, (p + 1):(p + q)])

  num_cats_mgm = rep(5L, p)   # total categories for mixedGM
  num_cats_bgms = rep(4L, p)  # max category index for bgms

  # --- mixedGM ---
  set.seed(42)
  mgm = mixedGM::mixed_sampler(
    x = x, y = y, num_categories = num_cats_mgm,
    n_warmup = 4000L, n_samples = 8000L,
    edge_selection = FALSE, verbose = FALSE,
    pseudolikelihood = "conditional"
  )
  mgm_est = mgm$final_parameters

  # --- bgms ---
  pq = p + q
  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  bgms_out = bgms:::test_mixed_mrf_sampler(
    x, y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, FALSE, "conditional",
    4000L, 8000L, 42L
  )
  bgms_est = extract_bgms_estimates(bgms_out$samples, p, q, num_cats_bgms)

  # --- correlations bgms vs mixedGM (r >= 0.95) ---
  kxx_ut = upper.tri(bgms_est$Kxx)
  expect_gt(
    cor(bgms_est$Kxx[kxx_ut], mgm_est$Kxx[kxx_ut]), 0.95,
    label = "Kxx bgms-mgm correlation (Wenchuan)"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kxy), as.vector(mgm_est$Kxy)), 0.95,
    label = "Kxy bgms-mgm correlation (Wenchuan)"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kyy), as.vector(mgm_est$Kyy)), 0.95,
    label = "Kyy bgms-mgm correlation (Wenchuan)"
  )
  expect_gt(
    cor(as.vector(bgms_est$mux), as.vector(mgm_est$mux)), 0.95,
    label = "mux bgms-mgm correlation (Wenchuan)"
  )

  # --- Kxx sign patterns agree ---
  expect_true(
    all(sign(bgms_est$Kxx[kxx_ut]) == sign(mgm_est$Kxx[kxx_ut])),
    label = "Kxx signs agree (Wenchuan)"
  )
})


# ===========================================================================
# Test 3: sampler does not crash on Kyy updates (c3 regression test)
# ===========================================================================
test_that("Kyy updates do not crash (c3 get_constants regression)", {
  # Minimal setup: p=2, q=2 binary triggers the c3 code path
  # where L(qm1,qm1)^2 must be included in the sum.
  set.seed(99)
  n = 50L; p = 2L; q = 2L
  x = matrix(sample(0:1, n * p, replace = TRUE), n, p)
  y = matrix(rnorm(n * q), n, q)

  pq = p + q
  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  # If c3 bug is present, this crashes with "chol(): decomposition failed"
  result = bgms:::test_mixed_mrf_sampler(
    x, y, as.integer(c(1L, 1L)),
    as.integer(c(1L, 1L)), c(0L, 0L),
    inc_prob, edge_ind, FALSE, "conditional",
    100L, 100L, 123L
  )

  expect_equal(nrow(result$samples), 100L)
  expect_true(all(is.finite(result$samples)))
})
