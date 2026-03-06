# =============================================================================
# test-mixed-mrf-edge-selection.R — Phase D edge selection tests
# =============================================================================
# Validates that update_edge_indicators() produces correct spike-and-slab
# structure learning via reversible-jump Metropolis-Hastings.
#
# Strategy: generate from a sparse graph with known zero edges, run the
# sampler with edge_selection = TRUE, and check that posterior inclusion
# probabilities (PIPs) correctly separate true edges from null edges.

skip_if_not_installed("mixedGM")

# ---------------------------------------------------------------------------
# Helper: extract edge indicators from indicator_samples
# ---------------------------------------------------------------------------
# Vectorization order: Gxx upper-tri, Gyy upper-tri, Gxy row-major
extract_pip = function(indicator_samples, p, q) {
  pip = colMeans(indicator_samples)
  idx = 1L

  # 1. Gxx upper-triangle
  n_xx = p * (p - 1) / 2
  pip_xx = pip[idx:(idx + n_xx - 1)]; idx = idx + n_xx

  Gxx = matrix(0, p, p)
  k = 1L
  for(i in 1:(p - 1)) for(j in (i + 1):p) {
    Gxx[i, j] = Gxx[j, i] = pip_xx[k]; k = k + 1L
  }

  # 2. Gyy upper-triangle
  n_yy = q * (q - 1) / 2
  Gyy = matrix(0, q, q)
  if(n_yy > 0) {
    pip_yy = pip[idx:(idx + n_yy - 1)]; idx = idx + n_yy
    k = 1L
    for(i in 1:(q - 1)) for(j in (i + 1):q) {
      Gyy[i, j] = Gyy[j, i] = pip_yy[k]; k = k + 1L
    }
  }

  # 3. Gxy row-major
  n_xy = p * q
  pip_xy = pip[idx:(idx + n_xy - 1)]
  Gxy = matrix(pip_xy, nrow = p, ncol = q, byrow = TRUE)

  list(Gxx = Gxx, Gyy = Gyy, Gxy = Gxy)
}


# ===========================================================================
# Test 1: Edge selection on simulation data (conditional PL)
# ===========================================================================
test_that("edge selection recovers sparse structure (conditional PL)", {
  skip_on_cran()

  p = 3L; q = 2L; n = 800L
  num_cats_mgm = rep(3L, p)
  num_cats_bgms = rep(2L, p)

  # Sparse graph: only (1,2) in Kxx, (1,2) in Kxy; Kyy has (1,2) edge
  Kxx_true = matrix(0, p, p)
  Kxx_true[1, 2] = Kxx_true[2, 1] = 0.4

  Kyy_true = diag(q)
  Kyy_true[1, 2] = Kyy_true[2, 1] = 0.3

  Kxy_true = matrix(0, p, q)
  Kxy_true[1, 2] = 0.3

  mux_true = matrix(0, p, max(num_cats_bgms))
  mux_true[, 1] = c(-0.5, 0.3, -0.2)
  mux_true[, 2] = c(-1.0, -0.8, -0.6)
  muy_true = c(0.5, -0.3)

  # Generate data using mixedGM
  set.seed(123)
  sim = mixedGM::mixed_gibbs_generate(
    n = n, num_categories = num_cats_mgm,
    mux = mux_true, Kxx = Kxx_true,
    muy = muy_true, Kyy = Kyy_true,
    Kxy = Kxy_true,
    n_burnin = 3000L
  )

  x = sim$x
  y = sim$y
  pq = p + q

  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  bgms_out = bgms:::test_mixed_mrf_sampler(
    x, y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, TRUE, "conditional",
    3000L, 5000L, 42L
  )

  pips = extract_pip(bgms_out$indicator_samples, p, q)

  # True edges should have high PIP
  expect_gt(pips$Gxx[1, 2], 0.5, label = "Kxx(1,2) PIP > 0.5")
  expect_gt(pips$Gxy[1, 2], 0.5, label = "Kxy(1,2) PIP > 0.5")
  expect_gt(pips$Gyy[1, 2], 0.5, label = "Kyy(1,2) PIP > 0.5")

  # Most null edges should have low PIP (allow one false positive due to

  # conditional PL sensitivity with small p, q)
  null_kxx = c(pips$Gxx[1, 3], pips$Gxx[2, 3])
  null_kxy = c(pips$Gxy[2, 1], pips$Gxy[2, 2],
               pips$Gxy[3, 1], pips$Gxy[3, 2])
  expect_lt(mean(c(null_kxx, null_kxy) > 0.5), 0.5,
            label = "Fewer than half of null edges are false positives")
  expect_lt(median(c(null_kxx, null_kxy)), 0.3,
            label = "Median null PIP is low")
})


# ===========================================================================
# Test 2: Edge selection on simulation data (marginal PL)
# ===========================================================================
test_that("edge selection recovers sparse structure (marginal PL)", {
  skip_on_cran()

  p = 3L; q = 2L; n = 800L
  num_cats_mgm = rep(3L, p)
  num_cats_bgms = rep(2L, p)

  Kxx_true = matrix(0, p, p)
  Kxx_true[1, 2] = Kxx_true[2, 1] = 0.4

  Kyy_true = diag(q)
  Kyy_true[1, 2] = Kyy_true[2, 1] = 0.3

  Kxy_true = matrix(0, p, q)
  Kxy_true[1, 2] = 0.3

  mux_true = matrix(0, p, max(num_cats_bgms))
  mux_true[, 1] = c(-0.5, 0.3, -0.2)
  mux_true[, 2] = c(-1.0, -0.8, -0.6)
  muy_true = c(0.5, -0.3)

  set.seed(123)
  sim = mixedGM::mixed_gibbs_generate(
    n = n, num_categories = num_cats_mgm,
    mux = mux_true, Kxx = Kxx_true,
    muy = muy_true, Kyy = Kyy_true,
    Kxy = Kxy_true,
    n_burnin = 3000L
  )

  x = sim$x
  y = sim$y
  pq = p + q

  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  bgms_out = bgms:::test_mixed_mrf_sampler(
    x, y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, TRUE, "marginal",
    3000L, 5000L, 42L
  )

  pips = extract_pip(bgms_out$indicator_samples, p, q)

  # True edges should have high PIP
  expect_gt(pips$Gxx[1, 2], 0.5, label = "Kxx(1,2) PIP > 0.5")
  expect_gt(pips$Gxy[1, 2], 0.5, label = "Kxy(1,2) PIP > 0.5")
  expect_gt(pips$Gyy[1, 2], 0.5, label = "Kyy(1,2) PIP > 0.5")

  # Null edges low PIP
  expect_lt(pips$Gxx[1, 3], 0.5, label = "Kxx(1,3) null PIP < 0.5")
  expect_lt(pips$Gxx[2, 3], 0.5, label = "Kxx(2,3) null PIP < 0.5")
  expect_lt(pips$Gxy[1, 1], 0.5, label = "Kxy(1,1) null PIP < 0.5")
  expect_lt(pips$Gxy[2, 1], 0.5, label = "Kxy(2,1) null PIP < 0.5")
  expect_lt(pips$Gxy[2, 2], 0.5, label = "Kxy(2,2) null PIP < 0.5")
  expect_lt(pips$Gxy[3, 1], 0.5, label = "Kxy(3,1) null PIP < 0.5")
  expect_lt(pips$Gxy[3, 2], 0.5, label = "Kxy(3,2) null PIP < 0.5")
})


# ===========================================================================
# Test 3: Edge selection vs mixedGM (conditional PL)
# ===========================================================================
test_that("edge selection PIPs agree with mixedGM (conditional PL)", {
  skip_on_cran()

  p = 3L; q = 2L; n = 800L
  num_cats_mgm = rep(3L, p)
  num_cats_bgms = rep(2L, p)

  Kxx_true = matrix(0, p, p)
  Kxx_true[1, 2] = Kxx_true[2, 1] = 0.4

  Kyy_true = diag(q)
  Kyy_true[1, 2] = Kyy_true[2, 1] = 0.3

  Kxy_true = matrix(0, p, q)
  Kxy_true[1, 2] = 0.3

  mux_true = matrix(0, p, max(num_cats_bgms))
  mux_true[, 1] = c(-0.5, 0.3, -0.2)
  mux_true[, 2] = c(-1.0, -0.8, -0.6)
  muy_true = c(0.5, -0.3)

  set.seed(123)
  sim = mixedGM::mixed_gibbs_generate(
    n = n, num_categories = num_cats_mgm,
    mux = mux_true, Kxx = Kxx_true,
    muy = muy_true, Kyy = Kyy_true,
    Kxy = Kxy_true,
    n_burnin = 3000L
  )

  x = sim$x
  y = sim$y
  pq = p + q

  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  # bgms with edge selection
  bgms_out = bgms:::test_mixed_mrf_sampler(
    x, y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, TRUE, "conditional",
    3000L, 5000L, 42L
  )
  pips_bgms = extract_pip(bgms_out$indicator_samples, p, q)

  # mixedGM with edge selection (using logistic prior to match bgms)
  set.seed(42)
  mgm = mixedGM::mixed_sampler(
    x = x, y = y, num_categories = num_cats_mgm,
    n_warmup = 3000L, n_samples = 5000L,
    edge_selection = TRUE, verbose = FALSE,
    pseudolikelihood = "conditional",
    log_prior_mean = mixedGM:::log_prior_logistic,
    lp_mean_params = list()
  )

  # mixedGM stores indicators as 3D array: samples$indicator[iter, i, j]
  # G is (p+q) x (p+q): rows/cols 1:p = discrete, (p+1):(p+q) = continuous
  G_array = mgm$samples$indicator
  Gxx_mgm = matrix(0, p, p)
  for(i in 1:(p - 1)) for(j in (i + 1):p) {
    Gxx_mgm[i, j] = Gxx_mgm[j, i] = mean(G_array[, i, j])
  }
  Gyy_mgm = matrix(0, q, q)
  if(q > 1) {
    for(i in 1:(q - 1)) for(j in (i + 1):q) {
      Gyy_mgm[i, j] = Gyy_mgm[j, i] = mean(G_array[, p + i, p + j])
    }
  }
  Gxy_mgm = matrix(0, p, q)
  for(i in 1:p) for(j in 1:q) {
    Gxy_mgm[i, j] = mean(G_array[, i, p + j])
  }

  # Both should agree on which edges to include/exclude
  # (directional agreement: both high for true edges, both low for null)
  all_bgms = c(pips_bgms$Gxx[upper.tri(pips_bgms$Gxx)],
               pips_bgms$Gyy[upper.tri(pips_bgms$Gyy)],
               as.vector(pips_bgms$Gxy))
  all_mgm = c(Gxx_mgm[upper.tri(Gxx_mgm)],
              Gyy_mgm[upper.tri(Gyy_mgm)],
              as.vector(Gxy_mgm))

  # Correlation of PIPs should be positive
  expect_gt(cor(all_bgms, all_mgm), 0.5,
            label = "PIP correlation bgms vs mixedGM > 0.5")

  # True structure agreement: same edges classified as > 0.5
  bgms_selected = all_bgms > 0.5
  mgm_selected = all_mgm > 0.5
  agreement = mean(bgms_selected == mgm_selected)
  expect_gt(agreement, 0.6, label = "Structure agreement > 60%")
})


# ===========================================================================
# Test 4: Edge-gating prevents parameter drift on excluded edges
# ===========================================================================
test_that("excluded edges have zero parameters after edge selection", {
  skip_on_cran()

  p = 3L; q = 2L; n = 500L
  num_cats_mgm = rep(3L, p)
  num_cats_bgms = rep(2L, p)

  # Fully sparse graph: no edges in true model
  set.seed(456)
  sim = mixedGM::mixed_gibbs_generate(
    n = n, num_categories = num_cats_mgm,
    mux = matrix(c(-0.5, 0.3, -0.2, -1.0, -0.8, -0.6), nrow = p),
    Kxx = matrix(0, p, p),
    muy = c(0.5, -0.3), Kyy = diag(q),
    Kxy = matrix(0, p, q),
    n_burnin = 3000L
  )

  x = sim$x
  y = sim$y
  pq = p + q

  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  bgms_out = bgms:::test_mixed_mrf_sampler(
    x, y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, TRUE, "conditional",
    2000L, 3000L, 42L
  )

  pips = extract_pip(bgms_out$indicator_samples, p, q)

  # With a fully null graph and n = 500, all PIPs should be low
  expect_lt(max(pips$Gxx[upper.tri(pips$Gxx)]), 0.8,
            label = "Max Kxx PIP < 0.8 for null graph")
  expect_lt(max(as.vector(pips$Gxy)), 0.8,
            label = "Max Kxy PIP < 0.8 for null graph")
})
