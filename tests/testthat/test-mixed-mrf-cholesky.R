# =============================================================================
# test-mixed-mrf-cholesky.R — Phase B+ Cholesky rank-1 update tests
# =============================================================================
# T28:  log_ggm_ratio_edge  matches brute-force log-likelihood difference
# T28b: log_ggm_ratio_diag  matches brute-force log-likelihood difference
# T29:  Cholesky / covariance / logdet fidelity after rank-2 edge update
# T29b: Cholesky / covariance / logdet fidelity after rank-1 diagonal update

# ------------------------------------------------------------------
# Helper: build parameter vector with nonzero Kyy (must be SPD)
# ------------------------------------------------------------------
make_cholesky_test = function(q = 3L) {
  set.seed(42)
  n = 20L; p = 2L
  num_cats = c(2L, 3L)
  is_ordinal = c(1L, 1L)
  baseline_cat = c(0L, 0L)

  x = matrix(0L, n, p)
  x[, 1] = sample(0:2, n, replace = TRUE)
  x[, 2] = sample(0:3, n, replace = TRUE)
  y = matrix(rnorm(n * q), n, q)

  pq = p + q
  inc_prob = matrix(0.5, pq, pq)
  edge_ind = matrix(1L, pq, pq)
  diag(edge_ind) = 0L

  # Param layout: main | Kxx | muy | Kyy (upper-tri) | Kxy
  num_main = sum(num_cats)     # 2 + 3 = 5
  num_kxx = p * (p - 1L) / 2L # 1
  num_kyy = q * (q + 1L) / 2L
  num_kxy = p * q
  total = num_main + num_kxx + q + num_kyy + num_kxy
  params = rep(0, total)

  # Thresholds
  params[1:2] = c(0.3, -0.2)
  params[3:5] = c(0.1, 0.4, -0.1)

  # Kxx(0,1)
  params[num_main + 1] = 0.2

  # muy
  muy_start = num_main + num_kxx
  params[muy_start + seq_len(q)] = rnorm(q, sd = 0.2)

  # Kyy: build SPD Q = L L^T with random L, then pack upper-tri
  L = matrix(0, q, q)
  diag(L) = runif(q, 1.0, 2.0)
  L[lower.tri(L)] = rnorm(q * (q - 1) / 2, sd = 0.3)
  Kyy = L %*% t(L)
  kyy_start = muy_start + q
  idx = 0L
  for(i in seq_len(q)) {
    for(j in i:q) {
      idx = idx + 1L
      params[kyy_start + idx] = Kyy[i, j]
    }
  }

  # Kxy
  kxy_start = kyy_start + num_kyy
  params[kxy_start + seq_len(num_kxy)] = rnorm(num_kxy, sd = 0.1)

  list(
    x = x, y = y,
    num_cats = num_cats,
    is_ordinal = as.integer(is_ordinal),
    baseline_cat = baseline_cat,
    inc_prob = inc_prob,
    edge_ind = edge_ind,
    params = params,
    q = q
  )
}


# ==============================================================================
# T28 + T28b: log-ratio agreement (off-diagonal and diagonal)
# ==============================================================================

test_that("log_ggm_ratio_edge matches brute-force (T28)", {
  d = make_cholesky_test(q = 3L)
  res = bgms:::test_mixed_mrf_cholesky(
    d$x, d$y, d$num_cats, d$is_ordinal, d$baseline_cat,
    d$inc_prob, d$edge_ind, d$params,
    seed = 1L, target_i = 0L, target_j = 1L
  )
  expect_equal(res$ratio_rank2, res$ratio_bruteforce, tolerance = 1e-10)
})

test_that("log_ggm_ratio_diag matches brute-force (T28b)", {
  d = make_cholesky_test(q = 3L)
  res = bgms:::test_mixed_mrf_cholesky(
    d$x, d$y, d$num_cats, d$is_ordinal, d$baseline_cat,
    d$inc_prob, d$edge_ind, d$params,
    seed = 1L, target_i = 0L, target_j = 1L
  )
  expect_equal(res$ratio_diag_rank1, res$ratio_diag_brute, tolerance = 1e-10)
})


# ==============================================================================
# T29 + T29b: Cholesky fidelity after rank-1/rank-2 updates
# ==============================================================================

test_that("Cholesky update after edge matches full recompute (T29)", {
  d = make_cholesky_test(q = 3L)
  res = bgms:::test_mixed_mrf_cholesky(
    d$x, d$y, d$num_cats, d$is_ordinal, d$baseline_cat,
    d$inc_prob, d$edge_ind, d$params,
    seed = 1L, target_i = 0L, target_j = 1L
  )
  expect_lt(res$chol_max_diff, 1e-10)
  expect_lt(res$cov_max_diff, 1e-10)
  expect_lt(res$logdet_diff, 1e-10)
})

test_that("Cholesky update after diagonal matches full recompute (T29b)", {
  d = make_cholesky_test(q = 3L)
  res = bgms:::test_mixed_mrf_cholesky(
    d$x, d$y, d$num_cats, d$is_ordinal, d$baseline_cat,
    d$inc_prob, d$edge_ind, d$params,
    seed = 1L, target_i = 0L, target_j = 1L
  )
  expect_lt(res$chol_diag_max_diff, 1e-10)
  expect_lt(res$cov_diag_max_diff, 1e-10)
  expect_lt(res$logdet_diag_diff, 1e-10)
})


# ==============================================================================
# Same tests with different target indices (i=1, j=2) and q=4
# ==============================================================================

test_that("log_ggm_ratio_edge matches brute-force (T28, q=4, ij=1,2)", {
  d = make_cholesky_test(q = 4L)
  res = bgms:::test_mixed_mrf_cholesky(
    d$x, d$y, d$num_cats, d$is_ordinal, d$baseline_cat,
    d$inc_prob, d$edge_ind, d$params,
    seed = 7L, target_i = 1L, target_j = 2L
  )
  expect_equal(res$ratio_rank2, res$ratio_bruteforce, tolerance = 1e-10)
  expect_equal(res$ratio_diag_rank1, res$ratio_diag_brute, tolerance = 1e-10)
})

test_that("Cholesky fidelity (T29/T29b, q=4, ij=1,2)", {
  d = make_cholesky_test(q = 4L)
  res = bgms:::test_mixed_mrf_cholesky(
    d$x, d$y, d$num_cats, d$is_ordinal, d$baseline_cat,
    d$inc_prob, d$edge_ind, d$params,
    seed = 7L, target_i = 1L, target_j = 2L
  )
  expect_lt(res$chol_max_diff, 1e-10)
  expect_lt(res$cov_max_diff, 1e-10)
  expect_lt(res$logdet_diff, 1e-10)
  expect_lt(res$chol_diag_max_diff, 1e-10)
  expect_lt(res$cov_diag_max_diff, 1e-10)
  expect_lt(res$logdet_diag_diff, 1e-10)
})


# ==============================================================================
# Test with q=2 (minimum for off-diagonal) — boundary case
# ==============================================================================

test_that("Rank-1 Cholesky correctness with q=2 boundary (T28-T29b)", {
  d = make_cholesky_test(q = 2L)
  res = bgms:::test_mixed_mrf_cholesky(
    d$x, d$y, d$num_cats, d$is_ordinal, d$baseline_cat,
    d$inc_prob, d$edge_ind, d$params,
    seed = 99L, target_i = 0L, target_j = 1L
  )
  expect_equal(res$ratio_rank2, res$ratio_bruteforce, tolerance = 1e-10)
  expect_equal(res$ratio_diag_rank1, res$ratio_diag_brute, tolerance = 1e-10)
  expect_lt(res$chol_max_diff, 1e-10)
  expect_lt(res$cov_max_diff, 1e-10)
  expect_lt(res$logdet_diff, 1e-10)
  expect_lt(res$chol_diag_max_diff, 1e-10)
  expect_lt(res$cov_diag_max_diff, 1e-10)
  expect_lt(res$logdet_diag_diff, 1e-10)
})
