# =============================================================================
# test-mixed-mrf-likelihoods.R — Phase B.1 likelihood tests
# =============================================================================
# Tests for log_conditional_omrf() and log_conditional_ggm() in the
# MixedMRFModel, validated against pure-R reference implementations.

# ------------------------------------------------------------------
# R reference: conditional OMRF pseudolikelihood for variable s
# ------------------------------------------------------------------
# Computes log f(x_s | x_{-s}, y) summed over n observations.
#
# @param x_int        n x p integer matrix (0-based categories)
# @param y            n x q continuous matrix
# @param Kxx          p x p symmetric, zero diagonal
# @param Kxy          p x q cross-interactions
# @param mux          p x max_cats  thresholds / BC coefficients
# @param num_cats     p-vector of category counts
# @param is_ordinal   p-vector (1 = ordinal, 0 = BC)
# @param baseline_cat p-vector of reference categories
# @param s            0-based variable index
#
# Returns: scalar log-likelihood.
ref_log_conditional_omrf = function(x_int, y, Kxx, Kxy, mux, num_cats,
                                    is_ordinal, baseline_cat, s) {
  n = nrow(x_int)
  s1 = s + 1L  # R indexing

  # Center BC observations (same as constructor)
  x = x_int
  for(j in seq_len(ncol(x))) {
    if(!is_ordinal[j]) x[, j] = x[, j] - baseline_cat[j]
  }
  x_dbl = matrix(as.double(x), nrow = n)

  # Rest score
  rest = x_dbl %*% Kxx[, s1] - x_dbl[, s1] * Kxx[s1, s1] +
    2.0 * y %*% t(Kxy[s1, , drop = FALSE])
  rest = as.numeric(rest)

  C_s = num_cats[s1]

  if(is_ordinal[s1]) {
    # Numerator: dot(x_s, rest) + sum_c count_c * mux(s, c)
    numer = sum(x_dbl[, s1] * rest)
    for(c in seq_len(C_s)) {
      count_c = sum(x[, s1] == c)
      numer = numer + count_c * mux[s1, c]
    }
    # Log-normalizer per person (log-sum-exp stabilized)
    log_Z = numeric(n)
    for(v in seq_len(n)) {
      terms = numeric(C_s + 1)
      terms[1] = 0  # reference category
      for(c in seq_len(C_s)) {
        terms[c + 1] = mux[s1, c] + c * rest[v]
      }
      mx = max(terms)
      log_Z[v] = mx + log(sum(exp(terms - mx)))
    }
    return(numer - sum(log_Z))
  } else {
    # Blume-Capel
    alpha = mux[s1, 1]
    beta = mux[s1, 2]
    ref = baseline_cat[s1]

    numer = sum(x_dbl[, s1] * rest) +
      alpha * sum(x[, s1]) + beta * sum(x[, s1]^2)

    log_Z = numeric(n)
    for(v in seq_len(n)) {
      cats = 0:C_s
      centered = cats - ref
      theta = alpha * centered + beta * centered^2
      terms = theta + centered * rest[v]
      mx = max(terms)
      log_Z[v] = mx + log(sum(exp(terms - mx)))
    }
    return(numer - sum(log_Z))
  }
}

# ------------------------------------------------------------------
# R reference: conditional GGM log-likelihood
# ------------------------------------------------------------------
# Computes log f(y | x) using conditional mean and precision Kyy.
#
# @param x_int        n x p integer matrix (0-based)
# @param y            n x q continuous matrix
# @param Kxy          p x q cross-interactions
# @param Kyy          q x q SPD precision matrix
# @param muy          q-vector of means
# @param is_ordinal   p-vector
# @param baseline_cat p-vector
#
# Returns: scalar log-likelihood.
ref_log_conditional_ggm = function(x_int, y, Kxy, Kyy, muy,
                                   is_ordinal, baseline_cat) {
  n = nrow(y)
  q = ncol(y)

  # Center BC observations
  x = x_int
  for(j in seq_len(ncol(x))) {
    if(!is_ordinal[j]) x[, j] = x[, j] - baseline_cat[j]
  }
  x_dbl = matrix(as.double(x), nrow = n)

  Kyy_inv = solve(Kyy)
  logdet_val = as.numeric(determinant(Kyy, logarithm = TRUE)$modulus)

  # Conditional mean: repmat(muy', n, 1) + 2 * x_dbl * Kxy * Kyy_inv
  cond_mean = matrix(rep(muy, each = n), nrow = n) +
    2.0 * x_dbl %*% Kxy %*% Kyy_inv

  D = y - cond_mean
  quad_sum = sum((D %*% Kyy) * D)

  ll = n / 2.0 * (-q * log(2 * pi) + logdet_val) - quad_sum / 2.0
  return(ll)
}


# ------------------------------------------------------------------
# Helper to build a zero parameter vector with correct structure
# ------------------------------------------------------------------
build_param_vec = function(p, q, num_cats, is_ordinal) {
  num_main = 0L
  for(s in seq_len(p)) {
    if(is_ordinal[s]) {
      num_main = num_main + num_cats[s]
    } else {
      num_main = num_main + 2L
    }
  }
  num_pairwise_xx = p * (p - 1L) / 2L
  num_kyy = q * (q + 1L) / 2L
  num_cross = p * q
  total = num_main + num_pairwise_xx + q + num_kyy + num_cross
  list(
    total = total,
    num_main = num_main,
    num_pairwise_xx = num_pairwise_xx,
    q = q,
    num_kyy = num_kyy,
    num_cross = num_cross,
    kyy_start = num_main + num_pairwise_xx + q  # 0-based start of Kyy block
  )
}


# ------------------------------------------------------------------
# Test setup shared across tests
# ------------------------------------------------------------------
make_test_data = function() {
  set.seed(42)
  n = 10L; p = 3L; q = 2L
  num_cats = c(2L, 3L, 2L)
  is_ordinal = c(1L, 1L, 0L)     # vars 1,2 ordinal; var 3 BC
  baseline_cat = c(0L, 0L, 1L)   # BC var has ref = 1

  x = matrix(0L, n, p)
  x[, 1] = sample(0:2, n, replace = TRUE)
  x[, 2] = sample(0:3, n, replace = TRUE)
  x[, 3] = sample(0:2, n, replace = TRUE)  # original BC obs in [0, C_s]
  y = matrix(rnorm(n * q), n, q)

  pq = p + q
  inc_prob = matrix(0.5, pq, pq)
  edge_ind = matrix(1L, pq, pq)
  diag(edge_ind) = 0L

  list(
    n = n, p = p, q = q,
    x = x, y = y,
    num_cats = num_cats,
    is_ordinal = as.integer(is_ordinal),
    baseline_cat = baseline_cat,
    inc_prob = inc_prob,
    edge_ind = edge_ind
  )
}


# ==============================================================================
# Tests
# ==============================================================================

test_that("log_conditional_omrf at zero parameters equals -n*log(C+1)", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)
  # Kyy = identity
  idx = info$kyy_start
  for(i in seq_len(d$q)) {
    for(j in i:d$q) {
      idx = idx + 1L
      if(i == j) params[idx] = 1.0
    }
  }

  res = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "conditional", params, 1L
  )

  # Ordinal vars: expected = -n * log(C_s + 1)
  expect_equal(res$omrf_ll[1], -d$n * log(d$num_cats[1] + 1), tolerance = 1e-10)
  expect_equal(res$omrf_ll[2], -d$n * log(d$num_cats[2] + 1), tolerance = 1e-10)

  # BC var: at zero params, all categories have theta = 0, rest = 0
  # so Z = num_cats + 1, numer = 0 => ll = -n * log(C+1)
  expect_equal(res$omrf_ll[3], -d$n * log(d$num_cats[3] + 1), tolerance = 1e-10)
})


test_that("log_conditional_ggm at zero parameters equals standard normal", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)
  idx = info$kyy_start
  for(i in seq_len(d$q)) {
    for(j in i:d$q) {
      idx = idx + 1L
      if(i == j) params[idx] = 1.0
    }
  }

  res = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "conditional", params, 1L
  )

  # At zero params: y ~ N(0, I), so ll = sum(dnorm(y, log = TRUE))
  expected = sum(dnorm(d$y, log = TRUE))
  expect_equal(res$ggm_ll, expected, tolerance = 1e-10)
})


test_that("log_conditional_omrf matches R reference (ordinal, nonzero params)", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)

  # Set some nonzero main effects
  params[1] = 0.3; params[2] = -0.2           # var 1 thresholds (C=2)
  params[3] = 0.1; params[4] = 0.4; params[5] = -0.1  # var 2 thresholds (C=3)
  params[6] = 0.5; params[7] = -0.3           # var 3 BC alpha, beta

  # Set some Kxx edges (upper-tri, row-major: (0,1), (0,2), (1,2))
  kxx_start = info$num_main
  params[kxx_start + 1] = 0.2    # Kxx(0,1)
  params[kxx_start + 2] = -0.1   # Kxx(0,2)
  params[kxx_start + 3] = 0.15   # Kxx(1,2)

  # muy
  muy_start = info$num_main + info$num_pairwise_xx
  params[muy_start + 1] = 0.1
  params[muy_start + 2] = -0.2

  # Kyy = SPD (2x2 with positive eigenvalues)
  kyy_start = info$kyy_start
  params[kyy_start + 1] = 2.0    # Kyy(0,0)
  params[kyy_start + 2] = 0.3    # Kyy(0,1)
  params[kyy_start + 3] = 1.5    # Kyy(1,1)

  # Kxy (p x q, row-major)
  kxy_start = kyy_start + info$num_kyy
  params[kxy_start + 1] = 0.1    # Kxy(0,0)
  params[kxy_start + 2] = -0.05  # Kxy(0,1)
  params[kxy_start + 3] = 0.08   # Kxy(1,0)
  params[kxy_start + 4] = 0.12   # Kxy(1,1)
  params[kxy_start + 5] = -0.07  # Kxy(2,0)
  params[kxy_start + 6] = 0.04   # Kxy(2,1)

  res = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "conditional", params, 1L
  )

  # Reconstruct parameter matrices for R reference
  mux = matrix(0, d$p, max(d$num_cats))
  mux[1, 1:2] = c(0.3, -0.2)
  mux[2, 1:3] = c(0.1, 0.4, -0.1)
  mux[3, 1:2] = c(0.5, -0.3)

  Kxx = matrix(0, d$p, d$p)
  Kxx[1, 2] = Kxx[2, 1] = 0.2
  Kxx[1, 3] = Kxx[3, 1] = -0.1
  Kxx[2, 3] = Kxx[3, 2] = 0.15

  Kxy = matrix(c(0.1, 0.08, -0.07, -0.05, 0.12, 0.04), nrow = d$p, ncol = d$q)
  muy = c(0.1, -0.2)

  # Check each ordinal variable against R reference
  for(s in 0:1) {
    expected = ref_log_conditional_omrf(
      d$x, d$y, Kxx, Kxy, mux, d$num_cats,
      d$is_ordinal, d$baseline_cat, s
    )
    expect_equal(res$omrf_ll[s + 1], expected, tolerance = 1e-8,
                 label = paste0("ordinal var ", s))
  }
})


test_that("log_conditional_omrf matches R reference (Blume-Capel, nonzero params)", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)

  # Same params as ordinal test
  params[1] = 0.3; params[2] = -0.2
  params[3] = 0.1; params[4] = 0.4; params[5] = -0.1
  params[6] = 0.5; params[7] = -0.3

  kxx_start = info$num_main
  params[kxx_start + 1] = 0.2
  params[kxx_start + 2] = -0.1
  params[kxx_start + 3] = 0.15

  muy_start = info$num_main + info$num_pairwise_xx
  params[muy_start + 1] = 0.1
  params[muy_start + 2] = -0.2

  kyy_start = info$kyy_start
  params[kyy_start + 1] = 2.0
  params[kyy_start + 2] = 0.3
  params[kyy_start + 3] = 1.5

  kxy_start = kyy_start + info$num_kyy
  params[kxy_start + 1] = 0.1
  params[kxy_start + 2] = -0.05
  params[kxy_start + 3] = 0.08
  params[kxy_start + 4] = 0.12
  params[kxy_start + 5] = -0.07
  params[kxy_start + 6] = 0.04

  res = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "conditional", params, 1L
  )

  mux = matrix(0, d$p, max(d$num_cats))
  mux[1, 1:2] = c(0.3, -0.2)
  mux[2, 1:3] = c(0.1, 0.4, -0.1)
  mux[3, 1:2] = c(0.5, -0.3)

  Kxx = matrix(0, d$p, d$p)
  Kxx[1, 2] = Kxx[2, 1] = 0.2
  Kxx[1, 3] = Kxx[3, 1] = -0.1
  Kxx[2, 3] = Kxx[3, 2] = 0.15

  Kxy = matrix(c(0.1, 0.08, -0.07, -0.05, 0.12, 0.04), nrow = d$p, ncol = d$q)

  # Test BC variable (s = 2, 0-based)
  expected = ref_log_conditional_omrf(
    d$x, d$y, Kxx, Kxy, mux, d$num_cats,
    d$is_ordinal, d$baseline_cat, s = 2
  )
  expect_equal(res$omrf_ll[3], expected, tolerance = 1e-8)
})


test_that("log_conditional_ggm matches R reference (nonzero params)", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)

  # Set same params
  params[1] = 0.3; params[2] = -0.2
  params[3] = 0.1; params[4] = 0.4; params[5] = -0.1
  params[6] = 0.5; params[7] = -0.3

  kxx_start = info$num_main
  params[kxx_start + 1] = 0.2
  params[kxx_start + 2] = -0.1
  params[kxx_start + 3] = 0.15

  muy_start = info$num_main + info$num_pairwise_xx
  params[muy_start + 1] = 0.1
  params[muy_start + 2] = -0.2

  kyy_start = info$kyy_start
  params[kyy_start + 1] = 2.0
  params[kyy_start + 2] = 0.3
  params[kyy_start + 3] = 1.5

  kxy_start = kyy_start + info$num_kyy
  params[kxy_start + 1] = 0.1
  params[kxy_start + 2] = -0.05
  params[kxy_start + 3] = 0.08
  params[kxy_start + 4] = 0.12
  params[kxy_start + 5] = -0.07
  params[kxy_start + 6] = 0.04

  res = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "conditional", params, 1L
  )

  Kyy = matrix(c(2.0, 0.3, 0.3, 1.5), 2, 2)
  Kxy = matrix(c(0.1, 0.08, -0.07, -0.05, 0.12, 0.04), nrow = d$p, ncol = d$q)
  muy = c(0.1, -0.2)

  expected = ref_log_conditional_ggm(
    d$x, d$y, Kxy, Kyy, muy, d$is_ordinal, d$baseline_cat
  )
  expect_equal(res$ggm_ll, expected, tolerance = 1e-8)
})


test_that("log_conditional_ggm invariant to Kxx changes", {
  # GGM likelihood depends on Kyy, Kxy, muy but not Kxx
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params1 = rep(0, info$total)
  params2 = rep(0, info$total)

  # Both: same Kyy = SPD, same Kxy, same muy
  kyy_start = info$kyy_start
  for(params in list(params1, params2)) {
    params[kyy_start + 1] = 2.0
    params[kyy_start + 2] = 0.3
    params[kyy_start + 3] = 1.5
  }
  # But different Kxx
  kxx_start = info$num_main
  params1[kxx_start + 1] = 0.5
  params2[kxx_start + 1] = -0.3

  # Need to re-assign since list() copies
  params1[kyy_start + 1] = 2.0
  params1[kyy_start + 2] = 0.3
  params1[kyy_start + 3] = 1.5
  params2[kyy_start + 1] = 2.0
  params2[kyy_start + 2] = 0.3
  params2[kyy_start + 3] = 1.5

  res1 = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "conditional", params1, 1L
  )
  res2 = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "conditional", params2, 1L
  )

  # GGM likelihood should be identical regardless of Kxx
  expect_equal(res1$ggm_ll, res2$ggm_ll, tolerance = 1e-12)
})


test_that("likelihood with p=1 ordinal, q=1 continuous works", {
  set.seed(99)
  n = 5L; p = 1L; q = 1L
  x = matrix(sample(0:1, n, replace = TRUE), n, p)
  y = matrix(rnorm(n), n, q)
  num_cats = 1L
  is_ordinal = 1L
  baseline_cat = 0L
  pq = p + q
  inc_prob = matrix(0.5, pq, pq)
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  # Params: 1 threshold + 0 Kxx + 1 muy + 1 Kyy + 1 Kxy = 4
  params = c(0.5,    # mux(0,0) = threshold
             0.2,    # muy
             1.5,    # Kyy(0,0)
             0.1)    # Kxy(0,0)

  res = bgms:::test_mixed_mrf_likelihoods(
    x, y, num_cats, is_ordinal, baseline_cat,
    inc_prob, edge_ind, FALSE, "conditional", params, 1L
  )

  # R reference
  mux = matrix(0.5, 1, 1)
  Kxx = matrix(0, 1, 1)
  Kxy = matrix(0.1, 1, 1)
  Kyy = matrix(1.5, 1, 1)
  muy = 0.2

  expected_omrf = ref_log_conditional_omrf(
    x, y, Kxx, Kxy, mux, num_cats, is_ordinal, baseline_cat, s = 0
  )
  expected_ggm = ref_log_conditional_ggm(
    x, y, Kxy, Kyy, muy, is_ordinal, baseline_cat
  )

  expect_equal(res$omrf_ll[1], expected_omrf, tolerance = 1e-10)
  expect_equal(res$ggm_ll, expected_ggm, tolerance = 1e-10)
})


test_that("all-Blume-Capel discrete variables work", {
  set.seed(77)
  n = 8L; p = 2L; q = 1L
  num_cats = c(3L, 4L)
  is_ordinal = c(0L, 0L)
  baseline_cat = c(1L, 2L)

  x = matrix(0L, n, p)
  x[, 1] = sample(0:3, n, replace = TRUE)
  x[, 2] = sample(0:4, n, replace = TRUE)
  y = matrix(rnorm(n), n, q)

  pq = p + q
  inc_prob = matrix(0.5, pq, pq)
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  # Params: 2 BC + 2 BC + 1 Kxx + 1 muy + 1 Kyy + 2 Kxy = 9
  # num_main = 4, num_pairwise_xx = 1, q = 1, num_kyy = 1, num_cross = 2
  params = c(0.3, -0.1,     # BC var 0: alpha, beta
             0.2, -0.05,    # BC var 1: alpha, beta
             0.15,           # Kxx(0,1)
             -0.1,           # muy
             2.0,            # Kyy
             0.08, -0.06)   # Kxy

  res = bgms:::test_mixed_mrf_likelihoods(
    x, y, num_cats, is_ordinal, baseline_cat,
    inc_prob, edge_ind, FALSE, "conditional", params, 1L
  )

  mux = matrix(0, 2, 4)
  mux[1, 1:2] = c(0.3, -0.1)
  mux[2, 1:2] = c(0.2, -0.05)

  Kxx = matrix(0, 2, 2)
  Kxx[1, 2] = Kxx[2, 1] = 0.15

  Kxy = matrix(c(0.08, -0.06), nrow = 2, ncol = 1)
  Kyy = matrix(2.0, 1, 1)
  muy = -0.1

  for(s in 0:1) {
    expected = ref_log_conditional_omrf(
      x, y, Kxx, Kxy, mux, num_cats, is_ordinal, baseline_cat, s
    )
    expect_equal(res$omrf_ll[s + 1], expected, tolerance = 1e-8,
                 label = paste0("BC var ", s))
  }

  expected_ggm = ref_log_conditional_ggm(
    x, y, Kxy, Kyy, muy, is_ordinal, baseline_cat
  )
  expect_equal(res$ggm_ll, expected_ggm, tolerance = 1e-8)
})
