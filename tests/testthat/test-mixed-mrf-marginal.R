# =============================================================================
# test-mixed-mrf-marginal.R — Phase C marginal pseudo-likelihood tests
# =============================================================================
# Tests for log_marginal_omrf() and marginal-mode sampler, validated against
# pure-R reference implementations and cross-package comparison with mixedGM.

# ------------------------------------------------------------------
# R reference: marginal OMRF pseudolikelihood for variable s
# ------------------------------------------------------------------
# Computes log f(x_s | x_{-s}) using Theta = Kxx + 2 Kxy Kyy^{-1} Kxy'
#
# @param x_int        n x p integer matrix (0-based categories)
# @param Theta        p x p marginal interaction matrix
# @param Kxy          p x q cross-interactions
# @param mux          p x max_cats thresholds / BC coefficients
# @param muy          q-vector of continuous means
# @param num_cats     p-vector of category counts
# @param is_ordinal   p-vector (1 = ordinal, 0 = BC)
# @param baseline_cat p-vector of reference categories
# @param s            0-based variable index
#
# Returns: scalar log-likelihood.
ref_log_marginal_omrf = function(x_int, Theta, Kxy, mux, muy, num_cats,
                                 is_ordinal, baseline_cat, s) {
  n = nrow(x_int)
  s1 = s + 1L

  # Center BC observations
  x = x_int
  for(j in seq_len(ncol(x))) {
    if(!is_ordinal[j]) x[, j] = x[, j] - baseline_cat[j]
  }
  x_dbl = matrix(as.double(x), nrow = n)

  theta_ss = Theta[s1, s1]

  # Rest score: Theta-based + Kxy*muy bias
  rest = x_dbl %*% Theta[, s1] - x_dbl[, s1] * theta_ss +
    2.0 * sum(Kxy[s1, ] * muy)
  rest = as.numeric(rest)

  C_s = num_cats[s1]

  if(is_ordinal[s1]) {
    # Numerator: dot(x_s, rest) + theta_ss * dot(x_s, x_s) + sum(count_c * mux)
    numer = sum(x_dbl[, s1] * rest) +
      theta_ss * sum(x_dbl[, s1]^2)
    for(c in seq_len(C_s)) {
      count_c = sum(x[, s1] == c)
      numer = numer + count_c * mux[s1, c]
    }

    # Denominator with col_offset = mux + c^2 * theta_ss
    log_Z = numeric(n)
    for(v in seq_len(n)) {
      terms = numeric(C_s + 1)
      terms[1] = 0  # reference category (c=0)
      for(c in seq_len(C_s)) {
        terms[c + 1] = mux[s1, c] + c^2 * theta_ss + c * rest[v]
      }
      mx = max(terms)
      log_Z[v] = mx + log(sum(exp(terms - mx)))
    }
    return(numer - sum(log_Z))
  } else {
    # Blume-Capel: alpha * sum(x) + beta * sum(x^2)
    alpha = mux[s1, 1]
    beta = mux[s1, 2]
    ref = baseline_cat[s1]

    numer = sum(x_dbl[, s1] * rest) +
      theta_ss * sum(x_dbl[, s1]^2) +
      alpha * sum(x[, s1]) + beta * sum(x[, s1]^2)

    # Effective beta includes theta_ss
    eff_beta = beta + theta_ss

    log_Z = numeric(n)
    for(v in seq_len(n)) {
      cats = 0:C_s
      centered = cats - ref
      theta = alpha * centered + eff_beta * centered^2
      terms = theta + centered * rest[v]
      mx = max(terms)
      log_Z[v] = mx + log(sum(exp(terms - mx)))
    }
    return(numer - sum(log_Z))
  }
}


# ------------------------------------------------------------------
# Test data setup (reuse from likelihood tests)
# ------------------------------------------------------------------
make_test_data = function() {
  set.seed(42)
  n = 10L; p = 3L; q = 2L
  num_cats = c(2L, 3L, 2L)
  is_ordinal = c(1L, 1L, 0L)
  baseline_cat = c(0L, 0L, 1L)

  x = matrix(0L, n, p)
  x[, 1] = sample(0:2, n, replace = TRUE)
  x[, 2] = sample(0:3, n, replace = TRUE)
  x[, 3] = sample(0:2, n, replace = TRUE)
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
    kyy_start = num_main + num_pairwise_xx + q
  )
}


# ==============================================================================
# Test: marginal OMRF at zero params equals conditional OMRF
# ==============================================================================

test_that("log_marginal_omrf at zero params equals log_conditional_omrf", {
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
    d$inc_prob, d$edge_ind, FALSE, "marginal", params, 1L
  )

  # With Kxy = 0 and Kyy = I, Theta = Kxx + 0 = Kxx and bias = 0,
  # so marginal OMRF = conditional OMRF
  for(s in seq_len(d$p)) {
    expect_equal(res$marg_omrf_ll[s], res$omrf_ll[s], tolerance = 1e-10,
                 label = paste0("var ", s - 1, " marginal == conditional at zero"))
  }
})


# ==============================================================================
# Test: marginal OMRF matches R reference (ordinal, nonzero params)
# ==============================================================================

test_that("log_marginal_omrf matches R reference (ordinal, nonzero params)", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)

  # Main effects
  params[1] = 0.3; params[2] = -0.2
  params[3] = 0.1; params[4] = 0.4; params[5] = -0.1
  params[6] = 0.5; params[7] = -0.3

  # Kxx edges
  kxx_start = info$num_main
  params[kxx_start + 1] = 0.2
  params[kxx_start + 2] = -0.1
  params[kxx_start + 3] = 0.15

  # muy
  muy_start = info$num_main + info$num_pairwise_xx
  params[muy_start + 1] = 0.1
  params[muy_start + 2] = -0.2

  # Kyy = SPD
  kyy_start = info$kyy_start
  params[kyy_start + 1] = 2.0
  params[kyy_start + 2] = 0.3
  params[kyy_start + 3] = 1.5

  # Kxy
  kxy_start = kyy_start + info$num_kyy
  params[kxy_start + 1] = 0.1
  params[kxy_start + 2] = -0.05
  params[kxy_start + 3] = 0.08
  params[kxy_start + 4] = 0.12
  params[kxy_start + 5] = -0.07
  params[kxy_start + 6] = 0.04

  res = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "marginal", params, 1L
  )

  # Build R reference matrices
  mux = matrix(0, d$p, max(d$num_cats))
  mux[1, 1:2] = c(0.3, -0.2)
  mux[2, 1:3] = c(0.1, 0.4, -0.1)
  mux[3, 1:2] = c(0.5, -0.3)

  Kxx = matrix(0, d$p, d$p)
  Kxx[1, 2] = Kxx[2, 1] = 0.2
  Kxx[1, 3] = Kxx[3, 1] = -0.1
  Kxx[2, 3] = Kxx[3, 2] = 0.15

  Kxy = matrix(c(0.1, 0.08, -0.07, -0.05, 0.12, 0.04), nrow = d$p, ncol = d$q)
  Kyy = matrix(c(2.0, 0.3, 0.3, 1.5), d$q, d$q)
  muy = c(0.1, -0.2)

  # Theta = Kxx + 2 Kxy Kyy^{-1} Kxy'
  Theta = Kxx + 2 * Kxy %*% solve(Kyy) %*% t(Kxy)

  # Check ordinal variables
  for(s in 0:1) {
    expected = ref_log_marginal_omrf(
      d$x, Theta, Kxy, mux, muy, d$num_cats,
      d$is_ordinal, d$baseline_cat, s
    )
    expect_equal(res$marg_omrf_ll[s + 1], expected, tolerance = 1e-8,
                 label = paste0("ordinal var ", s))
  }
})


# ==============================================================================
# Test: marginal OMRF matches R reference (Blume-Capel)
# ==============================================================================

test_that("log_marginal_omrf matches R reference (Blume-Capel, nonzero params)", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)

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
    d$inc_prob, d$edge_ind, FALSE, "marginal", params, 1L
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
  Kyy = matrix(c(2.0, 0.3, 0.3, 1.5), d$q, d$q)
  muy = c(0.1, -0.2)

  Theta = Kxx + 2 * Kxy %*% solve(Kyy) %*% t(Kxy)

  # Test BC variable (s = 2, 0-based)
  expected = ref_log_marginal_omrf(
    d$x, Theta, Kxy, mux, muy, d$num_cats,
    d$is_ordinal, d$baseline_cat, s = 2
  )
  expect_equal(res$marg_omrf_ll[3], expected, tolerance = 1e-8)
})


# ==============================================================================
# Test: marginal OMRF differs from conditional when Kxy != 0
# ==============================================================================

test_that("log_marginal_omrf differs from log_conditional_omrf with nonzero Kxy", {
  d = make_test_data()
  info = build_param_vec(d$p, d$q, d$num_cats, d$is_ordinal)
  params = rep(0, info$total)

  params[1] = 0.3; params[2] = -0.2
  params[3] = 0.1; params[4] = 0.4; params[5] = -0.1
  params[6] = 0.5; params[7] = -0.3

  kxx_start = info$num_main
  params[kxx_start + 1] = 0.2

  muy_start = info$num_main + info$num_pairwise_xx
  params[muy_start + 1] = 0.1

  kyy_start = info$kyy_start
  params[kyy_start + 1] = 2.0
  params[kyy_start + 2] = 0.3
  params[kyy_start + 3] = 1.5

  # Nonzero Kxy makes marginal != conditional
  kxy_start = kyy_start + info$num_kyy
  params[kxy_start + 1] = 0.3
  params[kxy_start + 2] = -0.2

  res = bgms:::test_mixed_mrf_likelihoods(
    d$x, d$y, d$num_cats, as.integer(d$is_ordinal), d$baseline_cat,
    d$inc_prob, d$edge_ind, FALSE, "marginal", params, 1L
  )

  # At least one variable should differ
  diffs = abs(res$marg_omrf_ll - res$omrf_ll)
  expect_true(max(diffs) > 0.01,
              label = "marginal and conditional OMRF differ with nonzero Kxy")
})


# ==============================================================================
# Helper: extract parameter matrices from bgms sample matrix
# ==============================================================================
# (Same as in test-mixed-mrf-sampling.R; duplicated here so marginal tests
# are self-contained.)
extract_bgms_estimates = function(samples, p, q, num_cats_bgms) {
  S = colMeans(samples)
  idx = 1L
  max_cat = max(num_cats_bgms)
  mux = matrix(0, p, max_cat)
  for(s in 1:p) {
    for(c in seq_len(num_cats_bgms[s])) {
      mux[s, c] = S[idx]; idx = idx + 1L
    }
  }
  Kxx = matrix(0, p, p)
  for(i in 1:(p - 1)) for(j in (i + 1):p) {
    Kxx[i, j] = Kxx[j, i] = S[idx]; idx = idx + 1L
  }
  muy = S[idx:(idx + q - 1)]; idx = idx + q
  Kyy = matrix(0, q, q)
  for(i in 1:q) for(j in i:q) {
    Kyy[i, j] = Kyy[j, i] = S[idx]; idx = idx + 1L
  }
  Kxy = matrix(S[idx:(idx + p * q - 1)], nrow = p, ncol = q, byrow = TRUE)
  list(mux = mux, Kxx = Kxx, muy = muy, Kyy = Kyy, Kxy = Kxy)
}


# ==============================================================================
# Test: marginal mode sampler agrees with mixedGM (simulation)
# ==============================================================================

test_that("marginal PL sampler agrees with mixedGM (ordinal simulation)", {
  skip_on_cran()
  skip_if_not_installed("mixedGM")

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

  # --- run mixedGM (marginal PL) ---
  set.seed(42)
  mgm = mixedGM::mixed_sampler(
    x = sim$x, y = sim$y, num_categories = rep(3L, p),
    n_warmup = 4000L, n_samples = 8000L,
    edge_selection = FALSE, verbose = FALSE,
    pseudolikelihood = "marginal"
  )
  mgm_est = mgm$final_parameters

  # --- run bgms (marginal PL) ---
  num_cats_bgms = rep(2L, p)
  pq = p + q
  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  bgms_out = bgms:::test_mixed_mrf_sampler(
    sim$x, sim$y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, FALSE, "marginal",
    4000L, 8000L, 42L
  )
  bgms_est = extract_bgms_estimates(bgms_out$samples, p, q, num_cats_bgms)

  # --- correlations bgms vs mixedGM (r >= 0.90) ---
  kxx_ut = upper.tri(bgms_est$Kxx)
  expect_gt(
    cor(bgms_est$Kxx[kxx_ut], mgm_est$Kxx[kxx_ut]), 0.90,
    label = "Kxx bgms-mgm correlation (marginal sim)"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kxy), as.vector(mgm_est$Kxy)), 0.90,
    label = "Kxy bgms-mgm correlation (marginal sim)"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kyy), as.vector(mgm_est$Kyy)), 0.90,
    label = "Kyy bgms-mgm correlation (marginal sim)"
  )
  expect_gt(
    cor(as.vector(bgms_est$mux), as.vector(mgm_est$mux)), 0.90,
    label = "mux bgms-mgm correlation (marginal sim)"
  )
})


# ==============================================================================
# Test: marginal mode sampler agrees with mixedGM on Wenchuan data
# ==============================================================================

test_that("marginal PL sampler agrees with mixedGM on Wenchuan data", {
  skip_on_cran()
  skip_if_not_installed("mixedGM")

  data(Wenchuan, package = "bgms")
  W = na.omit(Wenchuan)

  p = 4L; q = 2L
  x = as.matrix(W[, 1:p]) - 1L   # convert 1-5 to 0-4
  y = as.matrix(W[, (p + 1):(p + q)])

  num_cats_mgm = rep(5L, p)   # total categories for mixedGM
  num_cats_bgms = rep(4L, p)  # max category index for bgms

  # --- mixedGM (marginal PL) ---
  set.seed(42)
  mgm = mixedGM::mixed_sampler(
    x = x, y = y, num_categories = num_cats_mgm,
    n_warmup = 4000L, n_samples = 8000L,
    edge_selection = FALSE, verbose = FALSE,
    pseudolikelihood = "marginal"
  )
  mgm_est = mgm$final_parameters

  # --- bgms (marginal PL) ---
  pq = p + q
  inc_prob = matrix(0.5, pq, pq); diag(inc_prob) = 0
  edge_ind = matrix(1L, pq, pq); diag(edge_ind) = 0L

  bgms_out = bgms:::test_mixed_mrf_sampler(
    x, y, as.integer(num_cats_bgms),
    as.integer(rep(1L, p)), rep(0L, p),
    inc_prob, edge_ind, FALSE, "marginal",
    4000L, 8000L, 42L
  )
  bgms_est = extract_bgms_estimates(bgms_out$samples, p, q, num_cats_bgms)

  # --- correlations bgms vs mixedGM (r >= 0.90) ---
  kxx_ut = upper.tri(bgms_est$Kxx)
  expect_gt(
    cor(bgms_est$Kxx[kxx_ut], mgm_est$Kxx[kxx_ut]), 0.90,
    label = "Kxx bgms-mgm correlation (marginal Wenchuan)"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kxy), as.vector(mgm_est$Kxy)), 0.85,
    label = "Kxy bgms-mgm correlation (marginal Wenchuan)"
  )
  expect_gt(
    cor(as.vector(bgms_est$Kyy), as.vector(mgm_est$Kyy)), 0.90,
    label = "Kyy bgms-mgm correlation (marginal Wenchuan)"
  )
  expect_gt(
    cor(as.vector(bgms_est$mux), as.vector(mgm_est$mux)), 0.90,
    label = "mux bgms-mgm correlation (marginal Wenchuan)"
  )

  # --- Kxx sign patterns agree ---
  expect_true(
    all(sign(bgms_est$Kxx[kxx_ut]) == sign(mgm_est$Kxx[kxx_ut])),
    label = "Kxx signs agree (marginal Wenchuan)"
  )
})
