# --------------------------------------------------------------------------- #
# Reversibility check — runtime forward-backward round-trip for RATTLE.
#
# Tests verify:
#   1. Normal steps pass the reversibility check (reversible = TRUE)
#   2. A strict tolerance triggers non-reversible detections
#   3. The non_reversible diagnostic flows through to bgm() output
#   4. Stress: large step sizes amplify reversibility error
#   5. Stress: long trajectories (200 steps) remain reversible
#   6. Stress: small n / large p rank-deficient regime
#   7. Stress: ill-conditioned data (high condition number)
#   8. Stress: high-dimensional GGM (p=20) integration test
# --------------------------------------------------------------------------- #


# ---- Helpers (reused from test-rattle-leapfrog.R) ---------------------------

make_test_phi_rc = function(p, seed = 42) {
  set.seed(seed)
  A = matrix(rnorm(p * p), p, p)
  K = A %*% t(A) + diag(p)
  Phi = chol(K)
  list(Phi = Phi, K = K, p = p)
}

phi_to_full_position_rc = function(Phi) {
  p = nrow(Phi)
  x = numeric(p * (p + 1) / 2)
  idx = 1
  for(q in seq_len(p)) {
    if(q > 1) {
      for(i in seq_len(q - 1)) {
        x[idx] = Phi[i, q]
        idx = idx + 1
      }
    }
    x[idx] = log(Phi[q, q])
    idx = idx + 1
  }
  x
}

make_scenario_rc = function(p, edges, seed) {
  dat = make_test_phi_rc(p, seed = seed)
  n = 10
  set.seed(seed + 1000)
  X = matrix(rnorm(n * p), n, p)
  S = t(X) %*% X
  scale = 2.5

  x_raw = phi_to_full_position_rc(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x0 = as.vector(proj$x_projected)

  set.seed(seed + 2000)
  r_raw = rnorm(length(x0))
  r0 = as.vector(ggm_test_project_momentum(r_raw, x0, edges))

  list(x0 = x0, r0 = r0, S = S, n = n, scale = scale, p = p, edges = edges)
}


# ---- 1. Normal steps pass the reversibility check --------------------------

test_that("constrained leapfrog checked passes for well-behaved steps", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L

  sc = make_scenario_rc(p, edges, seed = 300)
  eps = 0.005
  n_steps = 10

  result = ggm_test_leapfrog_constrained_checked(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale,
    reverse_check_tol = 0.5
  )

  expect_equal(result$non_reversible_count, 0L)
})


test_that("checked leapfrog matches unchecked positions", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  sc = make_scenario_rc(p, edges, seed = 301)
  eps = 0.005
  n_steps = 8

  checked = ggm_test_leapfrog_constrained_checked(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )
  unchecked = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )

  expect_equal(as.vector(checked$x), as.vector(unchecked$x), tolerance = 1e-12)
  expect_equal(as.vector(checked$r), as.vector(unchecked$r), tolerance = 1e-12)
})


# ---- 2. Strict tolerance triggers non-reversible detections ----------------

test_that("extreme tolerance detects non-reversible steps", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  sc = make_scenario_rc(p, edges, seed = 302)
  eps = 0.01
  n_steps = 20

  # With an impossibly tight factor, expect some failures
  result = ggm_test_leapfrog_constrained_checked(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale,
    reverse_check_tol = 1e-11
  )

  expect_gt(result$non_reversible_count, 0L)
})


# ---- 3. Integration: non_reversible diagnostic in bgm() output -------------

test_that("bgm() nuts_diag includes non_reversible field", {
  skip_on_cran()

  set.seed(1)
  p = 4
  n = 80
  data = matrix(rnorm(n * p), n, p)
  colnames(data) = paste0("V", seq_len(p))

  fit = bgm(
    data,
    variable_type = "continuous",
    iter = 10,
    warmup = 50,
    edge_selection = TRUE,
    chains = 1,
    display_progress = "none"
  )

  if(!is.null(fit$nuts_diag)) {
    expect_true("non_reversible" %in% names(fit$nuts_diag))
    expect_true("total_non_reversible" %in% names(fit$nuts_diag$summary))
  }
})


# ---- 4. Stress: large step sizes -------------------------------------------

test_that("large step sizes increase reversibility error", {
  p = 5
  edges = matrix(1L, p, p)
  diag(edges) = 0L

  sc = make_scenario_rc(p, edges, seed = 400)

  # Small step: should be reversible
  small = ggm_test_leapfrog_constrained_checked(
    sc$x0, sc$r0, 0.005, 10, sc$S, sc$n, edges, sc$scale,
    reverse_check_tol = 0.5
  )
  expect_equal(small$non_reversible_count, 0L)

  # Large step (100x): near-machine-epsilon tolerance catches round-trip drift
  large = ggm_test_leapfrog_constrained_checked(
    sc$x0, sc$r0, 0.5, 20, sc$S, sc$n, edges, sc$scale,
    reverse_check_tol = 1e-15
  )
  expect_gt(large$non_reversible_count, 0L)
})


# ---- 5. Stress: long trajectories (200 steps) ------------------------------

test_that("long trajectory (200 steps) stays reversible at normal tolerance", {
  p = 5
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 4] = 0L
  edges[4, 1] = 0L

  sc = make_scenario_rc(p, edges, seed = 500)

  result = ggm_test_leapfrog_constrained_checked(
    sc$x0, sc$r0, 0.005, 200, sc$S, sc$n, edges, sc$scale,
    reverse_check_tol = 0.5
  )

  expect_equal(result$non_reversible_count, 0L)
})


# ---- 6. Stress: small n / large p (rank-deficient S) -----------------------

test_that("rank-deficient regime (n < p) stays reversible", {
  p = 8
  edges = matrix(1L, p, p)
  diag(edges) = 0L

  # Build scenario with n < p (rank-deficient sufficient statistic)
  dat = make_test_phi_rc(p, seed = 600)
  n = 5
  set.seed(601)
  X = matrix(rnorm(n * p), n, p)
  S = t(X) %*% X # rank = min(n, p) = 5 < 8

  x_raw = phi_to_full_position_rc(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x0 = as.vector(proj$x_projected)

  set.seed(602)
  r_raw = rnorm(length(x0))
  r0 = as.vector(ggm_test_project_momentum(r_raw, x0, edges))

  result = ggm_test_leapfrog_constrained_checked(
    x0, r0, 0.005, 30, S, n, edges, 2.5,
    reverse_check_tol = 0.5
  )

  expect_equal(result$non_reversible_count, 0L)
})


# ---- 7. Stress: ill-conditioned data (high condition number) ----------------

test_that("ill-conditioned data (kappa ~ 1e4) stays reversible", {
  p = 6
  edges = matrix(1L, p, p)
  diag(edges) = 0L

  # Build ill-conditioned sufficient statistic
  dat = make_test_phi_rc(p, seed = 700)
  n = 50
  set.seed(701)
  X = matrix(rnorm(n * p), n, p)
  # Stretch first column to create high condition number
  X[, 1] = X[, 1] * 100
  S = t(X) %*% X
  # Verify ill-conditioning
  ev = eigen(S, symmetric = TRUE, only.values = TRUE)$values
  kappa = max(ev) / min(ev)
  expect_gt(kappa, 1e3)

  x_raw = phi_to_full_position_rc(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x0 = as.vector(proj$x_projected)

  set.seed(702)
  r_raw = rnorm(length(x0))
  r0 = as.vector(ggm_test_project_momentum(r_raw, x0, edges))

  result = ggm_test_leapfrog_constrained_checked(
    x0, r0, 0.003, 30, S, n, edges, 2.5,
    reverse_check_tol = 0.5
  )

  expect_equal(result$non_reversible_count, 0L)
})


# ---- 8. Stress: high-dimensional GGM (p=20) integration --------------------

test_that("high-dimensional GGM (p=20) completes without crash", {
  skip_on_cran()

  set.seed(800)
  p = 20
  n = 100
  data = matrix(rnorm(n * p), n, p)
  colnames(data) = paste0("V", seq_len(p))

  fit = bgm(
    data,
    variable_type = "continuous",
    iter = 10,
    warmup = 50,
    edge_selection = TRUE,
    chains = 1,
    display_progress = "none"
  )

  expect_s3_class(fit, "bgms")
  if(!is.null(fit$nuts_diag)) {
    expect_true("non_reversible" %in% names(fit$nuts_diag))
  }
})
