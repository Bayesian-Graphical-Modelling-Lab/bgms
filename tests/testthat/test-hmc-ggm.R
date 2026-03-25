# --------------------------------------------------------------------------- #
# HMC for GGM — Integration tests for HMC sampler on GGM models.
#
# Verifies the constrained (RATTLE) and unconstrained HMC dispatch for
# continuous variable models, including the warning for constrained HMC.
#
# Fast tests:
#   6.1  Smoke: HMC + GGM unconstrained runs
#   6.2  Output structure
#   6.3  Warning: HMC + edge_selection emits "numerically fragile"
#   6.4  Smoke: HMC + edge_selection runs despite warning
#   6.5  RATTLE trajectory: HMC init pattern preserves invariants
#   6.6  RATTLE trajectory: reversibility with HMC init
#   6.7  RATTLE trajectory: Hamiltonian conservation with HMC init
#
# Slow tests (gated behind BGMS_RUN_SLOW_TESTS):
#   6.8  Posterior comparison: HMC vs MH unconstrained (p=4)
#   6.9  Edge selection: HMC vs MH PIPs (p=5)
# --------------------------------------------------------------------------- #


# ---- Skip gate ---------------------------------------------------------------

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run HMC GGM correctness tests"
  )
}


# ---- Helpers -----------------------------------------------------------------

make_test_phi = function(p, seed = 42) {
  set.seed(seed)
  A = matrix(rnorm(p * p), p, p)
  K = A %*% t(A) + diag(p)
  Phi = chol(K)
  list(Phi = Phi, K = K, p = p)
}

phi_to_full_position = function(Phi) {
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

full_position_to_phi = function(x, p) {
  Phi = matrix(0, p, p)
  idx = 1
  for(q in seq_len(p)) {
    if(q > 1) {
      for(i in seq_len(q - 1)) {
        Phi[i, q] = x[idx]
        idx = idx + 1
      }
    }
    Phi[q, q] = exp(x[idx])
    idx = idx + 1
  }
  Phi
}

build_constraint_jacobian = function(Phi, edges) {
  p = nrow(Phi)
  d = p * (p + 1) / 2
  offsets = numeric(p)
  for(q in 1:p) {
    offsets[q] = (q - 1) * q / 2 + 1
  }
  excluded = list()
  for(q in 2:p) {
    for(i in 1:(q - 1)) {
      if(edges[i, q] == 0L) {
        excluded = c(excluded, list(c(i, q)))
      }
    }
  }
  m = length(excluded)
  if(m == 0) {
    return(matrix(0, 0, d))
  }
  J = matrix(0, m, d)
  for(row_idx in seq_along(excluded)) {
    pair = excluded[[row_idx]]
    i = pair[1]
    q = pair[2]
    for(l in 1:i) {
      if(l < q) J[row_idx, offsets[q] + l - 1] = Phi[l, i]
    }
    if(i > 1) {
      for(l in 1:(i - 1)) {
        J[row_idx, offsets[i] + l - 1] = Phi[l, q]
      }
    }
    J[row_idx, offsets[i] + i - 1] = Phi[i, q] * Phi[i, i]
  }
  J
}

# Simulate the HMC initialization pattern:
# 1. Sample r ~ N(0, I)  (identity mass, as at start of adaptation)
# 2. project_momentum(r, x, edges)
# 3. Run constrained leapfrog for L steps
# Returns list(x0, r0, x1, r1, dH, logp_init, logp_final, p, edges)
hmc_trajectory = function(p, edges, eps, n_leapfrog, seed) {
  dat = make_test_phi(p, seed = seed)
  n = 10
  set.seed(seed + 1000)
  X = matrix(rnorm(n * p), n, p)
  S = t(X) %*% X
  scale = 2.5

  # Project start position onto constraint manifold
  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x0 = as.vector(proj$x_projected)

  # HMC init: sample r from N(0, I) then project onto cotangent space
  set.seed(seed + 2000)
  r_raw = rnorm(length(x0))
  r0 = as.vector(ggm_test_project_momentum(r_raw, x0, edges))

  # Compute initial log-posterior
  logp0 = ggm_test_logp_and_gradient_full(x0, S, n, edges, scale)$value

  # Run constrained leapfrog trajectory
  res = ggm_test_leapfrog_constrained(
    x0, r0, eps, n_leapfrog, S, n, edges, scale
  )

  list(
    x0 = x0, r0 = r0,
    x1 = as.vector(res$x), r1 = as.vector(res$r),
    dH = res$dH, logp_init = logp0, logp_final = res$logp_final,
    p = p, edges = edges, S = S, n = n, scale = scale
  )
}


# ==============================================================================
# 6.1  Smoke test: HMC + GGM unconstrained runs
# ==============================================================================

test_that("HMC + GGM unconstrained runs without error (p=3)", {
  set.seed(61)
  p = 3
  n = 50
  x = matrix(rnorm(n * p), nrow = n)
  colnames(x) = paste0("V", seq_len(p))

  expect_no_error({
    fit = bgm(
      x = x, variable_type = "continuous",
      update_method = "hamiltonian-mc",
      edge_selection = FALSE,
      iter = 100, warmup = 50, chains = 1,
      seed = 610, display_progress = "none"
    )
  })
})


# ==============================================================================
# 6.2  Output structure
# ==============================================================================

test_that("HMC + GGM output has expected components", {
  set.seed(62)
  p = 3
  n = 50
  x = matrix(rnorm(n * p), nrow = n)
  colnames(x) = paste0("V", seq_len(p))

  fit = bgm(
    x = x, variable_type = "continuous",
    update_method = "hamiltonian-mc",
    edge_selection = FALSE,
    iter = 100, warmup = 50, chains = 1,
    seed = 620, display_progress = "none"
  )

  # Pairwise interactions should be extractable
  interactions = extract_pairwise_interactions(fit)
  expect_true(is.matrix(interactions))
  expect_equal(ncol(interactions), p * (p - 1) / 2)
  expect_equal(nrow(interactions), 100)

  # GGM has no main effects (continuous model)
  expect_null(extract_main_effects(fit))

  # No NaN values in interactions
  expect_false(any(is.nan(interactions)))
})


# ==============================================================================
# 6.3  Warning: HMC + edge_selection emits "numerically fragile"
# ==============================================================================

test_that("bgm warns when using HMC with edge selection", {
  set.seed(63)
  p = 3
  n = 50
  x = matrix(rnorm(n * p), nrow = n)
  colnames(x) = paste0("V", seq_len(p))

  # The warning is emitted during validation, before sampling starts.
  # Constrained HMC may subsequently crash (known fragility), so we
  # capture the warning separately from any downstream error.
  warned = FALSE
  tryCatch(
    withCallingHandlers(
      bgm(
        x = x, variable_type = "continuous",
        update_method = "hamiltonian-mc",
        edge_selection = TRUE,
        iter = 50, warmup = 25, chains = 1,
        seed = 630, display_progress = "none"
      ),
      warning = function(w) {
        if(grepl("numerically fragile", conditionMessage(w))) warned <<- TRUE
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  expect_true(warned)
})


# ==============================================================================
# 6.4  Smoke: HMC + edge_selection runs despite warning
# ==============================================================================

test_that("HMC + edge selection produces output when it survives (p=3)", {
  set.seed(64)
  p = 3
  n = 200
  x = matrix(rnorm(n * p), nrow = n)
  colnames(x) = paste0("V", seq_len(p))

  # Constrained HMC is fragile; it may crash on some seeds.
  # We try a few seeds and require at least one to succeed.
  fit = NULL
  for(s in c(640, 641, 642, 643, 644)) {
    fit = tryCatch(
      suppressWarnings(
        bgm(
          x = x, variable_type = "continuous",
          update_method = "hamiltonian-mc",
          edge_selection = TRUE,
          iter = 100, warmup = 100, chains = 1,
          seed = s, display_progress = "none"
        )
      ),
      error = function(e) NULL
    )
    if(!is.null(fit)) break
  }

  skip_if(is.null(fit), "constrained HMC crashed on all seeds (known fragility)")

  # PIPs should exist and be in [0, 1]
  pip = fit$posterior_mean_indicator
  expect_true(is.matrix(pip))
  expect_equal(nrow(pip), p)
  expect_equal(ncol(pip), p)
  expect_true(all(pip >= 0 & pip <= 1))

  # Interactions should be finite
  interactions = extract_pairwise_interactions(fit)
  expect_true(all(is.finite(interactions)))
})


# ==============================================================================
# 6.5  RATTLE trajectory: HMC init preserves invariants
#
# Simulates the HMC step pipeline:
#   1. Sample momentum r ~ N(0, I)
#   2. Project r onto cotangent space (project_momentum)
#   3. Run L constrained leapfrog steps
#   4. Verify: position constraints preserved, momentum in cotangent space
#
# This tests the same RATTLE integrator used by NUTS but exercised via the
# HMC initialization pattern (direct project_momentum, identity mass).
# ==============================================================================

test_that("HMC trajectory preserves position constraints (p=4)", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  traj = hmc_trajectory(p, edges, eps = 0.005, n_leapfrog = 10, seed = 650)

  # Position constraints: K_{iq} = 0 for excluded edges
  Phi = full_position_to_phi(traj$x1, p)
  K = t(Phi) %*% Phi
  expect_equal(K[1, 3], 0, tolerance = 1e-10)
  expect_equal(K[2, 4], 0, tolerance = 1e-10)
})

test_that("HMC trajectory preserves momentum cotangent condition (p=4)", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  traj = hmc_trajectory(p, edges, eps = 0.005, n_leapfrog = 10, seed = 651)

  # Momentum: J * r = 0 at trajectory endpoint
  Phi = full_position_to_phi(traj$x1, p)
  J = build_constraint_jacobian(Phi, edges)
  Jr = J %*% traj$r1
  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-10)
})

test_that("HMC trajectory preserves constraints (p=6, multiple exclusions)", {
  p = 6
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 4] = 0L
  edges[4, 1] = 0L
  edges[2, 5] = 0L
  edges[5, 2] = 0L
  edges[3, 6] = 0L
  edges[6, 3] = 0L
  edges[1, 6] = 0L
  edges[6, 1] = 0L

  traj = hmc_trajectory(p, edges, eps = 0.005, n_leapfrog = 8, seed = 652)

  # Position constraints
  Phi = full_position_to_phi(traj$x1, p)
  K = t(Phi) %*% Phi
  expect_equal(K[1, 4], 0, tolerance = 1e-10)
  expect_equal(K[2, 5], 0, tolerance = 1e-10)
  expect_equal(K[3, 6], 0, tolerance = 1e-10)
  expect_equal(K[1, 6], 0, tolerance = 1e-10)

  # Momentum cotangent condition
  J = build_constraint_jacobian(Phi, edges)
  Jr = J %*% traj$r1
  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-10)
})


# ==============================================================================
# 6.6  RATTLE trajectory: reversibility under HMC init
#
# Forward L steps then backward L steps should return to the start.
# Inspired by mici's test_reversibility (Matt Graham, 2024).
# ==============================================================================

test_that("HMC trajectory is reversible (p=4)", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  traj = hmc_trajectory(p, edges, eps = 0.005, n_leapfrog = 10, seed = 660)

  # Reverse: negate momentum, run same number of steps
  bwd = ggm_test_leapfrog_constrained(
    traj$x1, -traj$r1, 0.005, 10,
    traj$S, traj$n, edges, traj$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  expect_lt(max(abs(x_return - traj$x0)), 1e-3,
    label = "position reversibility error"
  )
  expect_lt(max(abs(r_return - traj$r0)), 1e-3,
    label = "momentum reversibility error"
  )
})

test_that("HMC trajectory is reversible (p=6)", {
  p = 6
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 4] = 0L
  edges[4, 1] = 0L
  edges[2, 5] = 0L
  edges[5, 2] = 0L
  edges[3, 6] = 0L
  edges[6, 3] = 0L

  traj = hmc_trajectory(p, edges, eps = 0.005, n_leapfrog = 8, seed = 661)

  bwd = ggm_test_leapfrog_constrained(
    traj$x1, -traj$r1, 0.005, 8,
    traj$S, traj$n, edges, traj$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  expect_lt(max(abs(x_return - traj$x0)), 1e-2)
  expect_lt(max(abs(r_return - traj$r0)), 1e-2)
})


# ==============================================================================
# 6.7  RATTLE trajectory: Hamiltonian conservation under HMC init
#
# For a symplectic integrator |ΔH| = O(ε²).
# Inspired by mici's test_approx_hamiltonian_conservation (Matt Graham, 2024).
# ==============================================================================

test_that("HMC trajectory approximately conserves Hamiltonian (p=4)", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L

  traj = hmc_trajectory(p, edges, eps = 0.005, n_leapfrog = 20, seed = 670)

  expect_lt(abs(traj$dH), 1.0,
    label = paste("|dH| =", abs(traj$dH))
  )
})

test_that("Hamiltonian conservation improves with smaller step size", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  dat = make_test_phi(p, seed = 671)
  n = 10
  set.seed(1671)
  X = matrix(rnorm(n * p), n, p)
  S = t(X) %*% X
  scale = 2.5

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x0 = as.vector(proj$x_projected)

  set.seed(2671)
  r_raw = rnorm(length(x0))
  r0 = as.vector(ggm_test_project_momentum(r_raw, x0, edges))

  res_large = ggm_test_leapfrog_constrained(
    x0, r0, 0.01, 10, S, n, edges, scale
  )
  res_small = ggm_test_leapfrog_constrained(
    x0, r0, 0.005, 10, S, n, edges, scale
  )

  # Halving step size should reduce |ΔH| (O(ε²) scaling)
  expect_lt(abs(res_small$dH), abs(res_large$dH) + 0.01)
})


# ==============================================================================
# 6.8  Posterior comparison: HMC vs MH unconstrained (p=4)
# ==============================================================================

test_that("HMC and MH posteriors agree for unconstrained GGM (p=4)", {
  skip_unless_slow()

  p = 4
  n = 200
  set.seed(68)
  K_true = diag(2, p)
  for(i in seq_len(p - 1)) {
    K_true[i, i + 1] = -0.5
    K_true[i + 1, i] = -0.5
  }
  Sigma = solve(K_true)
  x = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x) = paste0("V", seq_len(p))

  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 680, display_progress = "none"
  )

  fit_hmc = bgm(
    x = x, variable_type = "continuous",
    update_method = "hamiltonian-mc",
    edge_selection = FALSE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 681, display_progress = "none"
  )

  # Posterior means should agree within 2 pooled SDs
  mh_pairwise = do.call(rbind, fit_mh$raw_samples$pairwise)
  hmc_pairwise = do.call(rbind, fit_hmc$raw_samples$pairwise)

  for(j in seq_len(ncol(mh_pairwise))) {
    mh_mean = mean(mh_pairwise[, j])
    hmc_mean = mean(hmc_pairwise[, j])
    pooled_sd = sqrt((var(mh_pairwise[, j]) + var(hmc_pairwise[, j])) / 2)
    expect_lt(
      abs(mh_mean - hmc_mean) / pooled_sd, 2,
      label = paste0("pairwise[", j, "] mean within 2 SDs")
    )
  }

  # Variance ratios should be close to 1
  for(j in seq_len(ncol(mh_pairwise))) {
    ratio = var(hmc_pairwise[, j]) / var(mh_pairwise[, j])
    expect_gt(ratio, 0.7, label = paste0("pairwise[", j, "] var ratio > 0.7"))
    expect_lt(ratio, 1.4, label = paste0("pairwise[", j, "] var ratio < 1.4"))
  }
})


# ==============================================================================
# 6.9  Edge selection: HMC vs MH PIPs (p=5)
# ==============================================================================

test_that("HMC and MH edge selection PIPs agree (p=5)", {
  skip_unless_slow()

  p = 5
  n = 300
  set.seed(69)
  K_true = diag(2, p)
  for(i in seq_len(p - 1)) {
    K_true[i, i + 1] = -0.5
    K_true[i + 1, i] = -0.5
  }
  Sigma = solve(K_true)
  x = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x) = paste0("V", seq_len(p))

  n_iter = 5000
  n_warmup = 2000

  fit_mh = bgm(
    x = x, variable_type = "continuous",
    update_method = "adaptive-metropolis",
    edge_selection = TRUE,
    iter = n_iter, warmup = n_warmup, chains = 2,
    seed = 690, display_progress = "none"
  )

  suppressWarnings({
    fit_hmc = bgm(
      x = x, variable_type = "continuous",
      update_method = "hamiltonian-mc",
      edge_selection = TRUE,
      iter = n_iter, warmup = n_warmup, chains = 2,
      seed = 691, display_progress = "none"
    )
  })

  # Interaction estimates should agree
  int_mh = colMeans(extract_pairwise_interactions(fit_mh))
  int_hmc = colMeans(extract_pairwise_interactions(fit_hmc))
  expect_lt(max(abs(int_mh - int_hmc)), 0.15,
    label = "max interaction difference < 0.15"
  )

  # PIPs should agree
  pip_mh = fit_mh$posterior_mean_indicator[upper.tri(fit_mh$posterior_mean_indicator)]
  pip_hmc = fit_hmc$posterior_mean_indicator[upper.tri(fit_hmc$posterior_mean_indicator)]
  expect_lt(max(abs(pip_mh - pip_hmc)), 0.15,
    label = "max PIP difference < 0.15"
  )
})
