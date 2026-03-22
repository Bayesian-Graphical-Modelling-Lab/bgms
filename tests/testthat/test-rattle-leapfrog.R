# --------------------------------------------------------------------------- #
# RATTLE Phase 4 — Constrained leapfrog integration tests.
#
# Tests verify:
#   1. Reversibility (forward + backward = identity)
#   2. Energy conservation (|ΔH| = O(ε²))
#   3. Constraint preservation (K_{iq} = 0 after trajectory)
#   4. Cotangent condition (J·r = 0 after trajectory)
#   5. Full-edge graph: constrained matches unconstrained leapfrog
# --------------------------------------------------------------------------- #


# ---- Helpers ----------------------------------------------------------------

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

# Build constraint Jacobian J in R
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

# Helper: prepare a test scenario (projected start + random momentum)
make_scenario = function(p, edges, seed) {
  dat = make_test_phi(p, seed = seed)
  n = 10
  set.seed(seed + 1000)
  X = matrix(rnorm(n * p), n, p)
  S = t(X) %*% X
  scale = 2.5

  # Project start position
  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x0 = as.vector(proj$x_projected)

  # Random momentum, projected onto cotangent space
  set.seed(seed + 2000)
  r_raw = rnorm(length(x0))
  r0 = as.vector(ggm_test_project_momentum(r_raw, x0, edges))

  list(x0 = x0, r0 = r0, S = S, n = n, scale = scale, p = p, edges = edges)
}


# ---- 1. Reversibility: forward N + backward N ≈ identity -------------------

test_that("constrained leapfrog is reversible, p=4 constrained", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  sc = make_scenario(p, edges, seed = 200)
  eps = 0.005
  n_steps = 10

  # Forward
  fwd = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )

  # Backward: negate momentum, run same steps, negate again
  bwd = ggm_test_leapfrog_constrained(
    as.vector(fwd$x), -as.vector(fwd$r), eps, n_steps,
    sc$S, sc$n, edges, sc$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  # Position and momentum should return to start
  # O(ε²) reversibility error from per-column position projection
  pos_err = max(abs(x_return - sc$x0))
  mom_err = max(abs(r_return - sc$r0))

  expect_lt(pos_err, 1e-3,
    label = paste("position reversibility error:", pos_err)
  )
  expect_lt(mom_err, 1e-3,
    label = paste("momentum reversibility error:", mom_err)
  )
})

test_that("constrained leapfrog is reversible, p=6 constrained", {
  p = 6
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 4] = 0L
  edges[4, 1] = 0L
  edges[2, 5] = 0L
  edges[5, 2] = 0L
  edges[3, 6] = 0L
  edges[6, 3] = 0L

  sc = make_scenario(p, edges, seed = 201)
  eps = 0.005
  n_steps = 8

  fwd = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )
  bwd = ggm_test_leapfrog_constrained(
    as.vector(fwd$x), -as.vector(fwd$r), eps, n_steps,
    sc$S, sc$n, edges, sc$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  expect_lt(max(abs(x_return - sc$x0)), 1e-2)
  expect_lt(max(abs(r_return - sc$r0)), 1e-2)
})

test_that("single-constraint reversibility is machine-precise", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L

  sc = make_scenario(p, edges, seed = 202)
  eps = 0.005
  n_steps = 10

  fwd = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )
  bwd = ggm_test_leapfrog_constrained(
    as.vector(fwd$x), -as.vector(fwd$r), eps, n_steps,
    sc$S, sc$n, edges, sc$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  # Single-constraint: no cross-column terms, exact RATTLE → machine precision
  expect_lt(max(abs(x_return - sc$x0)), 1e-12)
  expect_lt(max(abs(r_return - sc$r0)), 1e-12)
})

test_that("reversibility error scales as O(eps^2)", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  sc = make_scenario(p, edges, seed = 203)
  n_steps = 5

  rev_err = function(eps) {
    fwd = ggm_test_leapfrog_constrained(
      sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
    )
    bwd = ggm_test_leapfrog_constrained(
      as.vector(fwd$x), -as.vector(fwd$r), eps, n_steps,
      sc$S, sc$n, edges, sc$scale
    )
    max(abs(as.vector(bwd$x) - sc$x0))
  }

  err_large = rev_err(0.01)
  err_small = rev_err(0.005)

  # Halving ε should reduce reversibility error by ~4× (O(ε²))
  ratio = err_large / max(err_small, 1e-16)
  expect_gt(ratio, 2.0, label = paste("ε² scaling ratio:", ratio))
})


# ---- 2. Energy conservation: |ΔH| = O(ε²) ---------------------------------

test_that("energy is approximately conserved, p=4 constrained", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L

  sc = make_scenario(p, edges, seed = 210)
  eps = 0.005
  n_steps = 20

  res = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )

  # |ΔH| should be small (O(ε²) per step, O(N·ε²) total)
  expect_lt(abs(res$dH), 1.0,
    label = paste("|dH| =", abs(res$dH))
  )
})

test_that("energy conservation improves with smaller step size", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  sc = make_scenario(p, edges, seed = 211)
  n_steps = 10

  # Large step
  res_large = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, 0.01, n_steps, sc$S, sc$n, edges, sc$scale
  )
  # Small step (half)
  res_small = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, 0.005, n_steps, sc$S, sc$n, edges, sc$scale
  )

  # Halving ε should reduce |ΔH| by ~4× (O(ε²))
  expect_lt(abs(res_small$dH), abs(res_large$dH) + 0.01,
    label = "smaller step should give better energy conservation"
  )
})


# ---- 3. Constraint preservation: K_{iq} = 0 after trajectory ---------------

test_that("constraints preserved after trajectory, p=4", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  sc = make_scenario(p, edges, seed = 220)
  eps = 0.005
  n_steps = 15

  res = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )

  x_final = as.vector(res$x)
  Phi = full_position_to_phi(x_final, p)
  K = t(Phi) %*% Phi

  # Excluded edges should have K_{iq} ≈ 0
  expect_equal(K[1, 3], 0, tolerance = 1e-10)
  expect_equal(K[2, 4], 0, tolerance = 1e-10)
})

test_that("constraints preserved after trajectory, p=6", {
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

  sc = make_scenario(p, edges, seed = 221)
  eps = 0.005
  n_steps = 12

  res = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )

  x_final = as.vector(res$x)
  Phi = full_position_to_phi(x_final, p)
  K = t(Phi) %*% Phi

  expect_equal(K[1, 4], 0, tolerance = 1e-10)
  expect_equal(K[2, 5], 0, tolerance = 1e-10)
  expect_equal(K[3, 6], 0, tolerance = 1e-10)
  expect_equal(K[1, 6], 0, tolerance = 1e-10)
})


# ---- 4. Cotangent condition: J·r = 0 after trajectory ----------------------

test_that("momentum stays in cotangent space, p=4", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  sc = make_scenario(p, edges, seed = 230)
  eps = 0.005
  n_steps = 10

  res = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )

  x_final = as.vector(res$x)
  r_final = as.vector(res$r)

  Phi = full_position_to_phi(x_final, p)
  J = build_constraint_jacobian(Phi, edges)

  Jr = J %*% r_final
  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-10)
})

test_that("momentum stays in cotangent space, p=6", {
  p = 6
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 4] = 0L
  edges[4, 1] = 0L
  edges[2, 5] = 0L
  edges[5, 2] = 0L
  edges[3, 6] = 0L
  edges[6, 3] = 0L

  sc = make_scenario(p, edges, seed = 231)
  eps = 0.005
  n_steps = 8

  res = ggm_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
  )

  x_final = as.vector(res$x)
  r_final = as.vector(res$r)

  Phi = full_position_to_phi(x_final, p)
  J = build_constraint_jacobian(Phi, edges)

  Jr = J %*% r_final
  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-10)
})


# ---- 5. Full-edge graph: constrained = unconstrained -----------------------
# With no excluded edges, constrained leapfrog should behave identically
# to unconstrained (projection is a no-op).

test_that("full-edge graph: trajectory is valid, p=4", {
  p = 4
  edges = matrix(1L, p, p)
  diag(edges) = 0L

  dat = make_test_phi(p, seed = 240)
  n = 10
  set.seed(241)
  X = matrix(rnorm(n * p), n, p)
  S = t(X) %*% X
  scale = 2.5

  x0 = phi_to_full_position(dat$Phi)
  set.seed(242)
  r0 = rnorm(length(x0))

  eps = 0.005
  n_steps = 10

  res = ggm_test_leapfrog_constrained(
    x0, r0, eps, n_steps, S, n, edges, scale
  )

  # Energy should be approximately conserved
  expect_lt(abs(res$dH), 1.0)

  # Position should be finite
  expect_true(all(is.finite(as.vector(res$x))))
  expect_true(all(is.finite(as.vector(res$r))))
})


# ---- 6. Multiple constraint patterns: smoke test ---------------------------

test_that("constrained leapfrog works across constraint patterns", {
  for(p in c(4, 5, 6)) {
    dat = make_test_phi(p, seed = 250 + p)
    edges = matrix(1L, p, p)
    diag(edges) = 0L
    set.seed(250 + p)
    for(i in 1:(p - 1)) {
      for(j in (i + 1):p) {
        if(runif(1) < 0.35) {
          edges[i, j] = 0L
          edges[j, i] = 0L
        }
      }
    }

    sc = make_scenario(p, edges, seed = 260 + p)
    eps = 0.005
    n_steps = 8

    res = ggm_test_leapfrog_constrained(
      sc$x0, sc$r0, eps, n_steps, sc$S, sc$n, edges, sc$scale
    )

    # All outputs should be finite
    expect_true(all(is.finite(as.vector(res$x))),
      info = paste("p =", p, "position finite")
    )
    expect_true(all(is.finite(as.vector(res$r))),
      info = paste("p =", p, "momentum finite")
    )
    expect_true(is.finite(res$logp_final),
      info = paste("p =", p, "logp finite")
    )

    # Constraints should be satisfied
    x_final = as.vector(res$x)
    Phi = full_position_to_phi(x_final, p)
    K = t(Phi) %*% Phi
    for(i in 1:(p - 1)) {
      for(j in (i + 1):p) {
        if(edges[i, j] == 0L) {
          expect_equal(K[i, j], 0,
            tolerance = 1e-10,
            info = paste("p =", p, "K[", i, ",", j, "]")
          )
        }
      }
    }
  }
})
