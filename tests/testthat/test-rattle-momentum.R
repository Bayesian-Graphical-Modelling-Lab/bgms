# --------------------------------------------------------------------------- #
# RATTLE Phase 3 — Momentum projection tests.
#
# Tests verify:
#   1. J r = 0 after projection (cotangent space condition)
#   2. Idempotency (projecting twice = projecting once)
#   3. Full-edge graph: no-op (no constraints)
#   4. Projection removes only the normal component
#   5. Numerical Jacobian verification
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


# Build constraint Jacobian J in R (reference implementation)
# Returns m x d matrix where m = total excluded edges,
# d = p(p+1)/2 = length(x).
build_constraint_jacobian = function(Phi, edges) {
  p = nrow(Phi)
  d = p * (p + 1) / 2

  # Compute column offsets (1-indexed, matching C++ full_theta_offsets)
  # Column q (1-indexed) starts at position 1 + sum_{k=1}^{q-1} k
  # = 1 + q*(q-1)/2
  offsets = numeric(p)
  for(q in 1:p) {
    offsets[q] = 1 + (q - 1) * (q - 2) / 2 + (q - 1)
    # Simpler: offsets[q] = q*(q-1)/2 + 1
    offsets[q] = (q - 1) * q / 2 + 1
  }

  # Collect excluded edges
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
    i = pair[1] # 1-indexed R row
    q = pair[2] # 1-indexed R col

    # The offset for column q: entries are [x_{1,q}, ..., x_{q-1,q}, psi_q]
    # x_{l,q} is at position offsets[q] + (l - 1) for l = 1..q-1
    # psi_q is at position offsets[q] + (q - 1)

    # Type 1: dc/d(x_{l,q}) = Phi_{l,i} for l = 1..i (R 1-indexed)
    for(l in 1:i) {
      if(l < q) {
        J[row_idx, offsets[q] + l - 1] = Phi[l, i]
      }
    }

    # Type 2: dc/d(x_{l,i}) = Phi_{l,q} for l = 1..(i-1)
    if(i > 1) {
      for(l in 1:(i - 1)) {
        J[row_idx, offsets[i] + l - 1] = Phi[l, q]
      }
    }

    # Diagonal chain rule: dc/d(psi_i) = Phi_{i,q} * Phi_{i,i}
    # psi_i is the last entry in column i's block: offsets[i] + (i - 1)
    J[row_idx, offsets[i] + i - 1] = Phi[i, q] * Phi[i, i]
  }

  J
}


# ---- 1. J r = 0 after projection (cotangent space condition) ---------------

test_that("projected momentum satisfies J*r = 0, p=4", {
  p = 4
  dat = make_test_phi(p, seed = 100)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  # Project position first
  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  # Random momentum
  set.seed(100)
  r = rnorm(length(x))

  r_proj = as.vector(ggm_test_project_momentum(r, x, edges))

  # Build J in R and verify J * r_proj ≈ 0
  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)
  Jr = J %*% r_proj

  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-12)
})

test_that("projected momentum satisfies J*r = 0, p=6 many constraints", {
  p = 6
  dat = make_test_phi(p, seed = 101)
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
  edges[2, 3] = 0L
  edges[3, 2] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  set.seed(101)
  r = rnorm(length(x))
  r_proj = as.vector(ggm_test_project_momentum(r, x, edges))

  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)
  Jr = J %*% r_proj

  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-12)
})


# ---- 2. Idempotency -------------------------------------------------------

test_that("momentum projection is idempotent, p=4", {
  p = 4
  dat = make_test_phi(p, seed = 110)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  set.seed(110)
  r = rnorm(length(x))

  r1 = as.vector(ggm_test_project_momentum(r, x, edges))
  r2 = as.vector(ggm_test_project_momentum(r1, x, edges))

  expect_equal(r2, r1, tolerance = 1e-14)
})

test_that("momentum projection is idempotent, p=6", {
  p = 6
  dat = make_test_phi(p, seed = 111)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 4] = 0L
  edges[4, 1] = 0L
  edges[2, 5] = 0L
  edges[5, 2] = 0L
  edges[3, 6] = 0L
  edges[6, 3] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  set.seed(111)
  r = rnorm(length(x))

  r1 = as.vector(ggm_test_project_momentum(r, x, edges))
  r2 = as.vector(ggm_test_project_momentum(r1, x, edges))

  expect_equal(r2, r1, tolerance = 1e-14)
})


# ---- 3. Full-edge graph: no-op -------------------------------------------

test_that("full-edge graph: momentum unchanged, p=4", {
  p = 4
  dat = make_test_phi(p, seed = 120)
  edges = matrix(1L, p, p)
  diag(edges) = 0L

  x = phi_to_full_position(dat$Phi)

  set.seed(120)
  r = rnorm(length(x))

  r_proj = as.vector(ggm_test_project_momentum(r, x, edges))

  expect_equal(r_proj, r, tolerance = 1e-14)
})

test_that("full-edge graph: momentum unchanged, p=6", {
  p = 6
  dat = make_test_phi(p, seed = 121)
  edges = matrix(1L, p, p)
  diag(edges) = 0L

  x = phi_to_full_position(dat$Phi)

  set.seed(121)
  r = rnorm(length(x))

  r_proj = as.vector(ggm_test_project_momentum(r, x, edges))

  expect_equal(r_proj, r, tolerance = 1e-14)
})


# ---- 4. Projection removes only the normal component ----------------------
# After projecting, the change (r - r_proj) should lie in the row space of J
# (i.e., be a linear combination of J's rows).

test_that("projection change is in row space of J, p=4", {
  p = 4
  dat = make_test_phi(p, seed = 130)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  set.seed(130)
  r = rnorm(length(x))
  r_proj = as.vector(ggm_test_project_momentum(r, x, edges))

  delta_r = r - r_proj

  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)

  # delta_r should be in row space of J: delta_r = J^T lambda
  # Verify by checking that projecting delta_r onto null(J) gives 0
  # Equivalently, delta_r - J^T (J J^T)^{-1} J delta_r ≈ 0
  G = J %*% t(J)
  lambda = solve(G, J %*% delta_r)
  reconstructed = as.vector(t(J) %*% lambda)

  expect_equal(reconstructed, delta_r, tolerance = 1e-12)
})


# ---- 5. Numerical Jacobian verification -----------------------------------
# Build J numerically by perturbing x and checking constraint changes,
# then verify it matches the analytic J used by the C++ code.

test_that("numerical constraint Jacobian matches analytic, p=4", {
  p = 4
  dat = make_test_phi(p, seed = 140)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  # Constraint function: returns vector of K_{iq} for excluded edges
  constraint_fn = function(x_vec) {
    Phi = full_position_to_phi(x_vec, p)
    K = t(Phi) %*% Phi
    c_vals = numeric(0)
    for(q in 2:p) {
      for(i in 1:(q - 1)) {
        if(edges[i, q] == 0L) {
          c_vals = c(c_vals, K[i, q])
        }
      }
    }
    c_vals
  }

  # Numerical Jacobian
  eps = 1e-7
  c0 = constraint_fn(x)
  m = length(c0)
  d = length(x)
  J_num = matrix(0, m, d)
  for(k in seq_len(d)) {
    x_plus = x
    x_plus[k] = x_plus[k] + eps
    x_minus = x
    x_minus[k] = x_minus[k] - eps
    J_num[, k] = (constraint_fn(x_plus) - constraint_fn(x_minus)) / (2 * eps)
  }

  # Analytic Jacobian
  Phi = full_position_to_phi(x, p)
  J_ana = build_constraint_jacobian(Phi, edges)

  expect_equal(J_ana, J_num, tolerance = 1e-6)
})


# ---- 6. Momentum projection with varying constraint patterns ---------------

test_that("momentum projection correct across constraint patterns", {
  for(p in c(4, 5, 6, 8)) {
    dat = make_test_phi(p, seed = 150 + p)
    edges = matrix(1L, p, p)
    diag(edges) = 0L
    # Remove ~40% of edges
    set.seed(150 + p)
    for(i in 1:(p - 1)) {
      for(j in (i + 1):p) {
        if(runif(1) < 0.4) {
          edges[i, j] = 0L
          edges[j, i] = 0L
        }
      }
    }

    x_raw = phi_to_full_position(dat$Phi)
    proj = ggm_test_project_position(x_raw, edges)
    x = as.vector(proj$x_projected)

    set.seed(160 + p)
    r = rnorm(length(x))
    r_proj = as.vector(ggm_test_project_momentum(r, x, edges))

    # Verify J * r_proj = 0
    Phi = full_position_to_phi(x, p)
    J = build_constraint_jacobian(Phi, edges)
    if(nrow(J) > 0) {
      Jr = J %*% r_proj
      expect_equal(as.vector(Jr), rep(0, nrow(J)),
        tolerance = 1e-11,
        info = paste("p =", p)
      )
    }

    # Verify idempotency
    r_proj2 = as.vector(ggm_test_project_momentum(r_proj, x, edges))
    expect_equal(r_proj2, r_proj,
      tolerance = 1e-13,
      info = paste("p =", p, "idempotency")
    )
  }
})


# ---- 7. Already-projected momentum is unchanged ---------------------------

test_that("momentum already in cotangent space is unchanged", {
  p = 5
  dat = make_test_phi(p, seed = 170)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 4] = 0L
  edges[4, 1] = 0L
  edges[2, 5] = 0L
  edges[5, 2] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  # Create a momentum that is already in the cotangent space
  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)

  # Start with random r, project in R, then verify C++ agrees
  set.seed(170)
  r = rnorm(length(x))
  G = J %*% t(J)
  lambda = solve(G, J %*% r)
  r_tangent = r - as.vector(t(J) %*% lambda)

  # C++ projection should leave it unchanged
  r_cpp = as.vector(ggm_test_project_momentum(r_tangent, x, edges))

  expect_equal(r_cpp, r_tangent, tolerance = 1e-12)
})
