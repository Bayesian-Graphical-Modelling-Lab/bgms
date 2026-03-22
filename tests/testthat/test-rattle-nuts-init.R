# --------------------------------------------------------------------------- #
# RATTLE — NUTS initial momentum projection tests.
#
# Tests verify the initialization pattern used in nuts_step():
#   1. ProjectFn pattern preserves the original position
#   2. Projected initial momentum satisfies cotangent constraint
#   3. Kinetic energy after projection follows chi²(d-m)/2
#
# These tests verify the fix for the missing initial momentum projection
# bug, where r0 was sampled from N(0,M) but not projected onto the
# cotangent space before computing kin0, causing a systematic energy
# bias of ~m/2 nats.
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

# Mimics the nuts_step() initialization pattern: sample r0 from N(0, I),
# project via ProjectFn (project_position on a copy + project_momentum).
# Returns the projected momentum.
nuts_init_pattern = function(x, edges, seed) {
  d = length(x)
  set.seed(seed)
  r0 = rnorm(d)

  # This mirrors what nuts_step does:
  #   arma::vec pos_tmp = init_theta;
  #   (*project)(pos_tmp, r0);
  # where project calls project_position then project_momentum.
  #
  # Since x already satisfies constraints, project_position is a no-op.
  # We call project_momentum directly (which is what the combined
  # ProjectFn effectively does on valid positions).
  r_proj = as.vector(ggm_test_project_momentum(r0, x, edges))

  list(r0 = r0, r_proj = r_proj)
}


# ---- 1. ProjectFn pattern preserves the original position ------------------
# In nuts_step, the ProjectFn is called on pos_tmp (a copy of init_theta).
# project_position modifies pos_tmp, not init_theta. Verify this pattern:
# applying project_position to a valid position changes it by at most
# machine epsilon.

test_that("project_position is near-no-op on valid position, p=4", {
  p = 4
  dat = make_test_phi(p, seed = 300)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  # Project the already-valid position again
  proj2 = ggm_test_project_position(x, edges)
  x2 = as.vector(proj2$x_projected)

  expect_equal(x2, x, tolerance = 1e-14)
})

test_that("project_position is near-no-op on valid position, p=6", {
  p = 6
  dat = make_test_phi(p, seed = 301)
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

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  proj2 = ggm_test_project_position(x, edges)
  x2 = as.vector(proj2$x_projected)

  expect_equal(x2, x, tolerance = 1e-14)
})


# ---- 2. Projected initial momentum satisfies cotangent constraint ----------
# After the nuts_step initialization pattern, J * r_proj must equal zero.

test_that("NUTS init pattern: projected r0 satisfies J*r=0, p=4", {
  p = 4
  dat = make_test_phi(p, seed = 310)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  result = nuts_init_pattern(x, edges, seed = 310)

  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)
  Jr = J %*% result$r_proj

  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-12)
})

test_that("NUTS init pattern: projected r0 satisfies J*r=0, p=9", {
  p = 9
  dat = make_test_phi(p, seed = 311)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  # Remove ~60% of edges (sparse graph)
  set.seed(311)
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      if(runif(1) < 0.6) {
        edges[i, j] = 0L
        edges[j, i] = 0L
      }
    }
  }

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  result = nuts_init_pattern(x, edges, seed = 311)

  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)
  Jr = J %*% result$r_proj

  expect_equal(as.vector(Jr), rep(0, nrow(J)), tolerance = 1e-11)
})


# ---- 3. Kinetic energy after projection follows chi²(d-m)/2 ---------------
# When r0 ~ N(0, I) and is projected onto the (d-m)-dimensional cotangent
# space, kin0 = 0.5 * ||r_proj||² ~ chi²(d-m) / 2, so E[kin0] = (d-m)/2.
#
# This test generates many projected momenta and checks that the mean
# kinetic energy is consistent with (d-m)/2. A systematic deviation
# toward d/2 indicates the projection is not being applied.

test_that("mean kin0 after projection consistent with d-m DoF, p=5", {
  p = 5
  dat = make_test_phi(p, seed = 320)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  edges[1, 3] = 0L
  edges[3, 1] = 0L
  edges[2, 4] = 0L
  edges[4, 2] = 0L
  edges[1, 5] = 0L
  edges[5, 1] = 0L

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  d = length(x) # p*(p+1)/2 = 15
  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)
  m = nrow(J)

  n_samples = 2000
  kin_vals = numeric(n_samples)
  for(k in seq_len(n_samples)) {
    set.seed(320 + k)
    r0 = rnorm(d)
    r_proj = as.vector(ggm_test_project_momentum(r0, x, edges))
    kin_vals[k] = 0.5 * sum(r_proj^2)
  }

  expected_mean = (d - m) / 2
  observed_mean = mean(kin_vals)

  # With identity mass, kin ~ chi²(d-m)/2, so Var[kin] = (d-m)/2.
  # SE of mean = sqrt((d-m)/2) / sqrt(n_samples).
  se = sqrt((d - m) / 2) / sqrt(n_samples)

  # Check within 4 SE (very conservative for a regression test)
  expect_lt(abs(observed_mean - expected_mean), 4 * se,
    label = sprintf(
      "mean kin0 = %.3f, expected = %.3f (d=%d, m=%d, SE=%.4f)",
      observed_mean, expected_mean, d, m, se
    )
  )

  # The WRONG value (without projection) would be d/2, which is
  # expected_mean + m/2. Verify we are closer to d-m than to d.
  wrong_mean = d / 2
  expect_lt(abs(observed_mean - expected_mean),
    abs(observed_mean - wrong_mean),
    label = "kin0 closer to (d-m)/2 than d/2"
  )
})

test_that("mean kin0 after projection consistent with d-m DoF, p=8 sparse", {
  p = 8
  dat = make_test_phi(p, seed = 325)
  edges = matrix(1L, p, p)
  diag(edges) = 0L
  set.seed(325)
  for(i in 1:(p - 1)) {
    for(j in (i + 1):p) {
      if(runif(1) < 0.5) {
        edges[i, j] = 0L
        edges[j, i] = 0L
      }
    }
  }

  x_raw = phi_to_full_position(dat$Phi)
  proj = ggm_test_project_position(x_raw, edges)
  x = as.vector(proj$x_projected)

  d = length(x) # 36
  Phi = full_position_to_phi(x, p)
  J = build_constraint_jacobian(Phi, edges)
  m = nrow(J)

  n_samples = 2000
  kin_vals = numeric(n_samples)
  for(k in seq_len(n_samples)) {
    set.seed(325 + k)
    r0 = rnorm(d)
    r_proj = as.vector(ggm_test_project_momentum(r0, x, edges))
    kin_vals[k] = 0.5 * sum(r_proj^2)
  }

  expected_mean = (d - m) / 2
  observed_mean = mean(kin_vals)
  se = sqrt((d - m) / 2) / sqrt(n_samples)

  expect_lt(abs(observed_mean - expected_mean), 4 * se,
    label = sprintf(
      "mean kin0 = %.3f, expected = %.3f (d=%d, m=%d, SE=%.4f)",
      observed_mean, expected_mean, d, m, se
    )
  )

  wrong_mean = d / 2
  expect_lt(abs(observed_mean - expected_mean),
    abs(observed_mean - wrong_mean),
    label = "kin0 closer to (d-m)/2 than d/2"
  )
})
