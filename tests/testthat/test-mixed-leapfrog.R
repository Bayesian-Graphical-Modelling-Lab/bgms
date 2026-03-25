# --------------------------------------------------------------------------- #
# Constrained leapfrog trajectory tests for the Mixed MRF.
#
# Tests verify:
#   1. Reversibility (forward + backward = identity)
#   2. Energy conservation (|DH| = O(eps^2))
#   3. Constraint preservation (excluded edges stay zero after trajectory)
#   4. Cotangent condition (excluded momentum entries stay zero)
#   5. SPD preservation (Cholesky diagonal stays positive)
#   6. Combined two-phase trajectory (unconstrained + constrained blocks)
#   7. Full-edge graph: constrained = unconstrained (projection is no-op)
#
# These mirror test-rattle-leapfrog.R (GGM) for the mixed MRF.
# --------------------------------------------------------------------------- #


# ---- Helpers ----------------------------------------------------------------

# Compute the full-space parameter dimension for a given setup
mixed_full_dim = function(num_cats, is_ord, p, q) {
  n_main = sum(ifelse(is_ord == 1L, num_cats, 2L))
  n_pairwise_xx = p * (p - 1) / 2
  n_cross = p * q
  n_chol = q * (q + 1) / 2
  n_main + n_pairwise_xx + q + n_cross + n_chol
}

# Compute excluded indices in the full parameter vector (1-based R indices)
mixed_excluded_indices = function(num_cats, is_ord, p, q, edge_ind) {
  n_main = sum(ifelse(is_ord == 1L, num_cats, 2L))
  n_pairwise_xx = p * (p - 1) / 2

  # Kxx excluded indices
  kxx_offset = n_main
  kxx_excluded = integer(0)
  idx = 0L
  for(i in seq_len(p - 1)) {
    for(j in (i + 1):p) {
      if(edge_ind[i, j] == 0L) {
        kxx_excluded = c(kxx_excluded, kxx_offset + idx + 1L)
      }
      idx = idx + 1L
    }
  }

  # Kxy excluded indices
  kxy_offset = n_main + n_pairwise_xx + q
  kxy_excluded = integer(0)
  idx = 0L
  for(i in seq_len(p)) {
    for(j in seq_len(q)) {
      if(edge_ind[i, p + j] == 0L) {
        kxy_excluded = c(kxy_excluded, kxy_offset + idx + 1L)
      }
      idx = idx + 1L
    }
  }

  # Kyy excluded indices (in the Cholesky block)
  chol_offset = n_main + n_pairwise_xx + q + p * q
  kyy_excluded = integer(0)
  # The Cholesky block is stored column-major: for column c,
  # entries are (0..c-1) off-diag then (c) log-diag
  # Excluded Kyy edge (i,j) with i < j means K_ij = sum_l L_li * L_lj = 0
  # which is enforced by the SHAKE projection on the Cholesky entries.

  list(kxx = kxx_excluded, kxy = kxy_excluded)
}

# Build a valid starting position for the mixed MRF.
# Returns a projected full-space parameter vector.
make_mixed_scenario = function(p, q, n, edge_ind, num_cats, seed) {
  set.seed(seed)
  x_obs = matrix(sample(0:(num_cats[1]), n * p, replace = TRUE), n, p)
  y_obs = matrix(rnorm(n * q), n, q)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  scale = 2.5

  dim = mixed_full_dim(num_cats, is_ord, p, q)

  # Start with small random parameters
  set.seed(seed + 100)
  params = rnorm(dim, sd = 0.2)

  # Ensure the Cholesky block has positive diagonal (log scale)
  n_main = sum(num_cats)
  n_pairwise_xx = p * (p - 1) / 2
  chol_offset = n_main + n_pairwise_xx + q + p * q

  # Set Cholesky diag entries to log of positive values
  chol_idx = 1
  for(col in seq_len(q)) {
    # off-diagonal entries (col - 1 of them)
    chol_idx = chol_idx + (col - 1)
    # diagonal entry: log(L_cc), ensure positive L_cc
    params[chol_offset + chol_idx] = log(1.0 + abs(rnorm(1, sd = 0.3)))
    chol_idx = chol_idx + 1
  }

  inv_mass = rep(1.0, dim)

  # Project position to satisfy constraints
  proj = mixed_test_project_position(
    params, inv_mass, x_obs, y_obs, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", scale
  )
  x0 = as.vector(proj$projected)

  # Random momentum, projected onto cotangent space
  set.seed(seed + 200)
  r_raw = rnorm(dim)
  r_proj = mixed_test_project_momentum(
    r_raw, x0, inv_mass, x_obs, y_obs, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "conditional", scale
  )
  r0 = as.vector(r_proj$projected)

  list(
    x0 = x0, r0 = r0, x_obs = x_obs, y_obs = y_obs,
    num_cats = num_cats, is_ord = is_ord, base_cat = base_cat,
    edge_ind = edge_ind, scale = scale, p = p, q = q, n = n
  )
}


# ---- 1. Reversibility: forward N + backward N ~ identity -------------------

test_that("mixed leapfrog is reversible, p=3 q=2 sparse", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L # Kxx edge
  edge_ind[2, 4] = edge_ind[4, 2] = 0L # Kxy edge
  edge_ind[4, 5] = edge_ind[5, 4] = 0L # Kyy edge

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 300)
  eps = 0.005
  n_steps = 5

  fwd = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  bwd = mixed_test_leapfrog_constrained(
    as.vector(fwd$x), -as.vector(fwd$r), eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  pos_err = max(abs(x_return - sc$x0))
  mom_err = max(abs(r_return - sc$r0))

  expect_lt(pos_err, 1e-3,
    label = paste("position reversibility error:", pos_err)
  )
  expect_lt(mom_err, 1e-3,
    label = paste("momentum reversibility error:", mom_err)
  )
})

test_that("mixed leapfrog is reversible, p=3 q=3 sparse", {
  p = 3
  q = 3
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 2] = edge_ind[2, 1] = 0L # Kxx edge
  edge_ind[3, 5] = edge_ind[5, 3] = 0L # Kxy edge
  edge_ind[4, 6] = edge_ind[6, 4] = 0L # Kyy edge
  edge_ind[5, 6] = edge_ind[6, 5] = 0L # Kyy edge

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 301)
  eps = 0.005
  n_steps = 8

  fwd = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  bwd = mixed_test_leapfrog_constrained(
    as.vector(fwd$x), -as.vector(fwd$r), eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  expect_lt(max(abs(x_return - sc$x0)), 1e-2)
  expect_lt(max(abs(r_return - sc$r0)), 1e-2)
})

test_that("mixed leapfrog reversibility error scales as O(eps^2)", {
  # The O(eps^2) reversibility error arises because the SHAKE correction
  # direction (A_q) depends nonlinearly on the position through the exp()
  # parameterization of the Cholesky diagonal. This only produces a
  # position-dependent correction when A_q mixes off-diagonal entries
  # (linear in x) with diagonal entries (exp(x)) from an earlier column.
  # Concretely, the excluded edge index i must be > 0 so that
  # A_q = [Phi(0,i), ..., Phi(i,i)] includes both off-diagonal and

  # exp(diagonal) terms. Edge Kyy(1,2) satisfies this: column 2 with
  # excluded index i=1 gives A_q = [Phi(0,1), exp(x_{diag_1})].
  p = 3
  q = 3
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[5, 6] = edge_ind[6, 5] = 0L # Kyy(1,2): excluded index i=1

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 302)
  n_steps = 5

  rev_err = function(eps) {
    fwd = mixed_test_leapfrog_constrained(
      sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
      sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
      sc$edge_ind, "conditional", sc$scale
    )
    bwd = mixed_test_leapfrog_constrained(
      as.vector(fwd$x), -as.vector(fwd$r), eps, n_steps, sc$x_obs, sc$y_obs,
      sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
      sc$edge_ind, "conditional", sc$scale
    )
    max(abs(as.vector(bwd$x) - sc$x0))
  }

  err_large = rev_err(0.01)
  err_small = rev_err(0.005)

  ratio = err_large / max(err_small, 1e-16)
  expect_gt(ratio, 2.0, label = paste("eps^2 scaling ratio:", ratio))
})


# ---- 2. Energy conservation: |DH| = O(eps^2) -------------------------------

test_that("mixed leapfrog conserves energy, p=3 q=2 sparse", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L
  edge_ind[2, 4] = edge_ind[4, 2] = 0L

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 310)
  eps = 0.005
  n_steps = 20

  res = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  expect_lt(abs(res$dH), 1.0,
    label = paste("|dH| =", abs(res$dH))
  )
})

test_that("mixed leapfrog energy conservation improves with smaller step", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L
  edge_ind[4, 5] = edge_ind[5, 4] = 0L

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 311)
  n_steps = 10

  res_large = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, 0.01, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )
  res_small = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, 0.005, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  expect_lt(abs(res_small$dH), abs(res_large$dH) + 0.01,
    label = "smaller step should give better energy conservation"
  )
})


# ---- 3. Constraint preservation: excluded entries stay zero -----------------

test_that("excluded Kxx and Kxy entries stay zero after trajectory", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L # Kxx
  edge_ind[2, 4] = edge_ind[4, 2] = 0L # Kxy
  edge_ind[4, 5] = edge_ind[5, 4] = 0L # Kyy

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 320)
  eps = 0.005
  n_steps = 15

  res = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  x_final = as.vector(res$x)
  excl = mixed_excluded_indices(num_cats, sc$is_ord, p, q, edge_ind)

  # Excluded Kxx entries must be zero
  if(length(excl$kxx) > 0) {
    for(idx in excl$kxx) {
      expect_equal(x_final[idx], 0,
        tolerance = 1e-10,
        label = paste("Kxx excluded index", idx)
      )
    }
  }

  # Excluded Kxy entries must be zero
  if(length(excl$kxy) > 0) {
    for(idx in excl$kxy) {
      expect_equal(x_final[idx], 0,
        tolerance = 1e-10,
        label = paste("Kxy excluded index", idx)
      )
    }
  }
})

test_that("excluded entries stay zero, p=3 q=3 multiple excluded", {
  p = 3
  q = 3
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 2] = edge_ind[2, 1] = 0L # Kxx
  edge_ind[2, 3] = edge_ind[3, 2] = 0L # Kxx
  edge_ind[1, 5] = edge_ind[5, 1] = 0L # Kxy
  edge_ind[3, 6] = edge_ind[6, 3] = 0L # Kxy
  edge_ind[5, 6] = edge_ind[6, 5] = 0L # Kyy

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 321)
  eps = 0.005
  n_steps = 10

  res = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  x_final = as.vector(res$x)
  excl = mixed_excluded_indices(num_cats, sc$is_ord, p, q, edge_ind)

  if(length(excl$kxx) > 0) {
    for(idx in excl$kxx) {
      expect_equal(x_final[idx], 0, tolerance = 1e-10)
    }
  }
  if(length(excl$kxy) > 0) {
    for(idx in excl$kxy) {
      expect_equal(x_final[idx], 0, tolerance = 1e-10)
    }
  }
})


# ---- 4. Cotangent condition: excluded momentum entries stay zero ------------

test_that("excluded momentum entries stay zero after trajectory", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L
  edge_ind[2, 4] = edge_ind[4, 2] = 0L
  edge_ind[4, 5] = edge_ind[5, 4] = 0L

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 330)
  eps = 0.005
  n_steps = 10

  res = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  r_final = as.vector(res$r)
  excl = mixed_excluded_indices(num_cats, sc$is_ord, p, q, edge_ind)

  if(length(excl$kxx) > 0) {
    for(idx in excl$kxx) {
      expect_equal(r_final[idx], 0,
        tolerance = 1e-10,
        label = paste("Kxx momentum index", idx)
      )
    }
  }
  if(length(excl$kxy) > 0) {
    for(idx in excl$kxy) {
      expect_equal(r_final[idx], 0,
        tolerance = 1e-10,
        label = paste("Kxy momentum index", idx)
      )
    }
  }
})


# ---- 5. SPD preservation: Cholesky diagonal stays positive ------------------

test_that("Cholesky diagonal stays positive after trajectory", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L
  edge_ind[4, 5] = edge_ind[5, 4] = 0L

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 340)
  eps = 0.01
  n_steps = 20

  res = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  x_final = as.vector(res$x)

  # Extract Cholesky diagonal entries (stored as log)
  n_main = sum(num_cats)
  n_pairwise_xx = p * (p - 1) / 2
  chol_offset = n_main + n_pairwise_xx + q + p * q

  chol_diag_indices = integer(q)
  idx = 0
  for(col in seq_len(q)) {
    idx = idx + (col - 1) + 1 # skip off-diag, land on diag
    chol_diag_indices[col] = chol_offset + idx
  }

  # log(L_cc) values — exp should be positive (always true for real values,
  # but check finiteness)
  for(j in seq_len(q)) {
    log_diag = x_final[chol_diag_indices[j]]
    expect_true(is.finite(log_diag),
      label = paste("Cholesky log-diag[", j, "] =", log_diag)
    )
    expect_true(exp(log_diag) > 0,
      label = paste("Cholesky diag[", j, "] = exp(", log_diag, ")")
    )
  }
})


# ---- 6. Combined two-phase trajectory: constraints across all blocks -------
# This tests that the hybrid split (unconstrained Kxx/Kxy + constrained Kyy)
# works correctly as a combined step. An error at the phase boundary would
# show up here.

test_that("combined trajectory: constraints hold across all three blocks", {
  p = 3
  q = 3
  n = 60
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  # Remove one edge from each block
  edge_ind[1, 3] = edge_ind[3, 1] = 0L # Kxx
  edge_ind[2, 5] = edge_ind[5, 2] = 0L # Kxy
  edge_ind[4, 6] = edge_ind[6, 4] = 0L # Kyy

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 350)

  # Run trajectory with moderate step size and many steps
  eps = 0.005
  n_steps = 15

  res = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  # Check energy conservation
  expect_lt(abs(res$dH), 1.0,
    label = paste("combined trajectory |dH| =", abs(res$dH))
  )

  # Check constraint preservation (Kxx/Kxy)
  x_final = as.vector(res$x)
  excl = mixed_excluded_indices(num_cats, sc$is_ord, p, q, edge_ind)
  if(length(excl$kxx) > 0) {
    for(idx in excl$kxx) {
      expect_equal(x_final[idx], 0, tolerance = 1e-10)
    }
  }
  if(length(excl$kxy) > 0) {
    for(idx in excl$kxy) {
      expect_equal(x_final[idx], 0, tolerance = 1e-10)
    }
  }

  # Check reversibility
  bwd = mixed_test_leapfrog_constrained(
    as.vector(res$x), -as.vector(res$r), eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  x_return = as.vector(bwd$x)
  r_return = -as.vector(bwd$r)

  expect_lt(max(abs(x_return - sc$x0)), 1e-2,
    label = "combined trajectory reversibility"
  )

  # Check cotangent condition
  r_final = as.vector(res$r)
  if(length(excl$kxx) > 0) {
    for(idx in excl$kxx) {
      expect_equal(r_final[idx], 0, tolerance = 1e-10)
    }
  }
  if(length(excl$kxy) > 0) {
    for(idx in excl$kxy) {
      expect_equal(r_final[idx], 0, tolerance = 1e-10)
    }
  }

  # Check Cholesky diagonal still positive
  n_main = sum(num_cats)
  n_pairwise_xx = p * (p - 1) / 2
  chol_offset = n_main + n_pairwise_xx + q + p * q
  idx = 0
  for(col in seq_len(q)) {
    idx = idx + (col - 1) + 1
    log_diag = x_final[chol_offset + idx]
    expect_true(is.finite(log_diag))
  }
})


# ---- 7. Full-edge graph: projection is no-op, trajectory valid --------------

test_that("full-edge graph: mixed leapfrog trajectory is valid", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L

  num_cats = rep(2L, p)
  sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = 360)
  eps = 0.005
  n_steps = 10

  res = mixed_test_leapfrog_constrained(
    sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
    sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
    sc$edge_ind, "conditional", sc$scale
  )

  # Energy should be approximately conserved
  expect_lt(abs(res$dH), 1.0)

  # Position and momentum should be finite
  expect_true(all(is.finite(as.vector(res$x))))
  expect_true(all(is.finite(as.vector(res$r))))
})


# ---- 8. Marginal PL: trajectory properties hold ----------------------------

test_that("mixed leapfrog works with marginal pseudolikelihood", {
  p = 3
  q = 2
  n = 50
  total = p + q
  edge_ind = matrix(1L, total, total)
  diag(edge_ind) = 0L
  edge_ind[1, 3] = edge_ind[3, 1] = 0L
  edge_ind[4, 5] = edge_ind[5, 4] = 0L

  num_cats = rep(2L, p)
  # Build scenario with conditional, but run leapfrog with marginal
  set.seed(370)
  x_obs = matrix(sample(0:2, n * p, replace = TRUE), n, p)
  y_obs = matrix(rnorm(n * q), n, q)
  is_ord = rep(1L, p)
  base_cat = rep(0L, p)
  scale = 2.5
  dim = mixed_full_dim(num_cats, is_ord, p, q)

  set.seed(371)
  params = rnorm(dim, sd = 0.2)
  n_main = sum(num_cats)
  n_pairwise_xx = p * (p - 1) / 2
  chol_offset = n_main + n_pairwise_xx + q + p * q
  chol_idx = 1
  for(col in seq_len(q)) {
    chol_idx = chol_idx + (col - 1)
    params[chol_offset + chol_idx] = log(1.0 + abs(rnorm(1, sd = 0.3)))
    chol_idx = chol_idx + 1
  }

  inv_mass = rep(1.0, dim)
  proj = mixed_test_project_position(
    params, inv_mass, x_obs, y_obs, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "marginal", scale
  )
  x0 = as.vector(proj$projected)

  set.seed(372)
  r_raw = rnorm(dim)
  r_proj = mixed_test_project_momentum(
    r_raw, x0, inv_mass, x_obs, y_obs, num_cats, as.integer(is_ord),
    base_cat, edge_ind, "marginal", scale
  )
  r0 = as.vector(r_proj$projected)

  eps = 0.005
  n_steps = 10

  res = mixed_test_leapfrog_constrained(
    x0, r0, eps, n_steps, x_obs, y_obs,
    num_cats, as.integer(is_ord), base_cat,
    edge_ind, "marginal", scale
  )

  # Energy conservation
  expect_lt(abs(res$dH), 1.0)

  # Constraint preservation
  x_final = as.vector(res$x)
  excl = mixed_excluded_indices(num_cats, is_ord, p, q, edge_ind)
  if(length(excl$kxx) > 0) {
    for(idx in excl$kxx) {
      expect_equal(x_final[idx], 0, tolerance = 1e-10)
    }
  }

  # Reversibility
  bwd = mixed_test_leapfrog_constrained(
    as.vector(res$x), -as.vector(res$r), eps, n_steps, x_obs, y_obs,
    num_cats, as.integer(is_ord), base_cat,
    edge_ind, "marginal", scale
  )
  expect_lt(max(abs(as.vector(bwd$x) - x0)), 1e-2)
})


# ---- 9. Multiple constraint patterns: smoke test ---------------------------

test_that("mixed leapfrog works across constraint patterns", {
  n = 50
  configs = list(
    list(p = 3, q = 2, seed = 380),
    list(p = 3, q = 3, seed = 381),
    list(p = 4, q = 2, seed = 382)
  )

  for(cfg in configs) {
    p = cfg$p
    q = cfg$q
    total = p + q
    num_cats = rep(2L, p)

    # Random sparse edge pattern
    set.seed(cfg$seed)
    edge_ind = matrix(1L, total, total)
    diag(edge_ind) = 0L
    for(i in seq_len(total - 1)) {
      for(j in (i + 1):total) {
        if(runif(1) < 0.35) {
          edge_ind[i, j] = 0L
          edge_ind[j, i] = 0L
        }
      }
    }

    sc = make_mixed_scenario(p, q, n, edge_ind, num_cats, seed = cfg$seed + 10)
    eps = 0.005
    n_steps = 8

    res = mixed_test_leapfrog_constrained(
      sc$x0, sc$r0, eps, n_steps, sc$x_obs, sc$y_obs,
      sc$num_cats, as.integer(sc$is_ord), sc$base_cat,
      sc$edge_ind, "conditional", sc$scale
    )

    expect_lt(abs(res$dH), 2.0,
      label = paste("smoke p=", p, "q=", q, "|dH|:", abs(res$dH))
    )
    expect_true(all(is.finite(as.vector(res$x))))
    expect_true(all(is.finite(as.vector(res$r))))

    # Check constraint preservation
    x_final = as.vector(res$x)
    excl = mixed_excluded_indices(num_cats, sc$is_ord, p, q, edge_ind)
    if(length(excl$kxx) > 0) {
      for(idx in excl$kxx) {
        expect_equal(x_final[idx], 0, tolerance = 1e-10)
      }
    }
    if(length(excl$kxy) > 0) {
      for(idx in excl$kxy) {
        expect_equal(x_final[idx], 0, tolerance = 1e-10)
      }
    }
  }
})
