# --------------------------------------------------------------------------- #
# RATTLE Phase 1 — Position projection tests.
#
# Tests verify:
#   1. Full-position pack/unpack round-trips correctly
#   2. project_position enforces K_{iq} = 0 for excluded edges
#   3. Projection is idempotent (projecting twice = projecting once)
#   4. Projection preserves included edges (K_{iq} unchanged for ON edges)
#   5. Projection preserves diagonal entries
# --------------------------------------------------------------------------- #


# ---- Helper: build Phi from a PD precision matrix ---------------------------

make_test_phi = function(p, seed = 42) {
  set.seed(seed)
  # Random PD precision matrix
  A = matrix(rnorm(p * p), p, p)
  K = A %*% t(A) + diag(p)
  Phi = chol(K) # upper-triangular: K = t(Phi) %*% Phi
  list(Phi = Phi, K = K, p = p)
}


# ---- Helper: pack Phi into full-position vector ----------------------------

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


# ---- Helper: unpack full-position vector into Phi ---------------------------

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


# ---- 1. Pack/unpack round-trip ---------------------------------------------

test_that("get_full_position packs Phi correctly", {
  dat = make_test_phi(p = 4, seed = 1)
  x_r = phi_to_full_position(dat$Phi)
  x_cpp = as.vector(ggm_test_get_full_position(dat$Phi, matrix(1L, 4, 4)))
  expect_equal(x_cpp, x_r, tolerance = 1e-14)
})

test_that("full-position round-trip preserves Phi", {
  dat = make_test_phi(p = 5, seed = 2)
  x = phi_to_full_position(dat$Phi)
  Phi_roundtrip = full_position_to_phi(x, dat$p)
  expect_equal(Phi_roundtrip, dat$Phi, tolerance = 1e-14)
})


# ---- 2. Projection enforces zero constraints --------------------------------

test_that("project_position enforces K_{iq} = 0 for excluded edges (p=4)", {
  dat = make_test_phi(p = 4, seed = 10)
  x = phi_to_full_position(dat$Phi)

  # Exclude edges (1,3) and (2,4) — R 1-indexed
  edge_ind = matrix(1L, 4, 4)
  edge_ind[1, 3] = 0L
  edge_ind[3, 1] = 0L
  edge_ind[2, 4] = 0L
  edge_ind[4, 2] = 0L

  result = ggm_test_project_position(x, edge_ind)
  K = result$K

  expect_lt(abs(K[1, 3]), 1e-12)
  expect_lt(abs(K[3, 1]), 1e-12)
  expect_lt(abs(K[2, 4]), 1e-12)
  expect_lt(abs(K[4, 2]), 1e-12)
})

test_that("project_position enforces zeros for many excluded edges (p=6)", {
  dat = make_test_phi(p = 6, seed = 20)
  x = phi_to_full_position(dat$Phi)

  # Sparse graph: only nearest-neighbor edges
  edge_ind = diag(1L, 6)
  for(i in 1:5) {
    edge_ind[i, i + 1] = 1L
    edge_ind[i + 1, i] = 1L
  }

  result = ggm_test_project_position(x, edge_ind)
  K = result$K

  # All non-adjacent pairs should have K_{ij} = 0
  for(i in 1:6) {
    for(j in 1:6) {
      if(i != j && edge_ind[i, j] == 0L) {
        expect_lt(abs(K[i, j]), 1e-10,
          label = sprintf("K[%d,%d] should be zero", i, j)
        )
      }
    }
  }
})

test_that("project_position enforces zeros (p=8, dense exclusions)", {
  dat = make_test_phi(p = 8, seed = 30)
  x = phi_to_full_position(dat$Phi)

  # Remove ~60% of edges
  set.seed(99)
  edge_ind = diag(1L, 8)
  for(i in 1:7) {
    for(j in (i + 1):8) {
      if(runif(1) < 0.4) {
        edge_ind[i, j] = 1L
        edge_ind[j, i] = 1L
      }
    }
  }

  result = ggm_test_project_position(x, edge_ind)
  K = result$K

  for(i in 1:8) {
    for(j in 1:8) {
      if(i != j && edge_ind[i, j] == 0L) {
        expect_lt(abs(K[i, j]), 1e-10,
          label = sprintf("K[%d,%d] should be zero", i, j)
        )
      }
    }
  }
})


# ---- 3. Projection is idempotent -------------------------------------------

test_that("projecting twice gives same result as projecting once", {
  dat = make_test_phi(p = 5, seed = 40)
  x = phi_to_full_position(dat$Phi)

  edge_ind = diag(1L, 5)
  edge_ind[1, 2] = 1L
  edge_ind[2, 1] = 1L
  edge_ind[2, 3] = 1L
  edge_ind[3, 2] = 1L
  # Edges (1,3), (1,4), (1,5), (2,4), (2,5), (3,4), (3,5), (4,5) are excluded

  result1 = ggm_test_project_position(x, edge_ind)
  result2 = ggm_test_project_position(result1$x_projected, edge_ind)

  expect_equal(result2$x_projected, result1$x_projected, tolerance = 1e-12)
})


# ---- 4. Projection preserves included edges --------------------------------

test_that("included edge K entries are not zeroed by projection", {
  dat = make_test_phi(p = 5, seed = 50)
  x = phi_to_full_position(dat$Phi)

  # All edges on
  all_on = matrix(1L, 5, 5)
  result_all = ggm_test_project_position(x, all_on)

  # No projection should happen (no constraints)
  expect_equal(as.vector(result_all$x_projected), x, tolerance = 1e-14)
})


# ---- 5. Projection preserves diagonals -------------------------------------

test_that("diagonal entries are unchanged by projection", {
  dat = make_test_phi(p = 6, seed = 60)
  x = phi_to_full_position(dat$Phi)

  # Sparse graph
  edge_ind = diag(1L, 6)
  edge_ind[1, 2] = 1L
  edge_ind[2, 1] = 1L

  result = ggm_test_project_position(x, edge_ind)
  Phi_orig = dat$Phi
  Phi_proj = result$Phi

  for(q in 1:6) {
    expect_equal(Phi_proj[q, q], Phi_orig[q, q],
      tolerance = 1e-14,
      label = sprintf("Phi[%d,%d] diagonal preserved", q, q)
    )
  }
})


# ---- 6. Projection from a perturbed point ----------------------------------

test_that("projection recovers valid K from perturbed Cholesky entries", {
  dat = make_test_phi(p = 4, seed = 70)
  x = phi_to_full_position(dat$Phi)

  # Exclude edge (1,4) and (2,3)
  edge_ind = matrix(1L, 4, 4)
  edge_ind[1, 4] = 0L
  edge_ind[4, 1] = 0L
  edge_ind[2, 3] = 0L
  edge_ind[3, 2] = 0L

  # Perturb the position vector
  set.seed(71)
  x_perturbed = x + rnorm(length(x), sd = 0.1)

  result = ggm_test_project_position(x_perturbed, edge_ind)
  K = result$K

  expect_lt(abs(K[1, 4]), 1e-10)
  expect_lt(abs(K[2, 3]), 1e-10)

  # K should still be symmetric
  expect_equal(K, t(K), tolerance = 1e-14)
})
