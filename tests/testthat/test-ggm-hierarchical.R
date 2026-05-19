# --------------------------------------------------------------------------- #
# Phase 4a: hierarchical-spec MH integration smoke tests.
#
# Exercises the GGMModel between-edge MH path with graph_prior_spec_ set
# to GraphPriorSpec::Hierarchical: the joint MH ratio is multiplied by
# V(Γ_curr)/V(Γ_star) using the Phase 2 DEGORD sampler + Phase 3 V/RR
# estimator, and log_Z_NLO_curr is incremented on accept via Phase 1.
#
# These tests cover the wiring (chain runs, output shape is sensible,
# prior-family validation fires). Full SBC validation of the
# hierarchical-spec target is Phase 5.
# --------------------------------------------------------------------------- #


test_that("hierarchical-spec chain runs without crash and stays in valid state", {
  set.seed(7)
  p <- 5L
  n <- 100L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)

  out <- ggm_hierarchical_smoke_cpp(
    observations      = Y,
    inclusion_prob    = 0.5,
    interaction_scale = 1.0,
    diagonal_shape    = 1.0,
    diagonal_rate     = 1.0,
    delta             = 0.5,
    M_inner           = 100L,
    kappa             = 1.0,
    rho               = 0.5,
    n_sweeps          = 200L,
    seed              = 42L
  )

  expect_true(is.list(out))
  expect_true(is.matrix(out$final_edges))
  expect_equal(dim(out$final_edges), c(p, p))
  expect_true(isSymmetric(out$final_edges))
  # Edge counts in [0, p(p-1)/2] across the trajectory.
  max_edges <- p * (p - 1L) / 2L
  expect_true(all(out$n_edges_path >= 0L))
  expect_true(all(out$n_edges_path <= max_edges))
  # Some sweeps with non-trivial edge mass (chain isn't stuck at 0 or full).
  steady <- out$n_edges_path[101:200]
  expect_gt(length(unique(steady)), 1L)
})


test_that("hierarchical-spec chain is reproducible under fixed seed", {
  set.seed(13)
  p <- 4L; n <- 50L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  args <- list(
    observations = Y, inclusion_prob = 0.5,
    interaction_scale = 1.0,
    diagonal_shape = 1.0, diagonal_rate = 1.0,
    delta = 0.5, M_inner = 50L, kappa = 1.0, rho = 0.5,
    n_sweeps = 50L, seed = 99L
  )
  out_a <- do.call(ggm_hierarchical_smoke_cpp, args)
  out_b <- do.call(ggm_hierarchical_smoke_cpp, args)
  expect_identical(out_a$final_edges,  out_b$final_edges)
  expect_identical(out_a$n_edges_path, out_b$n_edges_path)
  # Different seed -> different trajectory (with overwhelming probability).
  args$seed <- 100L
  out_c <- do.call(ggm_hierarchical_smoke_cpp, args)
  expect_false(identical(out_a$n_edges_path, out_c$n_edges_path))
})


test_that("hierarchical-spec scales with delta sensibly", {
  # As δ increases, the |K|^δ tilt should push K further into the
  # interior of M+(Γ), making large connected graphs feasible. We don't
  # validate the SBC target here — only that the chain doesn't degenerate
  # to all-zero or all-one as δ varies. (Full SBC validation is Phase 5.)
  set.seed(91)
  p <- 4L; n <- 80L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  for (delta in c(0.0, 0.5, 1.0)) {
    out <- ggm_hierarchical_smoke_cpp(
      Y, 0.5, 1.0, 1.0, 1.0, delta, 50L, 1.0, 0.5, 100L, 17L
    )
    steady <- out$n_edges_path[51:100]
    expect_gt(length(unique(steady)), 1L,
              label = sprintf("steady is constant at delta=%g", delta))
  }
})
