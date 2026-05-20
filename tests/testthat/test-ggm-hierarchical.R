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
  p <- 4L
  n <- 50L
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


test_that("bgm() R API accepts graph_prior_spec = 'hierarchical' end-to-end", {
  set.seed(99)
  p <- 5L
  n <- 100L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  colnames(Y) <- paste0("V", seq_len(p))

  fit <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 0.5,
    graph_prior_spec = "hierarchical",
    z_ratio_tuning   = list(M_inner = 50L, kappa = 1.0, rho = 0.5),
    iter = 200L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 1L,
    display_progress = "none", verbose = FALSE
  )
  ind <- S7::prop(fit, "posterior_mean_indicator")
  expect_true(is.matrix(ind))
  expect_equal(dim(ind), c(p, p))
  expect_true(all(ind >= 0 & ind <= 1))
  expect_true(all(is.finite(ind)))
})


test_that("bgm() accepts update_method = 'nuts' + 'hierarchical' (2x2 API)", {
  # F1. The hierarchical-spec hook lives in the between-edge MH update;
  # NUTS only governs the within-model continuous block. The 2x2
  # cross-product (NUTS/AMH x joint/hierarchical) is supposed to be a clean
  # cross-product, so smoke-test the NUTS leg here. Mirrors the AMH smoke
  # right above; differs only in update_method.
  set.seed(99)
  p <- 5L
  n <- 100L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  colnames(Y) <- paste0("V", seq_len(p))

  fit <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 0.5,
    graph_prior_spec = "hierarchical",
    z_ratio_tuning   = list(M_inner = 50L, kappa = 1.0, rho = 0.5),
    iter = 200L, warmup = 50L,
    update_method = "nuts",
    chains = 1L, cores = 1L, seed = 1L,
    display_progress = "none", verbose = FALSE
  )
  ind <- S7::prop(fit, "posterior_mean_indicator")
  expect_true(is.matrix(ind))
  expect_equal(dim(ind), c(p, p))
  expect_true(all(ind >= 0 & ind <= 1))
  expect_true(all(is.finite(ind)))
  # Edge-count trajectory should be non-degenerate (chain isn't locked at
  # the empty or full graph for the whole post-warmup window).
  raw <- S7::prop(fit, "raw_samples")
  ind_chn <- raw$indicator[[1L]]
  n_edges_path <- rowSums(ind_chn)
  expect_gt(length(unique(n_edges_path)), 1L)
  max_edges <- p * (p - 1L) / 2L
  expect_true(all(n_edges_path >= 0L & n_edges_path <= max_edges))
})


test_that("bgm() exposes v_ratio_diagnostics under hierarchical (F2)", {
  # F2 smoke. Hierarchical fit must surface per-iteration sign(V_curr) and
  # log|V_curr| as a top-level v_ratio_diagnostics field. Joint-spec fits
  # must NOT have that field. Operational cell (q=5, δ=0.5, κ=1, ρ=0.5)
  # is well inside the sign-positive regime, so we additionally assert
  # that the vast majority of recorded signs are +1.
  set.seed(101)
  p <- 5L
  n <- 100L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  fit_hier <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 0.5,
    graph_prior_spec = "hierarchical",
    z_ratio_tuning   = list(M_inner = 50L, kappa = 1.0, rho = 0.5),
    iter = 200L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 1L,
    display_progress = "none", verbose = FALSE
  )
  expect_true(!is.null(fit_hier$v_ratio_diagnostics))
  expect_named(fit_hier$v_ratio_diagnostics, c("sign", "log_abs"))
  expect_length(fit_hier$v_ratio_diagnostics$sign, 1L)
  expect_length(fit_hier$v_ratio_diagnostics$log_abs, 1L)
  s <- fit_hier$v_ratio_diagnostics$sign[[1L]]
  la <- fit_hier$v_ratio_diagnostics$log_abs[[1L]]
  expect_length(s, 200L)
  expect_length(la, 200L)
  expect_true(all(s %in% c(-1L, 1L)))
  expect_true(all(is.finite(la)))
  # Operational cell: sign should be +1 nearly always.
  expect_gt(mean(s == 1L), 0.95)

  # Joint-spec fit should not surface the diagnostic.
  fit_joint <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 0.5,
    graph_prior_spec = "joint",
    iter = 100L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 1L,
    display_progress = "none", verbose = FALSE
  )
  expect_null(fit_joint$v_ratio_diagnostics)
})


test_that("bgms_posterior_mean() applies sign-correction (F3/F4)", {
  # F3 sign-corrected posterior mean. In the operational cell, sign === +1
  # so the sign-corrected mean must coincide bit-by-bit with the plain
  # posterior_mean_* fields. For a joint-spec fit the function falls
  # through to the plain fields unchanged.
  set.seed(102)
  p <- 5L
  n <- 100L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)

  # Hierarchical fit in the operational cell - all signs should be +1.
  fit_hier <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 0.5,
    graph_prior_spec = "hierarchical",
    z_ratio_tuning   = list(M_inner = 50L, kappa = 1.0, rho = 0.5),
    iter = 200L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 1L,
    display_progress = "none", verbose = FALSE
  )
  stopifnot(all(fit_hier$v_ratio_diagnostics$sign[[1L]] == 1L))
  pm_hier <- bgms_posterior_mean(fit_hier)
  expect_named(pm_hier, c("pairwise", "main", "indicator", "residual_variance"))
  # Bit-equal to the plain mean when sign === +1.
  expect_equal(pm_hier$pairwise,          fit_hier$posterior_mean_pairwise,
               tolerance = 1e-12)
  expect_equal(pm_hier$indicator,         fit_hier$posterior_mean_indicator,
               tolerance = 1e-12)
  expect_equal(pm_hier$residual_variance, fit_hier$posterior_mean_residual_variance,
               tolerance = 1e-12)
  expect_null(pm_hier$main)

  # Joint-spec fit: function falls through to the plain fields unchanged.
  fit_joint <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 0.5,
    graph_prior_spec = "joint",
    iter = 200L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 1L,
    display_progress = "none", verbose = FALSE
  )
  pm_joint <- bgms_posterior_mean(fit_joint)
  expect_equal(pm_joint$pairwise,          fit_joint$posterior_mean_pairwise)
  expect_equal(pm_joint$indicator,         fit_joint$posterior_mean_indicator)
  expect_equal(pm_joint$residual_variance, fit_joint$posterior_mean_residual_variance)
})


test_that("bgms_posterior_mean() applies sign correction when signs flip", {
  # Synthetic check: take a real hierarchical fit, flip half of the
  # recorded signs by hand, and confirm bgms_posterior_mean() returns a
  # different (sign-weighted) mean. This validates the weighting math
  # without needing to engineer a tuning that actually produces -1 signs.
  set.seed(103)
  p <- 4L
  n <- 80L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  fit <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 0.5,
    graph_prior_spec = "hierarchical",
    z_ratio_tuning   = list(M_inner = 50L, kappa = 1.0, rho = 0.5),
    iter = 200L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 7L,
    display_progress = "none", verbose = FALSE
  )
  # Force every other iteration to sign = -1. bgms_posterior_mean reads
  # `fit` via `$`, which routes through S7::prop for the real fit but
  # works on plain lists too - construct a list view here so we can
  # mutate the sign vector without fighting S7 immutability.
  n_iter <- length(fit$v_ratio_diagnostics$sign[[1L]])
  # 2/3 +1, 1/3 -1 - non-trivial weighting that avoids sum-zero degeneracy.
  fake_signs <- rep_len(c(1L, 1L, -1L), n_iter)
  fit_mut <- list(
    v_ratio_diagnostics              = list(sign    = list(fake_signs),
                                            log_abs = fit$v_ratio_diagnostics$log_abs),
    raw_samples                      = fit$raw_samples,
    posterior_mean_main              = fit$posterior_mean_main,
    posterior_mean_pairwise          = fit$posterior_mean_pairwise,
    posterior_mean_indicator         = fit$posterior_mean_indicator,
    posterior_mean_residual_variance = fit$posterior_mean_residual_variance
  )

  pm_plain <- bgms_posterior_mean(fit)
  pm_signed <- bgms_posterior_mean(fit_mut)

  # Hand-compute expected sign-weighted pairwise (precision scale, then
  # mapped to association scale).
  pairwise_mat <- do.call(rbind, fit_mut$raw_samples$pairwise)
  expected_pair_precision <- colSums(pairwise_mat * fake_signs) / sum(fake_signs)
  expected_pair_assoc <- matrix(0, p, p,
                                dimnames = dimnames(fit$posterior_mean_pairwise))
  expected_pair_assoc[lower.tri(expected_pair_assoc)] <- expected_pair_precision
  expected_pair_assoc <- expected_pair_assoc + t(expected_pair_assoc)
  expected_pair_assoc <- -0.5 * expected_pair_assoc

  expect_equal(pm_signed$pairwise, expected_pair_assoc, tolerance = 1e-12)
  # And the signed result must differ from the plain mean (signs were
  # half-flipped, so it almost surely differs by more than numerical noise).
  expect_gt(max(abs(pm_signed$pairwise - pm_plain$pairwise)), 1e-6)
})


test_that("bgm() with hierarchical errors helpfully for Cauchy slab", {
  set.seed(11)
  Y <- scale(matrix(rnorm(50 * 4L), 50, 4L), scale = FALSE)
  expect_error(
    bgm(Y, variable_type = "continuous",
        interaction_prior     = cauchy_prior(scale = 1),
        graph_prior_spec      = "hierarchical",
        iter = 50L, warmup = 25L,
        update_method = "adaptive-metropolis",
        chains = 1L, cores = 1L, seed = 1L,
        display_progress = "none", verbose = FALSE),
    regexp = "Normal slab"
  )
})


test_that("bgm() rejects hierarchical for non-continuous models", {
  set.seed(13)
  # Pure ordinal data should land in 'omrf' which has no precision matrix.
  X <- matrix(sample(0:3, 200 * 4L, replace = TRUE), 200, 4L)
  expect_error(
    bgm(X, variable_type = "ordinal",
        graph_prior_spec = "hierarchical",
        iter = 50L, warmup = 25L,
        update_method = "adaptive-metropolis",
        chains = 1L, cores = 1L, seed = 1L,
        display_progress = "none", verbose = FALSE),
    regexp = "continuous"
  )
})


test_that("z_ratio_tuning validation rejects out-of-range values", {
  set.seed(11)
  Y <- scale(matrix(rnorm(50 * 4L), 50, 4L), scale = FALSE)
  for (bad in list(
    list(M_inner = 0L, kappa = 1, rho = 0.5),
    list(M_inner = 100, kappa = -1, rho = 0.5),
    list(M_inner = 100, kappa = 1, rho = 0),
    list(M_inner = 100, kappa = 1, rho = 1)
  )) {
    expect_error(
      bgm(Y, variable_type = "continuous",
          interaction_prior     = normal_prior(scale = 1),
          precision_scale_prior = gamma_prior(shape = 1, rate = 1),
          graph_prior_spec = "hierarchical",
          z_ratio_tuning   = bad,
          iter = 50L, warmup = 25L,
          update_method = "adaptive-metropolis",
          chains = 1L, cores = 1L, seed = 1L,
          display_progress = "none", verbose = FALSE),
      regexp = "z_ratio_tuning"
    )
  }
})


test_that("hierarchical-spec stays finite at p = 20 where linear c underflows", {
  # F5 regression. At p = 20 with δ = 1, log_Z_NLO is on the order of -500
  # nats, so c = κ · exp(log_Z_NLO) flushes to 0 in double precision. The
  # pre-F5 linear V_at_Gamma_pi_degord then evaluates to NaN / Inf, the MH
  # ratio is auto-rejected, and the chain locks at its initial state. The
  # log-space V must keep ln_alpha finite and let the chain move.
  skip_if_not_installed("MASS")
  set.seed(2026)
  p <- 20L
  n <- 200L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  fit <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta = 1.0,
    graph_prior_spec = "hierarchical",
    z_ratio_tuning   = list(M_inner = 50L, kappa = 1.0, rho = 0.5),
    iter = 100L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 2026L,
    display_progress = "none", verbose = FALSE
  )
  ind <- S7::prop(fit, "posterior_mean_indicator")
  expect_true(is.matrix(ind))
  expect_equal(dim(ind), c(p, p))
  expect_true(all(is.finite(ind)))
  expect_true(all(ind >= 0 & ind <= 1))
  # The chain should NOT be stuck at its initial state (all-zero edges).
  # If F5 regressed and every proposal auto-rejected, every off-diagonal
  # ind entry would be exactly 0.
  off <- ind[upper.tri(ind)]
  expect_gt(sum(off > 0), 0L)
})


test_that("hierarchical-spec scales with delta sensibly", {
  # As δ increases, the |K|^δ tilt should push K further into the
  # interior of M+(Γ), making large connected graphs feasible. We don't
  # validate the SBC target here — only that the chain doesn't degenerate
  # to all-zero or all-one as δ varies. (Full SBC validation is Phase 5.)
  set.seed(91)
  p <- 4L
  n <- 80L
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
