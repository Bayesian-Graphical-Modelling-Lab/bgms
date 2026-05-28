# --------------------------------------------------------------------------- #
# Hierarchical-spec integration tests.
#
# `bgm(prior_factorization = "hierarchical")` is wired to the L-space
# Savage-Dickey between-edge MH step (closed-form Gaussian BF at α = 1,
# Gauss-Hermite quadrature at α > 1). End-to-end coverage via bgm(); the
# underlying SD primitive has its own unit tests + SBC sweeps under
# dev/sbc/.
# --------------------------------------------------------------------------- #



test_that("bgm() R API accepts prior_factorization = 'hierarchical' end-to-end", {
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
    prior_factorization = "hierarchical",
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
  # The SD between-edge step is orthogonal to the within-K sampler (AM or
  # NUTS), so the 2x2 cross-product should hold. Smoke-test the NUTS leg.
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
    prior_factorization = "hierarchical",
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


test_that("bgm() with Cauchy slab + hierarchical runs end-to-end", {
  # The Cauchy slab is supported via the scale-mixture-of-normals
  # representation (omega ~ IG(1/2, 1/2), K_ij | omega ~ N(0, sigma^2 omega))
  # with a per-edge omega slice-sampled by savage_dickey::slice_sample_cauchy_omega_active.
  # See math/savage_dickey/cauchy_omega.h and experiments/cauchy-slab/.
  set.seed(13)
  p <- 5L
  n <- 100L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  colnames(Y) <- paste0("V", seq_len(p))

  fit <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = cauchy_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    prior_factorization = "hierarchical",
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


test_that("Cauchy + hierarchical at alpha > 1 exercises the sinh primitive", {
  # alpha > 1 (diag prior shape) routes the SD between-step through the
  # sinh-midpoint quadrature primitive instead of the alpha = 1 closed
  # form. The Cauchy slab sigma^2 -> sigma^2 * omega_ij substitution
  # applies identically; this test guards against regressions in the
  # alpha > 1 dispatch under the Cauchy path.
  set.seed(7L)
  p <- 4L
  n <- 80L
  Y <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  colnames(Y) <- paste0("V", seq_len(p))

  fit <- bgm(
    Y, variable_type = "continuous",
    interaction_prior     = cauchy_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 3, rate = 1),
    prior_factorization = "hierarchical",
    iter = 200L, warmup = 50L,
    update_method = "adaptive-metropolis",
    chains = 1L, cores = 1L, seed = 1L,
    display_progress = "none", verbose = FALSE
  )
  ind <- S7::prop(fit, "posterior_mean_indicator")
  expect_true(all(is.finite(ind)))
  expect_true(all(ind >= 0 & ind <= 1))
  # Non-degenerate edge-count trajectory.
  raw <- S7::prop(fit, "raw_samples")
  ind_chn <- raw$indicator[[1L]]
  n_edges_path <- rowSums(ind_chn)
  expect_gt(length(unique(n_edges_path)), 1L)
})


test_that("Cauchy + hierarchical recovers a single true edge in q=4, n=120", {
  # End-to-end correctness: a clean q = 4 GGM with one strong edge
  # (K_12 = -0.5) should produce PIP near 1 on the true edge and PIPs
  # well below the prior 0.5 on the five null edges, even with a short
  # chain. This is a smoke test on the production routing through
  # cauchy_prior + hierarchical + the slice-sampled omega update.
  skip_if_not_installed("MASS")
  set.seed(42L)
  q <- 4L
  n <- 120L
  K_true <- diag(1, q)
  K_true[1, 2] <- K_true[2, 1] <- -0.5
  Y <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = solve(K_true))

  fit <- bgm(
    x = Y, variable_type = "continuous",
    interaction_prior     = cauchy_prior(scale = 0.5),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    prior_factorization   = "hierarchical",
    iter = 500L, warmup = 250L,
    update_method = "adaptive-metropolis",
    edge_selection = TRUE,
    chains = 1L, cores = 1L, seed = 11L,
    display_progress = "none", verbose = FALSE
  )
  pip <- S7::prop(fit, "posterior_mean_indicator")
  true_edge_pip <- pip[1, 2]
  null_pips <- c(pip[1, 3], pip[1, 4], pip[2, 3], pip[2, 4], pip[3, 4])

  expect_gt(true_edge_pip, 0.95)
  expect_true(all(null_pips < 0.4),
              info = sprintf("null PIPs: %s",
                              paste(sprintf("%.2f", null_pips), collapse = ", ")))
})


test_that("bgm() rejects hierarchical for non-continuous models", {
  set.seed(13)
  # Pure ordinal data should land in 'omrf' which has no precision matrix.
  X <- matrix(sample(0:3, 200 * 4L, replace = TRUE), 200, 4L)
  expect_error(
    bgm(X, variable_type = "ordinal",
        prior_factorization = "hierarchical",
        iter = 50L, warmup = 25L,
        update_method = "adaptive-metropolis",
        chains = 1L, cores = 1L, seed = 1L,
        display_progress = "none", verbose = FALSE),
    regexp = "continuous"
  )
})


