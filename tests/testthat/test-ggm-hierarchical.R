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


test_that("bgm() with hierarchical errors helpfully for Cauchy slab", {
  set.seed(11)
  Y <- scale(matrix(rnorm(50 * 4L), 50, 4L), scale = FALSE)
  expect_error(
    bgm(Y, variable_type = "continuous",
        interaction_prior     = cauchy_prior(scale = 1),
        prior_factorization      = "hierarchical",
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
        prior_factorization = "hierarchical",
        iter = 50L, warmup = 25L,
        update_method = "adaptive-metropolis",
        chains = 1L, cores = 1L, seed = 1L,
        display_progress = "none", verbose = FALSE),
    regexp = "continuous"
  )
})


