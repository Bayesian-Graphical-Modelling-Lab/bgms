# ==============================================================================
# Core bgm() Function Tests
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Reproducibility, correlation with sufficient statistics
#
# These tests verify core bgm() functionality:
#   - Reproducibility: identical seeds produce identical MCMC chains
#   - Sanity: posterior means correlate with classical sufficient statistics
#
# INTEGRATION NOTE: Many sampler configurations (HMC, adaptive-metropolis,
# Blume-Capel, missing data imputation, standardization) are tested via the
# parameterized fixture approach in test-methods.R. See:
#   - helper-fixtures.R: Cached fit functions (get_bgms_fit_hmc, etc.)
#   - test-methods.R: get_bgms_fixtures() loops over all configurations
#
# This file focuses on tests that require special setup or unique assertions.
# ==============================================================================

test_that("bgm is reproducible", {
  # Use cached fixture as fit1, run one fresh fit as fit2 with same params
  fit1 <- get_bgms_fit_ordinal()
  
  data("Wenchuan", package = "bgms")
  fit2 <- bgm(
    Wenchuan[1:50, 1:4],
    iter = 50, warmup = 100, chains = 2,
    seed = 12345,  # Same seed as fixture
    display_progress = "none"
  )

  testthat::expect_equal(fit1$raw_samples, fit2$raw_samples)
})

test_that("bgmCompare is reproducible", {
  # Use cached fixture as fit1, run one fresh fit as fit2 with same params
  fit1 <- get_bgmcompare_fit_ordinal()
  
  data("Wenchuan", package = "bgms")
  x <- Wenchuan[1:50, 1:4]
  group_ind <- rep(1:2, each = 25)
  fit2 <- bgmCompare(
    x = x, group_indicator = group_ind,
    iter = 50, warmup = 100, chains = 2,
    seed = 54321,  # Same seed as fixture
    display_progress = "none"
  )

  combine_chains <- function(fit) {
    pairs <- do.call(rbind, fit$raw_samples$pairwise)
    mains <- do.call(rbind, fit$raw_samples$main)
    cbind(mains, pairs)
  }

  testthat::expect_equal(combine_chains(fit1), combine_chains(fit2))
})

# ==============================================================================
# GGM Reproducibility and Structure Tests
# ==============================================================================

test_that("bgm GGM is reproducible", {
  fit1 <- get_bgms_fit_ggm()

  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) <- paste0("V", 1:4)

  fit2 <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = TRUE,
    iter = 50, warmup = 100, chains = 1,
    seed = 44442,  # Same seed as fixture
    display_progress = "none"
  )

  testthat::expect_equal(fit1$raw_samples$main, fit2$raw_samples$main)
  testthat::expect_equal(fit1$raw_samples$pairwise, fit2$raw_samples$pairwise)
})

test_that("bgm GGM output has correct dimensions", {
  fit <- get_bgms_fit_ggm()
  args <- extract_arguments(fit)
  p <- args$num_variables

  # main: p diagonal precision elements
  expect_equal(nrow(fit$posterior_summary_main), p)
  expect_equal(nrow(fit$posterior_mean_main), p)

  # pairwise: p*(p-1)/2 off-diagonal elements
  n_edges <- p * (p - 1) / 2
  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
  expect_equal(nrow(fit$posterior_mean_pairwise), p)
  expect_equal(ncol(fit$posterior_mean_pairwise), p)

  # indicators (edge selection = TRUE)
  expect_equal(nrow(fit$posterior_summary_indicator), n_edges)
  expect_equal(nrow(fit$posterior_mean_indicator), p)
  expect_equal(ncol(fit$posterior_mean_indicator), p)

  # raw samples
  expect_equal(ncol(fit$raw_samples$main[[1]]), p)
  expect_equal(ncol(fit$raw_samples$pairwise[[1]]), n_edges)
  expect_equal(nrow(fit$raw_samples$main[[1]]), args$iter)
})

test_that("bgm GGM without edge selection omits indicators", {
  fit <- get_bgms_fit_ggm_no_es()

  expect_s3_class(fit, "bgms")
  expect_null(fit$posterior_summary_indicator)
  expect_null(fit$posterior_mean_indicator)
})

test_that("bgm GGM posterior precision diagonals are positive", {
  fit <- get_bgms_fit_ggm_no_es()

  expect_true(all(fit$posterior_summary_main$mean > 0))
})

# ==============================================================================
# HMC Reproducibility Test
# ==============================================================================
# HMC sampler basic functionality is covered by get_bgms_fit_hmc fixture in
# test-methods.R. This test specifically verifies reproducibility with seeds.

test_that("bgm with HMC is reproducible", {
  # Use cached fixture as fit1, run one fresh fit as fit2 with same params
  fit1 <- get_bgms_fit_hmc()
  
  data("Wenchuan", package = "bgms")
  fit2 <- bgm(
    Wenchuan[1:50, 1:4],
    update_method = "hamiltonian-mc",
    iter = 25, warmup = 50, chains = 1,
    seed = 55555,  # Same seed as fixture
    display_progress = "none"
  )

  testthat::expect_equal(fit1$raw_samples, fit2$raw_samples)
})


# ==============================================================================
# Parameter Ordering Tests (p >= 4 required to detect row/column-major bugs)
# ==============================================================================
#
# For a symmetric p x p matrix, the off-diagonal elements can be vectorized in
# row-major or column-major upper-triangle order. These orderings are identical
# for p <= 3. At p = 4, the 3rd element differs: row-major gives (1,4) while
# column-major gives (2,3). Using p = 5 ensures robust detection.
#
# See helper-fixtures.R for check_summary_matrix_consistency() and
# check_extractor_matrix_consistency().
# ==============================================================================

test_that("bgm GGM output has correct parameter ordering", {
  skip_on_cran()

  # Precision matrix with distinctive values at swap positions (p=5):
  #   Row-major position 3 = (1,4), column-major position 3 = (2,3)
  #   Row-major position 7 = (2,5), column-major position 7 = (3,4)
  # Zero vs non-zero at these positions makes any swap detectable.
  p <- 5
  omega <- diag(p) * 2
  omega[1, 2] <- omega[2, 1] <-  0.6
  omega[1, 3] <- omega[3, 1] <- -0.4
  omega[1, 4] <- omega[4, 1] <-  0.0    # swaps with V2-V3 under wrong ordering
  omega[1, 5] <- omega[5, 1] <-  0.3
  omega[2, 3] <- omega[3, 2] <-  0.0    # swaps with V1-V4 under wrong ordering
  omega[2, 4] <- omega[4, 2] <- -0.5
  omega[2, 5] <- omega[5, 2] <-  0.0    # swaps with V3-V4 under wrong ordering
  omega[3, 4] <- omega[4, 3] <-  0.25
  omega[3, 5] <- omega[5, 3] <-  0.0
  omega[4, 5] <- omega[5, 4] <- -0.35

  n <- 1000
  x <- simulate_mrf(
    num_states = n, num_variables = p, pairwise = omega,
    variable_type = "continuous", seed = 42
  )
  colnames(x) <- paste0("V", 1:p)

  fit <- bgm(
    x, variable_type = "continuous",
    iter = 500, warmup = 500, chains = 1,
    edge_selection = FALSE, seed = 42,
    display_progress = "none"
  )

  # Summary names -> matrix positions (pairwise)
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_pairwise,
      fit$posterior_mean_pairwise
    )),
    info = "GGM pairwise summary names do not match matrix positions"
  )

  # Extractor column means -> matrix positions
  pw_means <- colMeans(extract_pairwise_interactions(fit))
  expect_true(
    all(check_extractor_matrix_consistency(pw_means, fit$posterior_mean_pairwise)),
    info = "GGM extract_pairwise_interactions() names do not match matrix positions"
  )

  # Truth-based swap-position checks:
  # V1-V4 (true = 0) should be near zero, not ~0.25 (V3-V4's value)
  expect_true(
    abs(fit$posterior_mean_pairwise["V1", "V4"]) < 0.15,
    info = sprintf("V1-V4 should be ~0 but is %.3f (possible swap with V2-V3)",
                   fit$posterior_mean_pairwise["V1", "V4"])
  )
  # V2-V3 (true = 0) should be near zero, not ~0.6 (V1-V2's value)
  expect_true(
    abs(fit$posterior_mean_pairwise["V2", "V3"]) < 0.15,
    info = sprintf("V2-V3 should be ~0 but is %.3f (possible swap with V1-V4)",
                   fit$posterior_mean_pairwise["V2", "V3"])
  )
  # V3-V4 (true = 0.25) should NOT be near zero
  expect_true(
    fit$posterior_mean_pairwise["V3", "V4"] > 0.1,
    info = sprintf("V3-V4 should be ~0.25 but is %.3f (possible swap with V2-V5)",
                   fit$posterior_mean_pairwise["V3", "V4"])
  )
  # V2-V4 (true = -0.5) should be strongly negative
  expect_true(
    fit$posterior_mean_pairwise["V2", "V4"] < -0.3,
    info = sprintf("V2-V4 should be ~-0.5 but is %.3f",
                   fit$posterior_mean_pairwise["V2", "V4"])
  )
})

test_that("bgm OMRF output has correct parameter ordering", {
  skip_on_cran()

  data("Wenchuan", package = "bgms")
  x <- na.omit(Wenchuan[, 1:5])  # p=5 to detect row/column-major bugs

  fit <- bgm(
    x, iter = 1000, warmup = 500, chains = 1,
    edge_selection = TRUE, seed = 42,
    display_progress = "none"
  )

  # Summary names -> matrix positions (pairwise)
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_pairwise,
      fit$posterior_mean_pairwise
    )),
    info = "OMRF pairwise summary names do not match matrix positions"
  )

  # Extractor column means -> matrix positions
  pw_means <- colMeans(extract_pairwise_interactions(fit))
  expect_true(
    all(check_extractor_matrix_consistency(pw_means, fit$posterior_mean_pairwise)),
    info = "OMRF extract_pairwise_interactions() names do not match matrix positions"
  )

  # Indicator summary names -> matrix positions
  expect_true(
    all(check_summary_matrix_consistency(
      fit$posterior_summary_indicator,
      fit$posterior_mean_indicator
    )),
    info = "OMRF indicator summary names do not match matrix positions"
  )
})
