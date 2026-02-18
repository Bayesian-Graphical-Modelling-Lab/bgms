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
