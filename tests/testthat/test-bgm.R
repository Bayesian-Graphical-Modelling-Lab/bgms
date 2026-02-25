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


# ==============================================================================
# GGM Expanded Test Suite (Part D)
# ==============================================================================
#
# Tests for GGM correctness, convergence, and edge detection.
# See dev/plans/ggm_cleanup.md Part D for the design rationale.
# ==============================================================================


# --- D.1: Multi-chain convergence ---------------------------------------------

test_that("bgm GGM multi-chain produces valid Rhat", {
  skip_on_cran()

  p <- 5
  omega <- diag(p) * 2
  omega[1, 2] <- omega[2, 1] <- 0.5
  omega[2, 3] <- omega[3, 2] <- -0.4
  omega[3, 4] <- omega[4, 3] <- 0.3

  x <- simulate_mrf(
    num_states = 200, num_variables = p, pairwise = omega,
    variable_type = "continuous", seed = 101
  )
  colnames(x) <- paste0("V", 1:p)

  fit <- bgm(
    x, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 1000, warmup = 500, chains = 2,
    seed = 202, display_progress = "none"
  )

  # All Rhat values should be below 1.1 for converged chains
  rhat_main <- fit$posterior_summary_main$Rhat
  rhat_pair <- fit$posterior_summary_pairwise$Rhat

  expect_true(
    all(rhat_main < 1.1),
    info = sprintf("Max main Rhat = %.3f (expected < 1.1)", max(rhat_main))
  )
  expect_true(
    all(rhat_pair < 1.1),
    info = sprintf("Max pairwise Rhat = %.3f (expected < 1.1)", max(rhat_pair))
  )
})


# --- D.2: Sufficient statistics / MLE convergence -----------------------------

test_that("bgm GGM posterior mean approaches MLE for large n", {
  # For large n without edge selection, the posterior mean should approach
  # the sample precision matrix, since the likelihood dominates the prior.
  p <- 4
  omega_true <- diag(p)
  omega_true[1, 2] <- omega_true[2, 1] <- 0.4
  omega_true[3, 4] <- omega_true[4, 3] <- -0.3

  n <- 500
  x <- simulate_mrf(
    num_states = n, num_variables = p, pairwise = omega_true,
    variable_type = "continuous", seed = 42
  )
  colnames(x) <- paste0("V", 1:p)

  # MLE of the precision matrix = solve of sample covariance (centered)
  x_centered <- scale(x, center = TRUE, scale = FALSE)
  S <- crossprod(x_centered) / n
  mle_precision <- solve(S)

  fit <- bgm(
    x, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 500, warmup = 500, chains = 1,
    seed = 43, display_progress = "none"
  )

  # Reconstruct posterior mean precision
  omega_hat <- fit$posterior_mean_pairwise
  diag(omega_hat) <- as.numeric(fit$posterior_mean_main)

  # Posterior mean should correlate highly with MLE (likelihood dominates)
  cor_offdiag <- cor(
    omega_hat[lower.tri(omega_hat)],
    mle_precision[lower.tri(mle_precision)]
  )
  cor_diag <- cor(diag(omega_hat), diag(mle_precision))

  expect_true(
    cor_offdiag > 0.95,
    info = sprintf("Off-diagonal cor with MLE = %.3f (expected > 0.95)", cor_offdiag)
  )
  expect_true(
    cor_diag > 0.95,
    info = sprintf("Diagonal cor with MLE = %.3f (expected > 0.95)", cor_diag)
  )
})


# --- D.3: Missing data handling -----------------------------------------------

test_that("bgm GGM with listwise deletion drops rows correctly", {
  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) <- paste0("V", 1:4)

  # Introduce NAs
  x[5, 2] <- NA
  x[10, 3] <- NA
  x[20, 1] <- NA

  fit <- bgm(
    x, variable_type = "continuous",
    na_action = "listwise",
    edge_selection = FALSE,
    iter = 50, warmup = 100, chains = 1,
    seed = 44, display_progress = "none"
  )

  # 3 rows removed -> n = 47
  expect_equal(extract_arguments(fit)$num_cases, 47L)
  expect_s3_class(fit, "bgms")
})

test_that("bgm GGM rejects na_action = 'impute'", {
  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) <- paste0("V", 1:4)
  x[5, 2] <- NA

  expect_error(
    bgm(x, variable_type = "continuous", na_action = "impute",
        iter = 10, warmup = 10, chains = 1, display_progress = "none"),
    "Imputation is not yet supported"
  )
})

test_that("bgm GGM column_means are stored in arguments", {
  fit <- get_bgms_fit_ggm()
  args <- extract_arguments(fit)

  expect_true(!is.null(args$column_means))
  expect_true(is.numeric(args$column_means))
  expect_equal(length(args$column_means), args$num_variables)
})


# --- D.4: Larger p (Cholesky stability) ---------------------------------------

test_that("bgm GGM with p = 15 produces valid output", {
  skip_on_cran()

  p <- 15
  # Diagonally dominant precision matrix (ensures positive definiteness)
  omega <- diag(p) * 3
  set.seed(77)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      if (runif(1) < 0.2) {  # ~20% non-zero edges
        val <- runif(1, -0.3, 0.3)
        omega[i, j] <- omega[j, i] <- val
      }
    }
  }

  x <- simulate_mrf(
    num_states = 300, num_variables = p, pairwise = omega,
    variable_type = "continuous", seed = 88
  )
  colnames(x) <- paste0("V", 1:p)

  fit <- bgm(
    x, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 500, warmup = 500, chains = 1,
    seed = 99, display_progress = "none"
  )

  # All precision diagonals should be positive
  expect_true(
    all(fit$posterior_summary_main$mean > 0),
    info = "Some diagonal precision elements are non-positive"
  )

  # All values should be finite
  expect_true(all(is.finite(fit$posterior_summary_main$mean)))
  expect_true(all(is.finite(fit$posterior_summary_pairwise$mean)))

  # Correct dimensions
  n_edges <- p * (p - 1) / 2
  expect_equal(nrow(fit$posterior_summary_main), p)
  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
})


# --- D.5: Edge detection power ------------------------------------------------

test_that("bgm GGM edge selection discriminates true edges", {
  skip_on_cran()

  p <- 6
  omega_true <- diag(p) * 2
  # 5 true edges (sparse graph)
  omega_true[1, 2] <- omega_true[2, 1] <-  0.6
  omega_true[2, 3] <- omega_true[3, 2] <- -0.5
  omega_true[3, 4] <- omega_true[4, 3] <-  0.4
  omega_true[4, 5] <- omega_true[5, 4] <- -0.5
  omega_true[5, 6] <- omega_true[6, 5] <-  0.3

  # True adjacency
  adj_true <- (omega_true != 0)
  diag(adj_true) <- FALSE
  true_edges <- adj_true[lower.tri(adj_true)]

  x <- simulate_mrf(
    num_states = 500, num_variables = p, pairwise = omega_true,
    variable_type = "continuous", seed = 321
  )
  colnames(x) <- paste0("V", 1:p)

  fit <- bgm(
    x, variable_type = "continuous",
    edge_selection = TRUE,
    iter = 3000, warmup = 500, chains = 2,
    seed = 654, display_progress = "none"
  )

  # Posterior inclusion probabilities
  pip <- fit$posterior_mean_indicator[lower.tri(fit$posterior_mean_indicator)]

  # True edges should have higher PIPs than non-edges on average
  mean_pip_true  <- mean(pip[true_edges])
  mean_pip_false <- mean(pip[!true_edges])

  expect_true(
    mean_pip_true > mean_pip_false,
    info = sprintf(
      "Mean PIP for true edges (%.3f) should exceed non-edges (%.3f)",
      mean_pip_true, mean_pip_false
    )
  )

  # True edges should mostly have PIP > 0.5
  expect_true(
    mean(pip[true_edges] > 0.5) >= 0.6,
    info = sprintf(
      "Only %.0f%% of true edges have PIP > 0.5 (expected >= 60%%)",
      100 * mean(pip[true_edges] > 0.5)
    )
  )

  # Non-edges should mostly have PIP < 0.5
  expect_true(
    mean(pip[!true_edges] < 0.5) >= 0.6,
    info = sprintf(
      "Only %.0f%% of non-edges have PIP < 0.5 (expected >= 60%%)",
      100 * mean(pip[!true_edges] < 0.5)
    )
  )
})


# --- D.7: Conditional regression check ----------------------------------------

test_that("bgm GGM implied regression matches OLS for large n", {
  skip_on_cran()

  p <- 5
  omega_true <- diag(p) * 2
  omega_true[1, 2] <- omega_true[2, 1] <-  0.5
  omega_true[2, 3] <- omega_true[3, 2] <- -0.4
  omega_true[3, 4] <- omega_true[4, 3] <-  0.3
  omega_true[1, 5] <- omega_true[5, 1] <-  0.2

  n <- 1000
  x <- simulate_mrf(
    num_states = n, num_variables = p, pairwise = omega_true,
    variable_type = "continuous", seed = 111
  )
  colnames(x) <- paste0("V", 1:p)

  fit <- bgm(
    x, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 2000, warmup = 500, chains = 2,
    seed = 222, display_progress = "none"
  )

  # Reconstruct posterior mean precision matrix
  omega_hat <- fit$posterior_mean_pairwise
  diag(omega_hat) <- as.numeric(fit$posterior_mean_main)

  # For each variable j, the implied regression coefficients are:
  #   beta_j = -omega_{j,-j} / omega_{jj}
  # This should match OLS coefficients from the data.
  column_means <- extract_arguments(fit)$column_means
  x_centered <- sweep(x, 2, column_means)

  for (j in seq_len(p)) {
    rest <- setdiff(seq_len(p), j)

    # Implied regression from precision matrix
    beta_implied <- -omega_hat[rest, j] / omega_hat[j, j]

    # OLS regression on centered data
    ols_fit <- lm(x_centered[, j] ~ x_centered[, rest] - 1)
    beta_ols <- coef(ols_fit)

    # For large n, these should agree closely
    expect_true(
      cor(beta_implied, beta_ols) > 0.95,
      info = sprintf(
        "Variable %d: cor(beta_implied, beta_ols) = %.3f (expected > 0.95)",
        j, cor(beta_implied, beta_ols)
      )
    )
  }
})
