# ==============================================================================
# Simulation-Recovery Tests (Correctness Tests)
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Self-consistency between estimation and simulation
#
# These tests verify that the estimation and simulation code are consistent:
#   1. Fit model on observed data → estimates A
#   2. Simulate new data from the fitted model
#   3. Refit model on simulated data → estimates B
#   4. Check: cor(A, B) should be high (model can reproduce its own structure)
#
# This approach:
#   - Does NOT require knowing "true" parameters
#   - Tests consistency between bgm()/bgmCompare() and simulate_mrf()/simulate.bgms()
#   - Detects bugs in likelihood, posterior, or simulation code
#
# These tests are computationally expensive and skipped on CRAN.
# ==============================================================================


# ------------------------------------------------------------------------------
# Helper Functions for Simulation-Recovery Tests
# ------------------------------------------------------------------------------

#' Run simulation-recovery test for a bgms fit
#'
#' @param fit A fitted bgms object
#' @param n_sim Number of observations to simulate (use >= 500 to avoid constant columns)
#' @param mcmc_args List of MCMC arguments for refitting
#' @param min_correlation Minimum acceptable correlation between estimates
#' @param seed Random seed for reproducibility
#'
#' @return List with correlation values and pass/fail status
run_simrec_test <- function(fit, n_sim = 350, mcmc_args = NULL,
                            min_correlation = 0.80, seed = 12345) {

  if (is.null(mcmc_args)) {
    mcmc_args <- list(iter = 1000, warmup = 1000, chains = 1,
                      display_progress = "none")
  }

  # Extract estimates from original fit
  original_pairwise <- colMeans(extract_pairwise_interactions(fit))
  original_main <- fit$posterior_summary_main$mean

  # Simulate data from the fitted model (use large n to avoid constant columns)
  set.seed(seed)
  simulated_data <- simulate(fit, nsim = n_sim, method = "posterior-mean",
                             seed = seed)

  # Validate: check for constant columns (would cause bgm to fail)
  # This can happen when the model predicts extreme probabilities for some categories
  col_vars <- apply(simulated_data, 2, function(x) length(unique(x)))
  if (any(col_vars < 2)) {
    # Return skipped result - model predictions are too extreme for this test
    return(list(
      cor_pairwise = NA_real_,
      cor_main = NA_real_,
      passed = NA,
      skipped = TRUE,
      reason = sprintf("Model produces degenerate predictions for variable(s): %s",
                       paste(which(col_vars < 2), collapse = ", "))
    ))
  }

  # Refit on simulated data
  args <- extract_arguments(fit)
  refit_args <- c(
    list(x = simulated_data, edge_selection = FALSE),
    mcmc_args
  )

  # Add variable_type if Blume-Capel
  if (any(args$variable_type == "blume-capel")) {
    refit_args$variable_type <- args$variable_type
    refit_args$baseline_category <- args$baseline_category
  }

  refit <- do.call(bgm, refit_args)

  # Extract estimates from refit
  refit_pairwise <- colMeans(extract_pairwise_interactions(refit))
  refit_main <- refit$posterior_summary_main$mean

  # Handle potential length mismatch in main parameters
  # (can happen when simulated data has fewer categories than original)
  n_main <- min(length(original_main), length(refit_main))
  original_main <- original_main[1:n_main]
  refit_main <- refit_main[1:n_main]

  # Calculate correlations
  # Use Pearson (not Spearman): with few parameters and many near zero,
  # Spearman rank correlation is dominated by noise in the ordering of
  # near-zero values, while Pearson correctly captures the linear agreement.
  cor_pairwise <- cor(original_pairwise, refit_pairwise)
  cor_main <- if (n_main >= 3) cor(original_main, refit_main) else NA_real_

  # If correlation is NA (zero variance or too few params), treat gracefully
  if (is.na(cor_pairwise)) cor_pairwise <- 0
  main_testable <- !is.na(cor_main)
  if (!main_testable) cor_main <- NA_real_

  list(
    cor_pairwise = cor_pairwise,
    cor_main = cor_main,
    passed = cor_pairwise >= min_correlation &&
             (!main_testable || cor_main >= min_correlation)
  )
}


#' Run simulation-recovery test for a GGM (continuous) bgms fit
#'
#' @param fit A fitted bgms object (GGM)
#' @param n_sim Number of observations to simulate
#' @param mcmc_args List of MCMC arguments for refitting
#' @param min_correlation Minimum acceptable correlation between estimates
#' @param seed Random seed for reproducibility
#'
#' @return List with correlation values and pass/fail status
run_simrec_test_ggm <- function(fit, n_sim = 500, mcmc_args = NULL,
                                min_correlation = 0.80, seed = 12345) {

  if(is.null(mcmc_args)) {
    mcmc_args <- list(iter = 1000, warmup = 1000, chains = 1,
                      display_progress = "none")
  }

  # Extract estimates from original fit
  # For GGM, pairwise contains full precision matrix (including diagonal)
  original_pairwise <- colMeans(extract_pairwise_interactions(fit))
  original_main <- diag(fit$posterior_mean_pairwise)

  # Simulate data from the fitted model
  set.seed(seed)
  simulated_data <- simulate(fit, nsim = n_sim, method = "posterior-mean",
                             seed = seed)

  # Refit on simulated data (must specify variable_type for GGM)
  refit_args <- c(
    list(x = simulated_data, variable_type = "continuous",
         edge_selection = FALSE),
    mcmc_args
  )

  refit <- do.call(bgm, refit_args)

  # Extract estimates from refit
  refit_pairwise <- colMeans(extract_pairwise_interactions(refit))
  refit_main <- diag(refit$posterior_mean_pairwise)

  # Calculate correlations
  cor_pairwise <- cor(original_pairwise, refit_pairwise)
  cor_main <- cor(original_main, refit_main)

  if(is.na(cor_pairwise)) cor_pairwise <- 0
  if(is.na(cor_main)) cor_main <- 0

  list(
    cor_pairwise = cor_pairwise,
    cor_main = cor_main,
    passed = cor_pairwise >= min_correlation && cor_main >= min_correlation
  )
}


#' Run simulation-recovery test for a bgmCompare fit
#'
#' @param fit A fitted bgmCompare object
#' @param n_per_group Number of observations per group to simulate (use >= 250)
#' @param mcmc_args List of MCMC arguments for refitting
#' @param min_correlation Minimum acceptable correlation
#' @param seed Random seed
#'
#' @return List with correlation values and pass/fail status
run_simrec_test_compare <- function(fit, n_per_group = 250, mcmc_args = NULL,
                                    min_correlation = 0.75, seed = 12345) {

  if (is.null(mcmc_args)) {
    mcmc_args <- list(iter = 1000, warmup = 1000, chains = 1,
                      display_progress = "none")
  }

  args <- extract_arguments(fit)
  n_groups <- args$num_groups

  # Extract baseline pairwise estimates
  original_pairwise <- colMeans(extract_pairwise_interactions(fit))

  # Simulate data for each group using group-specific parameters
  # For now, use baseline parameters (this is a simplification)
  interactions <- fit$posterior_mean_pairwise_baseline
  thresholds <- fit$posterior_mean_main_baseline

  set.seed(seed)
  simulated_datasets <- list()
  for (g in seq_len(n_groups)) {
    simulated_datasets[[g]] <- simulate_mrf(
      num_states = n_per_group,
      num_variables = args$num_variables,
      num_categories = args$num_categories,
      pairwise = interactions,
      main = thresholds,
      seed = seed + g
    )
    colnames(simulated_datasets[[g]]) <- args$data_columnnames
  }

  # Combine into single dataset with group indicator
  combined_data <- do.call(rbind, simulated_datasets)
  group_indicator <- rep(seq_len(n_groups), each = n_per_group)

  # Validate: check for constant columns (would cause bgmCompare to fail)
  col_vars <- apply(combined_data, 2, function(x) length(unique(x)))
  if (any(col_vars < 2)) {
    stop(sprintf("Simulated data has constant column(s): %s. Increase n_per_group or use different seed.",
                 paste(which(col_vars < 2), collapse = ", ")))
  }

  # Refit
  refit_args <- c(
    list(x = combined_data, group_indicator = group_indicator,
         difference_selection = FALSE),
    mcmc_args
  )

  refit <- do.call(bgmCompare, refit_args)

  # Extract estimates from refit
  refit_pairwise <- colMeans(extract_pairwise_interactions(refit))

  # Calculate correlation (handle zero variance edge case)
  cor_pairwise <- cor(original_pairwise, refit_pairwise)

  # If correlation is NA (zero variance), treat as failed
  if (is.na(cor_pairwise)) cor_pairwise <- 0

  list(
    cor_pairwise = cor_pairwise,
    passed = cor_pairwise >= min_correlation
  )
}


# ------------------------------------------------------------------------------
# bgm() Simulation-Recovery Tests
# ------------------------------------------------------------------------------
# These tests fit fresh models on larger datasets (matching data dimensions)
# rather than using the small session-cached fixtures.
# This takes longer but provides proper correctness validation.

test_that("bgm simulation-recovery: ordinal variables (NUTS)", {
  skip_on_cran()

  # Use full Wenchuan data
  data("Wenchuan", package = "bgms")
  x <- na.omit(Wenchuan[, 1:5])
  n_obs <- nrow(x)

  # Fit with adequate MCMC (1000 iter, 1000 warmup for proper convergence)
  fit <- bgm(x, iter = 1000, warmup = 1000, chains = 1,
             edge_selection = FALSE, seed = 11111,
             display_progress = "none")

  result <- run_simrec_test(
    fit,
    n_sim = n_obs,  # Match original sample size
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 11111
  )

  # Handle skipped case (model produces degenerate predictions)
  if (isTRUE(result$skipped)) {
    skip(result$reason)
  }

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
  expect_true(
    result$cor_main >= 0.80,
    info = sprintf("Main effects correlation = %.3f (expected >= 0.80)",
                   result$cor_main)
  )
})


test_that("bgm simulation-recovery: binary variables (NUTS)", {
  skip_on_cran()

  # Use full ADHD data
  data("ADHD", package = "bgms")
  x <- ADHD[, 2:6]
  n_obs <- nrow(x)

  fit <- bgm(x, iter = 1000, warmup = 1000, chains = 1,
             edge_selection = FALSE, seed = 22222,
             display_progress = "none")

  result <- run_simrec_test(
    fit,
    n_sim = n_obs,
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
                     display_progress = "none"),
    min_correlation = 0.75,
    seed = 22222
  )

  # Handle skipped case
  if (isTRUE(result$skipped)) {
    skip(result$reason)
  }

  expect_true(
    result$cor_pairwise >= 0.75,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.75)",
                   result$cor_pairwise)
  )
})


test_that("bgm simulation-recovery: Blume-Capel variables", {
  skip_on_cran()

  # Start from BC-simulated data to avoid category collapse:
  # Fitting BC on 5-category Wenchuan data then simulating produces only
  # 3 categories (the quadratic potential concentrates mass near baseline).
  # That makes the refit incomparable to the original. Using BC-simulated
  # data as the starting point keeps original and simulated data in the
  # same distributional regime.
  #
  # Uses a non-zero baseline_category to exercise the centering logic
  # in the OMRF C++ backend (observations_double_ centered around baseline).
  p <- 5
  n_obs <- 500
  pairwise <- matrix(0, p, p)
  pairwise[1, 2] <- pairwise[2, 1] <- 0.5
  pairwise[2, 3] <- pairwise[3, 2] <- 0.3
  pairwise[4, 5] <- pairwise[5, 4] <- -0.25

  main <- matrix(0, p, 2)
  main[, 1] <- c(-0.5, 0.0, 0.3, -0.2, 0.1)  # linear
  main[, 2] <- c(-0.3, -0.5, -0.4, -0.2, -0.6)  # quadratic

  x <- simulate_mrf(
    num_states = n_obs, num_variables = p,
    num_categories = rep(3, p),
    pairwise = pairwise, main = main,
    variable_type = "blume-capel",
    baseline_category = 1,
    seed = 33333
  )
  colnames(x) <- paste0("V", 1:p)

  fit <- bgm(x, iter = 5000, warmup = 1000, chains = 2,
             variable_type = "blume-capel", baseline_category = 1,
             edge_selection = FALSE, seed = 33333,
             display_progress = "none")

  result <- run_simrec_test(
    fit,
    n_sim = n_obs,
    mcmc_args = list(iter = 5000, warmup = 1000, chains = 2,
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 33333
  )

  # Handle skipped case
  if (isTRUE(result$skipped)) {
    skip(result$reason)
  }

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
})


test_that("bgm simulation-recovery: adaptive-metropolis", {
  skip_on_cran()

  # Use ADHD data with adaptive-metropolis sampler
  data("ADHD", package = "bgms")
  x <- ADHD[, 2:6]
  n_obs <- nrow(x)

  fit <- bgm(x, iter = 1000, warmup = 1000, chains = 1,
             update_method = "adaptive-metropolis",
             edge_selection = FALSE, seed = 44444,
             display_progress = "none")

  result <- run_simrec_test(
    fit,
    n_sim = n_obs,
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
                     update_method = "adaptive-metropolis",
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 44444
  )

  # Handle skipped case
  if (isTRUE(result$skipped)) {
    skip(result$reason)
  }

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
})


test_that("bgm simulation-recovery: GGM (continuous variables)", {
  skip_on_cran()

  # Generate continuous data from a known precision matrix
  p <- 5
  omega_true <- diag(p)
  omega_true[1, 2] <- omega_true[2, 1] <- 0.4
  omega_true[2, 3] <- omega_true[3, 2] <- 0.3
  omega_true[4, 5] <- omega_true[5, 4] <- -0.25
  omega_true[1, 3] <- omega_true[3, 1] <- 0.15

  n_obs <- 500
  x <- simulate_mrf(
    num_states = n_obs,
    num_variables = p,
    pairwise = omega_true,
    variable_type = "continuous",
    seed = 99999
  )
  colnames(x) <- paste0("V", 1:p)

  # Fit GGM
  fit <- bgm(x, variable_type = "continuous",
             iter = 1000, warmup = 1000, chains = 1,
             edge_selection = FALSE, seed = 99999,
             display_progress = "none")

  result <- run_simrec_test_ggm(
    fit,
    n_sim = n_obs,
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 99999
  )

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
  expect_true(
    result$cor_main >= 0.80,
    info = sprintf("Main effects (diagonal precision) correlation = %.3f (expected >= 0.80)",
                   result$cor_main)
  )
})


# ------------------------------------------------------------------------------
# bgmCompare() Simulation-Recovery Tests
# ------------------------------------------------------------------------------

test_that("bgmCompare simulation-recovery: ordinal variables", {
  skip_on_cran()

  # Use Boredom split into 2 groups
  data("Boredom", package = "bgms")
  x <- na.omit(Boredom[, 2:6])
  n_obs <- nrow(x)
  group_ind <- 1 * (Boredom[, 1] == "fr")

  fit <- bgmCompare(x, group_indicator = group_ind,
                    iter = 1000, warmup = 1000, chains = 1,
                    difference_selection = FALSE, seed = 55555,
                    display_progress = "none")

  result <- run_simrec_test_compare(
    fit,
    n_per_group = sum(group_ind),
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 55555
  )

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
})


test_that("bgmCompare simulation-recovery: binary variables", {
  skip_on_cran()

  # Use ADHD data with diagnosis group
  data("ADHD", package = "bgms")
  x <- ADHD[, 2:6]
  group_ind <- ADHD[, "group"]

  fit <- bgmCompare(x, group_indicator = group_ind,
                    iter = 1000, warmup = 1000, chains = 1,
                    difference_selection = FALSE, seed = 66666,
                    display_progress = "none")

  # Get group sizes for simulation
  n_per_group <- min(table(group_ind))

  result <- run_simrec_test_compare(
    fit,
    n_per_group = n_per_group,
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 66666
  )

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
})


test_that("bgmCompare simulation-recovery: adaptive-metropolis", {
  skip_on_cran()

  # Use ADHD data with adaptive-metropolis
  data("ADHD", package = "bgms")
  x <- ADHD[, 2:6]
  group_ind <- ADHD[, "group"]

  fit <- bgmCompare(x, group_indicator = group_ind,
                    iter = 1000, warmup = 1000, chains = 1,
                    update_method = "adaptive-metropolis",
                    difference_selection = FALSE, seed = 77777,
                    display_progress = "none")

  n_per_group <- min(table(group_ind))

  result <- run_simrec_test_compare(
    fit,
    n_per_group = n_per_group,
    mcmc_args = list(iter = 10000, warmup = 1000, chains = 1,
                     update_method = "adaptive-metropolis",
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 77777
  )

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
})


# ------------------------------------------------------------------------------
# Cross-Method Consistency Tests
# ------------------------------------------------------------------------------

test_that("NUTS and adaptive-metropolis produce consistent estimates",
{
  skip_on_cran()

  # Use larger dataset for meaningful comparison
  data("Wenchuan", package = "bgms")
  x <- na.omit(Wenchuan[, 1:5])

  # Fit with NUTS
  fit_nuts <- bgm(x, iter = 1000, warmup = 1000, chains = 1,
                  update_method = "nuts", edge_selection = FALSE,
                  seed = 88888, display_progress = "none")

  # Fit with adaptive-metropolis
  fit_am <- bgm(x, iter = 10000, warmup = 1000, chains = 1,
                update_method = "adaptive-metropolis", edge_selection = FALSE,
                seed = 88888, display_progress = "none")

  # Compare posterior means
  nuts_pairwise <- colMeans(extract_pairwise_interactions(fit_nuts))
  am_pairwise <- colMeans(extract_pairwise_interactions(fit_am))

  cor_val <- cor(nuts_pairwise, am_pairwise)

  expect_true(
    cor_val >= 0.80,
    info = sprintf("NUTS vs AM correlation = %.3f (expected >= 0.80)", cor_val)
  )
})


# ==============================================================================
# GGM Posterior Recovery Test
# ==============================================================================
#
# Verifies that the GGM sampler recovers known precision matrix parameters.

test_that("GGM posterior recovers parameters from simulated data", {

  n <- 1000
  p <- 10
  ne <- p * (p - 1) / 2

  # Fixed precision matrix (avoids BDgraph dependency)
  omega <- structure(c(6.240119, 0, 0, -0.370239, 0, 0, 0, 0, -1.622902,
              0, 0, 1.905013, 0, -0.194995, 0, 0, -2.468628, -0.557277, 0,
              0, 0, 0, 5.509142, -7.942389, 1.40081, 0, 0, -0.76775, 0, 0,
              -0.370239, -0.194995, -7.942389, 15.521405, -3.537489, 0, 4.60785,
              0, 3.278511, 0, 0, 0, 1.40081, -3.537489, 2.78257, 0, 0, 1.374641,
              0, -1.198092, 0, 0, 0, 0, 0, 1.350879, 0, 0.230677, -1.357952,
              0, 0, -2.468628, 0, 4.60785, 0, 0, 15.88698, 0, 1.20017, -1.973919,
              0, -0.557277, -0.76775, 0, 1.374641, 0.230677, 0, 7.007312, 1.597035,
              0, -1.622902, 0, 0, 3.278511, 0, -1.357952, 1.20017, 1.597035,
              13.378039, -4.769958, 0, 0, 0, 0, -1.198092, 0, -1.973919, 0,
              -4.769958, 5.536877), dim = c(10L, 10L))
  adj <- omega != 0
  diag(adj) <- 0
  covmat <- solve(omega)
  chol_cov <- chol(covmat)

  set.seed(43)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol_cov

  # Without edge selection
  fit_no_vs <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 3000, warmup = 500, chains = 2,
    display_progress = "none", seed = 42
  )

  expect_true(cor(fit_no_vs$posterior_summary_main$mean,     diag(omega)) > 0.9)
  expect_true(cor(fit_no_vs$posterior_summary_pairwise$mean, omega[lower.tri(omega)]) > 0.9)

  # With edge selection (Bernoulli prior)
  fit_vs <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = TRUE,
    iter = 5000, warmup = 500, chains = 2,
    display_progress = "none", seed = 42
  )

  expect_true(cor(fit_vs$posterior_summary_main$mean,     diag(omega)) > 0.9)
  expect_true(cor(fit_vs$posterior_summary_pairwise$mean, omega[lower.tri(omega)]) > 0.9)
  expect_true(cor(fit_vs$posterior_summary_indicator$mean, adj[lower.tri(adj)]) > 0.85)

  # With edge selection (SBM prior)
  fit_vs_sbm <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = TRUE,
    edge_prior = "Stochastic-Block",
    iter = 5000, warmup = 500, chains = 2,
    display_progress = "none", seed = 42
  )

  expect_true(cor(fit_vs_sbm$posterior_summary_main$mean,     diag(omega)) > 0.9)
  expect_true(cor(fit_vs_sbm$posterior_summary_pairwise$mean, omega[lower.tri(omega)]) > 0.9)
  expect_true(cor(fit_vs_sbm$posterior_summary_indicator$mean, adj[lower.tri(adj)]) > 0.85)

  # SBM-specific output
  expect_false(is.null(fit_vs_sbm$posterior_mean_coclustering_matrix))
  expect_equal(nrow(fit_vs_sbm$posterior_mean_coclustering_matrix), p)
  expect_equal(ncol(fit_vs_sbm$posterior_mean_coclustering_matrix), p)
  expect_false(is.null(fit_vs_sbm$posterior_num_blocks))
  expect_false(is.null(fit_vs_sbm$posterior_mode_allocations))
  expect_false(is.null(fit_vs_sbm$raw_samples$allocations))
})


# ==============================================================================
# Mixed MRF Simulation-Recovery
# ==============================================================================

#' Run simulation-recovery test for a mixed MRF bgms fit
#'
#' @param fit A fitted bgms object (mixed MRF)
#' @param n_sim Number of observations to simulate
#' @param mcmc_args List of MCMC arguments for refitting
#' @param min_correlation Minimum acceptable correlation between estimates
#' @param seed Random seed for reproducibility
#'
#' @return List with correlation values and pass/fail status
run_simrec_test_mixed <- function(fit, n_sim = 500, mcmc_args = NULL,
                                  min_correlation = 0.80, seed = 12345) {

  if(is.null(mcmc_args)) {
    mcmc_args <- list(iter = 1000, warmup = 1000, chains = 1,
                      display_progress = "none")
  }

  args <- extract_arguments(fit)

  # Extract pairwise estimates from original fit
  original_pairwise <- colMeans(extract_pairwise_interactions(fit))

  # Simulate data from the fitted model
  set.seed(seed)
  simulated_data <- simulate(fit, nsim = n_sim, method = "posterior-mean",
                             seed = seed)

  # Validate: check for constant discrete columns
  disc_idx <- args$discrete_indices
  col_vars <- apply(simulated_data[, disc_idx, drop = FALSE], 2,
                    function(x) length(unique(x)))
  if(any(col_vars < 2)) {
    return(list(
      cor_pairwise = NA_real_,
      passed = NA,
      skipped = TRUE,
      reason = sprintf("Degenerate discrete predictions for variable(s): %s",
                       paste(which(col_vars < 2), collapse = ", "))
    ))
  }

  # Refit on simulated data
  refit_args <- c(
    list(x = simulated_data,
         variable_type = args$variable_type,
         edge_selection = FALSE),
    mcmc_args
  )

  refit <- do.call(bgm, refit_args)

  # Extract estimates from refit
  refit_pairwise <- colMeans(extract_pairwise_interactions(refit))

  # Calculate correlation
  cor_pairwise <- cor(original_pairwise, refit_pairwise)
  if(is.na(cor_pairwise)) cor_pairwise <- 0

  list(
    cor_pairwise = cor_pairwise,
    passed = cor_pairwise >= min_correlation
  )
}


test_that("bgm simulation-recovery: mixed MRF (Boredom, ordinal + continuous)", {
  skip_on_cran()

  # Use Boredom data: 8 ordinal variables (1-7 scale), treat 3 as continuous.
  # Shift ordinal columns to 0-based to match bgm() expectations.
  data("Boredom", package = "bgms")
  x <- as.matrix(Boredom[, 2:6])
  x <- x - 1L  # shift from 1-7 to 0-6
  colnames(x) <- names(Boredom)[2:6]
  n_obs <- nrow(x)

  # Treat columns 2 and 4 as continuous, rest as ordinal
  vtype <- c("ordinal", "continuous", "ordinal", "continuous", "ordinal")

  fit <- bgm(x, variable_type = vtype,
             edge_selection = FALSE,
             iter = 1000, warmup = 1000, chains = 1,
             seed = 44321, display_progress = "none")

  result <- run_simrec_test_mixed(
    fit,
    n_sim = n_obs,
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
                     display_progress = "none"),
    min_correlation = 0.80,
    seed = 44321
  )

  if(isTRUE(result$skipped)) {
    skip(result$reason)
  }

  expect_true(
    result$cor_pairwise >= 0.80,
    info = sprintf("Mixed MRF pairwise correlation = %.3f (expected >= 0.80)",
                   result$cor_pairwise)
  )
})
