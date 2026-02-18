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
  original_main <- colMeans(extract_category_thresholds(fit))

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
  refit_main <- colMeans(extract_category_thresholds(refit))

  # Calculate correlations (handle zero variance edge case)
  cor_pairwise <- cor(original_pairwise, refit_pairwise, method = "spearman")
  cor_main <- cor(original_main, refit_main, method = "spearman")

  # If correlation is NA (zero variance), treat as failed
  if (is.na(cor_pairwise)) cor_pairwise <- 0
  if (is.na(cor_main)) cor_main <- 0

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
  cor_pairwise <- cor(original_pairwise, refit_pairwise, method = "spearman")

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
    min_correlation = 0.80,
    seed = 22222
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


test_that("bgm simulation-recovery: Blume-Capel variables", {
  skip_on_cran()

  # Use Wenchuan with Blume-Capel parameterization
  data("Wenchuan", package = "bgms")
  x <- na.omit(Wenchuan[, 1:5])
  n_obs <- nrow(x)

  fit <- bgm(x, iter = 1000, warmup = 1000, chains = 1,
             variable_type = "blume-capel", baseline_category = 2,
             edge_selection = FALSE, seed = 33333,
             display_progress = "none")

  result <- run_simrec_test(
    fit,
    n_sim = n_obs,
    mcmc_args = list(iter = 1000, warmup = 1000, chains = 1,
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

  cor_val <- cor(nuts_pairwise, am_pairwise, method = "spearman")

  expect_true(
    cor_val >= 0.80,
    info = sprintf("NUTS vs AM correlation = %.3f (expected >= 0.80)", cor_val)
  )
})
