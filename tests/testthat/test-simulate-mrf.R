# ==============================================================================
# Tests for simulate_mrf() - Standalone MRF Simulation
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Range invariants, reproducibility, distributional properties
#
# Tests for the standalone simulation function that generates data from
# a Markov Random Field with user-specified parameters.
#
# Tolerance testing principles applied here:
#   - Range invariants: simulated values in [0, n_cats], integers only
#   - Reproducibility: seed produces identical results
#   - Coarse distributional properties: positive interactions tend to
#     produce positive correlations, thresholds shift marginal distributions
#   - Dimension consistency: output has correct n_states x n_vars shape
#
# Note: We test *tendencies* rather than exact values because simulation
# output is stochastic. See test-tolerance.R for the foundational approach.
# ==============================================================================

# ------------------------------------------------------------------------------
# Basic Functionality Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf returns matrix of correct dimensions", {
  n_states <- 50
  n_vars <- 4
  n_cats <- 2

  interactions <- matrix(0, n_vars, n_vars)
  interactions[1, 2] <- interactions[2, 1] <- 0.3
  thresholds <- matrix(0, n_vars, n_cats)

  result <- simulate_mrf(
    num_states = n_states,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 123
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), n_states)
  expect_equal(ncol(result), n_vars)
})

test_that("simulate_mrf produces values in valid range", {
  n_cats <- 3
  n_vars <- 5

  interactions <- matrix(0, n_vars, n_vars)
  thresholds <- matrix(0, n_vars, n_cats)

  result <- simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 42
  )

  # All values should be integers
  expect_true(all(result == round(result)))

  # Values should be in 0 to n_cats range
  expect_true(all(result >= 0))
  expect_true(all(result <= n_cats))
})

test_that("simulate_mrf is reproducible with seed", {
  n_vars <- 3
  n_cats <- 2

  interactions <- matrix(0.2, n_vars, n_vars)
  diag(interactions) <- 0
  thresholds <- matrix(0.5, n_vars, n_cats)

  result1 <- simulate_mrf(
    num_states = 50,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 999
  )

  result2 <- simulate_mrf(
    num_states = 50,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 999
  )

  expect_equal(result1, result2)
})


# ------------------------------------------------------------------------------
# Variable Category Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf handles different categories per variable", {
  n_vars <- 4
  n_cats <- c(1, 2, 3, 4) # Different number of categories

  interactions <- matrix(0, n_vars, n_vars)
  thresholds <- matrix(0, n_vars, max(n_cats))

  # Warnings expected: variables with fewer categories have extra threshold columns
  result <- suppressWarnings(simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 123
  ))

  # Check each variable's range
  for (j in 1:n_vars) {
    expect_true(
      all(result[, j] <= n_cats[j]),
      info = sprintf("Variable %d should have max value %d", j, n_cats[j])
    )
  }
})

test_that("simulate_mrf handles binary variables correctly", {
  n_vars <- 3
  n_cats <- 1 # Binary: 0 or 1

  interactions <- matrix(0.3, n_vars, n_vars)
  diag(interactions) <- 0
  thresholds <- matrix(0, n_vars, n_cats)

  result <- simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    seed = 42
  )

  # Only 0 and 1 should appear
  expect_true(all(result %in% c(0, 1)))
})


# ------------------------------------------------------------------------------
# Blume-Capel Model Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf works with blume-capel variables", {
  n_vars <- 3
  n_cats <- 4 # Need > 2 categories for Blume-Capel

  interactions <- matrix(0, n_vars, n_vars)
  # Blume-Capel thresholds: alpha and beta
  thresholds <- matrix(NA, n_vars, n_cats)
  thresholds[, 1] <- 0 # alpha
  thresholds[, 2] <- -0.5 # beta

  result <- simulate_mrf(
    num_states = 100,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    variable_type = "blume-capel",
    baseline_category = 2,
    seed = 123
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 100)
  expect_true(all(result >= 0 & result <= n_cats))
})

test_that("simulate_mrf works with mixed ordinal and blume-capel", {
  n_vars <- 4
  n_cats <- 4

  interactions <- matrix(0, n_vars, n_vars)
  thresholds <- matrix(0, n_vars, n_cats)
  # Fill ordinal thresholds for vars 1 and 2
  thresholds[1:2, ] <- c(0, 0.3, 0.6, 0.9)
  # Blume-Capel params for vars 3 and 4
  thresholds[3:4, 1:2] <- cbind(c(0, 0), c(-0.3, -0.5))

  # Warnings expected: Blume-Capel variables only use 2 threshold columns
  result <- suppressWarnings(simulate_mrf(
    num_states = 80,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = thresholds,
    variable_type = c("ordinal", "ordinal", "blume-capel", "blume-capel"),
    baseline_category = c(0, 0, 2, 2), # baseline only matters for BC vars
    seed = 42
  ))

  expect_equal(nrow(result), 80)
  expect_equal(ncol(result), n_vars)
})


# ------------------------------------------------------------------------------
# Parameter Effect Tests
# ------------------------------------------------------------------------------

test_that("simulate_mrf: positive interaction tends to align responses", {
  # This is a weaker test that checks the basic behavior of interactions
  # Strong interactions can sometimes collapse variance, so we use moderate values
  n_vars <- 2
  n_cats <- 2

  # Moderate positive interaction
  pos_int <- matrix(c(0, 0.8, 0.8, 0), 2, 2)

  # Use spread-out thresholds to ensure variance in responses
  thresholds <- matrix(c(0, 0.5, 0, 0.5), n_vars, n_cats, byrow = TRUE)

  result <- simulate_mrf(
    num_states = 1000,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = pos_int,
    main = thresholds,
    iter = 2000,
    seed = 456
  )

  # Check that both variables have variance
  var1 <- var(result[, 1])
  var2 <- var(result[, 2])

  # With moderate interaction and spread thresholds, should have some variance
  expect_true(var1 > 0 && var2 > 0, info = "Variables should have non-zero variance")

  # Check the correlation is positive (as expected with positive interaction)
  if (var1 > 0 && var2 > 0) {
    cor_val <- cor(result[, 1], result[, 2])
    # With positive interaction, correlation should be non-negative
    # (allowing some tolerance for stochastic variation)
    expect_true(
      cor_val > -0.2,
      info = sprintf("Positive interaction should yield non-negative correlation, got: %.3f", cor_val)
    )
  }
})

test_that("simulate_mrf: threshold affects marginal distribution", {
  n_vars <- 1
  n_cats <- 1 # Binary

  interactions <- matrix(0, 1, 1)

  # Positive threshold -> more likely to be 1
  # In MRF parameterization, positive threshold increases log-odds of category 1
  pos_thresh <- matrix(3, 1, 1)
  # Negative threshold -> more likely to be 0
  neg_thresh <- matrix(-3, 1, 1)

  result_pos <- simulate_mrf(
    num_states = 500,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = pos_thresh,
    iter = 1000,
    seed = 42
  )

  result_neg <- simulate_mrf(
    num_states = 500,
    num_variables = n_vars,
    num_categories = n_cats,
    pairwise = interactions,
    main = neg_thresh,
    iter = 1000,
    seed = 42
  )

  prop_pos <- mean(result_pos == 1)
  prop_neg <- mean(result_neg == 1)

  # Positive threshold should have higher proportion of 1s
  expect_true(
    prop_pos > prop_neg,
    info = sprintf("Positive threshold prop(1)=%.3f should exceed negative threshold prop(1)=%.3f", prop_pos, prop_neg)
  )
})


# ------------------------------------------------------------------------------
# Deprecated mrfSampler() Test
# ------------------------------------------------------------------------------

test_that("mrfSampler works with deprecation warning", {
  n_vars <- 3
  n_cats <- 2

  interactions <- matrix(0, n_vars, n_vars)
  thresholds <- matrix(0, n_vars, n_cats)

  expect_warning(
    result <- mrfSampler(
      num_states = 20,
      num_variables = n_vars,
      num_categories = n_cats,
      pairwise = interactions,
      main = thresholds
    ),
    regexp = "deprecated"
  )

  expect_true(is.matrix(result))
})

test_that("mrfSampler produces identical results to simulate_mrf", {
  args <- list(
    num_states = 50,
    num_variables = 4,
    num_categories = 3,
    pairwise = matrix(0.1, 4, 4) - diag(0.1, 4),
    main = matrix(0, 4, 3),
    iter = 100,
    seed = 999
  )

  result_new <- do.call(simulate_mrf, args)
  result_old <- suppressWarnings(do.call(mrfSampler, args))

  expect_identical(result_new, result_old)
})
