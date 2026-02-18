# ==============================================================================
# Input Validation Tests
# ==============================================================================
#
# COMPLEMENTS: test-tolerance.R (tests the "failure" side of API contracts)
# PATTERN: Error conditions, boundary cases, invalid input handling
#
# While test-tolerance.R and other files test that VALID inputs produce
# outputs with correct properties, this file tests that INVALID inputs
# produce appropriate error messages.
#
# Together with tolerance tests, these form complete API contract testing:
#   - Tolerance tests: f(valid_input) -> output with correct properties
#   - Validation tests: f(invalid_input) -> informative error message
#
# See helper-fixtures.R for test data generators and testing philosophy.
# ==============================================================================

# ------------------------------------------------------------------------------
# bgm() Input Validation
# ------------------------------------------------------------------------------

test_that("bgm errors on non-matrix/data.frame input", {
  expect_error(bgm(x = 1:10), regexp = "matrix|data.frame|data frame")
  expect_error(bgm(x = list(a = 1, b = 2)), regexp = "matrix|data.frame|data frame")
})

test_that("bgm errors on data with too few variables", {
  bad_data <- matrix(c(0, 1, 0, 1), ncol = 1)
  expect_error(bgm(x = bad_data), regexp = "variable|column")
})

test_that("bgm errors on invalid iter values", {
  data <- generate_test_data(n = 20, p = 3)

  expect_error(bgm(x = data, iter = 0), regexp = "iter")
  expect_error(bgm(x = data, iter = -10), regexp = "iter")
  expect_error(bgm(x = data, iter = "100"), regexp = "iter|numeric")
})

test_that("bgm errors on invalid edge_prior", {
  data <- generate_test_data(n = 20, p = 3)

  expect_error(
    bgm(x = data, edge_selection = TRUE, edge_prior = "Invalid"),
    regexp = "should be one of|edge_prior"
  )
})

test_that("bgm errors on invalid na_action", {
  data <- generate_test_data(n = 20, p = 3)

  expect_error(
    bgm(x = data, na_action = "invalid_action"),
    regexp = "na_action"
  )
})

test_that("bgm errors on invalid update_method", {
  data <- generate_test_data(n = 20, p = 3)

  expect_error(
    bgm(x = data, update_method = "not-a-method"),
    regexp = "should be one of|update_method"
  )
})


# ------------------------------------------------------------------------------
# bgmCompare() Input Validation
# ------------------------------------------------------------------------------

test_that("bgmCompare errors on insufficient data", {
  # Too few groups
  data <- generate_test_data(n = 20, p = 3)
  group_ind <- rep(1, 20) # Only one group

  expect_error(
    bgmCompare(x = data, group_indicator = group_ind),
    regexp = "group"
  )
})

test_that("bgmCompare errors on mismatched group_indicator length", {
  data <- generate_test_data(n = 20, p = 3)
  group_ind <- rep(1:2, each = 5) # Only 10 elements for 20 rows

  expect_error(
    bgmCompare(x = data, group_indicator = group_ind),
    regexp = "group_indicator|length|match"
  )
})


# ------------------------------------------------------------------------------
# simulate_mrf() Input Validation
# ------------------------------------------------------------------------------

test_that("simulate_mrf errors on invalid num_states", {
  expect_error(
    simulate_mrf(
      num_states = 0,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 3, 3),
      main = matrix(0, 3, 2)
    ),
    regexp = "num_states"
  )

  expect_error(
    simulate_mrf(
      num_states = -5,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 3, 3),
      main = matrix(0, 3, 2)
    ),
    regexp = "num_states"
  )
})

test_that("simulate_mrf errors on non-symmetric pairwise", {
  non_sym <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), 3, 3)

  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = 3,
      num_categories = 2,
      pairwise = non_sym,
      main = matrix(0, 3, 2)
    ),
    regexp = "symmetric|pairwise"
  )
})

test_that("simulate_mrf errors on dimension mismatch", {
  # Interactions matrix wrong size
  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 4, 4), # Wrong: 4x4 for 3 variables
      main = matrix(0, 3, 2)
    ),
    regexp = "num_variables|dimension|size"
  )
})

test_that("simulate_mrf errors on missing thresholds", {
  expect_error(
    simulate_mrf(
      num_states = 10,
      num_variables = 3,
      num_categories = 2,
      pairwise = matrix(0, 3, 3),
      main = matrix(c(0, 0, NA, 0, 0, 0), 3, 2) # NA threshold
    ),
    regexp = "NA|threshold|missing"
  )
})


# ------------------------------------------------------------------------------
# predict.bgms() Input Validation
# ------------------------------------------------------------------------------

test_that("predict.bgms errors when newdata is missing", {
  fit <- get_bgms_fit()

  expect_error(predict(fit), regexp = "newdata")
})

test_that("predict.bgms errors on invalid variable names", {
  fit <- get_bgms_fit()
  args <- extract_arguments(fit)

  data("Wenchuan", package = "bgms")
  newdata <- Wenchuan[1:5, 1:args$num_variables]

  expect_error(
    predict(fit, newdata = newdata, variables = "NonexistentVar"),
    regexp = "not found|Variable"
  )
})

test_that("predict.bgms errors on out-of-range variable indices", {
  fit <- get_bgms_fit()
  args <- extract_arguments(fit)

  data("Wenchuan", package = "bgms")
  newdata <- Wenchuan[1:5, 1:args$num_variables]

  expect_error(
    predict(fit, newdata = newdata, variables = 999),
    regexp = "indices|between"
  )
})


# ------------------------------------------------------------------------------
# simulate.bgms() Input Validation
# ------------------------------------------------------------------------------

test_that("simulate.bgms errors on invalid seed", {
  fit <- get_bgms_fit()

  expect_error(simulate(fit, nsim = 10, seed = -1), regexp = "seed")
  expect_error(simulate(fit, nsim = 10, seed = "abc"), regexp = "seed")
})

test_that("simulate.bgms errors on invalid cores argument", {
  fit <- get_bgms_fit()

  expect_error(
    simulate(fit, nsim = 10, method = "posterior-sample", cores = 0),
    regexp = "cores"
  )
  expect_error(
    simulate(fit, nsim = 10, method = "posterior-sample", cores = "two"),
    regexp = "cores"
  )
})


# ------------------------------------------------------------------------------
# predict.bgmCompare() Input Validation
# ------------------------------------------------------------------------------

test_that("predict.bgmCompare errors when newdata is missing", {
  fit <- get_bgmcompare_fit()

  expect_error(predict(fit, group = 1), regexp = "newdata")
})

test_that("predict.bgmCompare errors when group is missing", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)

  newdata <- matrix(0L, nrow = 5, ncol = args$num_variables)

  expect_error(predict(fit, newdata = newdata), regexp = "group.*required")
})

test_that("predict.bgmCompare errors on invalid group argument", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)

  newdata <- matrix(0L, nrow = 5, ncol = args$num_variables)

  # Group out of range
  expect_error(
    predict(fit, newdata = newdata, group = 0),
    regexp = "group.*1"
  )
  expect_error(
    predict(fit, newdata = newdata, group = 999),
    regexp = "group"
  )
  # Invalid type
  expect_error(
    predict(fit, newdata = newdata, group = "a"),
    regexp = "group"
  )
})

test_that("predict.bgmCompare errors on invalid variable names", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)

  newdata <- matrix(0L, nrow = 5, ncol = args$num_variables)

  expect_error(
    predict(fit, newdata = newdata, group = 1, variables = "NonexistentVar"),
    regexp = "not found|Variable"
  )
})


# ------------------------------------------------------------------------------
# simulate.bgmCompare() Input Validation
# ------------------------------------------------------------------------------

test_that("simulate.bgmCompare errors when group is missing", {
  fit <- get_bgmcompare_fit()

  expect_error(simulate(fit, nsim = 10), regexp = "group.*required")
})

test_that("simulate.bgmCompare errors on invalid group argument", {
  fit <- get_bgmcompare_fit()

  # Group out of range
  expect_error(simulate(fit, nsim = 10, group = 0), regexp = "group.*1")
  expect_error(simulate(fit, nsim = 10, group = 999), regexp = "group")
  # Invalid type
  expect_error(simulate(fit, nsim = 10, group = "a"), regexp = "group")
})

test_that("simulate.bgmCompare errors on invalid seed", {
  fit <- get_bgmcompare_fit()

  expect_error(simulate(fit, nsim = 10, group = 1, seed = -1), regexp = "seed")
  expect_error(simulate(fit, nsim = 10, group = 1, seed = "abc"), regexp = "seed")
})


# ------------------------------------------------------------------------------
# Extractor Function Error Handling
# ------------------------------------------------------------------------------

test_that("extract_indicators errors correctly without edge selection", {
  data <- generate_test_data(n = 20, p = 3)
  args <- c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit <- do.call(bgm, args)

  expect_error(extract_indicators(fit), regexp = "edge_selection")
})

test_that("extract_posterior_inclusion_probabilities errors without edge selection", {
  data <- generate_test_data(n = 20, p = 3)
  args <- c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit <- do.call(bgm, args)

  expect_error(
    extract_posterior_inclusion_probabilities(fit),
    regexp = "edge_selection"
  )
})

test_that("extract_sbm errors with non-SBM prior", {
  fit <- get_bgms_fit()
  args <- extract_arguments(fit)

  # Skip if this fit actually used SBM
  if (identical(args$edge_prior, "Stochastic-Block")) {
    skip("Fit uses SBM prior")
  }

  expect_error(extract_sbm(fit), regexp = "Stochastic-Block")
})


# ------------------------------------------------------------------------------
# Edge Cases
# ------------------------------------------------------------------------------

test_that("bgm handles constant columns gracefully", {
  data <- generate_test_data(n = 20, p = 3)
  data[, 1] <- 0 # Make first column constant

  # This test verifies bgm doesn't crash unexpectedly on edge cases.
  # The function may error, warn, or succeed depending on implementation.
  result <- tryCatch(
    {
      args <- c(list(x = data), quick_mcmc_args())
      do.call(bgm, args)
    },
    error = function(e) list(type = "error", obj = e)
  )

  # Either it errors or it succeeds - both are acceptable behaviors
  if (is.list(result) && !is.null(result$type)) {
    # Got an error - bgm handled the edge case
    expect_true(TRUE, info = "bgm raised an error for constant column")
  } else {
    # If it succeeded, verify we got a valid bgms object
    expect_s3_class(result, "bgms")
  }
})

test_that("bgm handles all-NA columns", {
  data <- generate_test_data(n = 20, p = 3)
  data[, 1] <- NA # Make first column all NA

  expect_error(bgm(x = data), regexp = "NA|missing|incomplete")
})
