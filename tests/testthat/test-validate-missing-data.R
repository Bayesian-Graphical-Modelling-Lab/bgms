# ==============================================================================
# Unit tests for validate_missing_data()
# Phase A.5 of the R scaffolding refactor.
# ==============================================================================

# ==============================================================================
# 1. Listwise — no NAs
# ==============================================================================

test_that("listwise with no NAs returns data unchanged", {
  x <- matrix(1:12, nrow = 4, ncol = 3)
  result <- validate_missing_data(x, na_action = "listwise")
  expect_equal(result$x, x)
  expect_false(result$na_impute)
  expect_equal(result$missing_index, matrix(NA, 1, 1))
  expect_equal(result$n_removed, 0L)
})

# ==============================================================================
# 2. Listwise — some NAs removed
# ==============================================================================

test_that("listwise removes rows with NAs and messages", {
  x <- matrix(1:12, nrow = 4, ncol = 3)
  x[2, 1] <- NA
  x[4, 3] <- NA

  withr::with_options(list(bgms.verbose = TRUE), {
    expect_message(
      result <- validate_missing_data(x, na_action = "listwise"),
      "2 rows with missing values excluded"
    )
  })
  expect_equal(nrow(result$x), 2)
  expect_false(result$na_impute)
  expect_equal(result$n_removed, 2L)
})

test_that("listwise suppresses message when bgms.verbose is FALSE", {
  x <- matrix(1:12, nrow = 4, ncol = 3)
  x[2, 1] <- NA
  withr::with_options(list(bgms.verbose = FALSE), {
    expect_silent(
      result <- validate_missing_data(x, na_action = "listwise")
    )
  })
  expect_equal(nrow(result$x), 3)
})

# ==============================================================================
# 3. Listwise — all rows have NAs
# ==============================================================================

test_that("listwise errors when all rows have NAs", {
  x <- matrix(NA_real_, nrow = 3, ncol = 2)
  expect_error(
    validate_missing_data(x, na_action = "listwise"),
    "All rows in x contain at least one missing response"
  )
})

# ==============================================================================
# 4. Listwise — too few rows after removal
# ==============================================================================

test_that("listwise errors when < 2 rows remain", {
  x <- matrix(1:6, nrow = 3, ncol = 2)
  x[1, 1] <- NA
  x[2, 2] <- NA
  withr::with_options(list(bgms.verbose = FALSE), {
    expect_error(
      validate_missing_data(x, na_action = "listwise"),
      "less than two rows"
    )
  })
})

# ==============================================================================
# 5. Impute — no NAs
# ==============================================================================

test_that("impute with no NAs returns data unchanged", {
  x <- matrix(1:12, nrow = 4, ncol = 3)
  result <- validate_missing_data(x, na_action = "impute")
  expect_equal(result$x, x)
  expect_false(result$na_impute)
  expect_equal(result$missing_index, matrix(NA, 1, 1))
})

# ==============================================================================
# 6. Impute — with NAs
# ==============================================================================

test_that("impute fills NAs and builds missing_index", {
  set.seed(42)
  x <- matrix(c(1, 2, 3, 4, 5, NA, 7, 8, 9), nrow = 3, ncol = 3)
  result <- validate_missing_data(x, na_action = "impute")

  expect_true(result$na_impute)
  expect_false(anyNA(result$x))
  expect_equal(nrow(result$missing_index), 1)
  # 0-based indices: row 2 (0-based: 1), col 1 (0-based: 1)
  expect_equal(result$missing_index[1, 1], 2)  # row index (0-based)
  expect_equal(result$missing_index[1, 2], 1)  # col index (0-based)
})

test_that("impute builds correct missing_index for multiple NAs", {
  set.seed(123)
  x <- matrix(c(1, NA, 3, NA, 5, 6, 7, 8, NA), nrow = 3, ncol = 3)
  result <- validate_missing_data(x, na_action = "impute")

  expect_true(result$na_impute)
  expect_equal(nrow(result$missing_index), 3)
  expect_false(anyNA(result$x))
})

# ==============================================================================
# 7. GGM + impute guard
# ==============================================================================

test_that("GGM + impute errors", {
  x <- matrix(rnorm(12), nrow = 4, ncol = 3)
  expect_error(
    validate_missing_data(x, na_action = "impute", is_continuous = TRUE),
    "not yet supported for the Gaussian model"
  )
})

test_that("GGM + listwise works", {
  x <- matrix(rnorm(12), nrow = 4, ncol = 3)
  x[1, 2] <- NA
  withr::with_options(list(bgms.verbose = FALSE), {
    result <- validate_missing_data(x, na_action = "listwise", is_continuous = TRUE)
  })
  expect_equal(nrow(result$x), 3)
  expect_false(result$na_impute)
})

# ==============================================================================
# 8. Group vector filtering (bgmCompare path)
# ==============================================================================

test_that("listwise with group filters group vector", {
  x <- matrix(1:12, nrow = 4, ncol = 3)
  x[2, 1] <- NA
  group <- c(1, 1, 2, 2)

  withr::with_options(list(bgms.verbose = FALSE), {
    result <- validate_missing_data(x, na_action = "listwise", group = group)
  })
  expect_equal(length(result$group), nrow(result$x))
  expect_equal(result$group, c(1, 2, 2))
})

test_that("impute with group preserves group vector", {
  x <- matrix(1:12, nrow = 4, ncol = 3)
  group <- c(1, 1, 2, 2)
  result <- validate_missing_data(x, na_action = "impute", group = group)
  expect_equal(result$group, group)
})

test_that("group is absent from result when not provided", {
  x <- matrix(1:12, nrow = 4, ncol = 3)
  result <- validate_missing_data(x, na_action = "listwise")
  expect_null(result$group)
})
