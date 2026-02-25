# ==============================================================================
# Unit tests for validate_baseline_category()
# Phase A.2 of the scaffolding refactor
# ==============================================================================

# Helper: access the internal function
vbc = bgms:::validate_baseline_category

# --- All ordinal (no Blume-Capel) -------------------------------------------

test_that("all-ordinal returns zeros", {
  x = matrix(c(0, 1, 2, 1, 0, 2), nrow = 3, ncol = 2)
  res = vbc(
    baseline_category          = 999,
    baseline_category_provided = FALSE,
    x                          = x,
    variable_bool              = c(TRUE, TRUE)
  )
  expect_equal(res, c(0, 0))
})

test_that("all-ordinal ignores baseline_category value", {
  x = matrix(c(0, 1, 2, 1, 0, 2), nrow = 3, ncol = 2)
  res = vbc(
    baseline_category          = 42,
    baseline_category_provided = TRUE,
    x                          = x,
    variable_bool              = c(TRUE, TRUE)
  )
  expect_equal(res, c(0, 0))
})

# --- Blume-Capel: baseline_category not provided ----------------------------

test_that("BC variables without baseline_category errors", {
  x = matrix(c(0, 1, 2, 1, 0, 2), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = NULL,
      baseline_category_provided = FALSE,
      x                          = x,
      variable_bool              = c(FALSE, TRUE)
    ),
    "baseline_category is required"
  )
})

# --- Scalar baseline_category -----------------------------------------------

test_that("scalar baseline_category replicates to ncol(x)", {
  x = matrix(c(0, 1, 2, 0, 1, 2, 0, 1, 2), nrow = 3, ncol = 3)
  res = vbc(
    baseline_category          = 1,
    baseline_category_provided = TRUE,
    x                          = x,
    variable_bool              = c(FALSE, FALSE, FALSE)
  )
  expect_equal(res, c(1, 1, 1))
})

test_that("scalar non-integer baseline_category errors", {
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  # 1.3 is clearly non-integer
  expect_error(
    vbc(
      baseline_category          = 1.3,
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "integer value"
  )
  # 1.5 was NOT caught by the old code due to banker's rounding bug;
  # the fixed code (using abs()) now catches it correctly.
  expect_error(
    vbc(
      baseline_category          = 1.5,
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "integer value"
  )
})

test_that("scalar NA baseline_category errors", {
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = NA,
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "missing value"
  )
})

# --- Vector baseline_category ------------------------------------------------

test_that("vector baseline_category with correct length works", {
  x = matrix(c(0, 1, 2, 0, 1, 2, 0, 1, 2), nrow = 3, ncol = 3)
  res = vbc(
    baseline_category          = c(0, 1, 2),
    baseline_category_provided = TRUE,
    x                          = x,
    variable_bool              = c(FALSE, FALSE, FALSE)
  )
  expect_equal(res, c(0, 1, 2))
})

test_that("wrong-length vector baseline_category errors", {
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = c(0, 1, 2),
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "single integer or a vector"
  )
})

test_that("non-integer in vector baseline_category errors", {
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = c(1, 1.3),
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "need to be integer|needs to be an integer"
  )
  # 1.5 was missed by the old code; fixed with abs()
  expect_error(
    vbc(
      baseline_category          = c(1, 1.5),
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "need to be integer|needs to be an integer"
  )
})

test_that("NA in vector baseline_category errors", {
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = c(1, NA),
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "missing values"
  )
})

# --- Out of range -----------------------------------------------------------

test_that("baseline_category below minimum errors", {
  x = matrix(c(1, 2, 3, 1, 2, 3), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = 0,
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "within the range"
  )
})

test_that("baseline_category above maximum errors", {
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = 5,
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(FALSE, FALSE)
    ),
    "within the range"
  )
})

# --- Mixed ordinal + Blume-Capel --------------------------------------------

test_that("mixed ordinal+BC validates only BC variables", {
  x = matrix(c(0, 1, 2, 0, 1, 2, 0, 1, 2), nrow = 3, ncol = 3)
  # variable 1 = ordinal, variable 2+3 = BC
  res = vbc(
    baseline_category          = c(0, 1, 1),
    baseline_category_provided = TRUE,
    x                          = x,
    variable_bool              = c(TRUE, FALSE, FALSE)
  )
  expect_equal(res, c(0, 1, 1))
})

test_that("mixed: out-of-range on ordinal column is still caught", {
  # baseline_category for ordinal columns is checked against data range too
  x = matrix(c(0, 1, 2, 0, 1, 2), nrow = 3, ncol = 2)
  expect_error(
    vbc(
      baseline_category          = c(5, 1),
      baseline_category_provided = TRUE,
      x                          = x,
      variable_bool              = c(TRUE, FALSE)
    ),
    "within the range"
  )
})
