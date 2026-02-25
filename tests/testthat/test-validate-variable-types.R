# ==============================================================================
# Unit tests for validate_variable_types()
# Phase A.1 of the scaffolding refactor
# ==============================================================================

# Helper: access the internal function
vvt = bgms:::validate_variable_types

# --- Single-string inputs -----------------------------------------------------

test_that("single 'ordinal' replicates to all variables", {
  res = vvt("ordinal", num_variables = 5)
  expect_equal(res$variable_type, rep("ordinal", 5))
  expect_equal(res$variable_bool, rep(TRUE, 5))
  expect_false(res$is_continuous)
})

test_that("single 'blume-capel' replicates to all variables", {
  res = vvt("blume-capel", num_variables = 3)
  expect_equal(res$variable_type, rep("blume-capel", 3))
  expect_equal(res$variable_bool, rep(FALSE, 3))
  expect_false(res$is_continuous)
})

test_that("single 'continuous' works when allowed", {
  res = vvt("continuous", num_variables = 4, allow_continuous = TRUE)
  expect_equal(res$variable_type, rep("continuous", 4))
  expect_equal(res$variable_bool, rep(TRUE, 4))
  expect_true(res$is_continuous)
})

test_that("single 'continuous' is rejected when not allowed", {
  expect_error(
    vvt("continuous", num_variables = 4, allow_continuous = FALSE),
    "not of type continuous"
  )
})

test_that("single invalid type raises error with caller name", {
  expect_error(
    vvt("gaussian", num_variables = 3, caller = "bgm"),
    "bgm function supports"
  )
  expect_error(
    vvt("gaussian", num_variables = 3, caller = "bgmCompare"),
    "bgmCompare function supports"
  )
})

# --- Vector inputs ------------------------------------------------------------

test_that("vector of ordinal types", {
  res = vvt(c("ordinal", "ordinal", "ordinal"), num_variables = 3)
  expect_equal(res$variable_bool, c(TRUE, TRUE, TRUE))
  expect_false(res$is_continuous)
})

test_that("mixed ordinal and blume-capel vector", {
  res = vvt(c("ordinal", "blume-capel", "ordinal"), num_variables = 3)
  expect_equal(res$variable_bool, c(TRUE, FALSE, TRUE))
  expect_false(res$is_continuous)
})

test_that("vector of continuous types works when allowed", {
  res = vvt(c("continuous", "continuous"), num_variables = 2, allow_continuous = TRUE)
  expect_true(res$is_continuous)
  expect_equal(res$variable_bool, c(TRUE, TRUE))
})

test_that("mixed continuous + ordinal is rejected", {
  expect_error(
    vvt(c("continuous", "ordinal"), num_variables = 2, allow_continuous = TRUE),
    "all variables must be of type"
  )
})

test_that("vector of continuous rejected when not allowed", {
  expect_error(
    vvt(c("continuous", "continuous"), num_variables = 2, allow_continuous = FALSE),
    "not of type continuous"
  )
})

test_that("wrong-length vector raises error", {
  expect_error(
    vvt(c("ordinal", "ordinal"), num_variables = 5),
    "vector of character strings of length p"
  )
})

test_that("vector with invalid type raises error", {
  expect_error(
    vvt(c("ordinal", "gaussian"), num_variables = 2, caller = "bgm"),
    "not of type"
  )
})

# --- Default arguments --------------------------------------------------------

test_that("defaults: allow_continuous = TRUE, caller = 'bgm'", {
  # Continuous should work by default

  res = vvt("continuous", num_variables = 2)
  expect_true(res$is_continuous)

  # Error message should mention 'bgm' by default
  expect_error(vvt("invalid", num_variables = 2), "bgm function supports")
})
