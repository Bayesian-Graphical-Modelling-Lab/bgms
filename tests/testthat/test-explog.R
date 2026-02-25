# --------------------------------------------------------------------------- #
# Tests for the exp/log implementation used internally by bgms.
#
# The package can be compiled with either the standard C library exp/log or a
# custom IEEE 754 implementation (openlibm-derived). These tests verify that
# whichever is active produces correct results, agrees with R's built-in
# exp()/log(), and round-trips properly.
# --------------------------------------------------------------------------- #

test_that("get_explog_switch returns a valid mode string", {
  mode = get_explog_switch()
  expect_true(mode %in% c("custom", "standard"))
})


# ---- exp() correctness ----------------------------------------------------- #

test_that("internal exp matches R exp for typical values", {
  x = c(-100, -10, -1, -0.5, 0, 0.5, 1, 2, 10, 100, 500)
  # Tolerance allows for ulp-level differences between openlibm and platform exp
  expect_equal(rcpp_ieee754_exp(x), exp(x), tolerance = 1e-14)
})

test_that("internal exp handles special values", {
  expect_equal(rcpp_ieee754_exp(0), 1)
  expect_equal(rcpp_ieee754_exp(-Inf), 0)
  expect_equal(rcpp_ieee754_exp(Inf), Inf)
  expect_true(is.nan(rcpp_ieee754_exp(NaN)))
})

test_that("internal exp handles overflow/underflow gracefully", {
  # Should overflow to Inf
  expect_equal(rcpp_ieee754_exp(710), Inf)
  # Should underflow to 0 (not produce NaN or error)
  expect_equal(rcpp_ieee754_exp(-750), 0)
})


# ---- log() correctness ----------------------------------------------------- #

test_that("internal log matches R log for typical values", {
  x = c(0.001, 0.01, 0.1, 0.5, 1, 2, 10, 100, 1e10, 1e100)
  # Tolerance allows for ulp-level differences between openlibm and platform log
  expect_equal(rcpp_ieee754_log(x), log(x), tolerance = 1e-14)
})

test_that("internal log handles special values", {
  expect_equal(rcpp_ieee754_log(1), 0)
  expect_equal(rcpp_ieee754_log(Inf), Inf)
  expect_equal(rcpp_ieee754_log(0), -Inf)
  expect_true(is.nan(rcpp_ieee754_log(-1)))
  expect_true(is.nan(rcpp_ieee754_log(NaN)))
})


# ---- Round-trip ------------------------------------------------------------- #

test_that("log(exp(x)) round-trips for values in a safe range", {
  x = seq(-500, 500, length.out = 1001)
  expect_equal(rcpp_ieee754_log(rcpp_ieee754_exp(x)), x, tolerance = 1e-10)
})

test_that("exp(log(x)) round-trips for positive values", {
  x = c(1e-300, 1e-100, 0.01, 1, 100, 1e100, 1e300)
  expect_equal(rcpp_ieee754_exp(rcpp_ieee754_log(x)), x, tolerance = 1e-10)
})
