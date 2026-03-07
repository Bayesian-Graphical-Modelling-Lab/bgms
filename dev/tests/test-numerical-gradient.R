# dev/tests/test-numerical-gradient.R
#
# Numerical gradient check for MixedMRFModel::logp_and_gradient().
# Compares the analytical gradient against central finite differences.
# Tests both conditional and marginal pseudo-likelihood modes.

library(bgms)
library(testthat)

# --- Helper: central finite differences --------------------------------------
numerical_gradient = function(f, theta, eps = 1e-5) {
  grad = numeric(length(theta))
  for(i in seq_along(theta)) {
    theta_plus = theta
    theta_minus = theta
    theta_plus[i] = theta[i] + eps
    theta_minus[i] = theta[i] - eps
    grad[i] = (f(theta_plus) - f(theta_minus)) / (2 * eps)
  }
  grad
}

# --- Small test problem -------------------------------------------------------
set.seed(42)
n = 100
p = 2  # discrete variables, ordinal with 3 categories each
q = 2  # continuous variables

# Simulate discrete data (values 0, 1, 2)
x_discrete = matrix(sample(0:2, n * p, replace = TRUE), nrow = n, ncol = p)
storage.mode(x_discrete) = "integer"

# Simulate continuous data
x_continuous = matrix(rnorm(n * q), nrow = n, ncol = q)

num_categories = as.integer(rep(3L, p))
is_ordinal = as.integer(rep(1L, p))
baseline_category = as.integer(rep(0L, p))

# NUTS vector dimension: mux (2*3=6) + Kxx_ut (1) + muy (2) + Kxy (2*2=4) = 13
nuts_dim = sum(num_categories) + p * (p - 1) / 2 + q + p * q

make_input = function(pl) {
  list(
    discrete_observations   = x_discrete,
    continuous_observations = x_continuous,
    num_categories          = num_categories,
    is_ordinal_variable     = is_ordinal,
    baseline_category       = baseline_category,
    main_alpha              = 1.0,
    main_beta               = 1.0,
    pairwise_scale          = 2.5,
    pseudolikelihood        = pl
  )
}

# Start near zero with small random perturbations
theta0 = rnorm(nuts_dim, sd = 0.1)

# --- Test conditional PL ------------------------------------------------------
cat("Testing conditional pseudo-likelihood gradient...\n")
inp_cond = make_input("conditional")

logp_fn_cond = function(th) {
  bgms:::test_mixed_gradient(inp_cond, th)$logp
}

result_cond = bgms:::test_mixed_gradient(inp_cond, theta0)
grad_analytical = result_cond$gradient
grad_numerical = numerical_gradient(logp_fn_cond, theta0)

cat("  Max absolute difference:", max(abs(grad_analytical - grad_numerical)), "\n")
cat("  Max relative difference:",
    max(abs(grad_analytical - grad_numerical) /
      pmax(abs(grad_analytical), abs(grad_numerical), 1e-8)), "\n")

test_that("conditional PL gradient matches finite differences", {
  expect_equal(grad_analytical, grad_numerical, tolerance = 1e-4)
})

# --- Test marginal PL --------------------------------------------------------
cat("Testing marginal pseudo-likelihood gradient...\n")
inp_marg = make_input("marginal")

logp_fn_marg = function(th) {
  bgms:::test_mixed_gradient(inp_marg, th)$logp
}

result_marg = bgms:::test_mixed_gradient(inp_marg, theta0)
grad_analytical_m = result_marg$gradient
grad_numerical_m = numerical_gradient(logp_fn_marg, theta0)

cat("  Max absolute difference:", max(abs(grad_analytical_m - grad_numerical_m)), "\n")
cat("  Max relative difference:",
    max(abs(grad_analytical_m - grad_numerical_m) /
      pmax(abs(grad_analytical_m), abs(grad_numerical_m), 1e-8)), "\n")

test_that("marginal PL gradient matches finite differences", {
  expect_equal(grad_analytical_m, grad_numerical_m, tolerance = 1e-4)
})

# --- Test with Blume-Capel variable -------------------------------------------
cat("Testing with Blume-Capel variable...\n")
is_ordinal_bc = as.integer(c(1L, 0L))  # first ordinal, second Blume-Capel
# Blume-Capel has 2 main-effects (location + dispersion)
# NUTS dim: mux (3 + 2) + Kxx_ut (1) + muy (2) + Kxy (2*2=4) = 12
nuts_dim_bc = 3 + 2 + 1 + 2 + 4

inp_bc = make_input("conditional")
inp_bc$is_ordinal_variable = is_ordinal_bc

theta0_bc = rnorm(nuts_dim_bc, sd = 0.1)

logp_fn_bc = function(th) {
  bgms:::test_mixed_gradient(inp_bc, th)$logp
}

result_bc = bgms:::test_mixed_gradient(inp_bc, theta0_bc)
grad_analytical_bc = result_bc$gradient
grad_numerical_bc = numerical_gradient(logp_fn_bc, theta0_bc)

cat("  Max absolute difference:", max(abs(grad_analytical_bc - grad_numerical_bc)), "\n")
cat("  Max relative difference:",
    max(abs(grad_analytical_bc - grad_numerical_bc) /
      pmax(abs(grad_analytical_bc), abs(grad_numerical_bc), 1e-8)), "\n")

test_that("conditional PL gradient with Blume-Capel matches finite differences", {
  expect_equal(grad_analytical_bc, grad_numerical_bc, tolerance = 1e-4)
})

# --- Test marginal PL with Blume-Capel ----------------------------------------
cat("Testing marginal PL with Blume-Capel variable...\n")
inp_bc_marg = make_input("marginal")
inp_bc_marg$is_ordinal_variable = is_ordinal_bc

logp_fn_bc_marg = function(th) {
  bgms:::test_mixed_gradient(inp_bc_marg, th)$logp
}

result_bc_marg = bgms:::test_mixed_gradient(inp_bc_marg, theta0_bc)
grad_analytical_bc_m = result_bc_marg$gradient
grad_numerical_bc_m = numerical_gradient(logp_fn_bc_marg, theta0_bc)

cat("  Max absolute difference:", max(abs(grad_analytical_bc_m - grad_numerical_bc_m)), "\n")
cat("  Max relative difference:",
    max(abs(grad_analytical_bc_m - grad_numerical_bc_m) /
      pmax(abs(grad_analytical_bc_m), abs(grad_numerical_bc_m), 1e-8)), "\n")

test_that("marginal PL gradient with Blume-Capel matches finite differences", {
  expect_equal(grad_analytical_bc_m, grad_numerical_bc_m, tolerance = 1e-4)
})

cat("\nAll gradient checks passed!\n")
