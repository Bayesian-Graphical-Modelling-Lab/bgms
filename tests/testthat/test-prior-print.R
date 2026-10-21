# ==============================================================================
# Tests for the prior print methods
# ==============================================================================
#
# Kept apart from test-prior-interface.R, which skips on CRAN because of the
# fits further down it. These tests build a prior object and print it, with no
# sampler involved, so they run on CRAN.
# ==============================================================================


# Each print method must show the hyperparameters, not only the family label.
# One case per branch of each switch(), including the unrecognised-family
# fallbacks.

test_that("print.bgms_parameter_prior shows each family with its parameters", {
  expect_output(print(cauchy_prior(scale = 2.5)),
    "Parameter prior: Cauchy(0, 2.5)",
    fixed = TRUE
  )
  expect_output(print(normal_prior(scale = 0.5)),
    "Parameter prior: Normal(0, 0.5)",
    fixed = TRUE
  )
  expect_output(print(beta_prime_prior(alpha = 1.5, beta = 0.25)),
    "Parameter prior: Beta-prime(alpha = 1.5, beta = 0.25)",
    fixed = TRUE
  )

  # Fallback branch: a family the switch() does not name still prints.
  unknown = structure(
    list(family = "laplace", hyper.parameters = list(scale = 1)),
    class = "bgms_parameter_prior"
  )
  expect_output(print(unknown), "Parameter prior: laplace", fixed = TRUE)

  # print() returns its argument invisibly.
  expect_output(expect_invisible(print(cauchy_prior())))
})

test_that("print.bgms_scale_prior distinguishes the raw and standardized frames", {
  # Raw frame reports the rate it was given.
  expect_output(print(gamma_prior(shape = 2, rate = 0.5)),
    "Scale prior: Gamma(shape = 2, rate = 0.5)",
    fixed = TRUE
  )
  expect_output(print(exponential_prior(rate = 2)),
    "Scale prior: Exponential(rate = 2)",
    fixed = TRUE
  )

  # Standardized frame reports eta and labels it; the two frames give different
  # priors for the same number.
  expect_output(print(gamma_prior(shape = 2, eta = 0.5)),
    "Scale prior: Gamma(shape = 2, eta = 0.5, standardized frame)",
    fixed = TRUE
  )
  expect_output(print(exponential_prior(eta = 3)),
    "Scale prior: Exponential(eta = 3, standardized frame)",
    fixed = TRUE
  )

  unknown = structure(
    list(family = "half-cauchy", hyper.parameters = list(rate = 1)),
    class = "bgms_scale_prior"
  )
  expect_output(print(unknown), "Scale prior: half-cauchy", fixed = TRUE)

  # print() returns its argument invisibly.
  expect_output(expect_invisible(print(gamma_prior())))
})

test_that("print.bgms_indicator_prior shows each family with its parameters", {
  expect_output(print(bernoulli_prior(inclusion_probability = 0.25)),
    "Edge prior: Bernoulli(0.25)",
    fixed = TRUE
  )

  # Matrix branch: a variable-specific matrix has no single number to print,
  # so the method says what it is instead of formatting the first cell.
  expect_output(print(bernoulli_prior(inclusion_probability = matrix(0.3, 4, 4))),
    "Edge prior: Bernoulli (variable-specific inclusion probabilities)",
    fixed = TRUE
  )

  expect_output(print(beta_bernoulli_prior(alpha = 2, beta = 3)),
    "Edge prior: Beta-Bernoulli(alpha = 2, beta = 3)",
    fixed = TRUE
  )

  sbm = sbm_prior(
    alpha = 1, beta = 2, alpha_between = 3, beta_between = 4,
    dirichlet_alpha = 5, lambda = 6
  )
  expect_output(print(sbm), "Edge prior: Stochastic-Block", fixed = TRUE)
  expect_output(print(sbm), "Within:    Beta(1, 2)", fixed = TRUE)
  expect_output(print(sbm), "Between:   Beta(3, 4)", fixed = TRUE)
  expect_output(print(sbm), "Dirichlet: 5, Lambda: 6", fixed = TRUE)

  unknown = structure(
    list(family = "Erdos-Renyi", hyper.parameters = list()),
    class = "bgms_indicator_prior"
  )
  expect_output(print(unknown), "Edge prior: Erdos-Renyi", fixed = TRUE)

  # print() returns its argument invisibly.
  expect_output(expect_invisible(print(bernoulli_prior())))
})
