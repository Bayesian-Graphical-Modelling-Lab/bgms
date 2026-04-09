# ==============================================================================
# Tests for C++ polymorphic parameter prior classes
#
# 1. Correctness: logp and grad match R reference implementations
# 2. Numerical gradient: analytic grad matches finite differences
# 3. Scaled variants: logp(x, sf) matches logp at scaled parameters
# ==============================================================================


# ==============================================================================
# 1. Prior logp/grad correctness against R reference
# ==============================================================================

test_that("CauchyPrior logp matches R::dcauchy", {
  for(scale in c(0.5, 1.0, 2.5)) {
    for(x in c(-2.0, 0.0, 0.3, 1.5)) {
      res = test_parameter_prior("cauchy", x, scale = scale)
      expect_equal(res$logp, dcauchy(x, 0, scale, log = TRUE),
        tolerance = 1e-12,
        info = sprintf("cauchy logp at x=%.1f, scale=%.1f", x, scale)
      )
    }
  }
})

test_that("CauchyPrior grad matches analytic formula", {
  for(scale in c(0.5, 1.0, 2.5)) {
    for(x in c(-2.0, 0.0, 0.3, 1.5)) {
      res = test_parameter_prior("cauchy", x, scale = scale)
      expected = -2 * x / (scale^2 + x^2)
      expect_equal(res$grad, expected,
        tolerance = 1e-12,
        info = sprintf("cauchy grad at x=%.1f, scale=%.1f", x, scale)
      )
    }
  }
})

test_that("NormalPrior logp matches R::dnorm", {
  for(scale in c(0.5, 1.0, 2.5)) {
    for(x in c(-2.0, 0.0, 0.3, 1.5)) {
      res = test_parameter_prior("normal", x, scale = scale)
      expect_equal(res$logp, dnorm(x, 0, scale, log = TRUE),
        tolerance = 1e-12,
        info = sprintf("normal logp at x=%.1f, scale=%.1f", x, scale)
      )
    }
  }
})

test_that("NormalPrior grad matches analytic formula", {
  for(scale in c(0.5, 1.0, 2.5)) {
    for(x in c(-2.0, 0.1, 0.3, 1.5)) {
      res = test_parameter_prior("normal", x, scale = scale)
      expected = -x / scale^2
      expect_equal(res$grad, expected,
        tolerance = 1e-12,
        info = sprintf("normal grad at x=%.1f, scale=%.1f", x, scale)
      )
    }
  }
})

test_that("BetaPrimePrior logp matches logit-Beta formula", {
  for(alpha in c(0.5, 1.0, 2.0)) {
    for(beta in c(0.5, 1.0, 3.0)) {
      for(x in c(-2.0, 0.0, 0.5, 2.0)) {
        res = test_parameter_prior("beta-prime", x, alpha = alpha, beta = beta)
        expected = x * alpha - log1p(exp(x)) * (alpha + beta)
        expect_equal(res$logp, expected,
          tolerance = 1e-12,
          info = sprintf("beta-prime logp at x=%.1f, a=%.1f, b=%.1f", x, alpha, beta)
        )
      }
    }
  }
})

test_that("BetaPrimePrior grad matches analytic formula", {
  for(alpha in c(0.5, 1.0, 2.0)) {
    for(beta in c(0.5, 1.0, 3.0)) {
      for(x in c(-2.0, 0.0, 0.5, 2.0)) {
        res = test_parameter_prior("beta-prime", x, alpha = alpha, beta = beta)
        sigmoid = 1 / (1 + exp(-x))
        expected = alpha - (alpha + beta) * sigmoid
        expect_equal(res$grad, expected,
          tolerance = 1e-12,
          info = sprintf("beta-prime grad at x=%.1f, a=%.1f, b=%.1f", x, alpha, beta)
        )
      }
    }
  }
})

test_that("GammaScalePrior logp matches R::dgamma", {
  for(shape in c(1.0, 2.0, 0.5)) {
    for(rate in c(0.5, 1.0, 2.0)) {
      for(x in c(0.1, 0.5, 1.0, 3.0)) {
        res = test_scale_prior("gamma", x, shape = shape, rate = rate)
        expected = dgamma(x, shape = shape, rate = rate, log = TRUE)
        expect_equal(res$logp, expected,
          tolerance = 1e-10,
          info = sprintf("gamma logp at x=%.1f, shape=%.1f, rate=%.1f", x, shape, rate)
        )
      }
    }
  }
})

test_that("GammaScalePrior grad matches analytic formula", {
  for(shape in c(1.0, 2.0, 0.5)) {
    for(rate in c(0.5, 1.0, 2.0)) {
      for(x in c(0.1, 0.5, 1.0, 3.0)) {
        res = test_scale_prior("gamma", x, shape = shape, rate = rate)
        expected = (shape - 1) / x - rate
        expect_equal(res$grad, expected,
          tolerance = 1e-10,
          info = sprintf("gamma grad at x=%.1f, shape=%.1f, rate=%.1f", x, shape, rate)
        )
      }
    }
  }
})

test_that("ExponentialPrior (via create_scale_prior) matches Gamma(1, rate)", {
  for(rate in c(0.5, 1.0, 2.0)) {
    for(x in c(0.1, 0.5, 1.0, 3.0)) {
      res_exp = test_scale_prior("exponential", x, rate = rate)
      res_gamma = test_scale_prior("gamma", x, shape = 1.0, rate = rate)
      expect_equal(res_exp$logp, res_gamma$logp, tolerance = 1e-12)
      expect_equal(res_exp$grad, res_gamma$grad, tolerance = 1e-12)
    }
  }
})


# ==============================================================================
# 2. Scaled prior variants
# ==============================================================================

test_that("CauchyPrior logp_scaled matches Cauchy with scaled width", {
  scale = 1.0
  sf = 2.0
  for(x in c(-1.0, 0.0, 0.5, 2.0)) {
    res = test_parameter_prior("cauchy", x, scale = scale, scale_factor = sf)
    expected = dcauchy(x, 0, scale * sf, log = TRUE)
    expect_equal(res$logp_scaled, expected, tolerance = 1e-12)
  }
})

test_that("NormalPrior logp_scaled matches Normal with scaled sd", {
  scale = 1.0
  sf = 3.0
  for(x in c(-1.0, 0.0, 0.5, 2.0)) {
    res = test_parameter_prior("normal", x, scale = scale, scale_factor = sf)
    expected = dnorm(x, 0, scale * sf, log = TRUE)
    expect_equal(res$logp_scaled, expected, tolerance = 1e-12)
  }
})

test_that("CauchyPrior grad_scaled matches analytic at scaled width", {
  scale = 1.5
  sf = 2.0
  for(x in c(-1.0, 0.1, 0.5, 2.0)) {
    res = test_parameter_prior("cauchy", x, scale = scale, scale_factor = sf)
    s = scale * sf
    expected = -2 * x / (s^2 + x^2)
    expect_equal(res$grad_scaled, expected, tolerance = 1e-12)
  }
})

test_that("NormalPrior grad_scaled matches analytic at scaled sd", {
  scale = 0.5
  sf = 4.0
  for(x in c(-1.0, 0.1, 0.5, 2.0)) {
    res = test_parameter_prior("normal", x, scale = scale, scale_factor = sf)
    s = scale * sf
    expected = -x / s^2
    expect_equal(res$grad_scaled, expected, tolerance = 1e-12)
  }
})


# ==============================================================================
# 3. Numerical gradient verification for GGM with non-default priors
# ==============================================================================

make_edge_matrix = function(p, included_edges) {
  E = matrix(0L, nrow = p, ncol = p)
  if(length(included_edges) > 0) {
    for(k in seq_along(included_edges)) {
      ij = included_edges[[k]]
      E[ij[1], ij[2]] = 1L
      E[ij[2], ij[1]] = 1L
    }
  }
  E
}

theta_dim = function(edge_mat) {
  p = nrow(edge_mat)
  p + sum(edge_mat[upper.tri(edge_mat)] == 1L)
}

fd_gradient_prior = function(theta, suf_stat, n, edge_mat,
                             ipt, is, ia, ib, dpt, ds, dr,
                             eps = 1e-6) {
  g = numeric(length(theta))
  for(k in seq_along(theta)) {
    t_plus = theta
    t_minus = theta
    t_plus[k] = t_plus[k] + eps
    t_minus[k] = t_minus[k] - eps
    f_plus = ggm_test_logp_and_gradient_prior(
      t_plus, suf_stat, n, edge_mat,
      ipt, is, ia, ib, dpt, ds, dr
    )$value
    f_minus = ggm_test_logp_and_gradient_prior(
      t_minus, suf_stat, n, edge_mat,
      ipt, is, ia, ib, dpt, ds, dr
    )$value
    g[k] = (f_plus - f_minus) / (2 * eps)
  }
  g
}

check_gradient_prior = function(p, edge_mat, n = 200, seed = 1,
                                interaction_prior_type = "cauchy",
                                interaction_scale = 1.0,
                                interaction_alpha = 0.5,
                                interaction_beta = 0.5,
                                diagonal_prior_type = "gamma",
                                diagonal_shape = 1.0,
                                diagonal_rate = 1.0,
                                eps = 1e-6, tol = 1e-4) {
  set.seed(seed)
  X = matrix(rnorm(n * p), nrow = n, ncol = p)
  S = t(X) %*% X

  d = theta_dim(edge_mat)
  theta = rnorm(d, sd = 0.2)

  ag = ggm_test_logp_and_gradient_prior(
    theta, S, n, edge_mat,
    interaction_prior_type, interaction_scale,
    interaction_alpha, interaction_beta,
    diagonal_prior_type, diagonal_shape, diagonal_rate
  )
  fd = fd_gradient_prior(
    theta, S, n, edge_mat,
    interaction_prior_type, interaction_scale,
    interaction_alpha, interaction_beta,
    diagonal_prior_type, diagonal_shape, diagonal_rate,
    eps = eps
  )

  denom = pmax(abs(ag$gradient), abs(fd), 1)
  rel_err = abs(ag$gradient - fd) / denom
  max(rel_err)
}


test_that("GGM gradient correct with NormalPrior (full graph)", {
  p = 4
  E = matrix(1L, p, p)
  diag(E) = 0L
  err = check_gradient_prior(p, E,
    interaction_prior_type = "normal",
    interaction_scale = 2.0
  )
  expect_true(err < 1e-4, info = sprintf("max relative error: %g", err))
})

test_that("GGM gradient correct with BetaPrimePrior (full graph)", {
  p = 4
  E = matrix(1L, p, p)
  diag(E) = 0L
  err = check_gradient_prior(p, E,
    interaction_prior_type = "beta-prime",
    interaction_alpha = 1.0,
    interaction_beta = 1.0
  )
  expect_true(err < 1e-4, info = sprintf("max relative error: %g", err))
})

test_that("GGM gradient correct with Gamma(2, 0.5) diagonal prior", {
  p = 4
  E = matrix(1L, p, p)
  diag(E) = 0L
  err = check_gradient_prior(p, E,
    diagonal_prior_type = "gamma",
    diagonal_shape = 2.0,
    diagonal_rate = 0.5
  )
  expect_true(err < 1e-4, info = sprintf("max relative error: %g", err))
})

test_that("GGM gradient correct with NormalPrior (sparse graph)", {
  p = 4
  edges = list(c(1, 2), c(2, 3), c(3, 4))
  E = make_edge_matrix(p, edges)
  err = check_gradient_prior(p, E,
    interaction_prior_type = "normal",
    interaction_scale = 1.5
  )
  expect_true(err < 1e-4, info = sprintf("max relative error: %g", err))
})

test_that("GGM gradient correct with BetaPrimePrior + Gamma(2,1) diagonal (sparse)", {
  p = 4
  edges = list(c(1, 2), c(2, 3), c(3, 4))
  E = make_edge_matrix(p, edges)
  err = check_gradient_prior(p, E,
    interaction_prior_type = "beta-prime",
    interaction_alpha = 0.5,
    interaction_beta = 0.5,
    diagonal_prior_type = "gamma",
    diagonal_shape = 2.0,
    diagonal_rate = 1.0
  )
  expect_true(err < 1e-4, info = sprintf("max relative error: %g", err))
})


# ==============================================================================
# 4. Full-space gradient (RATTLE) with non-default priors
# ==============================================================================

fd_gradient_full_prior = function(x, suf_stat, n, edge_mat,
                                  ipt, is, ia, ib, dpt, ds, dr,
                                  eps = 1e-6) {
  g = numeric(length(x))
  for(k in seq_along(x)) {
    x_plus = x
    x_minus = x
    x_plus[k] = x_plus[k] + eps
    x_minus[k] = x_minus[k] - eps
    f_plus = ggm_test_logp_and_gradient_full_prior(
      x_plus, suf_stat, n, edge_mat,
      ipt, is, ia, ib, dpt, ds, dr
    )$value
    f_minus = ggm_test_logp_and_gradient_full_prior(
      x_minus, suf_stat, n, edge_mat,
      ipt, is, ia, ib, dpt, ds, dr
    )$value
    g[k] = (f_plus - f_minus) / (2 * eps)
  }
  g
}

test_that("GGM full gradient correct with NormalPrior + Gamma(2,1)", {
  p = 4
  edges = list(c(1, 2), c(2, 3), c(3, 4))
  E = make_edge_matrix(p, edges)
  full_dim = p * (p + 1) / 2

  set.seed(42)
  n = 200
  X = matrix(rnorm(n * p), nrow = n, ncol = p)
  S = t(X) %*% X
  x = rnorm(full_dim, sd = 0.2)

  ag = ggm_test_logp_and_gradient_full_prior(
    x, S, n, E,
    "normal", 1.5, 0.5, 0.5,
    "gamma", 2.0, 1.0
  )
  fd = fd_gradient_full_prior(
    x, S, n, E,
    "normal", 1.5, 0.5, 0.5,
    "gamma", 2.0, 1.0
  )

  denom = pmax(abs(ag$gradient), abs(fd), 1)
  rel_err = abs(ag$gradient - fd) / denom
  max_err = max(rel_err)
  expect_true(max_err < 1e-4, info = sprintf("max relative error: %g", max_err))
})


# ==============================================================================
# 5. Means prior and diagonal prior consumption tests
# ==============================================================================

test_that("precision_scale_prior affects GGM posterior", {
  set.seed(42)
  Y = as.data.frame(matrix(rnorm(200), nrow = 50, ncol = 4))

  fit_default = bgm(Y,
    variable_type = "continuous",
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    edge_selection = FALSE,
    iter = 100, warmup = 100, chains = 1,
    display_progress = "none"
  )

  fit_strong = bgm(Y,
    variable_type = "continuous",
    precision_scale_prior = gamma_prior(shape = 100, rate = 100),
    edge_selection = FALSE,
    iter = 100, warmup = 100, chains = 1,
    display_progress = "none"
  )

  # A Gamma(100, 100) prior strongly concentrates diagonal around 1,
  # while Gamma(1, 1) = Exp(1) is much more diffuse.
  # The posterior mean diagonals should differ.
  diag_default = fit_default$posterior_mean_residual_variance
  diag_strong = fit_strong$posterior_mean_residual_variance

  expect_false(
    isTRUE(all.equal(diag_default, diag_strong, tolerance = 0.01)),
    info = "Different diagonal priors should produce different posterior residual variances"
  )
})

test_that("means_prior affects mixed MRF posterior", {
  set.seed(42)
  data("Wenchuan", package = "bgms")
  dat = data.frame(Wenchuan[1:100, 1:3], V4 = rnorm(100) + 5)

  fit_default = bgm(dat,
    variable_type = c(rep("ordinal", 3), "continuous"),
    means_prior = normal_prior(scale = 1),
    edge_selection = FALSE,
    iter = 100, warmup = 100, chains = 1,
    display_progress = "none"
  )

  # Very tight prior centered at 0: should shrink the mean toward 0,
  # even though the data mean is ~5
  fit_tight = bgm(dat,
    variable_type = c(rep("ordinal", 3), "continuous"),
    means_prior = normal_prior(scale = 0.01),
    edge_selection = FALSE,
    iter = 100, warmup = 100, chains = 1,
    display_progress = "none"
  )

  # Extract continuous means from posterior
  means_default = fit_default$posterior_mean_main$continuous
  means_tight = fit_tight$posterior_mean_main$continuous

  # The tight prior should shrink the mean closer to 0
  expect_true(
    abs(means_tight[1]) < abs(means_default[1]),
    info = sprintf(
      "Tight means_prior should shrink mean toward 0: tight=%.3f, default=%.3f",
      means_tight[1], means_default[1]
    )
  )
})
