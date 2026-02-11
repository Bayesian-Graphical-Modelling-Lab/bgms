test_that("GGM runs and returns bgms object", {
  testthat::skip_on_cran()

  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) <- paste0("V", 1:4)

  fit <- bgm(
    x = x,
    variable_type = "continuous",
    iter = 200,
    warmup = 1000,
    chains = 1,
    display_progress = "none",
    seed = 123
  )

  expect_s3_class(fit, "bgms")
  expect_true(fit$arguments$is_continuous)
  expect_equal(fit$arguments$num_variables, 4)
  expect_equal(fit$arguments$iter, 200)
  expect_equal(fit$raw_samples$niter, 200)
  expect_equal(fit$raw_samples$nchains, 1)
})

test_that("GGM output has correct dimensions", {
  testthat::skip_on_cran()

  p <- 5
  set.seed(42)
  x <- matrix(rnorm(100 * p), nrow = 100, ncol = p)
  colnames(x) <- paste0("V", 1:p)

  fit <- bgm(
    x = x,
    variable_type = "continuous",
    edge_selection = TRUE,
    iter = 100,
    warmup = 1000,
    chains = 1,
    display_progress = "none",
    seed = 42
  )

  # main: p diagonal precision elements
  expect_equal(nrow(fit$posterior_summary_main), p)
  expect_equal(nrow(fit$posterior_mean_main), p)

  # pairwise: p*(p-1)/2 off-diagonal elements
  n_edges <- p * (p - 1) / 2
  expect_equal(nrow(fit$posterior_summary_pairwise), n_edges)
  expect_equal(nrow(fit$posterior_mean_pairwise), p)
  expect_equal(ncol(fit$posterior_mean_pairwise), p)

  # indicators
  expect_equal(nrow(fit$posterior_summary_indicator), n_edges)
  expect_equal(nrow(fit$posterior_mean_indicator), p)
  expect_equal(ncol(fit$posterior_mean_indicator), p)

  # raw samples
  expect_equal(ncol(fit$raw_samples$main[[1]]), p)
  expect_equal(ncol(fit$raw_samples$pairwise[[1]]), n_edges)
  expect_equal(nrow(fit$raw_samples$main[[1]]), 100)
})

test_that("GGM without edge selection works", {
  testthat::skip_on_cran()

  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) <- paste0("V", 1:4)

  fit <- bgm(
    x = x,
    variable_type = "continuous",
    edge_selection = FALSE,
    iter = 100,
    warmup = 1000,
    chains = 1,
    display_progress = "none",
    seed = 99
  )

  expect_s3_class(fit, "bgms")
  expect_null(fit$posterior_summary_indicator)
  expect_null(fit$posterior_mean_indicator)
})

test_that("GGM rejects NUTS and HMC", {
  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)

  expect_error(
    bgm(x = x, variable_type = "continuous", update_method = "nuts"),
    "only supports.*adaptive-metropolis"
  )
  expect_error(
    bgm(x = x, variable_type = "continuous", update_method = "hamiltonian-mc"),
    "only supports.*adaptive-metropolis"
  )
})

test_that("Mixed continuous and ordinal is rejected", {
  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)

  expect_error(
    bgm(x = x, variable_type = c("continuous", "ordinal", "ordinal", "ordinal")),
    "all variables must be of type"
  )
})

test_that("GGM is reproducible", {
  testthat::skip_on_cran()

  set.seed(42)
  x <- matrix(rnorm(200), nrow = 50, ncol = 4)
  colnames(x) <- paste0("V", 1:4)

  fit1 <- bgm(
    x = x, variable_type = "continuous",
    iter = 100, warmup = 1000, chains = 1,
    display_progress = "none", seed = 321
  )
  fit2 <- bgm(
    x = x, variable_type = "continuous",
    iter = 100, warmup = 1000, chains = 1,
    display_progress = "none", seed = 321
  )

  expect_equal(fit1$raw_samples$main, fit2$raw_samples$main)
  expect_equal(fit1$raw_samples$pairwise, fit2$raw_samples$pairwise)
})

test_that("GGM posterior precision diagonal means are positive", {
  testthat::skip_on_cran()

  set.seed(42)
  x <- matrix(rnorm(500), nrow = 100, ncol = 5)
  colnames(x) <- paste0("V", 1:5)

  fit <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 500, warmup = 1000, chains = 1,
    display_progress = "none", seed = 42
  )

  # Diagonal precision elements should be positive
  expect_true(all(fit$posterior_summary_main$mean > 0))
})


test_that("GGM posterior compare against simulated data", {
  testthat::skip_on_cran()

  n <- 1000
  p <- 10
  ne <- p * (p - 1) / 2
  # avoid a test dependency on BDgraph and a random graph structure by using a fixed precision matrix
  # set.seed(42)
  # adj <- matrix(0, p, p)
  # adj[lower.tri(adj)] <- runif(ne) < .5
  # adj <- adj + t(adj)
  # omega <- zapsmall(BDgraph::rgwish(1, adj))
  omega <- structure(c(6.240119, 0, 0, -0.370239, 0, 0, 0, 0, -1.622902,
              0, 0, 1.905013, 0, -0.194995, 0, 0, -2.468628, -0.557277, 0,
              0, 0, 0, 5.509142, -7.942389, 1.40081, 0, 0, -0.76775, 0, 0,
              -0.370239, -0.194995, -7.942389, 15.521405, -3.537489, 0, 4.60785,
              0, 3.278511, 0, 0, 0, 1.40081, -3.537489, 2.78257, 0, 0, 1.374641,
              0, -1.198092, 0, 0, 0, 0, 0, 1.350879, 0, 0.230677, -1.357952,
              0, 0, -2.468628, 0, 4.60785, 0, 0, 15.88698, 0, 1.20017, -1.973919,
              0, -0.557277, -0.76775, 0, 1.374641, 0.230677, 0, 7.007312, 1.597035,
              0, -1.622902, 0, 0, 3.278511, 0, -1.357952, 1.20017, 1.597035,
              13.378039, -4.769958, 0, 0, 0, 0, -1.198092, 0, -1.973919, 0,
              -4.769958, 5.536877), dim = c(10L, 10L))
  adj <- omega != 0
  diag(adj) <- 0
  covmat <- solve(omega)
  chol <- chol(covmat)

  set.seed(43)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol
  # cov(x) - covmat

  fit_no_vs <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 5000, warmup = 1000, chains = 2,
    display_progress = "none", seed = 42
  )

  expect_true(cor(fit_no_vs$posterior_summary_main$mean,     diag(omega)) > 0.9)
  expect_true(cor(fit_no_vs$posterior_summary_pairwise$mean, omega[upper.tri(omega)]) > 0.9)

  fit_vs <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = TRUE,
    iter = 10000, warmup = 2000, chains = 2,
    display_progress = "none", seed = 42
  )

  expect_true(cor(fit_vs$posterior_summary_main$mean,     diag(omega)) > 0.9)
  expect_true(cor(fit_vs$posterior_summary_pairwise$mean, omega[upper.tri(omega)]) > 0.9)
  expect_true(cor(fit_vs$posterior_summary_indicator$mean, adj[upper.tri(adj)]) > 0.85)

  # This test somehow fails? Can you fix this?
  fit_vs_sbm <- bgm(
    x = x, variable_type = "continuous",
    edge_selection = TRUE,
    edge_prior = "Stochastic-Block",
    iter = 5000, warmup = 1000, chains = 2,
    display_progress = "none", seed = 42
  )

  expect_true(cor(fit_vs_sbm$posterior_summary_main$mean,     diag(omega)) > 0.9)
  expect_true(cor(fit_vs_sbm$posterior_summary_pairwise$mean, omega[upper.tri(omega)]) > 0.9)
  expect_true(cor(fit_vs_sbm$posterior_summary_indicator$mean, adj[upper.tri(adj)]) > 0.85)

  # SBM-specific output
  expect_false(is.null(fit_vs_sbm$posterior_mean_coclustering_matrix))
  expect_equal(nrow(fit_vs_sbm$posterior_mean_coclustering_matrix), p)
  expect_equal(ncol(fit_vs_sbm$posterior_mean_coclustering_matrix), p)
  expect_false(is.null(fit_vs_sbm$posterior_num_blocks))
  expect_false(is.null(fit_vs_sbm$posterior_mode_allocations))
  expect_false(is.null(fit_vs_sbm$raw_samples$allocations))


})
