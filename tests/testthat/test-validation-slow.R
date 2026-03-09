# ==============================================================================
# Gated slow validation tests
# ==============================================================================
#
# Environment-gated wrappers around the key assertions from the validation
# scripts in dev/tests/validation/. These run in nightly CI (where
# BGMS_RUN_SLOW_TESTS=true) and are skipped during local devtools::test()
# and CRAN checks.
#
# Each test fits a small mixed MRF network with enough iterations to
# produce meaningful posterior estimates, then checks a quantitative
# recovery or agreement criterion.
#
# Source scripts: group1_parameter_recovery.R, group2_mh_vs_nuts.R,
#   group3_cond_vs_marg.R. Helpers in dev/tests/validation/helpers.R.
# ==============================================================================

# Load shared helpers once per file
helpers_path = file.path("dev", "tests", "validation", "helpers.R")
if(!file.exists(helpers_path)) {
  # When running from testthat, the working directory is tests/testthat/
  helpers_path = file.path("..", "..", "dev", "tests", "validation", "helpers.R")
}
helpers_available = file.exists(helpers_path)


# ==============================================================================
# Gate: skip unless BGMS_RUN_SLOW_TESTS=true
# ==============================================================================

skip_slow = function() {
  skip_if_not(
    isTRUE(as.logical(Sys.getenv("BGMS_RUN_SLOW_TESTS", "false"))),
    "Set BGMS_RUN_SLOW_TESTS=true to run slow validation tests"
  )
  skip_if(!helpers_available, "Validation helpers not found")
}


# ==============================================================================
# 1. Parameter recovery (adapted from group1)
# ==============================================================================

test_that("mixed MRF parameter recovery: cor > 0.8 (small network)", {
  skip_slow()
  source(helpers_path, local = TRUE)

  net = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
  dat = generate_data(net, n = 2000, source = "bgms", seed = 201)
  vt = c(rep("ordinal", 2), rep("continuous", 2))

  fit = bgm(dat,
    variable_type = vt,
    pseudolikelihood = "conditional", edge_selection = FALSE,
    iter = 10000, warmup = 5000, chains = 2, seed = 301
  )

  true_blocks = list(
    mux = net$mux, muy = net$muy,
    Kxx = net$Kxx, Kxy = net$Kxy, Kyy = net$Kyy
  )
  est_blocks = extract_bgms_blocks(fit, net)

  true_flat = flatten_params(true_blocks)
  est_flat = flatten_params(est_blocks)

  r = cor(true_flat, est_flat)
  expect_gt(r, 0.8,
    label = sprintf("recovery correlation (%.4f)", r)
  )
})


# ==============================================================================
# 2. Metropolis vs NUTS agreement (adapted from group2)
# ==============================================================================

test_that("MH vs NUTS posterior agreement: cor > 0.95", {
  skip_slow()
  source(helpers_path, local = TRUE)

  net = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
  dat = generate_data(net, n = 2000, source = "bgms", seed = 201)
  vt = c(rep("ordinal", 2), rep("continuous", 2))

  fit_mh = bgm(dat,
    variable_type = vt,
    pseudolikelihood = "conditional",
    update_method = "adaptive-metropolis",
    edge_selection = FALSE,
    iter = 15000, warmup = 10000, chains = 2, seed = 401
  )

  fit_nuts = bgm(dat,
    variable_type = vt,
    pseudolikelihood = "conditional",
    update_method = "nuts",
    edge_selection = FALSE,
    iter = 5000, warmup = 3000, chains = 2, seed = 402
  )

  est_mh = flatten_params(extract_bgms_blocks(fit_mh, net))
  est_nuts = flatten_params(extract_bgms_blocks(fit_nuts, net))

  r = cor(est_mh, est_nuts)
  expect_gt(r, 0.95,
    label = sprintf("MH vs NUTS agreement (%.4f)", r)
  )
})


# ==============================================================================
# 3. Conditional vs marginal PL agreement (adapted from group3)
# ==============================================================================

test_that("conditional vs marginal PL agreement: cor > 0.90", {
  skip_slow()
  source(helpers_path, local = TRUE)

  net = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
  dat = generate_data(net, n = 2000, source = "bgms", seed = 201)
  vt = c(rep("ordinal", 2), rep("continuous", 2))

  fit_cond = bgm(dat,
    variable_type = vt,
    pseudolikelihood = "conditional", edge_selection = FALSE,
    iter = 10000, warmup = 5000, chains = 2, seed = 501
  )

  fit_marg = bgm(dat,
    variable_type = vt,
    pseudolikelihood = "marginal", edge_selection = FALSE,
    iter = 10000, warmup = 5000, chains = 2, seed = 502
  )

  est_cond = flatten_params(extract_bgms_blocks(fit_cond, net))
  est_marg = flatten_params(extract_bgms_blocks(fit_marg, net))

  r = cor(est_cond, est_marg)
  expect_gt(r, 0.90,
    label = sprintf("cond vs marg PL agreement (%.4f)", r)
  )
})


# ==============================================================================
# 4. Estimate-simulate-re-estimate cycle (mixed MRF)
# ==============================================================================

test_that("estimate-simulate-re-estimate cycle: cor > 0.7 (mixed MRF)", {
  skip_slow()
  source(helpers_path, local = TRUE)

  net = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
  dat = generate_data(net, n = 2000, source = "bgms", seed = 201)
  vt = c(rep("ordinal", 2), rep("continuous", 2))

  fit1 = bgm(dat,
    variable_type = vt,
    edge_selection = FALSE, iter = 5000, warmup = 2000,
    chains = 1, seed = 601
  )

  simulated = simulate(fit1, nsim = 2000, seed = 701)

  fit2 = bgm(simulated,
    variable_type = vt,
    edge_selection = FALSE, iter = 5000, warmup = 2000,
    chains = 1, seed = 801
  )

  pw1 = as.vector(fit1$posterior_mean_pairwise)
  pw2 = as.vector(fit2$posterior_mean_pairwise)

  r = cor(pw1, pw2)
  expect_gt(r, 0.7,
    label = sprintf("cycle pairwise correlation (%.4f)", r)
  )
})
