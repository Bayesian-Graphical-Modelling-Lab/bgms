# ==============================================================================
# Phase A.3: MixedMRFModel Skeleton Tests
# ==============================================================================
#
# Verifies that the MixedMRFModel C++ class compiles, links, and that
# parameter_dimension() and get/set_vectorized_parameters() round-trip
# correctly. Tests cover:
#
#   1. All-ordinal variables (p ordinal discrete + q continuous)
#   2. All-Blume-Capel variables
#   3. Mixed ordinal + Blume-Capel variables
#   4. Both conditional and marginal pseudolikelihood modes
#   5. Clone round-trip
#   6. Dimension arithmetic (num_variables, num_pairwise, indicator length)
#
# ==============================================================================

# The helper-fixtures.R file loads bgms, but test_file() may skip it.
# Ensure the package is available:
if(!requireNamespace("bgms", quietly = TRUE)) {
  skip("bgms package not installed")
}
test_mixed_mrf_skeleton = bgms:::test_mixed_mrf_skeleton

# ---- Helper to build test inputs for a given scenario ----
make_mixed_mrf_inputs = function(n, p, q, num_categories, is_ordinal,
                                  baseline_category, edge_selection = FALSE,
                                  pseudolikelihood = "conditional",
                                  seed = 42L) {
  set.seed(seed)
  # Discrete observations: each column in 0..(num_categories[s])
  discrete_obs = matrix(0L, nrow = n, ncol = p)
  for(s in seq_len(p)) {
    discrete_obs[, s] = sample(0L:num_categories[s], n, replace = TRUE)
  }
  # Continuous observations
  continuous_obs = matrix(rnorm(n * q), nrow = n, ncol = q)

  total = p + q
  inclusion_prob = matrix(0.5, nrow = total, ncol = total)
  diag(inclusion_prob) = 0
  edge_ind = matrix(1L, nrow = total, ncol = total)
  diag(edge_ind) = 0L

  list(
    discrete_observations = discrete_obs,
    continuous_observations = continuous_obs,
    num_categories = as.integer(num_categories),
    is_ordinal_variable = as.integer(is_ordinal),
    baseline_category = as.integer(baseline_category),
    inclusion_probability = inclusion_prob,
    initial_edge_indicators = edge_ind,
    edge_selection = edge_selection,
    pseudolikelihood = pseudolikelihood,
    seed = seed
  )
}

run_skeleton_test = function(inputs) {
  do.call(test_mixed_mrf_skeleton, inputs)
}

# ---- Expected dimension calculator (mirrors C++ logic) ----
expected_full_dim = function(p, q, num_categories, is_ordinal) {
  num_main = sum(ifelse(is_ordinal == 1L, num_categories, 2L))
  num_pairwise_xx = p * (p - 1) / 2
  num_pairwise_yy_with_diag = q * (q + 1) / 2
  num_cross = p * q
  num_main + num_pairwise_xx + q + num_pairwise_yy_with_diag + num_cross
}

expected_num_pairwise = function(p, q) {
  p * (p - 1) / 2 + q * (q - 1) / 2 + p * q
}

expected_indicator_length = function(p, q) {
  expected_num_pairwise(p, q)
}

# ==============================================================================
# Test scenarios
# ==============================================================================

test_that("all-ordinal skeleton has correct dimensions", {
  p = 3L; q = 2L; n = 20L
  num_cats = c(2L, 3L, 2L)
  is_ord = c(1L, 1L, 1L)
  baseline = rep(0L, p)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  expect_equal(res$num_variables, p + q)
  expect_equal(res$parameter_dimension,
               expected_full_dim(p, q, num_cats, is_ord))
  expect_equal(res$full_parameter_dimension,
               expected_full_dim(p, q, num_cats, is_ord))
  expect_equal(res$params_length, res$full_parameter_dimension)
  expect_equal(res$num_pairwise, expected_num_pairwise(p, q))
  expect_equal(res$indicator_length, expected_indicator_length(p, q))
  expect_equal(res$edge_indicators_rows, p + q)
  expect_equal(res$edge_indicators_cols, p + q)
})

test_that("all-Blume-Capel skeleton has correct dimensions", {
  p = 3L; q = 2L; n = 20L
  num_cats = c(4L, 3L, 4L)
  is_ord = c(0L, 0L, 0L)
  baseline = c(2L, 1L, 2L)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  expect_equal(res$num_variables, p + q)
  # BC always contributes 2 per variable
  expect_equal(res$parameter_dimension,
               expected_full_dim(p, q, num_cats, is_ord))
  expect_equal(res$full_parameter_dimension,
               expected_full_dim(p, q, num_cats, is_ord))
})

test_that("mixed ordinal + BC skeleton has correct dimensions", {
  p = 4L; q = 3L; n = 25L
  num_cats = c(2L, 4L, 3L, 3L)
  is_ord = c(1L, 0L, 1L, 0L)
  baseline = c(0L, 2L, 0L, 1L)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  expect_equal(res$num_variables, p + q)
  expect_equal(res$parameter_dimension,
               expected_full_dim(p, q, num_cats, is_ord))
})

test_that("zero-parameter round-trip is exact", {
  p = 3L; q = 2L; n = 20L
  num_cats = c(2L, 3L, 2L)
  is_ord = c(1L, 1L, 1L)
  baseline = rep(0L, p)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  expect_equal(res$roundtrip_max_diff, 0.0)
})

test_that("non-trivial parameter round-trip is exact", {
  p = 3L; q = 2L; n = 20L
  num_cats = c(2L, 3L, 2L)
  is_ord = c(1L, 1L, 1L)
  baseline = rep(0L, p)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  expect_equal(res$nontrivial_roundtrip_max_diff, 0.0)
})

test_that("non-trivial round-trip with Blume-Capel variables", {
  p = 4L; q = 3L; n = 25L
  num_cats = c(2L, 4L, 3L, 3L)
  is_ord = c(1L, 0L, 1L, 0L)
  baseline = c(0L, 2L, 0L, 1L)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  expect_equal(res$nontrivial_roundtrip_max_diff, 0.0)
})

test_that("clone produces identical parameters", {
  p = 3L; q = 2L; n = 20L
  num_cats = c(2L, 3L, 2L)
  is_ord = c(1L, 1L, 1L)
  baseline = rep(0L, p)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  expect_equal(res$clone_max_diff, 0.0)
})

test_that("marginal pseudolikelihood mode works", {
  p = 3L; q = 2L; n = 20L
  num_cats = c(2L, 3L, 2L)
  is_ord = c(1L, 1L, 1L)
  baseline = rep(0L, p)

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline,
                               pseudolikelihood = "marginal")
  res = run_skeleton_test(inp)

  expect_equal(res$parameter_dimension,
               expected_full_dim(p, q, num_cats, is_ord))
  expect_equal(res$nontrivial_roundtrip_max_diff, 0.0)
})

test_that("capability queries return correct values", {
  p = 3L; q = 2L; n = 20L
  num_cats = c(2L, 3L, 2L)
  is_ord = c(1L, 1L, 1L)
  baseline = rep(0L, p)

  inp_no_es = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline,
                                     edge_selection = FALSE)
  res_no_es = run_skeleton_test(inp_no_es)
  expect_false(res_no_es$has_edge_selection)
  expect_true(res_no_es$has_adaptive_metropolis)
  expect_false(res_no_es$has_gradient)
  expect_false(res_no_es$has_missing_data)

  inp_es = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline,
                                  edge_selection = TRUE)
  res_es = run_skeleton_test(inp_es)
  expect_true(res_es$has_edge_selection)
})

test_that("single discrete + single continuous works", {
  p = 1L; q = 1L; n = 15L
  num_cats = 2L
  is_ord = 1L
  baseline = 0L

  inp = make_mixed_mrf_inputs(n, p, q, num_cats, is_ord, baseline)
  res = run_skeleton_test(inp)

  # p=1: 0 xx edges, q=1: 0 yy edges, 1 cross edge
  expect_equal(res$num_variables, 2L)
  expect_equal(res$num_pairwise, 1L)
  expect_equal(res$indicator_length, 1L)
  # mux: 2 thresholds, Kxx: 0, muy: 1, Kyy: 1 diag, Kxy: 1
  expect_equal(res$full_parameter_dimension, 2 + 0 + 1 + 1 + 1)
  expect_equal(res$nontrivial_roundtrip_max_diff, 0.0)
})
