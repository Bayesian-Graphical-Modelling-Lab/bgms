# Regression test for the OMRF residual-matrix invariant during stage-3b
# proposal-SD tuning (NUTS + edge-selection path).
#
# The model maintains residual_matrix_ == 2 * observations_double_ *
# pairwise_effects_ incrementally. tune_proposal_sd() moves the pairwise
# effects and must update the residual with the same factor of 2 that every
# other residual update uses. A missing factor (the historical bug) left the
# residual inconsistent after any accepted stage-3b move, biasing the
# between-model edge-proposal SD tuning. See test_omrf_residual_invariant().

test_that("tune_proposal_sd preserves the residual = 2*X*pairwise invariant", {
  set.seed(1)
  n = 80L
  p = 4L
  x = matrix(rbinom(n * p, 1, 0.5), n, p)
  storage.mode(x) = "integer"
  num_categories = rep(1L, p)
  pairwise = matrix(0.3, p, p)
  diag(pairwise) = 0

  res = test_omrf_residual_invariant(
    x, num_categories, pairwise,
    warmup = 400L, seed = 42L,
    enable_selection = TRUE, learn_sd = TRUE
  )

  # The harness must actually exercise stage-3b (otherwise the test is vacuous).
  expect_true(res$any_move_accepted)
  expect_gt(max(abs(res$pairwise_after - pairwise)), 0)

  # The invariant must hold to machine precision after tuning.
  expect_lt(res$max_abs_discrepancy, 1e-10)
})
