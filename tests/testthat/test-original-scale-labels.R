# ==============================================================================
# Ordinal category thresholds are reported using the original category values
# (the recode map, category_levels) rather than the internal rescored indices
# 1, 2, ... See ordinal_threshold_labels() and build_output naming.
# ==============================================================================

test_that("ordinal_threshold_labels uses original category values when available", {
  f = bgms:::ordinal_threshold_labels
  # Training {1,3,5} -> recoded {0,1,2}; thresholds for cats 1,2 are original 3,5.
  expect_equal(f(2, c(1, 3, 5)), c(3, 5))
  # 0-based contiguous: original == rescored.
  expect_equal(f(2, c(0, 1, 2)), c(1, 2))
  expect_equal(f(3, c(1, 2, 4, 8)), c(2, 4, 8))
})

test_that("ordinal_threshold_labels falls back to rescored indices without a map", {
  f = bgms:::ordinal_threshold_labels
  expect_equal(f(2, NULL), seq_len(2))
  # Defensive: a map whose length does not match K+1 falls back, too.
  expect_equal(f(2, c(1, 3)), seq_len(2))
})

test_that("summary labels ordinal thresholds in the original category scale", {
  set.seed(1)
  n = 60L
  x = cbind(sample(c(1, 3, 5), n, TRUE), sample(0:2, n, TRUE))
  colnames(x) = c("A", "B")
  fit = bgm(x,
    iter = 80, warmup = 80, chains = 1,
    display_progress = "none", update_method = "adaptive-metropolis"
  )
  rn = rownames(summary(fit)$main)
  # Variable A's thresholds carry the original values 3 and 5.
  expect_true(all(c("A (3)", "A (5)") %in% rn))
  # Variable B is 0-based, so labels are unchanged.
  expect_true(all(c("B (1)", "B (2)") %in% rn))
})

test_that("mixed-MRF discrete thresholds use the original category scale", {
  set.seed(1)
  n = 80L
  x = cbind(sample(c(1, 3, 5), n, TRUE), rnorm(n))
  colnames(x) = c("D", "C")
  fit = bgm(x,
    variable_type = c("ordinal", "continuous"),
    iter = 80, warmup = 80, chains = 1, display_progress = "none"
  )
  rn = rownames(summary(fit)$main)
  # Discrete variable D's thresholds carry the original values 3 and 5;
  # the continuous variable keeps its (mean) label.
  expect_true(all(c("D (3)", "D (5)") %in% rn))
  expect_true("C (mean)" %in% rn)
})
