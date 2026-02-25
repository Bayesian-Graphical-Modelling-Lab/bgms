# ==============================================================================
# Tests for collapse_categories_across_groups()
# Phase A.6b of the R scaffolding refactor.
# ==============================================================================


# ==============================================================================
# 1. No collapsing needed — all categories in all groups
# ==============================================================================

test_that("no collapsing when all categories appear in all groups", {
  x = matrix(c(
    0, 1, 2, 0, 1, 2,
    0, 1, 2, 0, 1, 2
  ), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 2)
  bc = c(0, 0)

  result = collapse_categories_across_groups(
    x, group, is_ordinal,
    num_categories, bc
  )
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(2, 2))
  expect_equal(result$baseline_category, bc)
})


# ==============================================================================
# 2. Basic collapsing — category missing from one group
# ==============================================================================

test_that("category missing from one group is collapsed", {
  # Variable 1: group 1 has {0,1,2}, group 2 has {0,2} → category 1 collapsed
  # After collapsing: 0→0, 1→0, 2→1
  x = matrix(c(
    0, 1, 2, 0, 2, 0,
    0, 1, 2, 0, 1, 2
  ), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 2)
  bc = c(0, 0)

  result = collapse_categories_across_groups(
    x, group, is_ordinal,
    num_categories, bc
  )
  expect_equal(result$x[, 1], c(0, 0, 1, 0, 1, 0))
  expect_equal(result$num_categories[1], 1)
  # Variable 2 unchanged (all categories in both groups)
  expect_equal(result$x[, 2], c(0, 1, 2, 0, 1, 2))
  expect_equal(result$num_categories[2], 2)
})

test_that("multiple categories collapsed from one variable", {
  # Variable: group 1 has {0,1,2,3}, group 2 has {0,3}
  # Categories 1 and 2 missing from group 2 → collapse
  # 0→0, 1→0, 2→0, 3→1
  x = matrix(c(0, 1, 2, 3, 0, 3, 0, 1, 0, 1, 0, 1), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 1, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(3, 1)
  bc = c(0, 0)

  result = collapse_categories_across_groups(
    x, group, is_ordinal,
    num_categories, bc
  )
  expect_equal(result$x[, 1], c(0, 0, 0, 1, 0, 1))
  expect_equal(result$num_categories[1], 1)
})


# ==============================================================================
# 3. BC variables pass through unchanged
# ==============================================================================

test_that("BC variables are not affected by group collapsing", {
  # Variable 1: BC, values 0,1,2,3
  # Variable 2: ordinal, values 0,1,2
  x = matrix(c(
    0, 1, 2, 3, 0, 1,
    0, 1, 2, 0, 2, 0
  ), nrow = 6, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2)
  is_ordinal = c(FALSE, TRUE)
  num_categories = c(3, 2)
  bc = c(1, 0)

  result = collapse_categories_across_groups(
    x, group, is_ordinal,
    num_categories, bc
  )
  # BC variable unchanged
  expect_equal(result$x[, 1], c(0, 1, 2, 3, 0, 1))
  expect_equal(result$num_categories[1], 3)
  expect_equal(result$baseline_category[1], 1)
  # Ordinal variable: group 1 has {0,1,2}, group 2 has {0,2} → 1 collapsed
  expect_equal(result$x[, 2], c(0, 0, 1, 0, 1, 0))
})


# ==============================================================================
# 4. More than two groups
# ==============================================================================

test_that("works with three groups", {
  # Variable: group 1 has {0,1,2}, group 2 has {0,1}, group 3 has {0,2}
  # Category 1 missing from group 3, category 2 missing from group 2
  # Only category 0 is in all groups → everything collapses to 0
  x = matrix(c(
    0, 1, 2, 0, 1, 0, 2, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 1, 0
  ), nrow = 9, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 1)
  bc = c(0, 0)

  # Variable 1: only 0 is in all groups → num_categories = 0 → error
  expect_error(
    collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    ),
    "Only one value was observed for variable 1"
  )
})

test_that("three groups with all categories present", {
  x = matrix(c(
    0, 1, 2, 0, 1, 2, 0, 1, 2,
    0, 1, 0, 0, 1, 0, 0, 1, 0
  ), nrow = 9, ncol = 2)
  group = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(2, 1)
  bc = c(0, 0)

  result = collapse_categories_across_groups(
    x, group, is_ordinal,
    num_categories, bc
  )
  expect_equal(result$x, x)
  expect_equal(result$num_categories, c(2, 1))
})


# ==============================================================================
# 5. Edge case: all-same-value variable after collapsing → error
# ==============================================================================

test_that("error when all categories collapse to one value", {
  # Group 1 has {0, 1}, group 2 has {2, 3} → no overlap → everything → 0
  x = matrix(c(
    0, 1, 2, 3,
    0, 1, 0, 1
  ), nrow = 4, ncol = 2)
  group = c(1, 1, 2, 2)
  is_ordinal = c(TRUE, TRUE)
  num_categories = c(3, 1)
  bc = c(0, 0)

  expect_error(
    collapse_categories_across_groups(
      x, group, is_ordinal,
      num_categories, bc
    ),
    "Only one value was observed for variable 1"
  )
})


# ==============================================================================
# 6. Integration: compare_reformat_data delegates correctly
# ==============================================================================

test_that("full pipeline works end-to-end with group collapsing", {
  # Variable 1: ordinal, values 1,2,3 across groups
  # Group 1 has {1,2,3}, group 2 has {1,3}
  # After ordinal recode: 0,1,2. Group 2 has {0,2} → category 1 missing
  # After collapse: 0→0, 1→0, 2→1
  x = matrix(c(
    1, 2, 3, 1, 2, 3, 1, 3, 1, 3,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 10, ncol = 2)
  group = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2)
  variable_bool = c(TRUE, TRUE)
  bc = c(0, 0)

  md = validate_missing_data(x,
    na_action = "listwise", is_continuous = FALSE,
    group = group
  )
  ord = reformat_ordinal_data(md$x,
    is_ordinal = variable_bool,
    baseline_category = bc
  )
  result = collapse_categories_across_groups(
    x = ord$x, group = md$group, is_ordinal = variable_bool,
    num_categories = ord$num_categories,
    baseline_category = ord$baseline_category
  )

  # Variable 1: 1→0, 2→0(collapsed), 3→1
  expect_equal(result$x[1:6, 1], c(0, 0, 1, 0, 0, 1))
  expect_equal(result$x[7:10, 1], c(0, 1, 0, 1))
  expect_equal(result$num_categories[1], 1)

  # Variable 2: 0,1 both in all groups → unchanged
  expect_equal(result$x[, 2], c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1))
  expect_equal(result$num_categories[2], 1)
})

test_that("full pipeline preserves BC variables through collapsing", {
  x = matrix(c(
    0, 1, 2, 3, 0, 1, 2, 3,
    0, 1, 2, 0, 0, 1, 2, 0
  ), nrow = 8, ncol = 2)
  group = c(1, 1, 1, 1, 2, 2, 2, 2)
  variable_bool = c(FALSE, TRUE)
  bc = c(1, 0)

  md = validate_missing_data(x,
    na_action = "listwise", is_continuous = FALSE,
    group = group
  )
  ord = reformat_ordinal_data(md$x,
    is_ordinal = variable_bool,
    baseline_category = bc
  )
  result = collapse_categories_across_groups(
    x = ord$x, group = md$group, is_ordinal = variable_bool,
    num_categories = ord$num_categories,
    baseline_category = ord$baseline_category
  )

  # BC variable: no collapsing, just standard ordinal recode
  expect_equal(result$x[, 1], c(0, 1, 2, 3, 0, 1, 2, 3))
  expect_equal(result$num_categories[1], 3)
  expect_equal(result$baseline_category[1], 1)

  # Ordinal variable: 0,1,2 all in both groups → unchanged
  expect_equal(result$x[, 2], c(0, 1, 2, 0, 0, 1, 2, 0))
  expect_equal(result$num_categories[2], 2)
})
