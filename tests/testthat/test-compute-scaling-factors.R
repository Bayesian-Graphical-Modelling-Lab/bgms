# ==============================================================================
# Unit tests for compute_scaling_factors()
# Phase A.8 of the R scaffolding refactor.
# ==============================================================================


# ==============================================================================
# 1. standardize = FALSE  →  ones matrix
# ==============================================================================

test_that("standardize = FALSE returns all-ones matrix with names", {
  res = compute_scaling_factors(
    num_variables     = 3,
    is_ordinal        = c(TRUE, TRUE, TRUE),
    num_categories    = c(3, 4, 2),
    baseline_category = c(0, 0, 0),
    standardize       = FALSE,
    varnames          = c("A", "B", "C")
  )
  expect_equal(dim(res), c(3, 3))
  expect_true(all(res == 1))
  expect_equal(rownames(res), c("A", "B", "C"))
  expect_equal(colnames(res), c("A", "B", "C"))
})

test_that("standardize = FALSE with NULL varnames gives default names", {
  res = compute_scaling_factors(
    num_variables     = 2,
    is_ordinal        = c(TRUE, TRUE),
    num_categories    = c(3, 4),
    baseline_category = c(0, 0),
    standardize       = FALSE,
    varnames          = NULL
  )
  expect_equal(rownames(res), c("Variable 1", "Variable 2"))
  expect_equal(colnames(res), c("Variable 1", "Variable 2"))
})


# ==============================================================================
# 2. Both ordinal  →  M_i * M_j
# ==============================================================================

test_that("two ordinal variables: scaling = M1 * M2", {
  # Variables with 4 and 3 categories (so max scores M=4 and M=3)
  res = compute_scaling_factors(
    num_variables     = 2,
    is_ordinal        = c(TRUE, TRUE),
    num_categories    = c(4, 3),
    baseline_category = c(0, 0),
    standardize       = TRUE,
    varnames          = c("V1", "V2")
  )
  expect_equal(res[1, 2], 4 * 3)
  expect_equal(res[2, 1], 4 * 3) # symmetric
  expect_equal(res[1, 1], 1) # diagonal stays 1
  expect_equal(res[2, 2], 1)
})

test_that("three ordinal variables: all pairs", {
  res = compute_scaling_factors(
    num_variables     = 3,
    is_ordinal        = c(TRUE, TRUE, TRUE),
    num_categories    = c(2, 5, 3),
    baseline_category = c(0, 0, 0),
    standardize       = TRUE,
    varnames          = c("A", "B", "C")
  )
  expect_equal(res[1, 2], 2 * 5)
  expect_equal(res[1, 3], 2 * 3)
  expect_equal(res[2, 3], 5 * 3)
  # Symmetry
  expect_equal(res[2, 1], res[1, 2])
  expect_equal(res[3, 1], res[1, 3])
  expect_equal(res[3, 2], res[2, 3])
})


# ==============================================================================
# 3. Both Blume-Capel  →  max |endpoints|
# ==============================================================================

test_that("two Blume-Capel: baseline at 0, range is (-0, M-0) = (0, M)", {
  # M1=4, b1=0 → endpoints (-0, 4); M2=3, b2=0 → endpoints (0, 3)
  # max(|outer|) = max(0*0, 0*3, 4*0, 4*3) = 12
  res = compute_scaling_factors(
    num_variables     = 2,
    is_ordinal        = c(FALSE, FALSE),
    num_categories    = c(4, 3),
    baseline_category = c(0, 0),
    standardize       = TRUE,
    varnames          = c("V1", "V2")
  )
  expect_equal(res[1, 2], 12)
  expect_equal(res[2, 1], 12)
})

test_that("two Blume-Capel: non-zero baselines", {
  # M1=4, b1=2 → endpoints (-2, 2); M2=5, b2=1 → endpoints (-1, 4)
  # outer: (-2)*(-1)=2, (-2)*4=-8, 2*(-1)=-2, 2*4=8
  # max(abs) = 8
  res = compute_scaling_factors(
    num_variables     = 2,
    is_ordinal        = c(FALSE, FALSE),
    num_categories    = c(4, 5),
    baseline_category = c(2, 1),
    standardize       = TRUE,
    varnames          = c("V1", "V2")
  )
  expect_equal(res[1, 2], 8)
})


# ==============================================================================
# 4. Mixed: one ordinal, one Blume-Capel
# ==============================================================================

test_that("ordinal + Blume-Capel: v1 ordinal, v2 BC", {
  # V1: ordinal, M=3, range (0, 3)
  # V2: BC, M=4, b=1, range (-1, 3)
  # outer: 0*(-1)=0, 0*3=0, 3*(-1)=-3, 3*3=9
  # max(abs) = 9
  res = compute_scaling_factors(
    num_variables     = 2,
    is_ordinal        = c(TRUE, FALSE),
    num_categories    = c(3, 4),
    baseline_category = c(0, 1),
    standardize       = TRUE,
    varnames          = c("Ord", "BC")
  )
  expect_equal(res[1, 2], 9)
  expect_equal(res[2, 1], 9)
})

test_that("Blume-Capel + ordinal: v1 BC, v2 ordinal", {
  # V1: BC, M=4, b=1, range (-1, 3)
  # V2: ordinal, M=3, range (0, 3)
  # Same result as above by symmetry of the formula
  res = compute_scaling_factors(
    num_variables     = 2,
    is_ordinal        = c(FALSE, TRUE),
    num_categories    = c(4, 3),
    baseline_category = c(1, 0),
    standardize       = TRUE,
    varnames          = c("BC", "Ord")
  )
  expect_equal(res[1, 2], 9)
  expect_equal(res[2, 1], 9)
})


# ==============================================================================
# 5. Mixed 3-variable scenario (all branch combos)
# ==============================================================================

test_that("3 variables: ordinal, BC, ordinal — covers all mixed pairs", {
  # V1: ordinal, M=3
  # V2: BC, M=5, b=2 → range (-2, 3)
  # V3: ordinal, M=2
  #
  # (1,2) ordinal+BC: endpoints1=(0,3), endpoints2=(-2,3)
  #   max(|0*(-2), 0*3, 3*(-2), 3*3|) = max(0, 0, 6, 9) = 9
  # (1,3) ordinal+ordinal: 3 * 2 = 6
  # (2,3) BC+ordinal: endpoints1=(-2,3), endpoints2=(0,2)
  #   max(|(-2)*0, (-2)*2, 3*0, 3*2|) = max(0, 4, 0, 6) = 6
  res = compute_scaling_factors(
    num_variables     = 3,
    is_ordinal        = c(TRUE, FALSE, TRUE),
    num_categories    = c(3, 5, 2),
    baseline_category = c(0, 2, 0),
    standardize       = TRUE,
    varnames          = c("Ord1", "BC", "Ord2")
  )
  expect_equal(res[1, 2], 9)
  expect_equal(res[1, 3], 6)
  expect_equal(res[2, 3], 6)
  # Symmetry
  expect_equal(res[2, 1], res[1, 2])
  expect_equal(res[3, 1], res[1, 3])
  expect_equal(res[3, 2], res[2, 3])
  # Diagonal
  expect_equal(unname(diag(res)), c(1, 1, 1))
})


# ==============================================================================
# 6. Matrix properties
# ==============================================================================

test_that("result is always symmetric", {
  res = compute_scaling_factors(
    num_variables     = 4,
    is_ordinal        = c(TRUE, FALSE, TRUE, FALSE),
    num_categories    = c(3, 4, 2, 5),
    baseline_category = c(0, 1, 0, 2),
    standardize       = TRUE,
    varnames          = paste0("V", 1:4)
  )
  expect_equal(res, t(res))
})

test_that("diagonal is always 1", {
  res = compute_scaling_factors(
    num_variables     = 4,
    is_ordinal        = c(TRUE, FALSE, TRUE, FALSE),
    num_categories    = c(3, 4, 2, 5),
    baseline_category = c(0, 1, 0, 2),
    standardize       = TRUE,
    varnames          = paste0("V", 1:4)
  )
  expect_equal(unname(diag(res)), rep(1, 4))
})

test_that("all scaling factors are positive", {
  res = compute_scaling_factors(
    num_variables     = 3,
    is_ordinal        = c(FALSE, FALSE, TRUE),
    num_categories    = c(4, 5, 3),
    baseline_category = c(2, 1, 0),
    standardize       = TRUE,
    varnames          = c("A", "B", "C")
  )
  expect_true(all(res > 0))
})


# ==============================================================================
# 7. Single variable (edge case)
# ==============================================================================

test_that("single variable returns 1x1 matrix", {
  res = compute_scaling_factors(
    num_variables     = 1,
    is_ordinal        = TRUE,
    num_categories    = 3,
    baseline_category = 0,
    standardize       = TRUE,
    varnames          = "Only"
  )
  expect_equal(dim(res), c(1, 1))
  expect_equal(res[1, 1], 1)
  expect_equal(rownames(res), "Only")
})
