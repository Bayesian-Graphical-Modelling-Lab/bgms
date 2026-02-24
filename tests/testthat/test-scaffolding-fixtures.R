# ==============================================================================
# Golden-Snapshot Fixture Verification Tests
# ==============================================================================
#
# Phase A-0 of the scaffolding refactor (dev/scaffolding/plan.md).
#
# These tests verify that the current validation and preprocessing functions
# produce outputs identical to the golden-snapshot fixtures captured by
# dev/generate_scaffolding_fixtures.R.
#
# During the refactor, these tests ensure that every extracted validator
# reproduces exactly the same intermediate outputs as the original monolithic
# functions. If any test fails, the refactored code has changed behavior.
#
# NOTE: The impute path (na_action = "impute") uses random sampling for
# starting values, so we set the seed before each reformat_data() call.
# The fixtures were generated with set.seed(42) at the top of the script,
# but individual impute calls depend on the accumulated RNG state. To handle
# this, impute fixtures store the expected OUTPUT (not re-derive it), so we
# compare check_model output (deterministic) field-by-field and for
# reformat_data we compare the deterministic fields (num_categories,
# baseline_category, na_impute) while skipping x and missing_index for
# impute cases (since those depend on RNG state).
# ==============================================================================

fixture_dir <- file.path("dev", "fixtures", "scaffolding")

# Skip all tests if fixtures haven't been generated yet
skip_if_no_fixtures <- function() {
  # When running via devtools::test(), the working directory is tests/testthat/
  # The fixtures live at the package root under dev/fixtures/scaffolding/
  pkg_root <- testthat::test_path("..", "..")
  fixture_dir <<- file.path(pkg_root, "dev", "fixtures", "scaffolding")
  manifest_path <- file.path(fixture_dir, "manifest.rds")
  if (!file.exists(manifest_path)) {
    skip("Scaffolding fixtures not found. Run: Rscript dev/generate_scaffolding_fixtures.R")
  }
}

# Helper: load a fixture by id
load_fixture <- function(id) {
  path <- file.path(fixture_dir, paste0(id, ".rds"))
  if (!file.exists(path)) {
    skip(paste("Fixture not found:", id))
  }
  readRDS(path)
}

# ==============================================================================
# check_model / check_compare_model verification
# ==============================================================================

test_that("check_model() output matches fixtures for bgm cases", {
  skip_if_no_fixtures()
  manifest <- readRDS(file.path(fixture_dir, "manifest.rds"))
  bgm_cases <- manifest$id[manifest$type == "bgm"]

  for (id in bgm_cases) {
    fixture <- load_fixture(id)
    input <- fixture$input

    result <- check_model(
      x = input$x,
      variable_type     = input$variable_type,
      baseline_category = input$baseline_category,
      edge_selection    = input$edge_selection %||% TRUE,
      edge_prior        = input$edge_prior %||% "Bernoulli",
      inclusion_probability = input$inclusion_probability %||% 0.5,
      beta_bernoulli_alpha  = input$beta_bernoulli_alpha %||% 1,
      beta_bernoulli_beta   = input$beta_bernoulli_beta %||% 1,
      beta_bernoulli_alpha_between = input$beta_bernoulli_alpha_between %||% 1,
      beta_bernoulli_beta_between  = input$beta_bernoulli_beta_between %||% 1,
      dirichlet_alpha = input$dirichlet_alpha %||% 1,
      lambda          = input$lambda %||% 1
    )

    expected <- fixture$check_model

    expect_identical(result$variable_bool, expected$variable_bool,
      label = paste(id, "- variable_bool"))
    expect_identical(result$baseline_category, expected$baseline_category,
      label = paste(id, "- baseline_category"))
    expect_identical(result$edge_selection, expected$edge_selection,
      label = paste(id, "- edge_selection"))
    expect_identical(result$edge_prior, expected$edge_prior,
      label = paste(id, "- edge_prior"))
    expect_equal(result$inclusion_probability, expected$inclusion_probability,
      label = paste(id, "- inclusion_probability"))
    expect_identical(result$is_continuous, expected$is_continuous,
      label = paste(id, "- is_continuous"))
  }
})

test_that("check_compare_model() output matches fixtures for compare cases", {
  skip_if_no_fixtures()
  manifest <- readRDS(file.path(fixture_dir, "manifest.rds"))
  compare_cases <- manifest$id[manifest$type == "compare"]

  for (id in compare_cases) {
    fixture <- load_fixture(id)
    input <- fixture$input

    result <- check_compare_model(
      x = input$x,
      y = NULL,
      group_indicator      = input$group_indicator,
      difference_selection = input$difference_selection %||% TRUE,
      variable_type        = input$variable_type,
      baseline_category    = input$baseline_category,
      difference_scale     = input$difference_scale %||% 1,
      difference_prior     = input$difference_prior %||% "Bernoulli",
      difference_probability = input$difference_probability %||% 0.5,
      beta_bernoulli_alpha = input$beta_bernoulli_alpha %||% 1,
      beta_bernoulli_beta  = input$beta_bernoulli_beta %||% 1
    )

    expected <- fixture$check_model

    expect_identical(result$variable_bool, expected$variable_bool,
      label = paste(id, "- variable_bool"))
    expect_identical(result$baseline_category, expected$baseline_category,
      label = paste(id, "- baseline_category"))
    expect_equal(result$x, expected$x,
      label = paste(id, "- x (combined data)"))
    expect_equal(result$group_indicator, expected$group_indicator,
      label = paste(id, "- group_indicator"))
    expect_identical(result$difference_prior, expected$difference_prior,
      label = paste(id, "- difference_prior"))
    expect_equal(result$inclusion_probability_difference,
      expected$inclusion_probability_difference,
      label = paste(id, "- inclusion_probability_difference"))
  }
})

# ==============================================================================
# reformat_data / compare_reformat_data verification
# ==============================================================================

test_that("reformat_data() output matches fixtures for bgm OMRF cases", {
  skip_if_no_fixtures()
  manifest <- readRDS(file.path(fixture_dir, "manifest.rds"))
  bgm_cases <- manifest$id[manifest$type == "bgm"]

  for (id in bgm_cases) {
    fixture <- load_fixture(id)

    # GGM cases don't call reformat_data
    if (is.null(fixture$reformat_data)) next

    input <- fixture$input
    cm <- fixture$check_model
    is_impute <- (input$na_action == "impute")

    result <- reformat_data(
      x = input$x,
      na_action     = input$na_action,
      variable_bool = cm$variable_bool,
      baseline_category = cm$baseline_category
    )

    expected <- fixture$reformat_data

    # These fields are always deterministic
    expect_equal(result$num_categories, expected$num_categories,
      label = paste(id, "- num_categories"))
    expect_equal(result$baseline_category, expected$baseline_category,
      label = paste(id, "- baseline_category"))
    expect_identical(result$na_impute, expected$na_impute,
      label = paste(id, "- na_impute"))

    # For listwise: data and missing_index are deterministic
    if (!is_impute) {
      expect_equal(result$x, expected$x,
        label = paste(id, "- x (recoded data)"))
      expect_equal(result$missing_index, expected$missing_index,
        label = paste(id, "- missing_index"))
    } else {
      # For impute: x and missing_index depend on RNG state.
      # Check structural properties instead.
      expect_equal(nrow(result$x), nrow(expected$x),
        label = paste(id, "- x nrow (impute)"))
      expect_equal(ncol(result$x), ncol(expected$x),
        label = paste(id, "- x ncol (impute)"))
      expect_equal(nrow(result$missing_index), nrow(expected$missing_index),
        label = paste(id, "- missing_index nrow (impute)"))
      expect_false(anyNA(result$x),
        label = paste(id, "- no NAs after impute"))
    }
  }
})

test_that("compare_reformat_data() output matches fixtures for compare cases", {
  skip_if_no_fixtures()
  manifest <- readRDS(file.path(fixture_dir, "manifest.rds"))
  compare_cases <- manifest$id[manifest$type == "compare"]

  for (id in compare_cases) {
    fixture <- load_fixture(id)
    input <- fixture$input
    cm <- fixture$check_model
    is_impute <- (input$na_action %||% "listwise") == "impute"

    result <- compare_reformat_data(
      x = cm$x,
      group = cm$group_indicator,
      na_action     = input$na_action %||% "listwise",
      variable_bool = cm$variable_bool,
      baseline_category = cm$baseline_category
    )

    expected <- fixture$reformat_data

    # These fields are deterministic
    expect_equal(result$num_categories, expected$num_categories,
      label = paste(id, "- num_categories"))
    expect_equal(result$baseline_category, expected$baseline_category,
      label = paste(id, "- baseline_category"))
    expect_identical(result$na_impute, expected$na_impute,
      label = paste(id, "- na_impute"))
    expect_equal(result$group, expected$group,
      label = paste(id, "- group"))

    if (!is_impute) {
      expect_equal(result$x, expected$x,
        label = paste(id, "- x (recoded data)"))
      expect_equal(result$missing_index, expected$missing_index,
        label = paste(id, "- missing_index"))
    } else {
      expect_equal(nrow(result$x), nrow(expected$x),
        label = paste(id, "- x nrow (impute)"))
      expect_equal(ncol(result$x), ncol(expected$x),
        label = paste(id, "- x ncol (impute)"))
      expect_equal(nrow(result$missing_index), nrow(expected$missing_index),
        label = paste(id, "- missing_index nrow (impute)"))
      expect_false(anyNA(result$x),
        label = paste(id, "- no NAs after impute"))
    }
  }
})

# ==============================================================================
# Structural sanity checks on the fixture set
# ==============================================================================

test_that("fixture manifest has expected number of cases", {
  skip_if_no_fixtures()
  manifest <- readRDS(file.path(fixture_dir, "manifest.rds"))
  expect_gte(nrow(manifest), 15)
})

test_that("fixture manifest covers both bgm and compare types", {
  skip_if_no_fixtures()
  manifest <- readRDS(file.path(fixture_dir, "manifest.rds"))
  expect_true("bgm" %in% manifest$type)
  expect_true("compare" %in% manifest$type)
})

test_that("all fixture files listed in manifest exist on disk", {
  skip_if_no_fixtures()
  manifest <- readRDS(file.path(fixture_dir, "manifest.rds"))
  for (id in manifest$id) {
    path <- file.path(fixture_dir, paste0(id, ".rds"))
    expect_true(file.exists(path), label = paste("File exists:", id))
  }
})
