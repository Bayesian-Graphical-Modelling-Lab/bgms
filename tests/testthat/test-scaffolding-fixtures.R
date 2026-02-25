# Null-coalescing operator for R < 4.4 compatibility
if(!exists("%||%", baseenv())) {
  `%||%` = function(x, y) if(is.null(x)) y else x
}

# ==============================================================================
# Golden-Snapshot Fixture Verification Tests
# ==============================================================================
#
# Phase A-0 of the scaffolding refactor (dev/scaffolding/plan.md).
#
# The original tests verified that check_model(), check_compare_model(),
# reformat_data(), and compare_reformat_data() reproduced golden-snapshot
# fixtures exactly. Those monolithic functions have now been deleted (Phase C.8)
# and their logic inlined into bgm_spec(). The golden-snapshot guardrails
# served their purpose during the transition.
#
# What remains: structural checks on the fixture set itself.
# ==============================================================================

fixture_dir = file.path(testthat::test_path("..", ".."), "dev", "fixtures", "scaffolding")

# Skip all tests if fixtures haven't been generated yet
skip_if_no_fixtures = function() {
  # When running via devtools::test(), the working directory is tests/testthat/
  # The fixtures live at the package root under dev/fixtures/scaffolding/
  pkg_root = testthat::test_path("..", "..")
  fixture_dir = file.path(pkg_root, "dev", "fixtures", "scaffolding")
  manifest_path = file.path(fixture_dir, "manifest.rds")
  if(!file.exists(manifest_path)) {
    skip("Scaffolding fixtures not found. Run: Rscript dev/generate_scaffolding_fixtures.R")
  }
}

# Helper: load a fixture by id
load_fixture = function(id) {
  path = file.path(fixture_dir, paste0(id, ".rds"))
  if(!file.exists(path)) {
    skip(paste("Fixture not found:", id))
  }
  readRDS(path)
}

# ==============================================================================
# Structural sanity checks on the fixture set
# ==============================================================================

test_that("fixture manifest has expected number of cases", {
  skip_if_no_fixtures()
  manifest = readRDS(file.path(fixture_dir, "manifest.rds"))
  expect_gte(nrow(manifest), 15)
})

test_that("fixture manifest covers both bgm and compare types", {
  skip_if_no_fixtures()
  manifest = readRDS(file.path(fixture_dir, "manifest.rds"))
  expect_true("bgm" %in% manifest$type)
  expect_true("compare" %in% manifest$type)
})

test_that("all fixture files listed in manifest exist on disk", {
  skip_if_no_fixtures()
  manifest = readRDS(file.path(fixture_dir, "manifest.rds"))
  for(id in manifest$id) {
    path = file.path(fixture_dir, paste0(id, ".rds"))
    expect_true(file.exists(path), label = paste("File exists:", id))
  }
})
