
# Debug helper:
# Run with this command to see progress + richer failure output:
# testthat::test_file("tests/testthat/test-snapshots.R", reporter = "progress")

test_that("bgms snapshots match expected outputs", {
  # ---------------------------------------------------------------------------
  # Purpose of this test
  # ---------------------------------------------------------------------------
  # This is a *snapshot / regression* test.
  #
  # Goal: detect unintended changes in the *returned object structure* or
  # printed representation over time (e.g., after refactors).
  #
  # Important: this is NOT a correctness test for the MCMC posterior.
  # Snapshot tests are inherently sensitive to small changes, so we:
  #   - skip on CRAN (platform/R-version differences can cause noise),
  #   - fix the RNG seed,
  #   - use a small, deterministic dataset subset.
  # ---------------------------------------------------------------------------
  
  # CRAN: snapshot tests are too brittle across platforms and often too slow.
  skip_on_cran()
  
  # Reproducibility: ensures stable draws for snapshot comparisons (within the
  # same R version / platform).
  set.seed(123)
  
  # Use a bundled dataset and fixed subset for speed + determinism.
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:5]
  
  # ---------------------------------------------------------------------------
  # Test specifications
  # ---------------------------------------------------------------------------
  # We snapshot multiple entry points using a shared harness.
  # Each spec defines how to call the function (fun + args) and how we label it
  # in the snapshot output.
  tests <- list(
    list(
      label     = "single_bgm",
      fun_label = "bgm",
      fun       = bgms::bgm,
      args      = list(
        x                = dat,
        iter             = 1000,
        warmup           = 1000,
        chains           = 2,
        edge_selection   = TRUE,
        edge_prior       = "Bernoulli",
        na_action        = "listwise",
        update_method    = "adaptive-metropolis",
        display_progress = "none"
      )
    ),
    list(
      label     = "compare_bgm",
      fun_label = "bgmCompare",
      fun       = bgms::bgmCompare,
      args      = list(
        x                    = dat,
        group_indicator      = rep(1:2, each = 20),
        iter                 = 1000,
        warmup               = 1000,
        chains               = 2,
        difference_selection = FALSE,
        na_action            = "listwise",
        update_method        = "adaptive-metropolis",
        display_progress     = "none"
      )
    )
  )
  
  # ---------------------------------------------------------------------------
  # Collect a snapshot-friendly representation
  # ---------------------------------------------------------------------------
  # We snapshot a structured list rather than the raw result object directly.
  # Reasons:
  #   - json2 snapshots are readable in diffs and stable to compare.
  #   - some bgms objects may not serialize cleanly as JSON.
  #   - capturing class + names catches common API regressions early.
  #
  # We also include a full textual dput() of the object for maximal coverage:
  # it will flag deep structural changes even if class/names stay the same.
  out <- list()
  
  for (spec in tests) {
    # ACT: run the function under test.
    result <- do.call(spec$fun, spec$args)
    
    # Arrange snapshot payload for this spec.
    out[[spec$label]] <- list(
      fun_label    = spec$fun_label,
      
      # Record class in a minimal form for stable diffs.
      result_class = unclass(class(result)),
      
      # Record top-level names to catch API changes (added/removed/renamed fields).
      result_names = names(result),
      
      # Full object captured as text so expect_snapshot_value(..., json2) can
      # store it. This is intentionally verbose: it’s the strongest guardrail
      # against subtle return-structure regressions.
      full_result_dput = paste(capture.output(dput(result)), collapse = "\n")
    )
  }
  
  # ASSERT: compare against stored snapshot in tests/testthat/_snaps/.
  # Update snapshots intentionally with:
  # testthat::snapshot_accept()
  expect_snapshot_value(out, style = "json2")
})
