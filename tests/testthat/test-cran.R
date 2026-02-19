## Debug helper:
## Run with this command to get more context when something fails:
## testthat::test_file("tests/testthat/test-cran.R", reporter = "progress")

test_that("bgms functions satisfy user-visible return contract", {
  # ---------------------------------------------------------------------------
  # Purpose of this test
  # ---------------------------------------------------------------------------
  # These tests are *API/contract tests*, not statistical correctness tests.
  #
  # We deliberately do NOT test whether the MCMC results are "right" (posterior
  # means, convergence, etc.) because:
  #   - MCMC is stochastic and output will vary across platforms / RNG versions.
  #   - CRAN checks should be stable and fast.
  #
  # Instead, we enforce a "user-visible return contract":
  #   1) The returned object contains expected named fields (stable public API).
  #   2) Key matrix outputs have the correct dimensions (p x p).
  #   3) Those matrices are not entirely NA (sanity check: something was produced).
  #
  # If these fail, downstream user code is likely to break even if the sampler
  # "runs".
  # ---------------------------------------------------------------------------
  
  # Make test reproducible and keep runtime reasonable.
  set.seed(123)
  
  # Use a small fixed subset of an included dataset so:
  #   - tests run quickly,
  #   - dimensional expectations are deterministic,
  #   - results are independent of external data availability.
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:5]
  p   <- ncol(dat)  # expected dimension for pairwise matrices
  
  # ---------------------------------------------------------------------------
  # Test specifications
  # ---------------------------------------------------------------------------
  # We test multiple entry points using the same assertion logic.
  # Each spec defines:
  #   - how to call the function (fun + args)
  #   - which fields must exist in the returned object
  #   - which fields must be p x p matrices
  tests <- list(
    list(
      label      = "single_bgm",
      fun_label  = "bgm",
      fun        = bgms::bgm,
      args       = list(
        x                = dat,
        iter             = 1000,
        warmup           = 1000,
        chains           = 2,
        edge_selection   = TRUE,
        edge_prior       = "Bernoulli",
        na_action        = "listwise",
        update_method    = "adaptive-metropolis",
        display_progress = "none"
      ),
      expected_fields = c(
        "posterior_summary_main",
        "posterior_summary_pairwise",
        "posterior_summary_indicator",
        "posterior_mean_main",
        "posterior_mean_pairwise",
        "posterior_mean_indicator",
        "raw_samples",
        "arguments"
      ),
      matrix_fields = c(
        "posterior_mean_pairwise",
        "posterior_mean_indicator"
      )
    ),
    list(
      label      = "compare_bgm",
      fun_label  = "bgmCompare",
      fun        = bgms::bgmCompare,
      args       = list(
        x                    = dat,
        group_indicator      = rep(1:2, each = 20),
        iter                 = 1000,
        warmup               = 1000,
        chains               = 2,
        difference_selection = FALSE,
        na_action            = "listwise",
        update_method        = "adaptive-metropolis",
        display_progress     = "none"
      ),
      expected_fields = c(
        "posterior_summary_main_baseline",
        "posterior_summary_pairwise_baseline",
        "posterior_summary_main_differences",
        "posterior_summary_pairwise_differences",
        "posterior_mean_main_baseline",
        "posterior_mean_pairwise_baseline",
        "raw_samples",
        "arguments"
      ),
      matrix_fields = c(
        "posterior_mean_pairwise_baseline"
      )
    )
  )
  
  # ---------------------------------------------------------------------------
  # Execute specs and assert contract
  # ---------------------------------------------------------------------------
  for (spec in tests) {
    
    # ACT: run the function under test with a controlled configuration.
    result <- do.call(spec$fun, spec$args)
    
    # ASSERT 1: required top-level fields exist.
    # This is the primary "public API doesn't change silently" check.
    missing <- setdiff(spec$expected_fields, names(result))
    expect_equal(
      length(missing), 0,
      info = sprintf(
        "[%s / %s] Missing fields: %s",
        spec$fun_label,
        spec$label,
        if (length(missing) == 0) "<none>" else paste(missing, collapse = ", ")
      )
    )
    
    # ASSERT 2: selected outputs are matrices with correct shape (p x p),
    # and contain at least some non-NA values.
    #
    # Note: we skip fields that are not present to avoid cascading errors
    # (but missing fields are already caught above).
    for (fld in spec$matrix_fields) {
      if (!fld %in% names(result)) next
      
      actual_dim <- if (!is.null(result[[fld]])) {
        paste(dim(result[[fld]]), collapse = "x")
      } else {
        "NULL"
      }
      
      expect_equal(
        dim(result[[fld]]),
        c(p, p),
        info = sprintf(
          "[%s / %s / %s] Wrong dim: expected %ix%i, got %s",
          spec$fun_label,
          spec$label,
          fld,
          p, p,
          actual_dim
        )
      )
      
      expect_false(
        all(is.na(result[[fld]])),
        info = sprintf(
          "[%s / %s / %s] All entries are NA",
          spec$fun_label,
          spec$label,
          fld
        )
      )
    }
  }
})