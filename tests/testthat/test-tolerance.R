## Debug helper:
## Run with this command to get more context when something fails:
## testthat::test_file("tests/testthat/test-tolerance.R", reporter = "progress")

test_that("bgms outputs are numerically sane (stochastic-robust)", {
  # ---------------------------------------------------------------------------
  # Purpose of this test
  # ---------------------------------------------------------------------------
  # These are *tolerance / sanity* tests (not correctness tests).
  # We avoid exact-value assertions because bgms is stochastic.
  #
  # We assert robust properties:
  #   - Range invariants (probabilities in [0,1], correlations in [-1,1])
  #   - Symmetry for pairwise matrices
  #   - Coarse aggregates within wide bounds
  # ---------------------------------------------------------------------------
  
  withr::local_seed(123)
  
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:5]
  p   <- ncol(dat)
  
  upper_vals <- function(M) M[upper.tri(M)]
  
  specs <- list(
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
      ),
      checks = list(
        # indicator sanity
        function(res, ctx) {
          fld <- "posterior_mean_indicator"
          M <- res[[fld]]
          
          actual_dim <- if (!is.null(dim(M))) paste(dim(M), collapse = "x") else "NULL"
          
          expect_true(is.matrix(M), info = sprintf("%s %s is not a matrix", ctx, fld))
          expect_equal(
            dim(M), c(p, p),
            info = sprintf("%s %s wrong dim: expected %ix%i, got %s", ctx, fld, p, p, actual_dim)
          )
          expect_false(all(is.na(M)), info = sprintf("%s %s is all NA", ctx, fld))
          
          expect_true(
            all(is.na(M) | (M >= 0 & M <= 1)),
            info = sprintf("%s %s has values outside [0,1]", ctx, fld)
          )
        },
        
        # pairwise sanity + symmetry
        function(res, ctx) {
          fld <- "posterior_mean_pairwise"
          M <- res[[fld]]
          
          actual_dim <- if (!is.null(dim(M))) paste(dim(M), collapse = "x") else "NULL"
          
          expect_true(is.matrix(M), info = sprintf("%s %s is not a matrix", ctx, fld))
          expect_equal(
            dim(M), c(p, p),
            info = sprintf("%s %s wrong dim: expected %ix%i, got %s", ctx, fld, p, p, actual_dim)
          )
          expect_false(all(is.na(M)), info = sprintf("%s %s is all NA", ctx, fld))
          
          expect_true(
            all(is.na(M) | (M >= -1 & M <= 1)),
            info = sprintf("%s %s has values outside [-1,1]", ctx, fld)
          )
          
          asym <- max(abs(M - t(M)), na.rm = TRUE)
          expect_true(
            is.finite(asym),
            info = sprintf("%s %s asymmetry not finite", ctx, fld)
          )
          expect_true(
            asym <= 1e-8,
            info = sprintf("%s %s asymmetry too large: %g", ctx, fld, asym)
          )
        },
        
        # coarse aggregate for pairwise (wide bounds; calibrate if you want tighter)
        function(res, ctx) {
          fld <- "posterior_mean_pairwise"
          M <- res[[fld]]
          vals <- abs(upper_vals(M))
          stat <- mean(vals, na.rm = TRUE)
          
          expect_true(
            is.finite(stat),
            info = sprintf("%s %s mean(|upper|) not finite", ctx, fld)
          )
          expect_true(
            stat >= 0.00,
            info = sprintf("%s %s mean(|upper|) too small: %0.3f", ctx, fld, stat)
          )
          expect_true(
            stat <= 0.80,
            info = sprintf("%s %s mean(|upper|) too large: %0.3f", ctx, fld, stat)
          )
        }
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
      ),
      checks = list(
        # baseline pairwise sanity + symmetry
        function(res, ctx) {
          fld <- "posterior_mean_pairwise_baseline"
          M <- res[[fld]]
          
          actual_dim <- if (!is.null(dim(M))) paste(dim(M), collapse = "x") else "NULL"
          
          expect_true(is.matrix(M), info = sprintf("%s %s is not a matrix", ctx, fld))
          expect_equal(
            dim(M), c(p, p),
            info = sprintf("%s %s wrong dim: expected %ix%i, got %s", ctx, fld, p, p, actual_dim)
          )
          expect_false(all(is.na(M)), info = sprintf("%s %s is all NA", ctx, fld))
          
          asym <- max(abs(M - t(M)), na.rm = TRUE)
          expect_true(
            is.finite(asym),
            info = sprintf("%s %s asymmetry not finite", ctx, fld)
          )
          expect_true(
            asym <= 1e-8,
            info = sprintf("%s %s asymmetry too large: %g", ctx, fld, asym)
          )
        },
        
        # coarse aggregate for baseline pairwise (wide bounds; calibrate if you want tighter)
        function(res, ctx) {
          fld <- "posterior_mean_pairwise_baseline"
          M <- res[[fld]]
          vals <- abs(upper_vals(M))
          stat <- mean(vals, na.rm = TRUE)
          
          expect_true(
            is.finite(stat),
            info = sprintf("%s %s mean(|upper|) not finite", ctx, fld)
          )
          expect_true(
            stat >= 0.00,
            info = sprintf("%s %s mean(|upper|) too small: %0.3f", ctx, fld, stat)
          )
          expect_true(
            stat <= 0.80,
            info = sprintf("%s %s mean(|upper|) too large: %0.3f", ctx, fld, stat)
          )
        }
      )
    )
  )
  
  for (spec in specs) {
    ctx <- sprintf("[%s / %s]", spec$fun_label, spec$label)
    res <- do.call(spec$fun, spec$args)
    for (chk in spec$checks) chk(res, ctx)
  }
})
