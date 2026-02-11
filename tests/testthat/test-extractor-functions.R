# ==============================================================================
# Contract Tests for Extractor Functions (Parameterized)
# ==============================================================================
#
# EXTENDS: test-tolerance.R (stochastic-robust testing approach)
# PATTERN: Range invariants, symmetry checks, dimension consistency
#
# This file uses parameterized testing (specs + loop) to reduce code repetition.
# Each extractor is tested across multiple fixture types with shared assertions.
#
# IMPORTANT: Changes to extractor function output structure may break easybgm!
# ==============================================================================

# ------------------------------------------------------------------------------
# Fixture Specifications
# ------------------------------------------------------------------------------
# Define all fixtures to test against, with their properties

get_all_fixtures <- function() {
  list(
    list(
      label = "bgms_binary",
      get_fit = get_bgms_fit,
      type = "bgms",
      var_type = "binary"
    ),
    list(
      label = "bgms_ordinal",
      get_fit = get_bgms_fit_ordinal,
      type = "bgms",
      var_type = "ordinal"
    ),
    list(
      label = "bgms_blumecapel",
      get_fit = get_bgms_fit_blumecapel,
      type = "bgms",
      var_type = "blume-capel"
    ),
    list(
      label = "bgmCompare_binary",
      get_fit = get_bgmcompare_fit,
      type = "bgmCompare",
      var_type = "binary"
    ),
    list(
      label = "bgmCompare_ordinal",
      get_fit = get_bgmcompare_fit_ordinal,
      type = "bgmCompare",
      var_type = "ordinal"
    )
  )
}


# ------------------------------------------------------------------------------
# extract_arguments() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_arguments returns complete argument list for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    # Basic structure
    expect_true(is.list(args), info = paste(ctx, "should be list"))
    expect_true(length(args) > 0, info = ctx)
    
    # Essential fields present
    expect_true("num_variables" %in% names(args), info = paste(ctx, "missing num_variables"))
    expect_true("num_cases" %in% names(args), info = paste(ctx, "missing num_cases"))
    expect_true("data_columnnames" %in% names(args), info = paste(ctx, "missing data_columnnames"))
    
    # Values are sensible
    expect_true(args$num_variables >= 1, info = ctx)
    expect_true(args$num_cases >= 1, info = ctx)
    
    # Type-specific fields
    if (spec$type == "bgms") {
      expect_true(is.logical(args$edge_selection), info = paste(ctx, "edge_selection should be logical"))
    } else {
      expect_true(args$num_groups >= 2, info = paste(ctx, "bgmCompare should have >= 2 groups"))
    }
  }
})

test_that("extract_arguments errors on non-bgms objects", {
  expect_error(extract_arguments(list()), class = "error")
  expect_error(extract_arguments(data.frame()), class = "error")
})


# ------------------------------------------------------------------------------
# extract_pairwise_interactions() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_pairwise_interactions returns valid matrix for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    pairwise <- extract_pairwise_interactions(fit)
    
    # Structure checks
    expect_true(is.matrix(pairwise), info = paste(ctx, "should be matrix"))
    
    p <- args$num_variables
    expected_cols <- p * (p - 1) / 2
    expect_equal(ncol(pairwise), expected_cols, 
                 info = paste(ctx, "wrong number of edge columns"))
    
    # Values finite
    expect_true(all(is.finite(pairwise)), info = paste(ctx, "should have finite values"))
    
    # Has column names
    expect_true(!is.null(colnames(pairwise)), info = paste(ctx, "should have column names"))
  }
})


# ------------------------------------------------------------------------------
# extract_category_thresholds() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_category_thresholds returns valid output for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    thresholds <- extract_category_thresholds(fit)
    
    # Structure checks
    expect_true(is.matrix(thresholds), info = paste(ctx, "should be matrix"))
    
    # Values finite where not NA
    vals <- thresholds[!is.na(thresholds)]
    expect_true(all(is.finite(vals)), info = paste(ctx, "non-NA values should be finite"))
  }
})


# ------------------------------------------------------------------------------
# extract_indicators() and extract_posterior_inclusion_probabilities() Tests
# ------------------------------------------------------------------------------
# These only apply to fits with edge_selection = TRUE

test_that("extract_indicators returns binary matrix for edge-selection fits", {
  # Only test fixtures with edge selection
  fixtures <- list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (!isTRUE(args$edge_selection)) {
      next
    }
    
    indicators <- extract_indicators(fit)
    
    # Structure
    expect_true(is.matrix(indicators), info = ctx)
    
    p <- args$num_variables
    expected_cols <- p * (p - 1) / 2
    expect_equal(ncol(indicators), expected_cols, info = paste(ctx, "wrong indicator columns"))
    
    # Binary values
    expect_true(all(indicators %in% c(0, 1)), 
                info = paste(ctx, "indicators should be 0 or 1"))
  }
})

test_that("extract_posterior_inclusion_probabilities returns symmetric PIP matrix", {
  # Only test fixtures with edge selection
  fixtures <- list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (!isTRUE(args$edge_selection)) {
      next
    }
    
    pip <- extract_posterior_inclusion_probabilities(fit)
    p <- args$num_variables
    
    # Structure
    expect_true(is.matrix(pip), info = ctx)
    expect_equal(dim(pip), c(p, p), info = paste(ctx, "should be p x p"))
    
    # Symmetry
    expect_true(is_symmetric(pip), info = paste(ctx, "should be symmetric"))
    
    # Range [0, 1]
    expect_true(values_in_range(pip, 0, 1), info = paste(ctx, "PIPs should be in [0,1]"))
    
    # Diagonal is zero (no self-loops)
    expect_true(all(diag(pip) == 0), info = paste(ctx, "diagonal should be 0"))
    
    # Has variable names
    expect_equal(colnames(pip), args$data_columnnames, info = ctx)
  }
})

test_that("extract_indicators errors when edge_selection = FALSE", {
  data <- generate_test_data(n = 20, p = 3)
  args <- c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit <- do.call(bgm, args)
  
  expect_error(extract_indicators(fit), regexp = "edge_selection")
})


# ------------------------------------------------------------------------------
# extract_rhat() and extract_ess() Tests (parameterized)
# ------------------------------------------------------------------------------

test_that("extract_rhat returns valid diagnostics for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    rhat <- extract_rhat(fit)
    
    expect_true(is.list(rhat), info = paste(ctx, "should be list"))
    
    if (spec$type == "bgms") {
      expect_true("pairwise" %in% names(rhat), info = paste(ctx, "missing pairwise"))
      expect_true(is.numeric(rhat$pairwise), info = ctx)
      expect_true(all(is.na(rhat$pairwise) | rhat$pairwise > 0), 
                  info = paste(ctx, "R-hat should be positive"))
    } else {
      expect_true("pairwise_baseline" %in% names(rhat), info = paste(ctx, "missing pairwise_baseline"))
      expect_true(is.numeric(rhat$pairwise_baseline), info = ctx)
      expect_true(all(is.na(rhat$pairwise_baseline) | rhat$pairwise_baseline > 0), info = ctx)
    }
  }
})

test_that("extract_ess returns valid diagnostics for all fit types", {
  fixtures <- get_all_fixtures()
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    ess <- extract_ess(fit)
    
    expect_true(is.list(ess), info = paste(ctx, "should be list"))
    
    if (spec$type == "bgms") {
      expect_true("pairwise" %in% names(ess), info = paste(ctx, "missing pairwise"))
      expect_true(is.numeric(ess$pairwise), info = ctx)
      expect_true(all(is.na(ess$pairwise) | ess$pairwise > 0), 
                  info = paste(ctx, "ESS should be positive"))
    } else {
      expect_true("pairwise_baseline" %in% names(ess), info = paste(ctx, "missing pairwise_baseline"))
      expect_true(is.numeric(ess$pairwise_baseline), info = ctx)
    }
  }
})

test_that("extract_rhat and extract_ess error on non-bgms objects", {
  expect_error(extract_rhat(list()), class = "error")
  expect_error(extract_rhat(data.frame()), class = "error")
  expect_error(extract_ess(list()), class = "error")
  expect_error(extract_ess(data.frame()), class = "error")
})


# ------------------------------------------------------------------------------
# extract_indicator_priors() Tests
# ------------------------------------------------------------------------------

test_that("extract_indicator_priors returns prior specification", {
  fit <- get_bgms_fit()
  args <- extract_arguments(fit)
  
  if (!isTRUE(args$edge_selection)) {
    skip("Fit object does not have edge_selection = TRUE")
  }
  
  priors <- extract_indicator_priors(fit)
  
  expect_type(priors, "list")
  expect_true("type" %in% names(priors))
  
  valid_types <- c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block")
  expect_true(priors$type %in% valid_types)
  
  # Type-specific checks
  if (priors$type == "Bernoulli") {
    expect_true("prior_inclusion_probability" %in% names(priors))
    pip <- priors$prior_inclusion_probability
    expect_true(all(pip >= 0 & pip <= 1))
  }
  
  if (priors$type == "Beta-Bernoulli") {
    expect_true(all(c("alpha", "beta") %in% names(priors)))
    expect_true(priors$alpha > 0 && priors$beta > 0)
  }
})

test_that("extract_indicator_priors errors when no selection performed", {
  data <- generate_test_data(n = 20, p = 3)
  args <- c(list(x = data, edge_selection = FALSE), quick_mcmc_args())
  fit <- do.call(bgm, args)
  
  expect_error(extract_indicator_priors(fit), regexp = "selection")
})


# ------------------------------------------------------------------------------
# bgmCompare-specific Tests
# ------------------------------------------------------------------------------

test_that("extract_group_params returns group-level parameters", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)
  
  group_params <- extract_group_params(fit)
  
  expect_type(group_params, "list")
  expect_true("main_effects_groups" %in% names(group_params))
  expect_true("pairwise_effects_groups" %in% names(group_params))
  
  # Dimensions match number of groups
  n_groups <- args$num_groups
  expect_equal(ncol(group_params$main_effects_groups), n_groups)
  expect_equal(ncol(group_params$pairwise_effects_groups), n_groups)
  
  # Values finite
  expect_true(all(is.finite(group_params$main_effects_groups)))
  expect_true(all(is.finite(group_params$pairwise_effects_groups)))
})


# ------------------------------------------------------------------------------
# Cross-Function Consistency Tests
# ------------------------------------------------------------------------------

test_that("extractor outputs are dimensionally consistent", {
  fixtures <- list(
    list(label = "bgms_binary", get_fit = get_bgms_fit)
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (!isTRUE(args$edge_selection)) {
      next
    }
    
    p <- args$num_variables
    n_edges <- p * (p - 1) / 2
    
    # All should agree on number of variables/edges
    pip <- extract_posterior_inclusion_probabilities(fit)
    expect_equal(nrow(pip), p, info = paste(ctx, "PIP rows"))
    
    indicators <- extract_indicators(fit)
    expect_equal(ncol(indicators), n_edges, info = paste(ctx, "indicator cols"))
    
    pairwise <- extract_pairwise_interactions(fit)
    expect_equal(ncol(pairwise), n_edges, info = paste(ctx, "pairwise cols"))
    
    thresholds <- extract_category_thresholds(fit)
    expect_equal(nrow(thresholds), p, info = paste(ctx, "threshold rows"))
  }
})


# ------------------------------------------------------------------------------
# Contract Tests for easybgm Integration
# ------------------------------------------------------------------------------

test_that("bgms fit contains all fields accessed by easybgm", {
  fixtures <- list(
    list(label = "bgms", get_fit = get_bgms_fit, type = "bgms"),
    list(label = "bgmCompare", get_fit = get_bgmcompare_fit, type = "bgmCompare")
  )
  
  for (spec in fixtures) {
    ctx <- sprintf("[%s]", spec$label)
    fit <- spec$get_fit()
    args <- extract_arguments(fit)
    
    if (spec$type == "bgms") {
      expect_true("posterior_summary_pairwise" %in% names(fit), info = ctx)
      expect_true(is.data.frame(fit$posterior_summary_pairwise), info = ctx)
      expect_true("Rhat" %in% names(fit$posterior_summary_pairwise), info = ctx)
      expect_true("n_eff" %in% names(fit$posterior_summary_pairwise), info = ctx)
      
      if (isTRUE(args$edge_selection)) {
        expect_true("posterior_summary_indicator" %in% names(fit), info = ctx)
        expect_true("n_eff" %in% names(fit$posterior_summary_indicator), info = ctx)
      }
    } else {
      expect_true("posterior_summary_pairwise_baseline" %in% names(fit), info = ctx)
      expect_true(is.data.frame(fit$posterior_summary_pairwise_baseline), info = ctx)
      expect_true("Rhat" %in% names(fit$posterior_summary_pairwise_baseline), info = ctx)
      expect_true("n_eff" %in% names(fit$posterior_summary_pairwise_baseline), info = ctx)
    }
  }
})


# ------------------------------------------------------------------------------
# extract_indicators.bgmCompare Tests
# ------------------------------------------------------------------------------

test_that("extract_indicators.bgmCompare returns indicator matrix for difference selection fits", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)
  
  if (!isTRUE(args$difference_selection)) {
    skip("Fit object does not have difference_selection = TRUE")
  }
  
  indicators <- extract_indicators(fit)
  
  # Structure
  expect_true(is.matrix(indicators), info = "should be matrix")
  
  # Binary values
  expect_true(all(indicators %in% c(0, 1)), 
              info = "indicators should be 0 or 1")
  
  # Has column names
  expect_true(!is.null(colnames(indicators)), info = "should have column names")
})

test_that("extract_indicators.bgmCompare errors when difference_selection = FALSE", {
  data <- generate_grouped_test_data(n_per_group = 15, p = 3)
  args <- c(
    list(x = data$x, group_indicator = data$group_indicator, difference_selection = FALSE),
    quick_mcmc_args()
  )
  fit <- do.call(bgmCompare, args)
  
  expect_error(extract_indicators(fit), regexp = "difference_selection")
})


# ------------------------------------------------------------------------------
# extract_posterior_inclusion_probabilities.bgmCompare Tests
# ------------------------------------------------------------------------------

test_that("extract_posterior_inclusion_probabilities.bgmCompare returns symmetric PIP matrix", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)
  
  if (!isTRUE(args$difference_selection)) {
    skip("Fit object does not have difference_selection = TRUE")
  }
  
  pip <- extract_posterior_inclusion_probabilities(fit)
  p <- args$num_variables
  
  # Structure
  expect_true(is.matrix(pip), info = "should be matrix")
  expect_equal(dim(pip), c(p, p), info = "should be p x p")
  
  # Symmetry
  expect_true(is_symmetric(pip), info = "should be symmetric")
  
  # Range [0, 1]
  expect_true(values_in_range(pip, 0, 1), info = "PIPs should be in [0,1]")
  
  # Has variable names
  expect_equal(colnames(pip), args$data_columnnames, info = "should have column names")
})

test_that("extract_posterior_inclusion_probabilities.bgmCompare errors when difference_selection = FALSE", {
  data <- generate_grouped_test_data(n_per_group = 15, p = 3)
  args <- c(
    list(x = data$x, group_indicator = data$group_indicator, difference_selection = FALSE),
    quick_mcmc_args()
  )
  fit <- do.call(bgmCompare, args)
  
  expect_error(extract_posterior_inclusion_probabilities(fit), regexp = "difference_selection")
})


# ------------------------------------------------------------------------------
# extract_indicator_priors.bgmCompare Tests
# ------------------------------------------------------------------------------

test_that("extract_indicator_priors.bgmCompare returns prior specification", {
  fit <- get_bgmcompare_fit()
  args <- extract_arguments(fit)
  
  if (!isTRUE(args$difference_selection)) {
    skip("Fit object does not have difference_selection = TRUE")
  }
  
  priors <- extract_indicator_priors(fit)
  
  # Returns the difference_prior from arguments
  expect_true(!is.null(priors), info = "should return prior specification")
})

test_that("extract_indicator_priors.bgmCompare errors when difference_selection = FALSE", {
  data <- generate_grouped_test_data(n_per_group = 15, p = 3)
  args <- c(
    list(x = data$x, group_indicator = data$group_indicator, difference_selection = FALSE),
    quick_mcmc_args()
  )
  fit <- do.call(bgmCompare, args)
  
  expect_error(extract_indicator_priors(fit), regexp = "selection")
})


# ------------------------------------------------------------------------------
# main_difference_selection Tests
# ------------------------------------------------------------------------------

test_that("bgmCompare with main_difference_selection = TRUE produces valid output", {
  fit <- get_bgmcompare_fit_main_selection()
  args <- extract_arguments(fit)
  
  # Verify main_difference_selection is TRUE in arguments

  expect_true(isTRUE(args$main_difference_selection), 
              info = "main_difference_selection should be TRUE")
  expect_true(isTRUE(args$difference_selection), 
              info = "difference_selection should be TRUE")
})

test_that("extract_indicators works with main_difference_selection = TRUE", {
  fit <- get_bgmcompare_fit_main_selection()
  args <- extract_arguments(fit)
  
  indicators <- extract_indicators(fit)
  
  # Structure
  expect_true(is.matrix(indicators), info = "should be matrix")
  
  # Binary values
  expect_true(all(indicators %in% c(0, 1)), 
              info = "indicators should be 0 or 1")
  
  # With main_difference_selection = TRUE, there should be more indicator columns
  # than just pairwise (includes main effect indicators)
  p <- args$num_variables
  n_pairwise <- p * (p - 1) / 2
  
  # Indicator dimensions should include main effects + pairwise
  # Exact count depends on number of categories, but should be > n_pairwise
  expect_true(ncol(indicators) >= n_pairwise, 
              info = "indicators should include at least pairwise effects")
})

test_that("extract_posterior_inclusion_probabilities works with main_difference_selection = TRUE", {
  fit <- get_bgmcompare_fit_main_selection()
  args <- extract_arguments(fit)
  
  pip <- extract_posterior_inclusion_probabilities(fit)
  p <- args$num_variables
  
  # Structure - should be p x p matrix
  expect_true(is.matrix(pip), info = "should be matrix")
  expect_equal(dim(pip), c(p, p), info = "should be p x p")
  
  # Symmetry
  expect_true(is_symmetric(pip), info = "should be symmetric")
  
  # Range [0, 1]
  expect_true(values_in_range(pip, 0, 1), info = "PIPs should be in [0,1]")
})

test_that("extract_group_params works with main_difference_selection = TRUE", {
  fit <- get_bgmcompare_fit_main_selection()
  args <- extract_arguments(fit)
  
  group_params <- extract_group_params(fit)
  
  expect_type(group_params, "list")
  expect_true("main_effects_groups" %in% names(group_params))
  expect_true("pairwise_effects_groups" %in% names(group_params))
  
  # Dimensions match number of groups
  n_groups <- args$num_groups
  expect_equal(ncol(group_params$main_effects_groups), n_groups)
  expect_equal(ncol(group_params$pairwise_effects_groups), n_groups)
  
  # Values finite
  expect_true(all(is.finite(group_params$main_effects_groups)))
  expect_true(all(is.finite(group_params$pairwise_effects_groups)))
})


# ------------------------------------------------------------------------------
# extract_sbm.bgms Tests (Stochastic Block Model)
# ------------------------------------------------------------------------------

test_that("extract_sbm.bgms returns SBM summaries for Stochastic-Block prior", {
  fit <- get_bgms_fit_sbm()
  args <- extract_arguments(fit)
  
  sbm <- extract_sbm(fit)
  
  # Structure
  expect_type(sbm, "list")
  
  # Required fields
  expect_true("posterior_num_blocks" %in% names(sbm), 
              info = "should have posterior_num_blocks")
  expect_true("posterior_mean_allocations" %in% names(sbm), 
              info = "should have posterior_mean_allocations")
  expect_true("posterior_mode_allocations" %in% names(sbm), 
              info = "should have posterior_mode_allocations")
  expect_true("posterior_mean_coclustering_matrix" %in% names(sbm), 
              info = "should have posterior_mean_coclustering_matrix")
  
  # Coclustering matrix should be symmetric
  ccm <- sbm$posterior_mean_coclustering_matrix
  expect_true(is.matrix(ccm), info = "coclustering matrix should be matrix")
  expect_true(is_symmetric(ccm), info = "coclustering matrix should be symmetric")
  
  # Values in [0, 1] for coclustering probabilities
  expect_true(values_in_range(ccm, 0, 1), 
              info = "coclustering probabilities should be in [0,1]")
})

test_that("extract_sbm.bgms errors for non-SBM prior", {
  fit <- get_bgms_fit()  # Uses default Bernoulli prior
  
  expect_error(extract_sbm(fit), regexp = "Stochastic-Block")
})


# ------------------------------------------------------------------------------
# extract_indicator_priors with Beta-Bernoulli Prior Tests
# ------------------------------------------------------------------------------

test_that("extract_indicator_priors returns Beta-Bernoulli parameters", {
  fit <- get_bgms_fit_beta_bernoulli()
  args <- extract_arguments(fit)
  
  priors <- extract_indicator_priors(fit)
  
  # Type check
  expect_type(priors, "list")
  expect_equal(priors$type, "Beta-Bernoulli")
  
  # Required parameters
  expect_true("alpha" %in% names(priors), info = "should have alpha parameter")
  expect_true("beta" %in% names(priors), info = "should have beta parameter")
  
  # Positive values

  expect_true(priors$alpha > 0, info = "alpha should be positive")
  expect_true(priors$beta > 0, info = "beta should be positive")
})

test_that("extract_indicator_priors returns Stochastic-Block parameters", {
  fit <- get_bgms_fit_sbm()
  
  priors <- extract_indicator_priors(fit)
  
  # Type check
  expect_type(priors, "list")
  expect_equal(priors$type, "Stochastic-Block")
  
  # Required parameters
  expect_true("beta_bernoulli_alpha" %in% names(priors), 
              info = "should have beta_bernoulli_alpha")
  expect_true("beta_bernoulli_beta" %in% names(priors), 
              info = "should have beta_bernoulli_beta")
  expect_true("dirichlet_alpha" %in% names(priors), 
              info = "should have dirichlet_alpha")
})

# ==============================================================================
# Legacy Format Support Tests
# ==============================================================================
# These tests verify backward compatibility with fit objects from older bgms versions.
# Legacy fixtures are stored in tests/testthat/fixtures/legacy/ (NOT shipped with package).
#
# To generate fixtures, run: Rscript dev/generate_legacy_fixtures.R
#
# Tests skip on CRAN since fixtures aren't available in installed package.
# ==============================================================================

# Helper to load legacy fixtures from tests/testthat/fixtures/legacy/
get_legacy_fit <- function(name) {
  skip_on_cran()  # Fixtures not shipped with package
  
  # Try relative path first (for devtools::test())
  path <- file.path("fixtures", "legacy", paste0(name, ".rds"))
  if (!file.exists(path)) {
    # Try from package root (for testthat::test_file())
    path <- file.path("tests", "testthat", "fixtures", "legacy", paste0(name, ".rds"))
  }
  if (!file.exists(path)) {
    skip(paste("Legacy fixture not found:", name, "- run dev/generate_legacy_fixtures.R"))
  }
  readRDS(path)
}

# Get all available legacy fixture versions
get_available_legacy_versions <- function() {
  skip_on_cran()
  
  # Try relative path first
  legacy_dir <- file.path("fixtures", "legacy")
  if (!dir.exists(legacy_dir)) {
    legacy_dir <- file.path("tests", "testthat", "fixtures", "legacy")
  }
  if (!dir.exists(legacy_dir)) {
    return(character(0))
  }
  
  files <- list.files(legacy_dir, pattern = "^fit_v.*\\.rds$")
  # Extract version from filename: fit_v0.1.4.2.rds → 0.1.4.2
  gsub("^fit_v(.*)\\.rds$", "\\1", files)
}

# Categorize versions by format era
categorize_version <- function(version) {
  v <- numeric_version(version)
  if (v < "0.1.4") {
    return("pre-0.1.4")      # Defunct: $gamma field → error
  } else if (v < "0.1.6") {
    return("0.1.4-0.1.5")    # Deprecated: $indicator at top level → warning
  } else {
    return("0.1.6+")         # Current: $raw_samples$indicator → no warning
  }
}

# ------------------------------------------------------------------------------
# Lifecycle Behavior Tests
# ------------------------------------------------------------------------------
# Verify correct lifecycle signals:
#   - pre-0.1.4: $gamma is DEFUNCT (error), but $interactions/$thresholds work (deprecated)
#   - 0.1.4-0.1.5: $indicator/$interactions/$thresholds are all deprecated (warning)
#   - 0.1.6+: current format, no warnings

test_that("pre-0.1.4 formats throw defunct errors for indicator extraction", {
  versions <- get_available_legacy_versions()
  pre_014_versions <- versions[vapply(versions, function(v) categorize_version(v) == "pre-0.1.4", logical(1))]
  skip_if(length(pre_014_versions) == 0, "No pre-0.1.4 fixtures available")
  
  for (version in pre_014_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    
    # extract_indicators uses $gamma → defunct
    expect_error(
      extract_indicators(fit),
      "defunct",
      info = paste("v", version, ": extract_indicators should error (defunct $gamma)")
    )
    
    # extract_posterior_inclusion_probabilities uses $gamma → defunct
    expect_error(
      extract_posterior_inclusion_probabilities(fit),
      "defunct",
      info = paste("v", version, ": extract_pip should error (defunct $gamma)")
    )
  }
})

test_that("pre-0.1.4 formats emit deprecation warnings for pairwise/thresholds", {
  versions <- get_available_legacy_versions()
  pre_014_versions <- versions[vapply(versions, function(v) categorize_version(v) == "pre-0.1.4", logical(1))]
  skip_if(length(pre_014_versions) == 0, "No pre-0.1.4 fixtures available")
  
  for (version in pre_014_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    
    # extract_pairwise_interactions uses $interactions → deprecated (not defunct)
    expect_warning(
      extract_pairwise_interactions(fit),
      "deprecated",
      info = paste("v", version, ": extract_pairwise should warn (deprecated $interactions)")
    )
    
    # extract_category_thresholds uses $thresholds → deprecated (not defunct)
    expect_warning(
      extract_category_thresholds(fit),
      "deprecated",
      info = paste("v", version, ": extract_thresholds should warn (deprecated $thresholds)")
    )
  }
})

test_that("0.1.4-0.1.5 formats emit deprecation warnings", {
  versions <- get_available_legacy_versions()
  deprecated_versions <- versions[vapply(versions, function(v) categorize_version(v) == "0.1.4-0.1.5", logical(1))]
  skip_if(length(deprecated_versions) == 0, "No 0.1.4-0.1.5 fixtures available")
  
  for (version in deprecated_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    
    # extract_indicators should warn about deprecated $indicator field
    expect_warning(
      extract_indicators(fit),
      "deprecated",
      info = paste("v", version, "should warn for extract_indicators")
    )
    
    # extract_posterior_inclusion_probabilities should also warn
    expect_warning(
      extract_posterior_inclusion_probabilities(fit),
      "deprecated",
      info = paste("v", version, "should warn for extract_posterior_inclusion_probabilities")
    )
    
    # extract_pairwise_interactions should warn about deprecated $interactions field
    expect_warning(
      extract_pairwise_interactions(fit),
      "deprecated",
      info = paste("v", version, "should warn for extract_pairwise_interactions")
    )
    
    # extract_category_thresholds should warn about deprecated $thresholds field
    expect_warning(
      extract_category_thresholds(fit),
      "deprecated",
      info = paste("v", version, "should warn for extract_category_thresholds")
    )
  }
})

test_that("0.1.6+ formats work without deprecation warnings", {
  versions <- get_available_legacy_versions()
  current_versions <- versions[vapply(versions, function(v) categorize_version(v) == "0.1.6+", logical(1))]
  skip_if(length(current_versions) == 0, "No 0.1.6+ fixtures available")
  
  for (version in current_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    
    # Should not produce any deprecation warnings
    expect_no_warning(extract_indicators(fit))
    expect_no_warning(extract_posterior_inclusion_probabilities(fit))
    expect_no_warning(extract_pairwise_interactions(fit))
    expect_no_warning(extract_category_thresholds(fit))
  }
})

# ------------------------------------------------------------------------------
# Functional Tests for Deprecated (0.1.4-0.1.5) Formats
# ------------------------------------------------------------------------------
# These versions still work (with warnings), so we test correctness

test_that("extract_indicators works with 0.1.4-0.1.5 formats", {
  versions <- get_available_legacy_versions()
  deprecated_versions <- versions[vapply(versions, function(v) categorize_version(v) == "0.1.4-0.1.5", logical(1))]
  skip_if(length(deprecated_versions) == 0, "No 0.1.4-0.1.5 fixtures available")
  
  for (version in deprecated_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    result <- suppressWarnings(extract_indicators(fit))
    
    expect_true(is.matrix(result), 
                info = paste("v", version, "should return matrix"))
    expect_true(nrow(result) > 0, 
                info = paste("v", version, "should have rows"))
    expect_true(ncol(result) > 0, 
                info = paste("v", version, "should have columns"))
    expect_true(all(result %in% c(0, 1)), 
                info = paste("v", version, "should have binary values"))
  }
})

test_that("extract_posterior_inclusion_probabilities works with 0.1.4-0.1.5 formats", {
  versions <- get_available_legacy_versions()
  deprecated_versions <- versions[vapply(versions, function(v) categorize_version(v) == "0.1.4-0.1.5", logical(1))]
  skip_if(length(deprecated_versions) == 0, "No 0.1.4-0.1.5 fixtures available")
  
  for (version in deprecated_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    result <- suppressWarnings(extract_posterior_inclusion_probabilities(fit))
    
    expect_true(is.matrix(result), 
                info = paste("v", version, "should return matrix"))
    expect_true(isSymmetric(result), 
                info = paste("v", version, "should be symmetric"))
    expect_true(all(result >= 0 & result <= 1), 
                info = paste("v", version, "should be in [0,1]"))
    expect_true(all(diag(result) == 0), 
                info = paste("v", version, "should have zero diagonal"))
  }
})

test_that("extract_pairwise_interactions works with pre-0.1.6 formats", {
  versions <- get_available_legacy_versions()
  # Both pre-0.1.4 and 0.1.4-0.1.5 have $interactions (deprecated but functional)
  pre_016_versions <- versions[vapply(versions, function(v) categorize_version(v) != "0.1.6+", logical(1))]
  skip_if(length(pre_016_versions) == 0, "No pre-0.1.6 fixtures available")
  
  for (version in pre_016_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    result <- suppressWarnings(extract_pairwise_interactions(fit))
    
    expect_true(is.matrix(result), 
                info = paste("v", version, "should return matrix"))
    expect_true(isSymmetric(result), 
                info = paste("v", version, "should be symmetric"))
    expect_true(all(diag(result) == 0), 
                info = paste("v", version, "should have zero diagonal"))
  }
})

test_that("extract_category_thresholds works with pre-0.1.6 formats", {
  versions <- get_available_legacy_versions()
  # Both pre-0.1.4 and 0.1.4-0.1.5 have $thresholds (deprecated but functional)
  pre_016_versions <- versions[vapply(versions, function(v) categorize_version(v) != "0.1.6+", logical(1))]
  skip_if(length(pre_016_versions) == 0, "No pre-0.1.6 fixtures available")
  
  for (version in pre_016_versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    result <- suppressWarnings(extract_category_thresholds(fit))
    
    expect_true(is.matrix(result), 
                info = paste("v", version, "should return matrix"))
    expect_true(nrow(result) > 0, 
                info = paste("v", version, "should have rows"))
  }
})

test_that("extract_arguments works with all legacy versions", {
  versions <- get_available_legacy_versions()
  skip_if(length(versions) == 0, "No legacy fixtures available")
  
  for (version in versions) {
    fit <- get_legacy_fit(paste0("fit_v", version))
    args <- extract_arguments(fit)
    
    expect_type(args, "list")
    expect_true(
      "no_variables" %in% names(args) || "num_variables" %in% names(args),
      info = paste("v", version, "should have variable count")
    )
    expect_true("data_columnnames" %in% names(args),
                info = paste("v", version, "should have column names"))
  }
})
