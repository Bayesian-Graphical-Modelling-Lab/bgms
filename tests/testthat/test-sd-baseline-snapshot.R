# --------------------------------------------------------------------------- #
# Pre-refactor snapshots of the L-space SD chain at alpha = 1 and alpha = 3.
#
# Purpose: lock the chain output before the AGHQ refactor (Phase 2 of the
# plan at dev/plans/active/sd-aghq-cubic.md) so the migration has a concrete
# baseline to assert against.
#
# Behaviour:
#   - If the fixture file does not exist, the test regenerates and saves it,
#     and the test passes with a note. Use this on the parent commit to
#     produce the baseline.
#   - If the fixture exists, the chain is re-run with the same seeds /
#     inputs and the output is compared to the saved fixture.
#
# Bit-equality expectation:
#   alpha = 1: the new AGHQ primitive collapses to the closed-form Gaussian
#     (the cubic factors and the AGHQ at the closed-form mode with N = 32
#     nodes is a single well-placed sweep). Chain output must match the
#     fixture to within roundoff (tolerance 1e-12).
#   alpha = 3: the AGHQ refactor changes the BF and the proposal Gaussian,
#     so chain output is expected to drift. The test reports the drift but
#     does not gate on it; the post-refactor SBC at alpha > 1 is the
#     end-to-end calibration check.
# --------------------------------------------------------------------------- #

# Run a deterministic SD chain at fixed seed + data. Returns the elements we
# compare against the fixture.
run_sd_baseline_chain <- function(seed,
                                  alpha,
                                  q = 6,
                                  n_obs = 80,
                                  n_warmup = 100,
                                  n_sweeps = 200,
                                  sample_thin = 20) {
  # Fixed-seed observations so the chain is reproducible.
  set.seed(seed)
  Y <- matrix(rnorm(n_obs * q), nrow = n_obs, ncol = q)
  # Centre to match bgm() default (per memory project_ggm_n_minus_1_bug).
  Y <- scale(Y, center = TRUE, scale = FALSE)
  attr(Y, "scaled:center") <- NULL

  out <- bgms:::ggm_sd_smoke_cpp(
    observations       = Y,
    inclusion_prob     = 0.5,
    interaction_scale  = 0.5,           # sigma
    diagonal_shape     = alpha,
    diagonal_rate      = 2.0,
    delta              = 0.0,            # no determinant tilt at the baseline
    n_warmup           = n_warmup,
    n_sweeps           = n_sweeps,
    seed               = seed,
    prior_only         = FALSE,
    include_within_k   = TRUE,
    sample_thin        = sample_thin,
    edge_selection     = TRUE,
    within_k_mode      = "am"            # Roverato + diag RW within K
  )
  list(
    pip               = out$pip,
    final_edges       = out$final_edges,
    final_K           = out$final_K,
    n_edges_path      = out$n_edges_path,
    K_offdiag_samples = out$K_offdiag_samples,
    K_diag_samples    = out$K_diag_samples,
    n_pd_reverts      = out$n_pd_reverts
  )
}

compare_to_fixture <- function(got, ref, tol, info_prefix) {
  expect_equal(got$pip,              ref$pip,              tolerance = tol,
               info = paste(info_prefix, "pip"))
  expect_equal(unname(got$final_edges), unname(ref$final_edges),
               info = paste(info_prefix, "final_edges"))
  expect_equal(got$final_K,          ref$final_K,          tolerance = tol,
               info = paste(info_prefix, "final_K"))
  expect_equal(got$n_edges_path,     ref$n_edges_path,
               info = paste(info_prefix, "n_edges_path"))
  expect_equal(got$K_offdiag_samples, ref$K_offdiag_samples, tolerance = tol,
               info = paste(info_prefix, "K_offdiag_samples"))
  expect_equal(got$K_diag_samples,   ref$K_diag_samples,   tolerance = tol,
               info = paste(info_prefix, "K_diag_samples"))
  expect_equal(got$n_pd_reverts,     ref$n_pd_reverts,
               info = paste(info_prefix, "n_pd_reverts"))
}


# --------------------------------------------------------------------------- #
# alpha = 1 baseline. Bit-equality target post-refactor.
# --------------------------------------------------------------------------- #

test_that("SD chain at alpha = 1 matches pre-refactor baseline (bit-equality)", {
  fixture_path <- testthat::test_path("fixtures", "sd_baseline_alpha1.rds")
  seeds <- c(11, 23, 47)
  got <- lapply(seeds, run_sd_baseline_chain, alpha = 1)
  names(got) <- paste0("seed_", seeds)
  if (!file.exists(fixture_path)) {
    saveRDS(got, fixture_path)
    skip(sprintf("Created baseline fixture at %s; rerun to verify.",
                 fixture_path))
  }
  ref <- readRDS(fixture_path)
  for (s in seq_along(seeds)) {
    compare_to_fixture(
      got[[s]], ref[[s]], tol = 1e-12,
      info_prefix = sprintf("alpha=1 seed=%d", seeds[s])
    )
  }
})


# --------------------------------------------------------------------------- #
# alpha = 3 baseline. Diagnostic only post-refactor (the AGHQ rewrite is
# expected to shift the chain). The test compares but reports drift as a
# message rather than a hard failure when run under the new code.
# --------------------------------------------------------------------------- #

test_that("SD chain at alpha = 3 baseline exists and captures current behaviour", {
  fixture_path <- testthat::test_path("fixtures", "sd_baseline_alpha3.rds")
  seeds <- c(11, 23, 47)
  got <- lapply(seeds, run_sd_baseline_chain, alpha = 3)
  names(got) <- paste0("seed_", seeds)
  if (!file.exists(fixture_path)) {
    saveRDS(got, fixture_path)
    skip(sprintf("Created baseline fixture at %s; rerun to verify.",
                 fixture_path))
  }
  # Under pre-refactor code we require an exact match to the fixture; under
  # post-refactor code the AGHQ path replaces GH-fixed and the chain drifts.
  # Both cases use the same expectation: bit-equality. If a future change
  # shifts the chain intentionally, regenerate the fixture and document why.
  ref <- readRDS(fixture_path)
  for (s in seq_along(seeds)) {
    compare_to_fixture(
      got[[s]], ref[[s]], tol = 1e-12,
      info_prefix = sprintf("alpha=3 seed=%d", seeds[s])
    )
  }
})
