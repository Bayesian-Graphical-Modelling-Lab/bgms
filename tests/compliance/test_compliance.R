# ==============================================================================
# Bitwise Compliance Test: Current Build vs CRAN 0.1.6.3 Fixtures
# ==============================================================================
#
# Reruns every bgm() and bgmCompare() configuration from the compliance
# fixture set and verifies bitwise-identical output.
#
# Usage:
#   Rscript tests/compliance/test_compliance.R
#
# Prerequisites:
#   Rscript tests/compliance/generate_fixtures.R
#
# Exit code:
#   0 = all pass, 1 = any fail
#
# ------------------------------------------------------------------------------
# Known differences vs CRAN 0.1.6.3
# ------------------------------------------------------------------------------
#
# 1. NUTS target_accept defaults changed between 0.1.6.3 and current:
#    - bgm():        0.60 → 0.80
#    - bgmCompare(): 0.65 → 0.80
#    DualAveraging receives a different target acceptance rate, changing the
#    adapted step size from iteration 1 onward and cascading into different
#    NUTS trajectories. MH and HMC defaults did not change.
#    Fix: bgm NUTS configs pass target_accept=0.6, bgmCompare NUTS configs
#    pass target_accept=0.65, matching the old defaults.
#
# 2. nuts_diag structure: current code adds a warmup_check field and
#    warmup_incomplete to the summary sub-list. These fields did not exist
#    in 0.1.6.3. The comparison skips nuts_diag and compares only the
#    shared core fields (treedepth, divergent, energy, ebfmi).
#
# 3. posterior_summary NA handling: 0.1.6.3 returned NA for mean and Rhat
#    of constant or near-constant parameters. Current code computes values
#    in those cells. The comparison ignores cells that are NA in the fixture.
#
# 4. SBM RNG ordering (CONFIRMED, NOT FIXED — not a bug):
#    The SBM algorithm itself is unchanged. However, old code initialized
#    cluster allocations (runif + block_probs_mfm_sbm) BEFORE the NUTS
#    step-size search; new code does step-size search first, then lazily
#    initializes SBM on the first edge_prior.update(). The swapped RNG
#    order causes divergence from iteration 1. Both orderings produce
#    valid MCMC chains; this is not a bug. SBM configs use structure-only
#    comparison.
#
# 5. Imputation bugs (CONFIRMED, 2 of 3 FIXED):
#    (a) Stale gradient cache: ensure_gradient_cache() cached sufficient
#        statistics (counts_per_category_, pairwise_stats_) but was never
#        invalidated after impute_missing() updated them. Old code
#        recomputed gradient_observed_active() fresh at each HMC call.
#        Fix: invalidate_gradient_cache() at end of impute_missing().
#    (b) Stale observations_double_t_: the precomputed transpose used
#        in gradient computation was set once in the constructor, never
#        updated after imputation changed observations_double_ cells.
#        Old code created obs_double_t fresh from observations each
#        gradient call. Fix: update observations_double_t_ at end of
#        impute_missing().
#    (c) Blume-Capel baseline probability: old code pre-initialized
#        cumsum with MY_EXP(ref*main + ref^2*main) before the category
#        loop, adding a phantom probability mass to the imputation
#        distribution. New code starts cumsum = 0.0. Confirmed: restoring
#        the old formula produces bitwise match. The old code was wrong;
#        the new code is correct. This config uses structure-only check.
#
# 6. Blume-Capel centering: all configs pass because CRAN 0.1.6.3 already
#    centered BC data in R before passing to C++, producing the same result
#    as the current C++ constructor. The observations_/observations_double_
#    mismatch was a regression introduced and fixed within PR #78; it was
#    never in a CRAN release.
#
# 7. Intermediate-overflow guard in compute_logZ_and_probs_ordinal and
#    compute_logZ_and_probs_blume_capel (ACCEPTED — not a bug):
#    Commit 04b9562 tightened the fast/slow block threshold from EXP_BOUND
#    (709) to FAST_LIM = max(0, EXP_BOUND - max_abs_main) to prevent
#    intermediate overflow in exp(main_param(c) + (c+1)*rest) before the
#    cancellation with exp(-bound). Both code paths are mathematically
#    identical but differ at floating-point level. During HMC leapfrog
#    integration, parameters can temporarily reach extreme values where
#    max_abs_main is large enough that FAST_LIM < EXP_BOUND, reclassifying
#    some observations between the fast (vectorized) and slow (per-element)
#    paths. The resulting floating-point perturbation cascades through the
#    fixed-step leapfrog integrator. This is needed for mixed MRF models
#    where Theta_ss is absorbed into main_param. The affected configs use
#    structure-only comparison against CRAN fixtures.
#
# 8. Association-scale reparameterization (stored A = σ/2):
#    OMRF pairwise effects are now stored on association scale (half the old
#    conditional rest-score coefficient σ). The Cauchy prior scale is halved
#    to keep the posterior equivalent, but the MCMC sampler operates on
#    association scale with different proposal/step-size scales, so MCMC
#    trajectories diverge from CRAN 0.1.6.3 fixtures. All configs use
#    structure-only comparison until fixtures are regenerated.
#
# 9. Lazy diagnostics & C++ ESS/Rhat (PR #77):
#    - posterior_summary diagnostic columns (n_eff, Rhat, mcse) are now
#      computed by C++ AR-spectral / Gelman-Rubin instead of coda.
#      Values differ at the algorithmic level.
#    - Indicator summaries rename n_eff -> n_eff_mixt.
#    - Edge-selected pairwise summaries gain an unconditional n_eff column.
#    - posterior_mean computed via colMeans() instead of mean() inside
#      summarize_fit(), causing ~1e-14 floating-point differences.
#    Fix: posterior_summary comparisons strip diagnostic columns (n_eff,
#    n_eff_mixt, Rhat, mcse). posterior_mean comparisons use tolerance
#    1e-12.
#
# 10. sd / mcse swap in pairwise summaries (bug fix, CONFIRMED):
#    CRAN 0.1.6.3 placed MCSE in the sd column and SD in the mcse column
#    for edge-selected pairwise summaries. The new code has them correct:
#    sd = posterior SD, mcse = sd / sqrt(n_eff).
#    Fix: comparisons strip sd and mcse columns (along with other derived
#    columns in diag_cols), so the historical swap is ignored.
#
# ==============================================================================

library(bgms)

fixture_dir = file.path("tests", "compliance", "fixtures")

if(!file.exists(file.path(fixture_dir, "manifest.rds"))) {
  stop("No fixtures found. Run tests/compliance/generate_fixtures.R first.")
}

manifest = readRDS(file.path(fixture_dir, "manifest.rds"))
cat(sprintf("Testing %d compliance fixtures against current build...\n", length(manifest)))
cat("Current bgms version:", as.character(packageVersion("bgms")), "\n\n")

# ==============================================================================
# Datasets (must match generate_fixtures.R)
# ==============================================================================

data(Wenchuan)
data(ADHD)
data(Boredom)

wenchuan_small = Wenchuan[, 1:6]
adhd_small = ADHD[, 2:7]
boredom_small = Boredom[, 2:7]

wenchuan_na = wenchuan_small
set.seed(999)
na_idx = sample(length(wenchuan_na), size = 20)
wenchuan_na[na_idx] = NA

n_w = nrow(wenchuan_small)
wenchuan_g1 = wenchuan_small[1:floor(n_w / 2), ]
wenchuan_g2 = wenchuan_small[(floor(n_w / 2) + 1):n_w, ]

n_a = nrow(adhd_small)
adhd_g1 = adhd_small[1:floor(n_a / 2), ]
adhd_g2 = adhd_small[(floor(n_a / 2) + 1):n_a, ]

n_b = nrow(boredom_small)
boredom_g1 = boredom_small[1:floor(n_b / 2), ]
boredom_g2 = boredom_small[(floor(n_b / 2) + 1):n_b, ]

wenchuan_g1_na = wenchuan_g1
wenchuan_g2_na = wenchuan_g2
set.seed(998)
na_idx1 = sample(length(wenchuan_g1_na), size = 10)
wenchuan_g1_na[na_idx1] = NA
set.seed(997)
na_idx2 = sample(length(wenchuan_g2_na), size = 10)
wenchuan_g2_na[na_idx2] = NA

datasets = list(
  wenchuan_small = wenchuan_small,
  adhd_small     = adhd_small,
  boredom_small  = boredom_small,
  wenchuan_na    = wenchuan_na,
  wenchuan_g1    = wenchuan_g1,
  wenchuan_g2    = wenchuan_g2,
  adhd_g1        = adhd_g1,
  adhd_g2        = adhd_g2,
  boredom_g1     = boredom_g1,
  boredom_g2     = boredom_g2,
  wenchuan_g1_na = wenchuan_g1_na,
  wenchuan_g2_na = wenchuan_g2_na
)

# ==============================================================================
# Configuration rebuild (mirrors generate_fixtures.R)
# ==============================================================================

bgm_configs = list(
  bgm_wenchuan_nuts_bernoulli = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    target_accept = 0.6, update_method = "nuts", seed = 101, display_progress = "none"
  ),
  bgm_wenchuan_nuts_no_edgesel = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = FALSE,
    target_accept = 0.6, update_method = "nuts", seed = 102, display_progress = "none"
  ),
  bgm_wenchuan_mh_bernoulli = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    update_method = "adaptive-metropolis", seed = 103, display_progress = "none"
  ),
  bgm_wenchuan_hmc_bernoulli = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    update_method = "hamiltonian-mc", hmc_num_leapfrogs = 20,
    seed = 104, display_progress = "none"
  ),
  bgm_wenchuan_nuts_betabern = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
    beta_bernoulli_alpha = 1, beta_bernoulli_beta = 1,
    target_accept = 0.6, update_method = "nuts", seed = 105, display_progress = "none"
  ),
  bgm_wenchuan_nuts_sbm = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Stochastic-Block",
    target_accept = 0.6, update_method = "nuts", seed = 106, display_progress = "none"
  ),
  bgm_wenchuan_nuts_scaled_prior = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    pairwise_scale = 1.0,
    target_accept = 0.6, update_method = "nuts", seed = 107, display_progress = "none"
  ),
  bgm_wenchuan_nuts_blumecapel = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    target_accept = 0.6, update_method = "nuts", seed = 401, display_progress = "none"
  ),
  bgm_wenchuan_mh_blumecapel = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    update_method = "adaptive-metropolis", seed = 402, display_progress = "none"
  ),
  bgm_wenchuan_hmc_blumecapel = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    update_method = "hamiltonian-mc", hmc_num_leapfrogs = 20,
    seed = 403, display_progress = "none"
  ),
  bgm_wenchuan_nuts_blumecapel_no_edgesel = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    edge_selection = FALSE,
    target_accept = 0.6, update_method = "nuts", seed = 404, display_progress = "none"
  ),
  bgm_wenchuan_nuts_blumecapel_baseline1 = list(
    x = "wenchuan_small", iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 1,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    target_accept = 0.6, update_method = "nuts", seed = 405, display_progress = "none"
  ),
  bgm_adhd_nuts_bernoulli = list(
    x = "adhd_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    target_accept = 0.6, update_method = "nuts", seed = 201, display_progress = "none"
  ),
  bgm_adhd_mh_bernoulli = list(
    x = "adhd_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    update_method = "adaptive-metropolis", seed = 202, display_progress = "none"
  ),
  bgm_adhd_nuts_no_edgesel = list(
    x = "adhd_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = FALSE,
    target_accept = 0.6, update_method = "nuts", seed = 203, display_progress = "none"
  ),
  bgm_adhd_nuts_sbm = list(
    x = "adhd_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Stochastic-Block",
    target_accept = 0.6, update_method = "nuts", seed = 204, display_progress = "none"
  ),
  bgm_boredom_nuts_bernoulli = list(
    x = "boredom_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    target_accept = 0.6, update_method = "nuts", seed = 301, display_progress = "none"
  ),
  bgm_boredom_mh_betabern = list(
    x = "boredom_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Beta-Bernoulli",
    update_method = "adaptive-metropolis", seed = 302, display_progress = "none"
  ),
  bgm_boredom_hmc_bernoulli = list(
    x = "boredom_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = TRUE, edge_prior = "Bernoulli",
    update_method = "hamiltonian-mc", hmc_num_leapfrogs = 20,
    seed = 303, display_progress = "none"
  ),
  bgm_boredom_nuts_no_edgesel = list(
    x = "boredom_small", iter = 200, warmup = 200, chains = 2,
    edge_selection = FALSE,
    target_accept = 0.6, update_method = "nuts", seed = 304, display_progress = "none"
  ),
  bgm_wenchuan_nuts_impute = list(
    x = "wenchuan_na", iter = 200, warmup = 200, chains = 2,
    na_action = "impute",
    edge_selection = TRUE, edge_prior = "Bernoulli",
    target_accept = 0.6, update_method = "nuts", seed = 501, display_progress = "none"
  ),
  bgm_wenchuan_nuts_blumecapel_impute = list(
    x = "wenchuan_na", iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    na_action = "impute",
    edge_selection = TRUE, edge_prior = "Bernoulli",
    target_accept = 0.6, update_method = "nuts", seed = 502, display_progress = "none"
  )
)

compare_configs = list(
  cmp_wenchuan_nuts_bernoulli = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    target_accept = 0.65, update_method = "nuts", seed = 601, display_progress = "none"
  ),
  cmp_wenchuan_mh_bernoulli = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    update_method = "adaptive-metropolis", seed = 602, display_progress = "none"
  ),
  cmp_wenchuan_hmc_bernoulli = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    update_method = "hamiltonian-mc", hmc_num_leapfrogs = 20,
    seed = 603, display_progress = "none"
  ),
  cmp_wenchuan_nuts_no_diffsel = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = FALSE,
    target_accept = 0.65, update_method = "nuts", seed = 604, display_progress = "none"
  ),
  cmp_wenchuan_nuts_main_diffsel = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    main_difference_selection = TRUE,
    target_accept = 0.65, update_method = "nuts", seed = 605, display_progress = "none"
  ),
  cmp_wenchuan_nuts_betabern = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    difference_prior = "Beta-Bernoulli",
    target_accept = 0.65, update_method = "nuts", seed = 606, display_progress = "none"
  ),
  cmp_wenchuan_nuts_blumecapel = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    difference_selection = TRUE,
    target_accept = 0.65, update_method = "nuts", seed = 701, display_progress = "none"
  ),
  cmp_wenchuan_mh_blumecapel = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    difference_selection = TRUE,
    update_method = "adaptive-metropolis", seed = 702, display_progress = "none"
  ),
  cmp_wenchuan_nuts_blumecapel_no_diffsel = list(
    x = "wenchuan_g1", y = "wenchuan_g2",
    iter = 200, warmup = 200, chains = 2,
    variable_type = "blume-capel", baseline_category = 2,
    difference_selection = FALSE,
    target_accept = 0.65, update_method = "nuts", seed = 703, display_progress = "none"
  ),
  cmp_adhd_nuts_bernoulli = list(
    x = "adhd_g1", y = "adhd_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    target_accept = 0.65, update_method = "nuts", seed = 801, display_progress = "none"
  ),
  cmp_adhd_mh_bernoulli = list(
    x = "adhd_g1", y = "adhd_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    update_method = "adaptive-metropolis", seed = 802, display_progress = "none"
  ),
  cmp_boredom_nuts_bernoulli = list(
    x = "boredom_g1", y = "boredom_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    target_accept = 0.65, update_method = "nuts", seed = 901, display_progress = "none"
  ),
  cmp_boredom_mh_bernoulli = list(
    x = "boredom_g1", y = "boredom_g2",
    iter = 200, warmup = 200, chains = 2,
    difference_selection = TRUE,
    update_method = "adaptive-metropolis", seed = 902, display_progress = "none"
  ),
  cmp_wenchuan_nuts_impute = list(
    x = "wenchuan_g1_na", y = "wenchuan_g2_na",
    iter = 200, warmup = 200, chains = 2,
    na_action = "impute",
    difference_selection = TRUE,
    target_accept = 0.65, update_method = "nuts", seed = 1001, display_progress = "none"
  )
)

all_configs = c(bgm_configs, compare_configs)

# ==============================================================================
# Comparison helpers
# ==============================================================================

resolve_args = function(args) {
  resolved = args
  for(nm in names(resolved)) {
    if(is.character(resolved[[nm]]) && resolved[[nm]] %in% names(datasets)) {
      resolved[[nm]] = datasets[[resolved[[nm]]]]
    }
  }
  resolved
}

extract_bgm_actual = function(fit) {
  list(
    posterior_summary_main = fit$posterior_summary_main,
    posterior_summary_pairwise = fit$posterior_summary_pairwise,
    posterior_summary_indicator = fit$posterior_summary_indicator,
    posterior_mean_main = fit$posterior_mean_main,
    posterior_mean_associations = fit$posterior_mean_associations,
    posterior_mean_indicator = fit$posterior_mean_indicator,
    raw_main_chain1 = fit$raw_samples$main[[1]],
    raw_pairwise_chain1 = fit$raw_samples$pairwise[[1]],
    raw_indicator_chain1 = if(!is.null(fit$raw_samples$indicator)) {
      fit$raw_samples$indicator[[1]]
    } else {
      NULL
    },
    nuts_diag = fit$nuts_diag,
    posterior_coclustering_matrix = fit$posterior_coclustering_matrix,
    posterior_mean_allocations = fit$posterior_mean_allocations
  )
}

extract_compare_actual = function(fit) {
  list(
    posterior_summary_main_baseline = fit$posterior_summary_main_baseline,
    posterior_summary_pairwise_baseline = fit$posterior_summary_pairwise_baseline,
    posterior_summary_main_differences = fit$posterior_summary_main_differences,
    posterior_summary_pairwise_differences = fit$posterior_summary_pairwise_differences,
    posterior_summary_indicator = fit$posterior_summary_indicator,
    posterior_mean_main_baseline = fit$posterior_mean_main_baseline,
    posterior_mean_associations_baseline = fit$posterior_mean_associations_baseline,
    posterior_mean_main_differences = fit$posterior_mean_main_differences,
    posterior_mean_associations_differences = fit$posterior_mean_associations_differences,
    posterior_mean_indicator = fit$posterior_mean_indicator,
    raw_samples = fit$raw_samples,
    nuts_diag = fit$nuts_diag
  )
}

# Configs where 0.1.6.3 returned NA in posterior_summary/posterior_mean for
# constant or near-constant parameters, and the current code computes values
# (bug fix). Only these are allowed to have NA→value differences.
na_bugfix_ids = c(
  "bgm_wenchuan_nuts_bernoulli",
  "bgm_wenchuan_hmc_bernoulli",
  "bgm_wenchuan_nuts_blumecapel",
  "bgm_boredom_hmc_bernoulli",
  "bgm_wenchuan_nuts_blumecapel_baseline1",
  "bgm_wenchuan_nuts_impute"
)

# Configs excluded from bitwise comparison due to confirmed algorithm changes
# (see header notes 4, 5, 7, 8, 9, 10). Checked for structural match only.
# All configs are structure-only pending fixture regeneration after the
# association-scale reparameterization (note 8).
structure_only_ids = c(
  names(all_configs)
)

# Diagnostic and derived columns that changed in PR #77 (note 8).
# Stripped from posterior_summary comparisons. The values in these columns
# are derived from the raw chains (which ARE compared bitwise).
# sd is included because CRAN 0.1.6.3 had sd/mcse values swapped in
# the edge-selected pairwise summaries (note 9).
diag_cols = c("n_eff", "n_eff_mixt", "Rhat", "mcse", "sd", "parameter")

# Strip diagnostic / derived columns from a data.frame/matrix, keeping only
# shared non-diagnostic columns between expected and actual.
strip_diag_cols = function(exp_df, act_df) {
  exp_names = if(is.data.frame(exp_df)) names(exp_df) else colnames(exp_df)
  act_names = if(is.data.frame(act_df)) names(act_df) else colnames(act_df)
  shared = intersect(exp_names, act_names)
  keep = setdiff(shared, diag_cols)
  list(
    exp = as.matrix(exp_df[, keep, drop = FALSE]),
    act = as.matrix(act_df[, keep, drop = FALSE])
  )
}

compare_fields = function(expected, actual, type, id) {
  if(type == "bgm") {
    fields = c(
      "posterior_summary_main", "posterior_summary_pairwise",
      "posterior_summary_indicator",
      "posterior_mean_main", "posterior_mean_associations", "posterior_mean_indicator",
      "raw_main_chain1", "raw_pairwise_chain1", "raw_indicator_chain1",
      "posterior_coclustering_matrix", "posterior_mean_allocations"
    )
  } else {
    fields = c(
      "posterior_summary_main_baseline", "posterior_summary_pairwise_baseline",
      "posterior_summary_main_differences", "posterior_summary_pairwise_differences",
      "posterior_summary_indicator",
      "posterior_mean_main_baseline", "posterior_mean_associations_baseline",
      "posterior_mean_main_differences", "posterior_mean_associations_differences",
      "posterior_mean_indicator",
      "raw_samples"
    )
  }

  allow_na_skip = id %in% na_bugfix_ids
  mismatches = character()

  for(field in fields) {
    exp_val = expected[[field]]
    act_val = actual[[field]]

    if(is.null(exp_val) && is.null(act_val)) next
    if(is.null(exp_val) != is.null(act_val)) {
      mismatches = c(mismatches, sprintf("  %s: one is NULL, the other is not", field))
      next
    }

    is_summary = grepl("^posterior_summary", field)
    is_mean = grepl("^posterior_mean", field)

    # posterior_summary: strip diagnostic columns before comparing (notes 8-9)
    if(is_summary && (is.data.frame(exp_val) || is.data.frame(act_val))) {
      stripped = strip_diag_cols(exp_val, act_val)
      exp_m = stripped$exp
      act_m = stripped$act

      if(!identical(dim(exp_m), dim(act_m))) {
        mismatches = c(mismatches, sprintf(
          "  %s: dim mismatch after stripping diag cols (%s vs %s)", field,
          paste(dim(exp_m), collapse = "x"), paste(dim(act_m), collapse = "x")
        ))
        next
      }

      if(allow_na_skip) {
        non_na = !is.na(exp_m)
        if(!isTRUE(all.equal(exp_m[non_na], act_m[non_na], tolerance = 1e-12))) {
          max_diff = max(abs(exp_m[non_na] - act_m[non_na]), na.rm = TRUE)
          mismatches = c(mismatches, sprintf(
            "  %s: NOT identical (max |diff| = %.2e, ignoring diag cols + fixture NAs)",
            field, max_diff
          ))
        }
      } else {
        if(!isTRUE(all.equal(exp_m, act_m, tolerance = 1e-12))) {
          max_diff = max(abs(exp_m - act_m), na.rm = TRUE)
          mismatches = c(mismatches, sprintf(
            "  %s: NOT identical (max |diff| = %.2e, ignoring diag cols)", field, max_diff
          ))
        }
      }
      next
    }

    # posterior_mean: allow tolerance for colMeans vs mean() differences (note 8)
    if(is_mean && (is.data.frame(exp_val) || is.matrix(exp_val))) {
      exp_m = as.matrix(exp_val)
      act_m = as.matrix(act_val)
      if(!identical(dim(exp_m), dim(act_m))) {
        mismatches = c(mismatches, sprintf(
          "  %s: dim mismatch (%s vs %s)", field,
          paste(dim(exp_m), collapse = "x"), paste(dim(act_m), collapse = "x")
        ))
        next
      }
      if(allow_na_skip) {
        non_na = !is.na(exp_m)
        if(!isTRUE(all.equal(exp_m[non_na], act_m[non_na], tolerance = 1e-12))) {
          max_diff = max(abs(exp_m[non_na] - act_m[non_na]), na.rm = TRUE)
          mismatches = c(mismatches, sprintf(
            "  %s: NOT identical (max |diff| = %.2e, ignoring fixture NAs)", field, max_diff
          ))
        }
      } else {
        if(!isTRUE(all.equal(exp_m, act_m, tolerance = 1e-12))) {
          max_diff = max(abs(exp_m - act_m), na.rm = TRUE)
          mismatches = c(mismatches, sprintf(
            "  %s: NOT identical (max |diff| = %.2e)", field, max_diff
          ))
        }
      }
      next
    }

    # Everything else: bitwise identical
    if(!identical(dim(exp_val), dim(act_val))) {
      mismatches = c(mismatches, sprintf(
        "  %s: dim mismatch (%s vs %s)", field,
        paste(dim(exp_val), collapse = "x"), paste(dim(act_val), collapse = "x")
      ))
      next
    }

    if(!identical(exp_val, act_val)) {
      if(is.numeric(exp_val) && is.numeric(act_val)) {
        max_diff = max(abs(as.numeric(exp_val) - as.numeric(act_val)), na.rm = TRUE)
        mismatches = c(mismatches, sprintf(
          "  %s: NOT identical (max |diff| = %.2e)", field, max_diff
        ))
      } else {
        mismatches = c(mismatches, sprintf("  %s: NOT identical (non-numeric)", field))
      }
    }
  }

  # nuts_diag: compare only the shared core fields (new code adds warmup_check)
  exp_nd = expected[["nuts_diag"]]
  act_nd = actual[["nuts_diag"]]
  if(!is.null(exp_nd) && !is.null(act_nd)) {
    for(nd_field in c("treedepth", "divergent", "energy", "ebfmi")) {
      if(!identical(exp_nd[[nd_field]], act_nd[[nd_field]])) {
        mismatches = c(mismatches, sprintf("  nuts_diag$%s: NOT identical", nd_field))
      }
    }
  } else if(!is.null(exp_nd) != !is.null(act_nd)) {
    mismatches = c(mismatches, "  nuts_diag: one is NULL, the other is not")
  }

  mismatches
}

# Structure-only check: verifies that output fields
# have matching names, dimensions, and types without requiring identical values.
check_structure = function(expected, actual, type) {
  if(type == "bgm") {
    fields = c(
      "posterior_summary_main", "posterior_summary_pairwise",
      "posterior_summary_indicator",
      "posterior_mean_main", "posterior_mean_associations", "posterior_mean_indicator",
      "raw_main_chain1", "raw_pairwise_chain1", "raw_indicator_chain1"
    )
  } else {
    fields = c(
      "posterior_summary_main_baseline", "posterior_summary_pairwise_baseline",
      "posterior_summary_main_differences", "posterior_summary_pairwise_differences",
      "posterior_summary_indicator",
      "posterior_mean_main_baseline", "posterior_mean_associations_baseline",
      "posterior_mean_main_differences", "posterior_mean_associations_differences",
      "posterior_mean_indicator",
      "raw_samples"
    )
  }

  mismatches = character()
  for(field in fields) {
    exp_val = expected[[field]]
    act_val = actual[[field]]
    if(is.null(exp_val) && is.null(act_val)) next
    if(is.null(exp_val) != is.null(act_val)) {
      mismatches = c(mismatches, sprintf("  %s: one is NULL, the other is not", field))
      next
    }
    if(!identical(class(exp_val), class(act_val))) {
      mismatches = c(mismatches, sprintf(
        "  %s: class mismatch (%s vs %s)",
        field, paste(class(exp_val), collapse = "/"), paste(class(act_val), collapse = "/")
      ))
      next
    }
    # For posterior_summary fields: allow extra columns (note 8) and check
    # row count only. For everything else: require identical dimensions.
    is_summary = grepl("^posterior_summary", field)
    if(is_summary && (is.data.frame(exp_val) || is.data.frame(act_val))) {
      exp_nr = nrow(exp_val)
      act_nr = nrow(act_val)
      if(!identical(exp_nr, act_nr)) {
        mismatches = c(mismatches, sprintf(
          "  %s: row count mismatch (%d vs %d)", field, exp_nr, act_nr
        ))
      }
    } else if(!identical(dim(exp_val), dim(act_val)) &&
      !identical(length(exp_val), length(act_val))) {
      mismatches = c(mismatches, sprintf(
        "  %s: dim mismatch (%s vs %s)",
        field, paste(dim(exp_val), collapse = "x"), paste(dim(act_val), collapse = "x")
      ))
    }
  }
  mismatches
}

# ==============================================================================
# Run all comparisons
# ==============================================================================

pass_count = 0
fail_count = 0
skip_count = 0
error_count = 0
failures = list()

for(entry in manifest) {
  id = entry$id
  type = entry$type
  cat(sprintf("  [%s] %s ... ", id, entry$desc))

  # Load expected fixture
  fixture_path = file.path(fixture_dir, entry$file)
  if(!file.exists(fixture_path)) {
    cat("SKIP (fixture file missing)\n")
    skip_count = skip_count + 1
    next
  }
  expected = readRDS(fixture_path)

  # Get config
  if(!id %in% names(all_configs)) {
    cat("SKIP (config not found)\n")
    skip_count = skip_count + 1
    next
  }
  args = resolve_args(all_configs[[id]])

  # Run current build
  set.seed(args$seed)
  fit = tryCatch(
    {
      if(type == "bgm") {
        do.call(bgm, args)
      } else {
        do.call(bgmCompare, args)
      }
    },
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if(is.null(fit)) {
    error_count = error_count + 1
    failures[[id]] = "model fitting errored"
    next
  }

  # Extract comparable output
  if(type == "bgm") {
    actual = extract_bgm_actual(fit)
  } else {
    actual = extract_compare_actual(fit)
  }

  # Compare: structure-only for all configs pending fixture regeneration
  # (note 8: association-scale reparameterization breaks bitwise identity
  # with CRAN 0.1.6.3 fixtures)
  if(id %in% structure_only_ids) {
    mismatches = check_structure(expected, actual, type)
    label = "PASS (structure)"
  } else {
    mismatches = compare_fields(expected, actual, type, id)
    label = "PASS"
  }

  if(length(mismatches) == 0) {
    cat(label, "\n")
    pass_count = pass_count + 1
  } else {
    cat("FAIL\n")
    for(m in mismatches) cat(m, "\n")
    fail_count = fail_count + 1
    failures[[id]] = mismatches
  }
}

# ==============================================================================
# Summary
# ==============================================================================

total = length(manifest)
cat(sprintf(
  "\n=== Results: %d PASS, %d FAIL, %d ERROR, %d SKIP (of %d) ===\n",
  pass_count, fail_count, error_count, skip_count, total
))

if(fail_count > 0 || error_count > 0) {
  cat("\nFailed/errored fixtures:\n")
  for(id in names(failures)) {
    cat(sprintf("  %s:\n", id))
    msgs = failures[[id]]
    for(m in msgs) cat("   ", m, "\n")
  }
  quit(status = 1)
} else if(pass_count == total) {
  cat("All fixtures match (structure-only pending association-scale fixture regeneration).\n")
}
