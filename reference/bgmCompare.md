# Bayesian Estimation and Variable Selection for Group Differences in Markov Random Fields

The `bgmCompare` function estimates group differences in category
threshold parameters (main effects) and pairwise interactions (pairwise
effects) of a Markov Random Field (MRF) for binary and ordinal
variables. Groups can be defined either by supplying two separate
datasets (`x` and `y`) or by a group membership vector. Optionally,
Bayesian variable selection can be applied to identify differences
across groups.

## Usage

``` r
bgmCompare(
  x,
  y,
  group_indicator,
  difference_selection = TRUE,
  main_difference_selection = FALSE,
  variable_type = "ordinal",
  baseline_category,
  difference_scale = 1,
  difference_prior = c("Bernoulli", "Beta-Bernoulli"),
  difference_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  pairwise_scale = 2.5,
  main_alpha = 0.5,
  main_beta = 0.5,
  iter = 1000,
  warmup = 1000,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis", "hamiltonian-mc"),
  target_accept,
  hmc_num_leapfrogs = 100,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  standardize = FALSE,
  verbose = getOption("bgms.verbose", TRUE),
  main_difference_model,
  reference_category,
  main_difference_scale,
  pairwise_difference_scale,
  pairwise_difference_prior,
  main_difference_prior,
  pairwise_difference_probability,
  main_difference_probability,
  pairwise_beta_bernoulli_alpha,
  pairwise_beta_bernoulli_beta,
  main_beta_bernoulli_alpha,
  main_beta_bernoulli_beta,
  interaction_scale,
  threshold_alpha,
  threshold_beta,
  burnin,
  save
)
```

## Arguments

- x:

  A data frame or matrix of binary and ordinal responses for Group 1.
  Variables should be coded as nonnegative integers starting at 0. For
  ordinal variables, unused categories are collapsed; for Blume–Capel
  variables, all categories are retained.

- y:

  Optional data frame or matrix for Group 2 (two-group designs). Must
  have the same variables (columns) as `x`.

- group_indicator:

  Optional integer vector of group memberships for rows of `x`
  (multi-group designs). Ignored if `y` is supplied.

- difference_selection:

  Logical. If `TRUE`, spike-and-slab priors are applied to difference
  parameters. Default: `TRUE`.

- main_difference_selection:

  Logical. If `TRUE`, apply spike-and-slab selection to main effect
  (threshold) differences. If `FALSE`, main effect differences are
  always included (no selection). Since main effects are often nuisance
  parameters and their selection can interfere with pairwise selection
  under the Beta-Bernoulli prior, the default is `FALSE`. Only used when
  `difference_selection = TRUE`.

- variable_type:

  Character vector specifying type of each variable: `"ordinal"`
  (default) or `"blume-capel"`.

- baseline_category:

  Integer or vector giving the baseline category for Blume–Capel
  variables.

- difference_scale:

  Double. Scale of the Cauchy prior for difference parameters. Default:
  `1`.

- difference_prior:

  Character. Prior for difference inclusion: `"Bernoulli"` or
  `"Beta-Bernoulli"`. Default: `"Bernoulli"`.

- difference_probability:

  Numeric. Prior inclusion probability for differences (Bernoulli
  prior). Default: `0.5`.

- beta_bernoulli_alpha, beta_bernoulli_beta:

  Doubles. Shape parameters of the Beta prior for inclusion
  probabilities in the Beta–Bernoulli model. Defaults: `1`.

- pairwise_scale:

  Double. Scale of the Cauchy prior for baseline pairwise interactions.
  Default: `2.5`.

- main_alpha, main_beta:

  Doubles. Shape parameters of the beta-prime prior for baseline
  threshold parameters. Defaults: `0.5`.

- iter:

  Integer. Number of post–warmup iterations per chain. Default: `1e3`.

- warmup:

  Integer. Number of warmup iterations before sampling. Default: `1e3`.

- na_action:

  Character. How to handle missing data: `"listwise"` (drop rows) or
  `"impute"` (impute within Gibbs). Default: `"listwise"`.

- update_method:

  Character. Sampling algorithm: `"adaptive-metropolis"`,
  `"hamiltonian-mc"`, or `"nuts"`. Default: `"nuts"`.

- target_accept:

  Numeric between 0 and 1. Target acceptance rate. Defaults: 0.44
  (Metropolis), 0.65 (HMC), 0.80 (NUTS).

- hmc_num_leapfrogs:

  Integer. Leapfrog steps for HMC. Default: `100`.

- nuts_max_depth:

  Integer. Maximum tree depth for NUTS. Default: `10`.

- learn_mass_matrix:

  Logical. If `TRUE`, adapts a diagonal mass matrix during warmup
  (HMC/NUTS only). Default: `TRUE`.

- chains:

  Integer. Number of parallel chains. Default: `4`.

- cores:

  Integer. Number of CPU cores. Default:
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).

- display_progress:

  Character. Controls progress reporting: `"per-chain"`, `"total"`, or
  `"none"`. Default: `"per-chain"`.

- seed:

  Optional integer. Random seed for reproducibility.

- standardize:

  Logical. If `TRUE`, the Cauchy prior scale for each pairwise
  interaction (both baseline and difference) is adjusted based on the
  range of response scores. Without standardization, pairs with more
  response categories experience less shrinkage because their naturally
  smaller interaction effects make a fixed prior relatively wide.
  Standardization equalizes relative shrinkage across all pairs, with
  `pairwise_scale` itself applying to the unit interval (binary) case.
  See
  [`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  for details on the adjustment. Default: `FALSE`.

- verbose:

  Logical. If `TRUE`, prints informational messages during data
  processing (e.g., missing data handling, variable recoding). Defaults
  to `getOption("bgms.verbose", TRUE)`. Set
  `options(bgms.verbose = FALSE)` to suppress messages globally.

- main_difference_model, reference_category, pairwise_difference_scale,
  main_difference_scale, pairwise_difference_prior,
  main_difference_prior, pairwise_difference_probability,
  main_difference_probability, pairwise_beta_bernoulli_alpha,
  pairwise_beta_bernoulli_beta, main_beta_bernoulli_alpha,
  main_beta_bernoulli_beta, interaction_scale, threshold_alpha,
  threshold_beta, burnin, save:

  \`r lifecycle::badge("deprecated")\` Deprecated arguments as of
  \*\*bgms 0.1.6.0\*\*. Use \`difference_scale\`, \`difference_prior\`,
  \`difference_probability\`, \`beta_bernoulli_alpha\`,
  \`beta_bernoulli_beta\`, \`baseline_category\`, \`pairwise_scale\`,
  and \`warmup\` instead.

## Value

A list of class `"bgmCompare"` containing posterior summaries, posterior
mean matrices, and raw MCMC samples:

- `posterior_summary_main_baseline`,
  `posterior_summary_pairwise_baseline`: summaries of baseline
  thresholds and pairwise interactions.

- `posterior_summary_main_differences`,
  `posterior_summary_pairwise_differences`: summaries of group
  differences in thresholds and pairwise interactions.

- `posterior_summary_indicator`: summaries of inclusion indicators (if
  `difference_selection = TRUE`).

- `posterior_mean_main_baseline`, `posterior_mean_pairwise_baseline`:
  posterior mean matrices (legacy style).

- `raw_samples`: list of raw draws per chain for main, pairwise, and
  indicator parameters.

- `arguments`: list of function call arguments and metadata.

The [`summary()`](https://rdrr.io/r/base/summary.html) method prints
formatted summaries, and [`coef()`](https://rdrr.io/r/stats/coef.html)
extracts posterior means.

NUTS diagnostics (tree depth, divergences, energy, E-BFMI) are included
in `fit$nuts_diag` if `update_method = "nuts"`.

## Details

This function extends the ordinal MRF framework Marsman et al. (2025) to
multiple groups. The basic idea of modeling, analyzing, and testing
group differences in MRFs was introduced in Marsman et al. (2025) ,
where two–group comparisons were conducted using adaptive Metropolis
sampling. The present implementation generalizes that approach to more
than two groups and supports additional samplers (HMC and NUTS) with
staged warmup adaptation.

Key components of the model:

## Pairwise Interactions

For variables \\i\\ and \\j\\, the group-specific interaction is
represented as: \$\$\theta\_{ij}^{(g)} = \phi\_{ij} +
\delta\_{ij}^{(g)},\$\$ where \\\phi\_{ij}\\ is the baseline effect and
\\\delta\_{ij}^{(g)}\\ are group differences constrained to sum to zero.

## Ordinal Variables

**Regular ordinal variables**: category thresholds are decomposed into a
baseline plus group differences for each category.

**Blume–Capel variables**: category thresholds are quadratic in the
category index, with both the linear and quadratic terms split into a
baseline plus group differences.

## Variable Selection

When `difference_selection = TRUE`, spike-and-slab priors are applied to
difference parameters:

- **Bernoulli**: fixed prior inclusion probability.

- **Beta–Bernoulli**: inclusion probability given a Beta prior.

## Sampling Algorithms and Warmup

Parameters are updated within a Gibbs framework, using the same sampling
algorithms and staged warmup scheme described in
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md):

- **Adaptive Metropolis–Hastings**: componentwise random–walk proposals
  with Robbins–Monro adaptation of proposal SDs.

- **Hamiltonian Monte Carlo (HMC)**: joint updates with fixed leapfrog
  trajectories; step size and optionally the mass matrix are adapted
  during warmup.

- **No–U–Turn Sampler (NUTS)**: an adaptive HMC variant with dynamic
  trajectory lengths; warmup uses the same staged adaptation schedule as
  HMC.

For details on the staged adaptation schedule (fast–slow–fast phases),
see
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md).
In addition, when `difference_selection = TRUE`, updates of inclusion
indicators are delayed until late warmup. In HMC/NUTS, this appends two
extra phases (Stage-3b and Stage-3c), so that the total number of warmup
iterations exceeds the user-specified `warmup`.

After warmup, adaptation is disabled: step size and mass matrix are
fixed at their learned values, and proposal SDs remain constant.

## References

Marsman M, van den Bergh D, Haslbeck JMB (2025). “Bayesian analysis of
the ordinal Markov random field.” *Psychometrika*, **90**(1), 146–182.
[doi:10.1017/psy.2024.4](https://doi.org/10.1017/psy.2024.4) .  
  
Marsman M, Waldorp LJ, Sekulovski N, Haslbeck JMB (2025). “Bayes factor
tests for group differences in ordinal and binary graphical models.”
*Psychometrika*, **90**(5), 1809–1842.
[doi:10.1017/psy.2025.10060](https://doi.org/10.1017/psy.2025.10060) .

## See also

[`vignette("comparison", package = "bgms")`](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/comparison.md)
for a worked example.

## Examples

``` r
# \dontrun{
# Run bgmCompare on subset of the Boredom dataset
x = Boredom[Boredom$language == "fr", 2:6]
y = Boredom[Boredom$language != "fr", 2:6]

fit <- bgmCompare(x, y)
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2000 (2.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 47/2000 (2.4%)
#> Chain 3 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 53/2000 (2.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 48/2000 (2.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 198/8000 (2.5%)
#> Elapsed: 2s | ETA: 1m 18s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 106/2000 (5.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 194/2000 (9.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 126/2000 (6.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 526/8000 (6.6%)
#> Elapsed: 5s | ETA: 1m 11s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/2000 (10.0%)
#> Chain 2 (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 202/2000 (10.1%)
#> Chain 3 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 315/2000 (15.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 231/2000 (11.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 948/8000 (11.8%)
#> Elapsed: 6s | ETA: 45s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2000 (15.0%)
#> Chain 2 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 304/2000 (15.2%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 438/2000 (21.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 328/2000 (16.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1370/8000 (17.1%)
#> Elapsed: 7s | ETA: 34s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2000 (20.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 399/2000 (20.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 531/2000 (26.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 430/2000 (21.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1760/8000 (22.0%)
#> Elapsed: 7s | ETA: 25s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2000 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 499/2000 (24.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 641/2000 (32.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 522/2000 (26.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2162/8000 (27.0%)
#> Elapsed: 8s | ETA: 22s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 600/2000 (30.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 602/2000 (30.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 761/2000 (38.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 628/2000 (31.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2591/8000 (32.4%)
#> Elapsed: 8s | ETA: 17s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/2000 (35.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 708/2000 (35.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 857/2000 (42.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 742/2000 (37.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3007/8000 (37.6%)
#> Elapsed: 9s | ETA: 15s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 809/2000 (40.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 898/2000 (44.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 846/2000 (42.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3353/8000 (41.9%)
#> Elapsed: 10s | ETA: 14s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2000 (45.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 909/2000 (45.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 1001/2000 (50.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 924/2000 (46.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3734/8000 (46.7%)
#> Elapsed: 11s | ETA: 13s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 950/2000 (47.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 970/2000 (48.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1109/2000 (55.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 1004/2000 (50.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 4033/8000 (50.4%)
#> Elapsed: 11s | ETA: 11s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1050/2000 (52.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 1065/2000 (53.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1214/2000 (60.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1108/2000 (55.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 4437/8000 (55.5%)
#> Elapsed: 12s | ETA: 10s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2000 (57.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1161/2000 (58.1%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1312/2000 (65.6%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1208/2000 (60.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 4831/8000 (60.4%)
#> Elapsed: 13s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1250/2000 (62.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1263/2000 (63.1%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1416/2000 (70.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1310/2000 (65.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 5239/8000 (65.5%)
#> Elapsed: 13s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1350/2000 (67.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1363/2000 (68.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1512/2000 (75.6%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1408/2000 (70.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 5633/8000 (70.4%)
#> Elapsed: 14s | ETA: 6s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1450/2000 (72.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1464/2000 (73.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1607/2000 (80.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1505/2000 (75.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 6026/8000 (75.3%)
#> Elapsed: 15s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1550/2000 (77.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1562/2000 (78.1%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 1706/2000 (85.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1609/2000 (80.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 6427/8000 (80.3%)
#> Elapsed: 15s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2000 (82.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 1656/2000 (82.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 1806/2000 (90.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 1708/2000 (85.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 6820/8000 (85.2%)
#> Elapsed: 16s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/2000 (87.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1753/2000 (87.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 1907/2000 (95.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 1809/2000 (90.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 7219/8000 (90.2%)
#> Elapsed: 17s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1850/2000 (92.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 1858/2000 (92.9%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 1922/2000 (96.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 7630/8000 (95.4%)
#> Elapsed: 17s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1950/2000 (97.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 1960/2000 (98.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7910/8000 (98.9%)
#> Elapsed: 18s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 18s | ETA: 0s

# Posterior inclusion probabilities
summary(fit)$indicator
#>                            parameter    mean        sd        mcse
#> 1                  loose_ends (main) 1.00000 0.0000000          NA
#> 2    loose_ends-entertain (pairwise) 0.01650 0.1273882 0.002074373
#> 3   loose_ends-repetitive (pairwise) 0.04600 0.2094851 0.004625482
#> 4  loose_ends-stimulation (pairwise) 0.22700 0.4188926 0.018746283
#> 5    loose_ends-motivated (pairwise) 0.02725 0.1628110 0.003274486
#> 6                   entertain (main) 1.00000 0.0000000          NA
#> 7    entertain-repetitive (pairwise) 0.02425 0.1538244 0.002842817
#> 8   entertain-stimulation (pairwise) 0.19025 0.3924983 0.016798292
#> 9     entertain-motivated (pairwise) 0.04750 0.2127058 0.005123297
#> 10                 repetitive (main) 1.00000 0.0000000          NA
#> 11 repetitive-stimulation (pairwise) 0.02625 0.1598779 0.003121224
#> 12   repetitive-motivated (pairwise) 0.02300 0.1499033 0.002588036
#> 13                stimulation (main) 1.00000 0.0000000          NA
#> 14  stimulation-motivated (pairwise) 0.02175 0.1458662 0.003734278
#> 15                  motivated (main) 1.00000 0.0000000          NA
#>    n0->0 n0->1 n1->0 n1->1     n_eff     Rhat
#> 1      0     0     0  3999        NA       NA
#> 2   3870    63    63     3 3771.2438 1.003821
#> 3   3696   119   119    65 2051.1250 1.007849
#> 4   2936   155   156   752  499.3155 1.000207
#> 5   3809    81    81    28 2472.1876 1.006916
#> 6      0     0     0  3999        NA       NA
#> 7   3822    80    80    17 2927.8739 1.002743
#> 8   3090   148   148   613  545.9406 1.013966
#> 9   3700   109   109    81 1723.6916 1.011465
#> 10     0     0     0  3999        NA       NA
#> 11  3813    81    81    24 2623.7766 1.003714
#> 12  3825    82    82    10 3354.9173 1.001688
#> 13     0     0     0  3999        NA       NA
#> 14  3865    47    47    40 1525.7938 1.031318
#> 15     0     0     0  3999        NA       NA

# Bayesian model averaged main effects for the groups
coef(fit)$main_effects_groups
#>                      group1     group2
#> loose_ends(c1)  -0.95025506 -0.9131525
#> loose_ends(c2)  -2.74190680 -2.2417457
#> loose_ends(c3)  -4.00888563 -3.5538812
#> loose_ends(c4)  -5.31360791 -4.8356803
#> loose_ends(c5)  -7.61366103 -7.4277155
#> loose_ends(c6)  -9.84945159 -9.9581958
#> entertain(c1)   -0.74939401 -1.0389573
#> entertain(c2)   -2.20053482 -2.2787430
#> entertain(c3)   -3.99885872 -3.6908604
#> entertain(c4)   -5.06958045 -5.1797478
#> entertain(c5)   -7.04809994 -6.9801153
#> entertain(c6)   -9.71088014 -9.4708179
#> repetitive(c1)  -0.05089123 -0.2843362
#> repetitive(c2)  -0.49986516 -0.9195860
#> repetitive(c3)  -1.03128548 -1.1357746
#> repetitive(c4)  -1.96644007 -1.7319975
#> repetitive(c5)  -3.56419074 -2.9698960
#> repetitive(c6)  -5.28894869 -4.6859081
#> stimulation(c1) -0.35028051 -0.8660676
#> stimulation(c2) -1.76197416 -1.8651591
#> stimulation(c3) -2.44743966 -2.6910173
#> stimulation(c4) -3.42830954 -3.8713596
#> stimulation(c5) -5.06380697 -5.3137620
#> stimulation(c6) -6.73175409 -7.4268522
#> motivated(c1)   -0.45995233 -0.7025423
#> motivated(c2)   -1.74171483 -1.8709442
#> motivated(c3)   -3.42290624 -3.1548319
#> motivated(c4)   -5.05139285 -4.5817679
#> motivated(c5)   -6.63647459 -6.6852194
#> motivated(c6)   -9.31041096 -8.8970538

# Bayesian model averaged pairwise effects for the groups
coef(fit)$pairwise_effects_groups
#>                            group1     group2
#> loose_ends-entertain   0.16957163 0.16970185
#> loose_ends-repetitive  0.05641909 0.05763802
#> loose_ends-stimulation 0.12225965 0.13369412
#> loose_ends-motivated   0.14026518 0.13967640
#> entertain-repetitive   0.06419560 0.06460285
#> entertain-stimulation  0.10448950 0.11332850
#> entertain-motivated    0.08449783 0.08613167
#> repetitive-stimulation 0.05632405 0.05682491
#> repetitive-motivated   0.13468029 0.13508042
#> stimulation-motivated  0.10767714 0.10810825
# }
```
