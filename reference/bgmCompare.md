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
  learn_mass_matrix = FALSE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
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

  Logical. If `TRUE`, adapt the mass matrix during warmup (HMC/NUTS
  only). Default: `FALSE`.

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
group differences in MRFs was introduced in Marsman et al. (in press) ,
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
the ordinal Markov random field.” *Psychometrika*, **90**, 146–-182.  
  
Marsman M, Waldorp LJ, Sekulovski N, Haslbeck JMB (in press). “Bayes
factor tests for group differences in ordinal and binary graphical
models.” *Psychometrika*.

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
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2200 (2.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 46/2200 (2.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 45/2200 (2.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 47/2200 (2.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 188/8800 (2.1%)
#> Elapsed: 10s | ETA: 7m 38s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2200 (4.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 101/2200 (4.6%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 98/2200 (4.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 102/2200 (4.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 401/8800 (4.6%)
#> Elapsed: 26s | ETA: 9m 4s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 150/2200 (6.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 152/2200 (6.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 147/2200 (6.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 147/2200 (6.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 596/8800 (6.8%)
#> Elapsed: 40s | ETA: 9m 10s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/2200 (9.1%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 203/2200 (9.2%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 194/2200 (8.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 202/2200 (9.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 799/8800 (9.1%)
#> Elapsed: 55s | ETA: 9m 10s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/2200 (11.4%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 253/2200 (11.5%)
#> Chain 3 (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 238/2200 (10.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/2200 (11.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 991/8800 (11.3%)
#> Elapsed: 1m 7s | ETA: 8m 47s
#> Chain 1 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2200 (13.6%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2200 (13.6%)
#> Chain 3 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 282/2200 (12.8%)
#> Chain 4 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 292/2200 (13.3%)
#> Total   (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1174/8800 (13.3%)
#> Elapsed: 1m 19s | ETA: 8m 33s
#> Chain 1 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 350/2200 (15.9%)
#> Chain 2 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 354/2200 (16.1%)
#> Chain 3 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 338/2200 (15.4%)
#> Chain 4 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 352/2200 (16.0%)
#> Total   (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1394/8800 (15.8%)
#> Elapsed: 1m 34s | ETA: 8m 19s
#> Chain 1 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2200 (18.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 404/2200 (18.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 383/2200 (17.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2200 (18.2%)
#> Total   (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1587/8800 (18.0%)
#> Elapsed: 1m 45s | ETA: 7m 57s
#> Chain 1 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/2200 (20.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 453/2200 (20.6%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 433/2200 (19.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/2200 (20.5%)
#> Total   (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1786/8800 (20.3%)
#> Elapsed: 1m 58s | ETA: 7m 43s
#> Chain 1 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2200 (22.7%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 496/2200 (22.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 483/2200 (22.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 492/2200 (22.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1971/8800 (22.4%)
#> Elapsed: 2m 9s | ETA: 7m 26s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 550/2200 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 546/2200 (24.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 538/2200 (24.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 545/2200 (24.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2179/8800 (24.8%)
#> Elapsed: 2m 21s | ETA: 7m 8s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 600/2200 (27.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 591/2200 (26.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 581/2200 (26.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 584/2200 (26.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2356/8800 (26.8%)
#> Elapsed: 2m 32s | ETA: 6m 55s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 650/2200 (29.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 645/2200 (29.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 628/2200 (28.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 631/2200 (28.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2554/8800 (29.0%)
#> Elapsed: 2m 45s | ETA: 6m 43s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/2200 (31.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 691/2200 (31.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 670/2200 (30.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 682/2200 (31.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2743/8800 (31.2%)
#> Elapsed: 2m 56s | ETA: 6m 28s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 750/2200 (34.1%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 741/2200 (33.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 718/2200 (32.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 722/2200 (32.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2931/8800 (33.3%)
#> Elapsed: 3m 7s | ETA: 6m 14s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2200 (36.4%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 788/2200 (35.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 768/2200 (34.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 775/2200 (35.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3131/8800 (35.6%)
#> Elapsed: 3m 19s | ETA: 6m
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/2200 (38.6%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 846/2200 (38.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 819/2200 (37.2%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 835/2200 (38.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3350/8800 (38.1%)
#> Elapsed: 3m 33s | ETA: 5m 46s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2200 (40.9%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 899/2200 (40.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 870/2200 (39.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 889/2200 (40.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 3558/8800 (40.4%)
#> Elapsed: 3m 46s | ETA: 5m 32s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 950/2200 (43.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 950/2200 (43.2%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 917/2200 (41.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 937/2200 (42.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3754/8800 (42.7%)
#> Elapsed: 3m 58s | ETA: 5m 19s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1000/2200 (45.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 996/2200 (45.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 959/2200 (43.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 988/2200 (44.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3943/8800 (44.8%)
#> Elapsed: 4m 10s | ETA: 5m 7s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1050/2200 (47.7%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1046/2200 (47.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1018/2200 (46.3%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1040/2200 (47.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4154/8800 (47.2%)
#> Elapsed: 4m 24s | ETA: 4m 55s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/2200 (50.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1095/2200 (49.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1073/2200 (48.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1091/2200 (49.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4359/8800 (49.5%)
#> Elapsed: 4m 37s | ETA: 4m 42s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2200 (52.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1143/2200 (52.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1098/2200 (49.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1130/2200 (51.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4521/8800 (51.4%)
#> Elapsed: 4m 44s | ETA: 4m 28s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1200/2200 (54.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 1182/2200 (53.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1137/2200 (51.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1183/2200 (53.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 4702/8800 (53.4%)
#> Elapsed: 4m 50s | ETA: 4m 12s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1250/2200 (56.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1232/2200 (56.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1198/2200 (54.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1226/2200 (55.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 4906/8800 (55.8%)
#> Elapsed: 4m 57s | ETA: 3m 55s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1300/2200 (59.1%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1275/2200 (58.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1244/2200 (56.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1274/2200 (57.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 5093/8800 (57.9%)
#> Elapsed: 5m 4s | ETA: 3m 41s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1350/2200 (61.4%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1326/2200 (60.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1306/2200 (59.4%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1339/2200 (60.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 5321/8800 (60.5%)
#> Elapsed: 5m 13s | ETA: 3m 24s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1400/2200 (63.6%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1376/2200 (62.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1354/2200 (61.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1400/2200 (63.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 5530/8800 (62.8%)
#> Elapsed: 5m 21s | ETA: 3m 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1450/2200 (65.9%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1420/2200 (64.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1397/2200 (63.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1447/2200 (65.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5714/8800 (64.9%)
#> Elapsed: 5m 27s | ETA: 2m 56s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1500/2200 (68.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1463/2200 (66.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1445/2200 (65.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1491/2200 (67.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5899/8800 (67.0%)
#> Elapsed: 5m 35s | ETA: 2m 44s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1550/2200 (70.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1507/2200 (68.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1496/2200 (68.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1539/2200 (70.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6092/8800 (69.2%)
#> Elapsed: 5m 43s | ETA: 2m 32s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1600/2200 (72.7%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1560/2200 (70.9%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1540/2200 (70.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1591/2200 (72.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6291/8800 (71.5%)
#> Elapsed: 5m 50s | ETA: 2m 19s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2200 (75.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1610/2200 (73.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1585/2200 (72.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1641/2200 (74.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 6486/8800 (73.7%)
#> Elapsed: 5m 58s | ETA: 2m 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/2200 (77.3%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1668/2200 (75.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1635/2200 (74.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1694/2200 (77.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 6697/8800 (76.1%)
#> Elapsed: 6m 6s | ETA: 1m 54s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/2200 (79.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1709/2200 (77.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1672/2200 (76.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1736/2200 (78.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 6867/8800 (78.0%)
#> Elapsed: 6m 12s | ETA: 1m 44s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/2200 (81.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1756/2200 (79.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1715/2200 (78.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1787/2200 (81.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 7058/8800 (80.2%)
#> Elapsed: 6m 20s | ETA: 1m 33s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1850/2200 (84.1%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1796/2200 (81.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1762/2200 (80.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1846/2200 (83.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7254/8800 (82.4%)
#> Elapsed: 6m 26s | ETA: 1m 22s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2200 (86.4%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1854/2200 (84.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 1826/2200 (83.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1906/2200 (86.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 7486/8800 (85.1%)
#> Elapsed: 6m 36s | ETA: 1m 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1950/2200 (88.6%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1907/2200 (86.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 1877/2200 (85.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1950/2200 (88.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7684/8800 (87.3%)
#> Elapsed: 6m 44s | ETA: 59s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 2000/2200 (90.9%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1951/2200 (88.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1919/2200 (87.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 2000/2200 (90.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7870/8800 (89.4%)
#> Elapsed: 6m 51s | ETA: 49s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2050/2200 (93.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 1997/2200 (90.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1964/2200 (89.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2044/2200 (92.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8055/8800 (91.5%)
#> Elapsed: 6m 58s | ETA: 39s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2100/2200 (95.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2055/2200 (93.4%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2022/2200 (91.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2091/2200 (95.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8268/8800 (94.0%)
#> Elapsed: 7m 6s | ETA: 27s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 2150/2200 (97.7%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2097/2200 (95.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2065/2200 (93.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2135/2200 (97.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 8447/8800 (96.0%)
#> Elapsed: 7m 13s | ETA: 18s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 2149/2200 (97.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2129/2200 (96.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2186/2200 (99.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 8664/8800 (98.5%)
#> Elapsed: 7m 21s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8800/8800 (100.0%)
#> Elapsed: 7m 29s | ETA: 0s
#> NUTS Diagnostics Summary:
#>   Total divergences:         0 
#>   Max tree depth hits:       4 
#>   Min E-BFMI across chains:  1.469 
#> Note: 0.10% of transitions hit the maximum tree depth (4 of 4000).
#> Check efficiency metrics such as effective sample size (ESS) to ensure
#> sufficient exploration of the posterior.

# Posterior inclusion probabilities
summary(fit)$indicator
#>                            parameter    mean         sd         mcse
#> 1                  loose_ends (main) 0.00625 0.07880950 0.0024824187
#> 2    loose_ends-entertain (pairwise) 0.02225 0.14749555 0.0023059720
#> 3   loose_ends-repetitive (pairwise) 0.48275 0.49970235 0.0186102123
#> 4  loose_ends-stimulation (pairwise) 0.04575 0.20894243 0.0033923711
#> 5    loose_ends-motivated (pairwise) 0.01925 0.13740247 0.0023374976
#> 6                   entertain (main) 0.00175 0.04179638 0.0007619244
#> 7    entertain-repetitive (pairwise) 0.12525 0.33100217 0.0083343098
#> 8   entertain-stimulation (pairwise) 0.07075 0.25640678 0.0043579995
#> 9     entertain-motivated (pairwise) 0.11400 0.31781126 0.0073160717
#> 10                 repetitive (main) 0.02700 0.16208331 0.0086387280
#> 11 repetitive-stimulation (pairwise) 0.02475 0.15536228 0.0026019808
#> 12   repetitive-motivated (pairwise) 0.42850 0.49486134 0.0216635435
#> 13                stimulation (main) 0.00375 0.06112232 0.0012663693
#> 14  stimulation-motivated (pairwise) 0.01575 0.12450678 0.0020334255
#> 15                  motivated (main) 0.00325 0.05691606 0.0012337458
#>    n0->0 n0->1 n1->0 n1->1     n_eff      Rhat
#> 1   3964    10    10    15 1007.8760 1.0392609
#> 2   3822    88    88     1 4091.1913 1.0075779
#> 3   1763   305   305  1626  720.9749 1.0064673
#> 4   3646   170   170    13 3793.5566 1.0007457
#> 5   3852    70    70     7 3455.3074 1.0087667
#> 6   3986     6     6     1 3009.2180 1.0587785
#> 7   3251   247   248   253 1577.3294 1.0039605
#> 8   3472   244   244    39 3461.6621 1.0011649
#> 9   3284   259   259   197 1887.0468 1.0107878
#> 10  3874    17    17    91  352.0278 1.0239003
#> 11  3809    91    91     8 3565.1925 0.9998473
#> 12  2060   226   226  1487  521.8048 1.0030641
#> 13  3973    11    11     4 2329.5865 1.0247554
#> 14  3876    60    60     3 3749.1211 1.0093500
#> 15  3977     9     9     4 2128.2283 1.0398287

# Bayesian model averaged main effects for the groups
coef(fit)$main_effects_groups
#>                      group1      group2
#> loose_ends(c1)   -0.9359677  -0.9361159
#> loose_ends(c2)   -2.5175513  -2.5152350
#> loose_ends(c3)   -3.7967203  -3.7945120
#> loose_ends(c4)   -5.0808053  -5.0787380
#> loose_ends(c5)   -7.6032434  -7.6036143
#> loose_ends(c6)  -10.1017243 -10.1040494
#> entertain(c1)    -0.8679665  -0.8687533
#> entertain(c2)    -2.2408504  -2.2410082
#> entertain(c3)    -3.8071858  -3.8066712
#> entertain(c4)    -5.1542741  -5.1547468
#> entertain(c5)    -7.0196039  -7.0194251
#> entertain(c6)    -9.5430678  -9.5424886
#> repetitive(c1)   -0.1352162  -0.1425647
#> repetitive(c2)   -0.6558459  -0.6688758
#> repetitive(c3)   -1.0827212  -1.0866392
#> repetitive(c4)   -1.8624076  -1.8577802
#> repetitive(c5)   -3.2567125  -3.2438692
#> repetitive(c6)   -5.0336933  -5.0199061
#> stimulation(c1)  -0.5423604  -0.5443373
#> stimulation(c2)  -1.7853398  -1.7856367
#> stimulation(c3)  -2.5267431  -2.5278867
#> stimulation(c4)  -3.6092398  -3.6107213
#> stimulation(c5)  -5.0871853  -5.0878656
#> stimulation(c6)  -7.0499181  -7.0523056
#> motivated(c1)    -0.5549315  -0.5560584
#> motivated(c2)    -1.8034165  -1.8043094
#> motivated(c3)    -3.2901070  -3.2895281
#> motivated(c4)    -4.7852350  -4.7838308
#> motivated(c5)    -6.7792060  -6.7799159
#> motivated(c6)    -9.1506426  -9.1496620

# Bayesian model averaged pairwise effects for the groups
coef(fit)$pairwise_effects_groups
#>                            group1     group2
#> loose_ends-entertain   0.16890271 0.16914940
#> loose_ends-repetitive  0.04797785 0.06836349
#> loose_ends-stimulation 0.12622236 0.12726009
#> loose_ends-motivated   0.14088843 0.14054627
#> entertain-repetitive   0.06446456 0.06859517
#> entertain-stimulation  0.10799696 0.10989246
#> entertain-motivated    0.08437295 0.08818267
#> repetitive-stimulation 0.05590889 0.05630306
#> repetitive-motivated   0.12859022 0.14577087
#> stimulation-motivated  0.10751862 0.10764585
# }
```
