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
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2200 (2.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 55/2200 (2.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 54/2200 (2.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 52/2200 (2.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 211/8800 (2.4%)
#> Elapsed: 12s | ETA: 8m 8s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2200 (4.5%)
#> Chain 2 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 121/2200 (5.5%)
#> Chain 3 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 183/2200 (8.3%)
#> Chain 4 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 137/2200 (6.2%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 541/8800 (6.1%)
#> Elapsed: 26s | ETA: 6m 36s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 150/2200 (6.8%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 189/2200 (8.6%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 255/2200 (11.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 203/2200 (9.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 797/8800 (9.1%)
#> Elapsed: 27s | ETA: 4m 31s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/2200 (9.1%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 251/2200 (11.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 311/2200 (14.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 264/2200 (12.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1026/8800 (11.7%)
#> Elapsed: 28s | ETA: 3m 32s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/2200 (11.4%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 288/2200 (13.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 363/2200 (16.5%)
#> Chain 4 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 294/2200 (13.4%)
#> Total   (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1195/8800 (13.6%)
#> Elapsed: 29s | ETA: 3m 4s
#> Chain 1 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2200 (13.6%)
#> Chain 2 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 352/2200 (16.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 442/2200 (20.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 361/2200 (16.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1455/8800 (16.5%)
#> Elapsed: 30s | ETA: 2m 31s
#> Chain 1 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 350/2200 (15.9%)
#> Chain 2 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 409/2200 (18.6%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 491/2200 (22.3%)
#> Chain 4 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 409/2200 (18.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1659/8800 (18.9%)
#> Elapsed: 31s | ETA: 2m 13s
#> Chain 1 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2200 (18.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 458/2200 (20.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 545/2200 (24.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 461/2200 (21.0%)
#> Total   (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1864/8800 (21.2%)
#> Elapsed: 32s | ETA: 1m 59s
#> Chain 1 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/2200 (20.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 493/2200 (22.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 592/2200 (26.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 493/2200 (22.4%)
#> Total   (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2028/8800 (23.0%)
#> Elapsed: 32s | ETA: 1m 46s
#> Chain 1 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2200 (22.7%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 573/2200 (26.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 671/2200 (30.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 571/2200 (26.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2315/8800 (26.3%)
#> Elapsed: 33s | ETA: 1m 32s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 550/2200 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 638/2200 (29.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 741/2200 (33.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 644/2200 (29.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2573/8800 (29.2%)
#> Elapsed: 34s | ETA: 1m 22s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 600/2200 (27.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 694/2200 (31.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 798/2200 (36.3%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 703/2200 (32.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2795/8800 (31.8%)
#> Elapsed: 35s | ETA: 1m 15s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 650/2200 (29.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 747/2200 (34.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 851/2200 (38.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 762/2200 (34.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3010/8800 (34.2%)
#> Elapsed: 36s | ETA: 1m 9s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/2200 (31.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 811/2200 (36.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 885/2200 (40.2%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 820/2200 (37.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3216/8800 (36.5%)
#> Elapsed: 36s | ETA: 1m 2s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 750/2200 (34.1%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 860/2200 (39.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 916/2200 (41.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 868/2200 (39.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3394/8800 (38.6%)
#> Elapsed: 37s | ETA: 59s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2200 (36.4%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 895/2200 (40.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 952/2200 (43.3%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 908/2200 (41.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 3555/8800 (40.4%)
#> Elapsed: 38s | ETA: 56s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/2200 (38.6%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 925/2200 (42.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 993/2200 (45.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 939/2200 (42.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3707/8800 (42.1%)
#> Elapsed: 38s | ETA: 52s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2200 (40.9%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 968/2200 (44.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1036/2200 (47.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 989/2200 (45.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3893/8800 (44.2%)
#> Elapsed: 39s | ETA: 49s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 950/2200 (43.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1022/2200 (46.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1079/2200 (49.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1032/2200 (46.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4083/8800 (46.4%)
#> Elapsed: 40s | ETA: 46s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1000/2200 (45.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1062/2200 (48.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1138/2200 (51.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1073/2200 (48.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 4273/8800 (48.6%)
#> Elapsed: 41s | ETA: 43s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1050/2200 (47.7%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1136/2200 (51.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1243/2200 (56.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1146/2200 (52.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4575/8800 (52.0%)
#> Elapsed: 42s | ETA: 39s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/2200 (50.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1238/2200 (56.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1342/2200 (61.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1245/2200 (56.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 4925/8800 (56.0%)
#> Elapsed: 44s | ETA: 35s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2200 (52.3%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1288/2200 (58.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1390/2200 (63.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1292/2200 (58.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 5120/8800 (58.2%)
#> Elapsed: 44s | ETA: 32s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1200/2200 (54.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1336/2200 (60.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1439/2200 (65.4%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1338/2200 (60.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 5313/8800 (60.4%)
#> Elapsed: 45s | ETA: 30s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1250/2200 (56.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1388/2200 (63.1%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1491/2200 (67.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1389/2200 (63.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 5518/8800 (62.7%)
#> Elapsed: 45s | ETA: 27s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1300/2200 (59.1%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1441/2200 (65.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1543/2200 (70.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1443/2200 (65.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 5727/8800 (65.1%)
#> Elapsed: 46s | ETA: 25s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1350/2200 (61.4%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1490/2200 (67.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1596/2200 (72.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1491/2200 (67.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5927/8800 (67.4%)
#> Elapsed: 47s | ETA: 23s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1400/2200 (63.6%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1548/2200 (70.4%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2200 (75.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1539/2200 (70.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6137/8800 (69.7%)
#> Elapsed: 47s | ETA: 20s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1450/2200 (65.9%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1601/2200 (72.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1699/2200 (77.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1587/2200 (72.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6337/8800 (72.0%)
#> Elapsed: 48s | ETA: 19s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1500/2200 (68.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2200 (75.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1757/2200 (79.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1638/2200 (74.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6545/8800 (74.4%)
#> Elapsed: 49s | ETA: 17s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1550/2200 (70.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1703/2200 (77.4%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1810/2200 (82.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1688/2200 (76.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6751/8800 (76.7%)
#> Elapsed: 49s | ETA: 15s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1600/2200 (72.7%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/2200 (79.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1863/2200 (84.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1736/2200 (78.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6949/8800 (79.0%)
#> Elapsed: 50s | ETA: 13s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2200 (75.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1809/2200 (82.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1926/2200 (87.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1789/2200 (81.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7174/8800 (81.5%)
#> Elapsed: 51s | ETA: 12s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/2200 (77.3%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1865/2200 (84.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 1983/2200 (90.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1844/2200 (83.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7392/8800 (84.0%)
#> Elapsed: 51s | ETA: 10s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/2200 (79.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1922/2200 (87.4%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2043/2200 (92.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1899/2200 (86.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7614/8800 (86.5%)
#> Elapsed: 52s | ETA: 8s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/2200 (81.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1976/2200 (89.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2099/2200 (95.4%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1952/2200 (88.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7827/8800 (88.9%)
#> Elapsed: 53s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1850/2200 (84.1%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2031/2200 (92.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 2149/2200 (97.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 2007/2200 (91.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8037/8800 (91.3%)
#> Elapsed: 53s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2200 (86.4%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2082/2200 (94.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2062/2200 (93.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 8244/8800 (93.7%)
#> Elapsed: 54s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 2000/2200 (90.9%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 2161/2200 (98.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2188/2200 (99.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8549/8800 (97.1%)
#> Elapsed: 55s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2100/2200 (95.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8700/8800 (98.9%)
#> Elapsed: 56s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8800/8800 (100.0%)
#> Elapsed: 56s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8800/8800 (100.0%)
#> Elapsed: 56s | ETA: 0s
#> NUTS Diagnostics Summary:
#>   Total divergences:         1 
#>   Max tree depth hits:       0 
#>   Min E-BFMI across chains:  0.957 
#> Note: 0.025% of transitions ended with a divergence (1 of 4000).
#> Check R-hat and effective sample size (ESS) to ensure the chains are
#> mixing well.

# Posterior inclusion probabilities
summary(fit)$indicator
#>                            parameter    mean        sd        mcse
#> 1                  loose_ends (main) 1.00000 0.0000000          NA
#> 2    loose_ends-entertain (pairwise) 0.01825 0.1338542 0.002196714
#> 3   loose_ends-repetitive (pairwise) 0.04000 0.1959592 0.004929257
#> 4  loose_ends-stimulation (pairwise) 0.25550 0.4361419 0.019977919
#> 5    loose_ends-motivated (pairwise) 0.02875 0.1671031 0.003371249
#> 6                   entertain (main) 1.00000 0.0000000          NA
#> 7    entertain-repetitive (pairwise) 0.04075 0.1977105 0.005986537
#> 8   entertain-stimulation (pairwise) 0.18875 0.3913099 0.017643556
#> 9     entertain-motivated (pairwise) 0.04925 0.2163896 0.005708486
#> 10                 repetitive (main) 1.00000 0.0000000          NA
#> 11 repetitive-stimulation (pairwise) 0.02075 0.1425463 0.002796884
#> 12   repetitive-motivated (pairwise) 0.02625 0.1598779 0.003252704
#> 13                stimulation (main) 1.00000 0.0000000          NA
#> 14  stimulation-motivated (pairwise) 0.02125 0.1442166 0.002610804
#> 15                  motivated (main) 1.00000 0.0000000          NA
#>    n0->0 n0->1 n1->0 n1->1     n_eff     Rhat
#> 1      0     0     0  3999        NA       NA
#> 2   3857    69    69     4 3712.9295 1.009502
#> 3   3752    87    87    73 1580.4044 1.013567
#> 4   2815   162   162   860  476.6012 1.012846
#> 5   3799    85    85    30 2456.8987 1.007793
#> 6      0     0     0  3999        NA       NA
#> 7   3769    67    67    96 1090.7070 1.063158
#> 8   3111   134   134   620  491.8916 1.014424
#> 9   3703    99    99    98 1436.9138 1.001459
#> 10     0     0     0  3999        NA       NA
#> 11  3852    64    64    19 2597.5439 1.055136
#> 12  3817    77    77    28 2415.9493 1.005520
#> 13     0     0     0  3999        NA       NA
#> 14  3842    72    72    13 3051.2807 1.009465
#> 15     0     0     0  3999        NA       NA

# Bayesian model averaged main effects for the groups
coef(fit)$main_effects_groups
#>                     group1     group2
#> loose_ends(c1)  -0.9459476 -0.9051717
#> loose_ends(c2)  -2.7291997 -2.2334680
#> loose_ends(c3)  -3.9951805 -3.5433756
#> loose_ends(c4)  -5.2946544 -4.8180857
#> loose_ends(c5)  -7.5877125 -7.4060686
#> loose_ends(c6)  -9.8089630 -9.9354416
#> entertain(c1)   -0.7467576 -1.0389311
#> entertain(c2)   -2.1940643 -2.2774013
#> entertain(c3)   -3.9875284 -3.6816013
#> entertain(c4)   -5.0540803 -5.1660209
#> entertain(c5)   -7.0295576 -6.9623738
#> entertain(c6)   -9.6763147 -9.4396278
#> repetitive(c1)  -0.0461532 -0.2804461
#> repetitive(c2)  -0.4975123 -0.9130021
#> repetitive(c3)  -1.0281447 -1.1327559
#> repetitive(c4)  -1.9590802 -1.7263721
#> repetitive(c5)  -3.5529729 -2.9710856
#> repetitive(c6)  -5.2795611 -4.6836230
#> stimulation(c1) -0.3517567 -0.8553277
#> stimulation(c2) -1.7551303 -1.8495732
#> stimulation(c3) -2.4342393 -2.6756204
#> stimulation(c4) -3.4086345 -3.8512182
#> stimulation(c5) -5.0291214 -5.2880739
#> stimulation(c6) -6.6905115 -7.3868133
#> motivated(c1)   -0.4607087 -0.7064547
#> motivated(c2)   -1.7434850 -1.8718262
#> motivated(c3)   -3.4267307 -3.1608480
#> motivated(c4)   -5.0560447 -4.5834851
#> motivated(c5)   -6.6354589 -6.6888996
#> motivated(c6)   -9.3130441 -8.9030977

# Bayesian model averaged pairwise effects for the groups
coef(fit)$pairwise_effects_groups
#>                            group1     group2
#> loose_ends-entertain   0.16893198 0.16901540
#> loose_ends-repetitive  0.05625349 0.05726649
#> loose_ends-stimulation 0.12112503 0.13377721
#> loose_ends-motivated   0.14030417 0.13973978
#> entertain-repetitive   0.06398351 0.06486809
#> entertain-stimulation  0.10363000 0.11255487
#> entertain-motivated    0.08449238 0.08593339
#> repetitive-stimulation 0.05584406 0.05622727
#> repetitive-motivated   0.13493574 0.13551281
#> stimulation-motivated  0.10735979 0.10764248
# }
```
