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
  difference_prior = bernoulli_prior(0.5),
  difference_probability,
  interaction_prior = cauchy_prior(scale = 1),
  threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
  iter = 2000,
  warmup = 2000,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis"),
  target_accept,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  standardize = FALSE,
  verbose = getOption("bgms.verbose", TRUE),
  progress_callback = NULL,
  pairwise_scale,
  main_alpha,
  main_beta,
  beta_bernoulli_alpha,
  beta_bernoulli_beta,
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

  An indicator prior specification object for difference selection,
  created by one of:

  - [`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md):
    Fixed inclusion probability (default).

  - [`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md):
    Beta-distributed inclusion.

  - [`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md):
    Stochastic Block Model.

  Legacy character strings `"Bernoulli"` and `"Beta-Bernoulli"` are
  still accepted but deprecated. Default: `bernoulli_prior(0.5)`.

- difference_probability:

  **\[deprecated\]** Numeric. Use
  `difference_prior = bernoulli_prior(probability)` instead. Default:
  `0.5`.

- interaction_prior:

  A prior specification object for baseline pairwise interaction
  parameters, created by one of the prior constructor functions:

  - [`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md):
    Cauchy(0, scale) prior (default).

  - [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md):
    Normal(0, scale) prior.

  When supplied, overrides `pairwise_scale`. Default:
  `cauchy_prior(scale = 1)`.

- threshold_prior:

  A prior specification object for threshold (main effect) parameters,
  created by one of the prior constructor functions:

  - [`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md):
    Beta-prime prior (default).

  - [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md):
    Normal(0, scale) prior.

  When supplied, overrides `main_alpha` and `main_beta`. Default:
  `beta_prime_prior(alpha = 0.5, beta = 0.5)`.

- iter:

  Integer. Number of post–warmup iterations per chain. Default: `2e3`.

- warmup:

  Integer. Number of warmup iterations before sampling. Default: `2e3`.

- na_action:

  Character. How to handle missing data: `"listwise"` (drop rows) or
  `"impute"` (impute within Gibbs). Default: `"listwise"`.

- update_method:

  Character. Sampling algorithm: `"adaptive-metropolis"` or `"nuts"`.
  Default: `"nuts"`.

- target_accept:

  Numeric between 0 and 1. Target acceptance rate. Defaults: 0.44
  (Metropolis), 0.80 (NUTS).

- nuts_max_depth:

  Integer. Maximum tree depth for NUTS. Default: `10`.

- learn_mass_matrix:

  Logical. If `TRUE`, adapts a diagonal mass matrix during warmup (NUTS
  only). Default: `TRUE`.

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
  the `interaction_prior` (e.g. `cauchy_prior(scale)`) scale itself
  applying to the unit interval (binary) case. See
  [`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  for details on the adjustment. Default: `FALSE`.

- verbose:

  Logical. If `TRUE`, prints informational messages during data
  processing (e.g., missing data handling, variable recoding). Defaults
  to `getOption("bgms.verbose", TRUE)`. Set
  `options(bgms.verbose = FALSE)` to suppress messages globally.

- progress_callback:

  An optional R function with signature `function(completed, total)`
  that is called at regular intervals during sampling, where `completed`
  is the number of iterations completed across all chains and `total` is
  the total number of iterations. Useful for external front-ends (e.g.,
  JASP) that supply their own progress reporting. When `NULL` (the
  default), no callback is invoked.

- pairwise_scale:

  **\[deprecated\]** Double. Scale of the Cauchy prior for baseline
  pairwise interactions. Use `interaction_prior = cauchy_prior(scale)`
  instead.

- main_alpha, main_beta:

  **\[deprecated\]** Doubles. Shape parameters of the beta-prime prior
  for baseline threshold parameters. Use
  `threshold_prior = beta_prime_prior(alpha, beta)` instead.

- beta_bernoulli_alpha, beta_bernoulli_beta:

  Doubles. Shape parameters of the Beta prior for inclusion
  probabilities in the Beta–Bernoulli model. Defaults: `1`.

- main_difference_model, reference_category, pairwise_difference_scale,
  main_difference_scale, pairwise_difference_prior,
  main_difference_prior, pairwise_difference_probability,
  main_difference_probability, pairwise_beta_bernoulli_alpha,
  pairwise_beta_bernoulli_beta, main_beta_bernoulli_alpha,
  main_beta_bernoulli_beta, interaction_scale, threshold_alpha,
  threshold_beta, burnin, save:

  **\[deprecated\]** Deprecated arguments as of **bgms 0.1.6.0**. Use
  `difference_scale`, `difference_prior`, `difference_probability`,
  `beta_bernoulli_alpha`, `beta_bernoulli_beta`, `baseline_category`,
  `interaction_prior`, `threshold_prior`, and `warmup` instead.

## Value

An S7 object of class `bgmCompare` supporting list-style `$` / `[[`
access for backward compatibility, containing posterior summaries,
posterior mean matrices, and raw MCMC samples. Some
`posterior_summary_*` fields are computed lazily on first access:

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

Group-specific parameters are decomposed into a shared baseline plus
group differences that sum to zero. Difference selection uses
spike-and-slab priors (Bernoulli or Beta-Bernoulli). Parameters are
sampled with NUTS (default) or adaptive Metropolis–Hastings, using the
same multi-stage warmup schedule as
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md).

For full details on model specification, prior choices, and output
interpretation, see the package website at
<https://bayesian-graphical-modelling-lab.github.io/bgms-docs/>.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

[`vignette("comparison", package = "bgms")`](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/comparison.md)
for a worked example.

Other model-fitting:
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

## Examples

``` r
# \dontrun{
# Run bgmCompare on subset of the Boredom dataset
x = Boredom[Boredom$language == "fr", 2:6]
y = Boredom[Boredom$language != "fr", 2:6]

fit = bgmCompare(x, y, chains = 2)
#> Chain 1 (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/4000 (1.2%)
#> Chain 2 (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 48/4000 (1.2%)
#> Total   (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 98/8000 (1.2%)
#> Elapsed: 4s | ETA: 5m 22s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/4000 (2.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 98/4000 (2.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 198/8000 (2.5%)
#> Elapsed: 10s | ETA: 6m 34s
#> Chain 1 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 150/4000 (3.8%)
#> Chain 2 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 113/4000 (2.8%)
#> Total   (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 263/8000 (3.3%)
#> Elapsed: 11s | ETA: 5m 23s
#> Chain 1 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/4000 (6.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 153/4000 (3.8%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 403/8000 (5.0%)
#> Elapsed: 12s | ETA: 3m 46s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/4000 (7.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 254/4000 (6.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 554/8000 (6.9%)
#> Elapsed: 14s | ETA: 3m 8s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/4000 (10.0%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 345/4000 (8.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 745/8000 (9.3%)
#> Elapsed: 15s | ETA: 2m 26s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/4000 (12.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 456/4000 (11.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 956/8000 (11.9%)
#> Elapsed: 16s | ETA: 1m 57s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 600/4000 (15.0%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 536/4000 (13.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1136/8000 (14.2%)
#> Elapsed: 17s | ETA: 1m 42s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/4000 (17.5%)
#> Chain 2 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 626/4000 (15.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1326/8000 (16.6%)
#> Elapsed: 17s | ETA: 1m 25s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/4000 (20.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 714/4000 (17.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1514/8000 (18.9%)
#> Elapsed: 18s | ETA: 1m 17s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/4000 (22.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 833/4000 (20.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1733/8000 (21.7%)
#> Elapsed: 19s | ETA: 1m 8s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1000/4000 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 924/4000 (23.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1924/8000 (24.1%)
#> Elapsed: 20s | ETA: 1m 3s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/4000 (27.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1019/4000 (25.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2119/8000 (26.5%)
#> Elapsed: 20s | ETA: 56s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1200/4000 (30.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1119/4000 (28.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2319/8000 (29.0%)
#> Elapsed: 21s | ETA: 51s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1300/4000 (32.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1213/4000 (30.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2513/8000 (31.4%)
#> Elapsed: 21s | ETA: 46s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1400/4000 (35.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1324/4000 (33.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2724/8000 (34.1%)
#> Elapsed: 22s | ETA: 43s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1500/4000 (37.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1437/4000 (35.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2937/8000 (36.7%)
#> Elapsed: 23s | ETA: 40s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1600/4000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1540/4000 (38.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3140/8000 (39.2%)
#> Elapsed: 23s | ETA: 36s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/4000 (42.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1665/4000 (41.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3365/8000 (42.1%)
#> Elapsed: 27s | ETA: 37s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/4000 (43.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 1723/4000 (43.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3473/8000 (43.4%)
#> Elapsed: 27s | ETA: 35s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/4000 (45.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1774/4000 (44.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3574/8000 (44.7%)
#> Elapsed: 28s | ETA: 35s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1850/4000 (46.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1825/4000 (45.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 3675/8000 (45.9%)
#> Elapsed: 29s | ETA: 34s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/4000 (47.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1879/4000 (47.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3779/8000 (47.2%)
#> Elapsed: 30s | ETA: 34s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/4000 (50.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1946/4000 (48.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3946/8000 (49.3%)
#> Elapsed: 30s | ETA: 31s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2100/4000 (52.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 2028/4000 (50.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4128/8000 (51.6%)
#> Elapsed: 31s | ETA: 29s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/4000 (55.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 2125/4000 (53.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4325/8000 (54.1%)
#> Elapsed: 31s | ETA: 26s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2300/4000 (57.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 2226/4000 (55.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4526/8000 (56.6%)
#> Elapsed: 32s | ETA: 25s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2400/4000 (60.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 2330/4000 (58.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4730/8000 (59.1%)
#> Elapsed: 32s | ETA: 22s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2500/4000 (62.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 2430/4000 (60.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4930/8000 (61.6%)
#> Elapsed: 33s | ETA: 21s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2600/4000 (65.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 2526/4000 (63.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5126/8000 (64.1%)
#> Elapsed: 34s | ETA: 19s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2700/4000 (67.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 2615/4000 (65.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5315/8000 (66.4%)
#> Elapsed: 34s | ETA: 17s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2800/4000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 2711/4000 (67.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5511/8000 (68.9%)
#> Elapsed: 35s | ETA: 16s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2900/4000 (72.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 2808/4000 (70.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5708/8000 (71.4%)
#> Elapsed: 35s | ETA: 14s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3000/4000 (75.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 2904/4000 (72.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5904/8000 (73.8%)
#> Elapsed: 36s | ETA: 13s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3100/4000 (77.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2995/4000 (74.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 6095/8000 (76.2%)
#> Elapsed: 36s | ETA: 11s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3200/4000 (80.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3093/4000 (77.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 6293/8000 (78.7%)
#> Elapsed: 37s | ETA: 10s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3300/4000 (82.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3186/4000 (79.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 6486/8000 (81.1%)
#> Elapsed: 37s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3400/4000 (85.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3288/4000 (82.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 6688/8000 (83.6%)
#> Elapsed: 38s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3500/4000 (87.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3385/4000 (84.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 6885/8000 (86.1%)
#> Elapsed: 39s | ETA: 6s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3600/4000 (90.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3477/4000 (86.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 7077/8000 (88.5%)
#> Elapsed: 39s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3700/4000 (92.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3577/4000 (89.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 7277/8000 (91.0%)
#> Elapsed: 40s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3800/4000 (95.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3674/4000 (91.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 7474/8000 (93.4%)
#> Elapsed: 40s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3900/4000 (97.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3758/4000 (94.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 7658/8000 (95.7%)
#> Elapsed: 41s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 3849/4000 (96.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 7849/8000 (98.1%)
#> Elapsed: 41s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 42s | ETA: 0s

# Posterior inclusion probabilities
summary(fit)$indicator
#>                            parameter    mean        mcse        sd
#> 1                  loose_ends (main) 1.00000          NA 0.0000000
#> 2    loose_ends-entertain (pairwise) 0.02075 0.002465560 0.1425463
#> 3   loose_ends-repetitive (pairwise) 0.03400 0.004292323 0.1812291
#> 4  loose_ends-stimulation (pairwise) 0.18525 0.015385589 0.3885002
#> 5    loose_ends-motivated (pairwise) 0.02700 0.003277419 0.1620833
#> 6                   entertain (main) 1.00000          NA 0.0000000
#> 7    entertain-repetitive (pairwise) 0.02550 0.003236806 0.1576380
#> 8   entertain-stimulation (pairwise) 0.19600 0.014945503 0.3969685
#> 9     entertain-motivated (pairwise) 0.05275 0.005091077 0.2235340
#> 10                 repetitive (main) 1.00000          NA 0.0000000
#> 11 repetitive-stimulation (pairwise) 0.03250 0.003959743 0.1773239
#> 12   repetitive-motivated (pairwise) 0.02575 0.003076864 0.1583886
#> 13                stimulation (main) 1.00000          NA 0.0000000
#> 14  stimulation-motivated (pairwise) 0.02150 0.002717689 0.1450440
#> 15                  motivated (main) 1.00000          NA 0.0000000
#>    n0->0 n0->1 n1->0 n1->1 n_eff_mixt      Rhat
#> 1      0     0     0  3999         NA        NA
#> 2   3842    74    74     9  3342.5691 0.9998241
#> 3   3782    81    81    55  1782.6716 1.0015384
#> 4   3092   166   166   575   637.6088 1.0056350
#> 5   3812    80    79    28  2445.7532 1.0025602
#> 6      0     0     0  3999         NA        NA
#> 7   3823    74    74    28  2371.8576 1.0153935
#> 8   3026   189   189   595   705.4903 1.0011247
#> 9   3658   130   130    81  1927.8257 1.0125021
#> 10     0     0     0  3999         NA        NA
#> 11  3785    84    84    46  2005.3968 1.0335602
#> 12  3817    79    80    23  2649.9102 1.0021137
#> 13     0     0     0  3999         NA        NA
#> 14  3843    70    70    16  2848.3920 1.0008539
#> 15     0     0     0  3999         NA        NA

# Bayesian model averaged main effects for the groups
coef(fit)$main_effects_groups
#>                      group1     group2
#> loose_ends(c1)  -0.95016614 -0.9117283
#> loose_ends(c2)  -2.74702132 -2.2364094
#> loose_ends(c3)  -4.01588140 -3.5409158
#> loose_ends(c4)  -5.32208706 -4.8121972
#> loose_ends(c5)  -7.62241777 -7.3964166
#> loose_ends(c6)  -9.86066047 -9.9215537
#> entertain(c1)   -0.74624512 -1.0395753
#> entertain(c2)   -2.19242809 -2.2775450
#> entertain(c3)   -3.98731883 -3.6870683
#> entertain(c4)   -5.05364349 -5.1653601
#> entertain(c5)   -7.02064761 -6.9615941
#> entertain(c6)   -9.68859337 -9.4460012
#> repetitive(c1)  -0.04906668 -0.2836083
#> repetitive(c2)  -0.50124957 -0.9170959
#> repetitive(c3)  -1.03147910 -1.1356898
#> repetitive(c4)  -1.96462050 -1.7341976
#> repetitive(c5)  -3.56189321 -2.9714787
#> repetitive(c6)  -5.29598764 -4.6905399
#> stimulation(c1) -0.34748453 -0.8538983
#> stimulation(c2) -1.75528308 -1.8505253
#> stimulation(c3) -2.44253926 -2.6726572
#> stimulation(c4) -3.42234052 -3.8450436
#> stimulation(c5) -5.05128231 -5.2768116
#> stimulation(c6) -6.72244716 -7.3737635
#> motivated(c1)   -0.46094605 -0.7073385
#> motivated(c2)   -1.74051862 -1.8692412
#> motivated(c3)   -3.41909357 -3.1581994
#> motivated(c4)   -5.04804954 -4.5797147
#> motivated(c5)   -6.63256734 -6.6827406
#> motivated(c6)   -9.29664271 -8.8925731

# Bayesian model averaged pairwise effects for the groups
coef(fit)$pairwise_effects_groups
#>                            group1     group2
#> loose_ends-entertain   0.16906685 0.16925516
#> loose_ends-repetitive  0.05688985 0.05778803
#> loose_ends-stimulation 0.12297675 0.13176343
#> loose_ends-motivated   0.14052529 0.13992543
#> entertain-repetitive   0.06415804 0.06472679
#> entertain-stimulation  0.10422423 0.11309280
#> entertain-motivated    0.08383659 0.08549094
#> repetitive-stimulation 0.05581212 0.05657889
#> repetitive-motivated   0.13477296 0.13528928
#> stimulation-motivated  0.10738071 0.10765510
# }
```
