# Extract Posterior Mean Log-Odds

Retrieves the posterior mean pairwise interaction matrix for discrete
variables from a model fitted with
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md).
These are the log-odds parameters of the discrete (Markov random field)
block. GGM models have no discrete variables and return `NULL`.

## Usage

``` r
extract_log_odds(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class `bgms` (from
  [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)).

## Value

A named numeric matrix of posterior mean log-odds interactions, or
`NULL` for GGM models.

- OMRF:

  A symmetric matrix with zero diagonal and one row and column per
  variable.

- Mixed MRF:

  A symmetric matrix with zero diagonal and one row and column per
  discrete variable.

- GGM:

  `NULL` (invisibly).

## Details

Extract Posterior Mean Log-Odds (Pairwise Interactions)

## See also

[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_precision()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_precision.md)

Other extractors:
[`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md),
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md),
[`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md),
[`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md),
[`extract_main_effects()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_main_effects.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_partial_correlations()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_partial_correlations.md),
[`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md),
[`extract_precision()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_precision.md),
[`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)

## Examples

``` r
# \donttest{
fit = bgm(x = Wenchuan[, 1:3])
#> 2 rows with missing values excluded (n = 360 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2000 (2.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 48/2000 (2.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 49/2000 (2.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 38/2000 (1.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 185/8000 (2.3%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 86/2000 (4.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 89/2000 (4.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 77/2000 (3.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 352/8000 (4.4%)
#> Elapsed: 1s | ETA: 22s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2000 (15.0%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 162/2000 (8.1%)
#> Chain 3 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 260/2000 (13.0%)
#> Chain 4 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 162/2000 (8.1%)
#> Total   (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 884/8000 (11.1%)
#> Elapsed: 3s | ETA: 24s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2000 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 251/2000 (12.6%)
#> Chain 3 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 324/2000 (16.2%)
#> Chain 4 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 252/2000 (12.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1327/8000 (16.6%)
#> Elapsed: 3s | ETA: 15s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 750/2000 (37.5%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 257/2000 (12.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 553/2000 (27.7%)
#> Chain 4 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 259/2000 (13.0%)
#> Total   (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1819/8000 (22.7%)
#> Elapsed: 4s | ETA: 14s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/2000 (42.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 578/2000 (28.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 808/2000 (40.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 713/2000 (35.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2949/8000 (36.9%)
#> Elapsed: 5s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1000/2000 (50.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 849/2000 (42.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 804/2000 (40.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3453/8000 (43.2%)
#> Elapsed: 6s | ETA: 8s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1200/2000 (60.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 1006/2000 (50.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1044/2000 (52.2%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 811/2000 (40.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 4061/8000 (50.8%)
#> Elapsed: 6s | ETA: 6s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1400/2000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1195/2000 (59.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1227/2000 (61.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 957/2000 (47.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4779/8000 (59.7%)
#> Elapsed: 7s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1600/2000 (80.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1379/2000 (69.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1414/2000 (70.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1170/2000 (58.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5563/8000 (69.5%)
#> Elapsed: 8s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/2000 (90.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1587/2000 (79.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1598/2000 (79.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1387/2000 (69.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6372/8000 (79.7%)
#> Elapsed: 8s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1783/2000 (89.1%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1781/2000 (89.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1580/2000 (79.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7144/8000 (89.3%)
#> Elapsed: 9s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 10s | ETA: 0s
extract_log_odds(fit)
#>           intrusion    dreams     flash
#> intrusion 0.0000000 0.6960187 0.4257605
#> dreams    0.6960187 0.0000000 0.5604572
#> flash     0.4257605 0.5604572 0.0000000
# }
```
