# Extract Posterior Mean Precision Matrix

Retrieves the posterior mean precision matrix from a model fitted with
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md).
For GGM models this is the full precision matrix. For mixed MRF models
this is the precision matrix of the continuous (Gaussian) block. OMRF
models have no precision matrix and return `NULL`.

For mixed MRF models the precision matrix is reconstructed from the
internal association-scale parameterization.

## Usage

``` r
extract_precision(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class `bgms` (from
  [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)).

## Value

A named numeric matrix containing the posterior mean precision matrix,
or `NULL` for OMRF models.

- GGM:

  A symmetric matrix with one row and column per variable.

- Mixed MRF:

  A symmetric matrix with one row and column per continuous variable.

- OMRF:

  `NULL` (invisibly).

## Details

Extract Posterior Mean Precision Matrix

## See also

[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md),
[`coef.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgms.md),
[`extract_partial_correlations()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_partial_correlations.md)

Other extractors:
[`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md),
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md),
[`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md),
[`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md),
[`extract_log_odds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_log_odds.md),
[`extract_main_effects()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_main_effects.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_partial_correlations()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_partial_correlations.md),
[`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md),
[`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)

## Examples

``` r
# \donttest{
fit = bgm(
  x = Wenchuan[, 1:3],
  variable_type = rep("continuous", 3)
)
#> 2 rows with missing values excluded (n = 360 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 0s | ETA: 0s
extract_precision(fit)
#>            intrusion     dreams      flash
#> intrusion  1.7031353 -0.8230416 -0.4882264
#> dreams    -0.8230416  1.6166990 -0.6701382
#> flash     -0.4882264 -0.6701382  1.5437290
# }
```
