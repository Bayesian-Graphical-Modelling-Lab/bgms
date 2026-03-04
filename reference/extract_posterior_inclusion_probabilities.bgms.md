# Extract Posterior Inclusion Probabilities from bgms Objects

Computes the posterior probability of edge inclusion from a model fitted
with \[bgm()\] using edge selection.

## Usage

``` r
# S3 method for class 'bgms'
extract_posterior_inclusion_probabilities(bgms_object)
```

## Arguments

- bgms_object:

  An object of class \`bgms\` returned by \[bgm()\].

## Value

A symmetric p x p matrix of posterior inclusion probabilities, with
variable names as row and column names.

## See also

\[bgm()\], \[summary.bgms()\], \[coef.bgms()\]

Other extractors:
[`extract_arguments.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.bgmCompare.md),
[`extract_arguments.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.bgms.md),
[`extract_category_thresholds.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.bgmCompare.md),
[`extract_category_thresholds.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.bgms.md),
[`extract_ess.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.bgmCompare.md),
[`extract_ess.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.bgms.md),
[`extract_group_params.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.bgmCompare.md),
[`extract_indicator_priors.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.bgmCompare.md),
[`extract_indicator_priors.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.bgms.md),
[`extract_indicators.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.bgmCompare.md),
[`extract_indicators.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.bgms.md),
[`extract_pairwise_interactions.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.bgmCompare.md),
[`extract_pairwise_interactions.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.bgms.md),
[`extract_posterior_inclusion_probabilities.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.bgmCompare.md),
[`extract_rhat.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.bgmCompare.md),
[`extract_rhat.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.bgms.md),
[`extract_sbm.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.bgms.md)
