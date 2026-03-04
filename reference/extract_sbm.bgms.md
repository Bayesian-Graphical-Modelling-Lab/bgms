# Extract Stochastic Block Model Summaries from bgms Objects

Retrieves posterior summaries from a model fitted with \[bgm()\] using
the Stochastic Block prior on edge inclusion.

## Usage

``` r
# S3 method for class 'bgms'
extract_sbm(bgms_object)
```

## Arguments

- bgms_object:

  An object of class \`bgms\` returned by \[bgm()\].

## Value

A list with \`posterior_num_blocks\`, \`posterior_mean_allocations\`,
\`posterior_mode_allocations\`, and
\`posterior_mean_coclustering_matrix\`.

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
[`extract_posterior_inclusion_probabilities.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.bgms.md),
[`extract_rhat.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.bgmCompare.md),
[`extract_rhat.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.bgms.md)
