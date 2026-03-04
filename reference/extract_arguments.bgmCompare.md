# Extract Model Arguments from bgmCompare Objects

Retrieves the arguments used when fitting a model with \[bgmCompare()\].

## Usage

``` r
# S3 method for class 'bgmCompare'
extract_arguments(bgms_object)
```

## Arguments

- bgms_object:

  An object of class \`bgmCompare\` returned by \[bgmCompare()\].

## Value

A named list containing all arguments passed to \[bgmCompare()\],
including data dimensions, prior settings, and MCMC configuration.

## See also

\[bgmCompare()\], \[summary.bgmCompare()\], \[coef.bgmCompare()\]

Other extractors:
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
[`extract_rhat.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.bgms.md),
[`extract_sbm.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.bgms.md)
