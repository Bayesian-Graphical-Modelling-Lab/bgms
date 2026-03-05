# Extract Model Arguments

Retrieves the arguments used when fitting a model with \[bgm()\] or
\[bgmCompare()\].

## Usage

``` r
extract_arguments(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class \`bgms\` (from \[bgm()\]) or
  \`bgmCompare\` (from \[bgmCompare()\]).

## Value

A named list containing all arguments passed to the fitting function,
including data dimensions, prior settings, and MCMC configuration.

## See also

\[bgm()\], \[bgmCompare()\], \[summary.bgms()\],
\[summary.bgmCompare()\]

Other extractors:
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md),
[`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md),
[`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md),
[`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)
