# Extract R-hat Convergence Diagnostics

Retrieves R-hat convergence diagnostics for all parameters from a model
fitted with \[bgm()\] or \[bgmCompare()\].

## Usage

``` r
extract_rhat(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class \`bgms\` (from \[bgm()\]) or
  \`bgmCompare\` (from \[bgmCompare()\]).

## Value

A named list with R-hat values for each parameter type present in the
model (e.g., \`main\`, \`pairwise\`, \`indicator\`).

## See also

\[bgm()\], \[bgmCompare()\], \[extract_ess()\]

Other extractors:
[`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md),
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md),
[`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md),
[`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)
