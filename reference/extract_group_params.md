# Extract Group-Specific Parameters

Computes group-specific parameter estimates by combining baseline
parameters and group differences from a model fitted with
\[bgmCompare()\].

## Usage

``` r
extract_group_params(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class \`bgmCompare\` (from \[bgmCompare()\]).

## Value

A list with elements \`main_effects_groups\` (main effects per group)
and \`pairwise_effects_groups\` (pairwise effects per group).

## See also

\[bgmCompare()\], \[extract_pairwise_interactions()\],
\[extract_category_thresholds()\]

Other extractors:
[`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md),
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md),
[`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md),
[`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)
