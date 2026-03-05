# Extract Indicator Samples

Retrieves posterior samples of inclusion indicators from a model fitted
with \[bgm()\] (edge inclusion indicators) or \[bgmCompare()\]
(difference indicators).

## Usage

``` r
extract_indicators(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class \`bgms\` (from \[bgm()\]) or
  \`bgmCompare\` (from \[bgmCompare()\]).

## Value

A matrix with one row per post-warmup iteration and one column per
indicator, containing binary (0/1) samples.

- bgms:

  One column per edge. Requires \`edge_selection = TRUE\`.

- bgmCompare:

  Columns for main-effect and pairwise difference indicators. Requires
  \`difference_selection = TRUE\`.

## See also

\[bgm()\], \[bgmCompare()\],
\[extract_posterior_inclusion_probabilities()\]

Other extractors:
[`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md),
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md),
[`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md),
[`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)
