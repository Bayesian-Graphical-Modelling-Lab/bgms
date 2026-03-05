# Extract Indicator Prior Structure

Retrieves the prior specification used for inclusion indicators in a
model fitted with \[bgm()\] (edge indicators) or \[bgmCompare()\]
(difference indicators).

## Usage

``` r
extract_indicator_priors(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class \`bgms\` (from \[bgm()\]) or
  \`bgmCompare\` (from \[bgmCompare()\]).

## Value

A named list describing the prior structure, including the prior type
and any hyperparameters.

- bgms:

  Requires \`edge_selection = TRUE\`. Returns a list with the prior type
  (\`"Bernoulli"\`, \`"Beta-Bernoulli"\`, or \`"Stochastic-Block"\`) and
  associated hyperparameters.

- bgmCompare:

  Requires \`difference_selection = TRUE\`. Returns the difference prior
  specification.

## See also

\[bgm()\], \[bgmCompare()\], \[extract_indicators()\]

Other extractors:
[`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md),
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md),
[`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md),
[`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)
