# Extract Posterior Inclusion Probabilities

Computes posterior inclusion probabilities from a model fitted with
\[bgm()\] (edge inclusion) or \[bgmCompare()\] (difference inclusion).

## Usage

``` r
extract_posterior_inclusion_probabilities(bgms_object)
```

## Arguments

- bgms_object:

  A fitted model object of class \`bgms\` (from \[bgm()\]) or
  \`bgmCompare\` (from \[bgmCompare()\]).

## Value

A symmetric p x p matrix of posterior inclusion probabilities, with
variable names as row and column names.

- bgms:

  Off-diagonal entries are edge inclusion probabilities. Requires
  \`edge_selection = TRUE\`.

- bgmCompare:

  Diagonal entries are main-effect inclusion probabilities; off-diagonal
  entries are pairwise difference inclusion probabilities. Requires
  \`difference_selection = TRUE\`.

## See also

\[bgm()\], \[bgmCompare()\], \[extract_indicators()\]

Other extractors:
[`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md),
[`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md),
[`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md),
[`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md),
[`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md),
[`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md),
[`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md),
[`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md),
[`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)
