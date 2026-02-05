# Changelog

## bgms 0.1.6.3

### New features

- added
  [`simulate_mrf()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate_mrf.md)
  function for standalone MRF data simulation with user-specified
  parameters
- added
  [`simulate.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgms.md)
  S3 method to generate new observations from a fitted model using
  estimated parameters. Supports parallel processing via `cores`
  argument when using `method = "posterior-sample"` with optional
  progress bar
- added
  [`predict.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgms.md)
  S3 method to compute conditional probability distributions P(X_j \|
  X\_{-j}) for one or more variables given observed data
- both methods support using posterior mean parameters
  (`method = "posterior-mean"`) or averaging over posterior draws
  (`method = "posterior-sample"`) for full uncertainty propagation
- `baseline_category` is now stored in the fitted object’s arguments for
  use with Blume-Capel variables in simulation and prediction

### Bug fixes

- fixed mass matrix adaptation for NUTS/HMC: inverse mass matrix now
  correctly uses variance (following STAN convention) instead of
  precision, substantially improving sampling efficiency.
- fixed step size heuristic to re-run after each mass matrix update
  during warmup, ensuring step size is appropriate for the current mass
  matrix.
- fixed step size heuristic to resample momentum on each iteration
  (matching STAN’s init_stepsize behavior), improving step size
  selection.
- fixed energy diagnostic computation in NUTS to use actual accepted
  trajectory momentum instead of a random sample, making E-BFMI
  diagnostic meaningful.
- fixed Blume-Capel interaction term computation to use centered scores
  `(c - ref)` instead of raw category `c` in the pseudolikelihood
  denominator. This bug caused severe NUTS performance degradation (100%
  max tree depth hits) when using non-zero baseline categories.

### Deprecated

- [`mrfSampler()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/mrfSampler.md)
  is deprecated in favor of
  [`simulate_mrf()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate_mrf.md)
  for consistency with S3 method naming conventions

## bgms 0.1.6.2

CRAN release: 2026-01-20

### New features

- added option to separately specify beta priors for the within- and
  between-cluster probabilities for the SBM prior.

### Other changes

- reparameterized the Blume-capel model to use (score-baseline) instead
  of score.
- implemented a new way to compute the denominators and probabilities.
  This made their computation both faster and more stable.
- refactored c++ code for better maintainability.
- removed the prepared_data field from bgm objects.

### Bug fixes

- fixed numerical problems with Blume-Capel variables using HMC and
  NUTS.
- fixed a reporting bug where category thresholds for ordinal variables
  with a single category were incorrectly expanded to two parameters,
  resulting in spurious NA values.

## bgms 0.1.6.1

CRAN release: 2025-10-04

### Other changes

- added extractor function for joint SBM output
- cleaned up documentation, and c++ files
- changed length of warmup phase I in warmup scheduler HMC / NUTS (15% →
  7.5%)

### Bug fixes

- fixed a problem with warmup scheduling for adaptive-metropolis in
  bgmCompare()
- fixed stability problems with parallel sampling for bgm()
- fixed spurious output errors printing to console after user interrupt.

## bgms 0.1.6.0

CRAN release: 2025-09-27

### New features

- added NUTS and HMC options for sampling
  [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  and
  [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  models
- added support for running multiple chains in parallel
- added user interrupt handling for parallel sampling
- added Markov chain diagnostics (effective sample size and R-hat) for
  sampled parameters
- added [`summary()`](https://rdrr.io/r/base/summary.html),
  [`print()`](https://rdrr.io/r/base/print.html), and
  [`coef()`](https://rdrr.io/r/stats/coef.html) methods for fitted
  objects
- MCMC sampling in
  [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  and
  [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  is now reproducible when a `seed` argument is specified

### Other changes

- improved progress bar for parallel sampling
- [`summary()`](https://rdrr.io/r/base/summary.html) now integrates the
  functionality of the old `summary_SBM()`
- removed options for modeling main differences; main differences are
  now always estimated or selected, equivalent to the previous
  `main_difference_model = "collapse"` setting

### Bug fixes

- fixed an out-of-bounds error in
  [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  when handling missing data
- fixed a bug in the SBM prior computation

### Deprecated

- In
  [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md),
  the following arguments are deprecated:
  - `interaction_scale` → use `pairwise_scale`
  - `burnin` → use `warmup`
  - `save` → no longer needed (all outputs are returned by default)
  - `threshold_alpha`, `threshold_beta` → use `main_alpha`, `main_beta`
- In
  [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md),
  arguments related to difference models are deprecated:
  - `main_difference_model` (removed without replacement)
  - `reference_category` → use `baseline_category`
  - `pairwise_difference_*`, `main_difference_*` → use unified
    `difference_*` arguments
  - `pairwise_beta_bernoulli_*`, `main_beta_bernoulli_*` → use unified
    `beta_bernoulli_*` arguments
  - `interaction_scale` → use `pairwise_scale`
  - `threshold_alpha`, `threshold_beta` → use `main_alpha`, `main_beta`
  - `burnin` → use `warmup`
  - `save` → no longer needed
- Deprecated extractor functions:
  - [`extract_edge_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extractor_functions.md)
    → use
    [`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extractor_functions.md)
  - [`extract_pairwise_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extractor_functions.md)
    → use
    [`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extractor_functions.md)
- Deprecated object fields:
  - `$gamma` (pre-0.1.4) and `$indicator` (0.1.4–0.1.5) → replaced by
    `$raw_samples$indicator`
  - `$main_effects` (pre-0.1.4) and `$posterior_mean_main` (0.1.4–0.1.5)
    → replaced by `$raw_samples$main` (raw samples) and
    `$posterior_summary_main` (summaries)

## bgms 0.1.5.0 (GitHub only)

### New features

- The bgmCompare function now allows for network comparison for two or
  more groups.
- The new summary_sbm function can be used to summarize the output from
  the bgm function with the “Stochastic-Block” prior.
- Two new data sets are included in the package: ADHD and Boredom.

### Other changes

- The bgm function with the “Stochastic-Block” prior can now also return
  the sampled allocations and block probabilities, and sample and return
  the number of blocks.
- The underlying R and c++ functions received a massive update to
  improve their efficiency and maintainance.
- Repository moved to the Bayesian Graphical Modelling Lab organization.
- Included custom c++ implementations for exp and log on Windows.

### Bug fixes

- Fixed a bug in the bgmCompare function with selecting group
  differences of blume-capel parameters. Parameter differences that were
  not selected and should be fixed to zero were still updated.
- Fixed a bug in the bgmCompare function with handling the samples of
  blume-capel parameters. Output was not properly stored.
- Fixed a bug in the bgmCompare function with handling threshold
  estimation when missing categories and main_model = “Free”. The
  sufficient statistics and number of categories were not computed
  correctly.
- Partially fixed a bug in which the bgms package is slower on Windows
  than on Linux or MacOS. This is because the computation of exp and log
  using the gcc compiler for Windows is really slow. With a custom c++
  implementation, the speed is now closer to the speed achieved on Linux
  and MacOS.

## bgms 0.1.4.2

CRAN release: 2024-12-05

### Bug fixes

- fixed a bug with adjusting the variance of the proposal distributions
- fixed a bug with recoding data under the “collapse” condition

### Other changes

- when `selection = TRUE`, the burnin phase now runs `2 * burnin`
  iterations instead of `1 * burnin`. This ensures the chain starts with
  well-calibrated parameter values
- changed the maximum standard deviation of the adaptive proposal from
  20 back to 2

## bgms 0.1.4.1

CRAN release: 2024-11-12

This is a minor release that adds some documentation and output bug
fixes.

## bgms 0.1.4

CRAN release: 2024-10-20

### New features

- Comparing the category threshold and pairwise interaction parameters
  in two independent samples with bgmCompare().
- The Stochastic Block model is a new prior option for the network
  structure in bgm().

### Other changes

- Exported extractor functions to extract results from bgm objects in a
  safe way.
- Changed the maximum standard deviation of the adaptive proposal from 2
  to 20.
- Some small bug fixes.

## bgms 0.1.3

CRAN release: 2024-02-25

### New features

- Added support for Bayesian estimation without edge selection to bgm().
- Added support for simulating data from a (mixed) binary, ordinal, and
  Blume-Capel MRF to mrfSampler()
- Added support for analyzing (mixed) binary, ordinal, and Blume-Capel
  variables to bgm()

### User level changes

- Removed support of optimization based functions, mple(), mppe(), and
  bgm.em()
- Removed support for the Unit-Information prior from bgm()
- Removed support to do non-adaptive Metropolis from bgm()
- Reduced file size when saving raw MCMC samples

## bgms 0.1.2

CRAN release: 2023-10-13

This is a minor release that adds some bug fixes.

## bgms 0.1.1

CRAN release: 2023-09-01

This is a minor release adding some new features and fixing some minor
bugs.

### New features

- Missing data imputation for the bgm function. See the `na.action`
  option.
- Prior distributions for the network structure in the bgm function. See
  the `edge_prior` option.
- Adaptive Metropolis as an alternative to the current random walk
  Metropolis algorithm in the bgm function. See the `adaptive` option.

### User level changes

- Changed the default specification of the interaction prior from
  UnitInfo to Cauchy. See the `interaction_prior` option.
- Changed the default threshold hyperparameter specification from 1.0 to
  0.5. See the `threshold_alpha` and `threshold_beta` options.
- Analysis output now uses the column names of the data.
