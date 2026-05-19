# Package index

## Model fitting

Fit Bayesian graphical models and compare groups.

- [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  : Bayesian Estimation or Edge Selection for Markov Random Fields
- [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  : Bayesian Estimation and Variable Selection for Group Differences in
  Markov Random Fields

## Posterior methods

Inspect and summarize fitted models.

- [`print(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgms.md)
  :

  Print method for `bgms` objects

- [`summary(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)
  :

  Summary method for `bgms` objects

- [`coef(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgms.md)
  : Extract Coefficients from a bgms Object

- [`print(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgmCompare.md)
  :

  Print method for `bgmCompare` objects

- [`summary(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgmCompare.md)
  :

  Summary method for `bgmCompare` objects

- [`coef(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgmCompare.md)
  : Extract Coefficients from a bgmCompare Object

## Simulation and prediction

Generate observations or predict new ones. `simulate_*` functions draw
data given known parameters (forward simulation, exact draws);
`sample_*` functions run an MCMC sampler to draw from a target
distribution (e.g. the prior on the precision matrix).

- [`simulate_mrf()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate_mrf.md)
  : Simulate Observations from a Markov Random Field
- [`simulate(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgms.md)
  : Simulate Data from a Fitted bgms Model
- [`simulate(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgmCompare.md)
  : Simulate Data from a Fitted bgmCompare Model
- [`predict(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgms.md)
  : Predict Conditional Probabilities from a Fitted bgms Model
- [`predict(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgmCompare.md)
  : Predict Conditional Probabilities from a Fitted bgmCompare Model
- [`sample_ggm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sample_ggm_prior.md)
  : Sample from the GGM (Partial-Association) Prior

## Prior constructors

Specify priors for model parameters, precision scale, and edge
inclusion.

- [`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md)
  : Cauchy Prior for Model Parameters
- [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md)
  : Normal Prior for Model Parameters
- [`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md)
  : Beta-Prime Prior for Model Parameters
- [`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md)
  : Gamma Prior for Scale Parameters
- [`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md)
  : Exponential Prior for Scale Parameters
- [`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md)
  : Bernoulli Prior for Inclusion Indicators
- [`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md)
  : Beta-Bernoulli Prior for Inclusion Indicators
- [`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)
  : Stochastic Block Model Prior for Inclusion Indicators

## Extractors

Extract specific components from fitted model objects.

- [`extract_arguments()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.md)
  : Extract Model Arguments
- [`extract_category_thresholds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.md)
  **\[deprecated\]** : Extract Category Threshold Estimates
- [`extract_ess()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.md)
  : Extract Effective Sample Size
- [`extract_group_params()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.md)
  : Extract Group-Specific Parameters
- [`extract_indicator_priors()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.md)
  : Extract Indicator Prior Structure
- [`extract_indicators()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.md)
  : Extract Indicator Samples
- [`extract_log_odds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_log_odds.md)
  : Extract Posterior Mean Log-Odds
- [`extract_main_effects()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_main_effects.md)
  : Extract Main Effect Estimates
- [`extract_pairwise_interactions()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.md)
  : Extract Pairwise Interaction Samples
- [`extract_partial_correlations()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_partial_correlations.md)
  : Extract Posterior Mean Partial Correlations
- [`extract_posterior_inclusion_probabilities()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.md)
  : Extract Posterior Inclusion Probabilities
- [`extract_precision()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_precision.md)
  : Extract Posterior Mean Precision Matrix
- [`extract_rhat()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.md)
  : Extract R-hat Convergence Diagnostics
- [`extract_sbm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.md)
  : Extract Stochastic Block Model Summaries

## Legacy

Deprecated functions retained for backwards compatibility.

- [`mrfSampler()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/mrfSampler.md)
  **\[deprecated\]** : Sample observations from the ordinal MRF

## Datasets

- [`ADHD`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/ADHD.md)
  : ADHD Symptom Checklist for Children Aged 6–8 Years
- [`Boredom`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/Boredom.md)
  : Short Boredom Proneness Scale Responses
- [`Wenchuan`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/Wenchuan.md)
  : PTSD Symptoms in Wenchuan Earthquake Survivors Who Lost a Child

## Package

- [`bgms`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgms-package.md)
  [`bgms-package`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgms-package.md)
  : bgms: Bayesian Analysis of Graphical Models
