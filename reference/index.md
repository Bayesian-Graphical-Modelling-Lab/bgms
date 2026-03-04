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
  : Print method for \`bgms\` objects
- [`summary(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)
  : Summary method for \`bgms\` objects
- [`coef(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgms.md)
  : Extract Coefficients from a bgms Object
- [`print(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgmCompare.md)
  : Print method for \`bgmCompare\` objects
- [`summary(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgmCompare.md)
  : Summary method for \`bgmCompare\` objects
- [`coef(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgmCompare.md)
  : Extract Coefficients from a bgmCompare Object

## Simulation and prediction

Simulate from Markov random fields and predict new observations.

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

## Extractors

Extract specific components from fitted model objects.

- [`extract_arguments(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.bgmCompare.md)
  : Extract Model Arguments from bgmCompare Objects
- [`extract_arguments(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_arguments.bgms.md)
  : Extract Model Arguments from bgms Objects
- [`extract_category_thresholds(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.bgmCompare.md)
  : Extract Category Threshold Samples from bgmCompare Objects
- [`extract_category_thresholds(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_category_thresholds.bgms.md)
  : Extract Category Threshold Estimates from bgms Objects
- [`extract_ess(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.bgmCompare.md)
  : Extract Effective Sample Size from bgmCompare Objects
- [`extract_ess(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_ess.bgms.md)
  : Extract Effective Sample Size from bgms Objects
- [`extract_group_params(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_group_params.bgmCompare.md)
  : Extract Group-Specific Parameters from bgmCompare Objects
- [`extract_indicator_priors(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.bgmCompare.md)
  : Extract Indicator Prior Structure from bgmCompare Objects
- [`extract_indicator_priors(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicator_priors.bgms.md)
  : Extract Indicator Prior Structure from bgms Objects
- [`extract_indicators(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.bgmCompare.md)
  : Extract Difference Indicator Samples from bgmCompare Objects
- [`extract_indicators(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_indicators.bgms.md)
  : Extract Edge Indicator Samples from bgms Objects
- [`extract_pairwise_interactions(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.bgmCompare.md)
  : Extract Pairwise Interaction Samples from bgmCompare Objects
- [`extract_pairwise_interactions(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_pairwise_interactions.bgms.md)
  : Extract Pairwise Interaction Samples from bgms Objects
- [`extract_posterior_inclusion_probabilities(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.bgmCompare.md)
  : Extract Posterior Inclusion Probabilities from bgmCompare Objects
- [`extract_posterior_inclusion_probabilities(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_posterior_inclusion_probabilities.bgms.md)
  : Extract Posterior Inclusion Probabilities from bgms Objects
- [`extract_rhat(`*`<bgmCompare>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.bgmCompare.md)
  : Extract R-hat Diagnostics from bgmCompare Objects
- [`extract_rhat(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_rhat.bgms.md)
  : Extract R-hat Diagnostics from bgms Objects
- [`extract_sbm(`*`<bgms>`*`)`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_sbm.bgms.md)
  : Extract Stochastic Block Model Summaries from bgms Objects

## Legacy

Deprecated functions retained for backwards compatibility.

- [`mrfSampler()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/mrfSampler.md)
  : Sample observations from the ordinal MRF

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
  : bgms: Bayesian Analysis of Networks of Binary and/or Ordinal
  Variables
