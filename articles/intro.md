# Getting Started with bgms

## Introduction

The **bgms** package implements Bayesian methods for analyzing graphical
models. It supports three variable types:

- **ordinal** (including binary) — Markov random field (MRF) models,
- **Blume–Capel** — ordinal MRF with a reference category,
- **continuous** — Gaussian graphical models (GGM).

The package estimates main effects and pairwise interactions, with
optional Bayesian edge selection via spike-and-slab priors. It provides
two main entry points:

- [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  for one-sample designs (single network),
- [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  for independent-sample designs (group comparisons).

This vignette walks through the basic workflow for ordinal data. For
continuous data, set `variable_type = "continuous"` in
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
to fit a Gaussian graphical model.

## Wenchuan dataset

The dataset `Wenchuan` contains responses from survivors of the 2008
Wenchuan earthquake on posttraumatic stress items. Here, we analyze a
subset of the first five items as a demonstration.

``` r

library(bgms)

# Analyse a subset of the Wenchuan dataset
?Wenchuan
data = Wenchuan[, 1:5]
head(data)
#>      intrusion dreams flash upset physior
#> [1,]         2      2     2     2       3
#> [2,]         2      2     2     3       3
#> [3,]         2      4     4     4       3
#> [4,]         2      1     2     2       1
#> [5,]         2      2     2     2       2
#> [6,]         4      3     2     2       2
```

## Fitting a model

The main entry point is
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
for single-group models and
[`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
for multiple-group comparisons.

``` r

fit = bgm(data, seed = 1234)
```

## Posterior summaries

``` r

summary(fit)
#> Posterior summaries from Bayesian estimation:
#> 
#> Category thresholds: 
#>                 mean  mcse    sd    n_eff  Rhat
#> intrusion (2)  0.476 0.004 0.229 3220.312 1.000
#> intrusion (3) -1.905 0.009 0.331 1216.778 1.000
#> intrusion (4) -4.850 0.015 0.541 1294.639 1.002
#> intrusion (5) -9.519 0.026 0.873 1104.327 1.003
#> dreams (2)    -0.592 0.003 0.190 3504.797 1.000
#> dreams (3)    -3.789 0.007 0.342 2270.566 1.002
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.315 0.000 0.032 4340.561            1.002
#> intrusion-flash   0.170 0.000 0.031 4139.371            1.001
#> intrusion-upset   0.095 0.002 0.037  276.610    245.619 1.001
#> intrusion-physior 0.101 0.001 0.031  619.137    596.378 1.000
#> dreams-flash      0.250 0.000 0.031 4700.667            1.000
#> dreams-upset      0.115 0.001 0.028 2135.159   1854.237 1.000
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  3999           
#> intrusion-flash   1.000       0.000     0     0     0  3999           
#> intrusion-upset   0.935 0.022 0.246   244    14    14  3727    119.477
#> intrusion-physior 0.978 0.011 0.146    80     7     7  3905    171.553
#> dreams-flash      1.000       0.000     0     0     0  3999           
#> dreams-upset      0.999 0.001 0.032     2     2     2  3993           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.008
#> intrusion-physior 1.003
#> dreams-flash           
#> dreams-upset      1.292
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;
#> `summary(fit)$indicator` still contains all computed values.
#> 
#> Use `summary(fit)$<component>` to access full results.
#> Use `extract_log_odds(fit)` for log odds ratios.
#> See the `easybgm` package for other summary and plotting tools.
```

You can also access posterior means or inclusion probabilities directly:

``` r

coef(fit)
#> $main
#>              cat (1)   cat (2)   cat (3)    cat (4)
#> intrusion  0.4761194 -1.905055 -4.849646  -9.519353
#> dreams    -0.5923305 -3.788923 -7.115211 -11.545634
#> flash     -0.1017229 -2.566101 -5.371361  -9.679987
#> upset      0.4158459 -1.305727 -3.372244  -7.035173
#> physior   -0.6128987 -3.169453 -6.222240 -10.569988
#> 
#> $pairwise
#>           intrusion      dreams       flash       upset     physior
#> intrusion 0.0000000 0.314801118 0.169977835 0.095469103 0.100587283
#> dreams    0.3148011 0.000000000 0.249816272 0.114758607 0.002218409
#> flash     0.1699778 0.249816272 0.000000000 0.004178344 0.152727316
#> upset     0.0954691 0.114758607 0.004178344 0.000000000 0.355014255
#> physior   0.1005873 0.002218409 0.152727316 0.355014255 0.000000000
#> 
#> $indicator
#>           intrusion  dreams   flash   upset physior
#> intrusion   0.00000 1.00000 1.00000 0.93550 0.97825
#> dreams      1.00000 0.00000 1.00000 0.99900 0.05525
#> flash       1.00000 1.00000 0.00000 0.08725 1.00000
#> upset       0.93550 0.99900 0.08725 0.00000 1.00000
#> physior     0.97825 0.05525 1.00000 1.00000 0.00000
```

## Network plot

To visualize the network structure, we threshold the posterior inclusion
probabilities at 0.5 and plot the resulting adjacency matrix.

``` r

library(qgraph)

median_probability_network = coef(fit)$pairwise
median_probability_network[coef(fit)$indicator < 0.5] = 0.0

qgraph(median_probability_network,
  theme = "TeamFortress",
  maximum = 1,
  fade = FALSE,
  color = c("#f0ae0e"), vsize = 10, repulsion = .9,
  label.cex = 1, label.scale = "FALSE",
  labels = colnames(data)
)
```

![](intro_files/figure-html/unnamed-chunk-7-1.png)

## Continuous data (GGM)

For continuous variables,
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
fits a Gaussian graphical model when `variable_type = "continuous"`. The
workflow is the same:

``` r

fit_ggm = bgm(continuous_data, variable_type = "continuous", seed = 1234)
summary(fit_ggm)
```

The pairwise effects are unstandardized partial associations derived
from the off-diagonal entries of the precision matrix (stored as
`-0.5 * K_ij`). They are not standardized partial correlations. Use
[`extract_partial_correlations()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_partial_correlations.md)
to convert them to actual partial correlations, or
[`extract_precision()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_precision.md)
to obtain the precision matrix. Missing values can be imputed during
sampling with `na_action = "impute"`.

## Next steps

- For comparing groups, see
  [`?bgmCompare`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  or the *Model Comparison* vignette.
- For diagnostics and convergence checks, see the *Diagnostics*
  vignette.
