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
#> intrusion (1)  0.484 0.006 0.227 1257.007 1.008
#> intrusion (2) -1.882 0.012 0.332  791.224 1.013
#> intrusion (3) -4.809 0.022 0.544  613.411 1.012
#> intrusion (4) -9.449 0.035 0.876  612.156 1.015
#> dreams (1)    -0.585 0.005 0.188 1508.300 1.004
#> dreams (2)    -3.767 0.009 0.327 1306.396 1.004
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.315 0.001 0.032 2024.905            1.008
#> intrusion-flash   0.168 0.001 0.031 1884.688            1.001
#> intrusion-upset   0.098 0.003 0.036  249.474    188.357 1.001
#> intrusion-physior 0.095 0.004 0.038  251.613    118.044 1.050
#> dreams-flash      0.249 0.001 0.030 1666.601            1.000
#> dreams-upset      0.113 0.001 0.027 1462.892            1.001
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  1999           
#> intrusion-flash   1.000       0.000     0     0     0  1999           
#> intrusion-upset   0.953 0.023 0.212    87     7     7  1898      81.32
#> intrusion-physior 0.922 0.033 0.268   147     9     9  1834     64.597
#> dreams-flash      1.000       0.000     0     0     0  1999           
#> dreams-upset      1.000       0.000     0     0     0  1999           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset       1
#> intrusion-physior 1.157
#> dreams-flash           
#> dreams-upset           
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
#> intrusion  0.4835406 -1.882415 -4.808577  -9.449140
#> dreams    -0.5849534 -3.767083 -7.079554 -11.498049
#> flash     -0.1043367 -2.569951 -5.374069  -9.695435
#> upset      0.4191383 -1.313351 -3.388606  -7.073257
#> physior   -0.6116235 -3.159539 -6.193801 -10.527694
#> 
#> $pairwise
#>            intrusion      dreams       flash       upset     physior
#> intrusion 0.00000000 0.315073338 0.168341600 0.097807766 0.095269549
#> dreams    0.31507334 0.000000000 0.248862725 0.113270460 0.002065929
#> flash     0.16834160 0.248862725 0.000000000 0.006060051 0.154306636
#> upset     0.09780777 0.113270460 0.006060051 0.000000000 0.355636928
#> physior   0.09526955 0.002065929 0.154306636 0.355636928 0.000000000
#> 
#> $indicator
#>           intrusion dreams  flash  upset physior
#> intrusion     0.000 1.0000 1.0000 0.9530  0.9220
#> dreams        1.000 0.0000 1.0000 1.0000  0.0545
#> flash         1.000 1.0000 0.0000 0.1195  1.0000
#> upset         0.953 1.0000 0.1195 0.0000  1.0000
#> physior       0.922 0.0545 1.0000 1.0000  0.0000
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

The pairwise effects are partial correlations (off-diagonal entries of
the standardized precision matrix). Missing values can be imputed during
sampling with `na_action = "impute"`.

## Next steps

- For comparing groups, see
  [`?bgmCompare`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  or the *Model Comparison* vignette.
- For diagnostics and convergence checks, see the *Diagnostics*
  vignette.
