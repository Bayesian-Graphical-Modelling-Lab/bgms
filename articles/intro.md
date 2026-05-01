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
#> intrusion (1)  0.481 0.005 0.232 2430.679 1.001
#> intrusion (2) -1.897 0.009 0.340 1427.735 1.001
#> intrusion (3) -4.834 0.017 0.557 1095.271 1.002
#> intrusion (4) -9.493 0.027 0.905 1109.597 1.001
#> dreams (1)    -0.589 0.003 0.194 3559.828 1.000
#> dreams (2)    -3.786 0.007 0.352 2485.547 1.000
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.316 0.001 0.034 3460.024            1.000
#> intrusion-flash   0.170 0.000 0.031 3848.193            1.001
#> intrusion-upset   0.101 0.002 0.035  559.970    325.012 1.027
#> intrusion-physior 0.091 0.003 0.039  173.147    155.484 1.001
#> dreams-flash      0.249 0.000 0.030 4348.861            1.001
#> dreams-upset      0.112 0.001 0.028 1776.515    776.053 1.006
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  3999           
#> intrusion-flash   1.000       0.000     0     0     0  3999           
#> intrusion-upset   0.960 0.017 0.197   151    10    11  3827     134.32
#> intrusion-physior 0.902 0.031 0.297   374    16    16  3593     93.033
#> dreams-flash      1.000       0.000     0     0     0  3999           
#> dreams-upset      0.993 0.007 0.083    26     2     2  3969           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.156
#> intrusion-physior 1.001
#> dreams-flash           
#> dreams-upset      1.297
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
#> intrusion  0.4807884 -1.897073 -4.833674  -9.493188
#> dreams    -0.5893051 -3.786156 -7.111724 -11.534887
#> flash     -0.1057963 -2.566706 -5.374301  -9.687611
#> upset      0.4067246 -1.327807 -3.409470  -7.090275
#> physior   -0.6018852 -3.139148 -6.163577 -10.474474
#> 
#> $pairwise
#>           intrusion      dreams       flash       upset     physior
#> intrusion 0.0000000 0.315856256 0.169889151 0.101432767 0.091231397
#> dreams    0.3158563 0.000000000 0.248946836 0.112403788 0.004325468
#> flash     0.1698892 0.248946836 0.000000000 0.003855482 0.154300291
#> upset     0.1014328 0.112403788 0.003855482 0.000000000 0.354919341
#> physior   0.0912314 0.004325468 0.154300291 0.354919341 0.000000000
#> 
#> $indicator
#>           intrusion dreams flash  upset physior
#> intrusion    0.0000  1.000 1.000 0.9595  0.9025
#> dreams       1.0000  0.000 1.000 0.9930  0.0890
#> flash        1.0000  1.000 0.000 0.0730  1.0000
#> upset        0.9595  0.993 0.073 0.0000  1.0000
#> physior      0.9025  0.089 1.000 1.0000  0.0000
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
