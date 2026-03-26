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
#>                 mean  mcse    sd   n_eff  Rhat
#> intrusion (1)  0.481 0.008 0.244 849.010 1.002
#> intrusion (2) -1.887 0.013 0.353 685.881 1.005
#> intrusion (3) -4.812 0.024 0.572 560.736 1.006
#> intrusion (4) -9.467 0.037 0.906 605.623 1.005
#> dreams (1)    -0.599 0.006 0.193 896.581 1.001
#> dreams (2)    -3.796 0.014 0.361 697.001 1.004
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.314 0.001 0.034 1165.909            1.000
#> intrusion-flash   0.169 0.001 0.032 1433.416            1.001
#> intrusion-upset   0.097 0.003 0.037  182.327    141.315 1.036
#> intrusion-physior 0.096 0.003 0.037  152.058    122.914 1.009
#> dreams-flash      0.249 0.001 0.031 1438.002            1.000
#> dreams-upset      0.114 0.002 0.030  666.101     223.74 1.022
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  1999           
#> intrusion-flash   1.000       0.000     0     0     0  1999           
#> intrusion-upset   0.944 0.028 0.231   106     7     7  1879     67.887
#> intrusion-physior 0.935 0.032 0.246   122     7     7  1863      59.74
#> dreams-flash      1.000       0.000     0     0     0  1999           
#> dreams-upset      0.989 0.015 0.104    21     1     1  1976           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset    1.14
#> intrusion-physior 1.033
#> dreams-flash           
#> dreams-upset      1.301
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
#> intrusion  0.4807844 -1.887276 -4.812235  -9.466882
#> dreams    -0.5990120 -3.796343 -7.120031 -11.554333
#> flash     -0.1010139 -2.563203 -5.356205  -9.657865
#> upset      0.4130636 -1.311834 -3.387882  -7.056063
#> physior   -0.6078418 -3.157839 -6.199070 -10.534263
#> 
#> $pairwise
#>            intrusion      dreams       flash       upset     physior
#> intrusion 0.00000000 0.314016539 0.169250494 0.097478982 0.095610613
#> dreams    0.31401654 0.000000000 0.249094590 0.113843953 0.004991406
#> flash     0.16925049 0.249094590 0.000000000 0.004239581 0.152961036
#> upset     0.09747898 0.113843953 0.004239581 0.000000000 0.355081270
#> physior   0.09561061 0.004991406 0.152961036 0.355081270 0.000000000
#> 
#> $indicator
#>           intrusion dreams  flash  upset physior
#> intrusion    0.0000  1.000 1.0000 0.9435  0.9355
#> dreams       1.0000  0.000 1.0000 0.9890  0.1050
#> flash        1.0000  1.000 0.0000 0.0835  1.0000
#> upset        0.9435  0.989 0.0835 0.0000  1.0000
#> physior      0.9355  0.105 1.0000 1.0000  0.0000
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
