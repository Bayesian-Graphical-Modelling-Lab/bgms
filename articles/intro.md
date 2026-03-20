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
#> intrusion (1)  0.505 0.008 0.238 949.566 1.007
#> intrusion (2) -1.834 0.017 0.363 468.832 1.041
#> intrusion (3) -4.710 0.032 0.599 342.754 1.051
#> intrusion (4) -9.302 0.052 0.966 342.180 1.048
#> dreams (1)    -0.596 0.006 0.194 926.562 1.002
#> dreams (2)    -3.794 0.013 0.346 735.954 1.010
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean    sd  mcse    n_eff  Rhat
#> intrusion-dreams  0.315 0.001 0.034 1049.013 1.003
#> intrusion-flash   0.171 0.001 0.033 1332.705 1.000
#> intrusion-upset   0.088 0.046 0.007   42.611 1.102
#> intrusion-physior 0.093 0.044 0.005   75.344 1.018
#> dreams-flash      0.249 0.001 0.030 1587.148 1.003
#> dreams-upset      0.118 0.001 0.029  891.284 1.044
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean    sd  mcse n0->0 n0->1 n1->0 n1->1  n_eff
#> intrusion-dreams  1.000 0.000           0     0     0  1999       
#> intrusion-flash   1.000 0.000           0     0     0  1999       
#> intrusion-upset   0.840 0.367 0.067   313     8     8  1670 30.137
#> intrusion-physior 0.877 0.328 0.045   234    11    11  1743 52.513
#> dreams-flash      1.000 0.000           0     0     0  1999       
#> dreams-upset      1.000 0.000           0     0     0  1999       
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.459
#> intrusion-physior 1.008
#> dreams-flash           
#> dreams-upset           
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant (all 0 or all 1) across all iterations, so sd/mcse/n_eff/Rhat
#> are undefined; `summary(fit)$indicator` still contains the NA values.
#> 
#> Use `summary(fit)$<component>` to access full results.
#> See the `easybgm` package for other summary and plotting tools.
```

You can also access posterior means or inclusion probabilities directly:

``` r
coef(fit)
#> $main
#>              cat (1)   cat (2)   cat (3)    cat (4)
#> intrusion  0.5049730 -1.834465 -4.709887  -9.301568
#> dreams    -0.5960347 -3.793921 -7.134361 -11.576948
#> flash     -0.1109266 -2.589919 -5.407930  -9.743585
#> upset      0.4277309 -1.282217 -3.335943  -6.987220
#> physior   -0.6060884 -3.147660 -6.183574 -10.507612
#> 
#> $pairwise
#>            intrusion      dreams       flash       upset     physior
#> intrusion 0.00000000 0.315425796 0.171408964 0.088392559 0.092684081
#> dreams    0.31542580 0.000000000 0.249067188 0.117866731 0.001358329
#> flash     0.17140896 0.249067188 0.000000000 0.003699253 0.156461423
#> upset     0.08839256 0.117866731 0.003699253 0.000000000 0.356633286
#> physior   0.09268408 0.001358329 0.156461423 0.356633286 0.000000000
#> 
#> $indicator
#>           intrusion dreams  flash  upset physior
#> intrusion    0.0000 1.0000 1.0000 0.8395  0.8775
#> dreams       1.0000 0.0000 1.0000 1.0000  0.0315
#> flash        1.0000 1.0000 0.0000 0.0735  1.0000
#> upset        0.8395 1.0000 0.0735 0.0000  1.0000
#> physior      0.8775 0.0315 1.0000 1.0000  0.0000
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
