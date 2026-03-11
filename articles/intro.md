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
#> intrusion (1)  0.486 0.008 0.234  955.016 1.000
#> intrusion (2) -1.876 0.016 0.347  446.256 1.005
#> intrusion (3) -4.792 0.029 0.567  373.846 1.003
#> intrusion (4) -9.417 0.048 0.910  351.887 1.005
#> dreams (1)    -0.599 0.006 0.195 1082.710 1.001
#> dreams (2)    -3.786 0.013 0.363  800.508 1.001
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean    sd  mcse    n_eff  Rhat
#> intrusion-dreams  0.630 0.002 0.068 1093.158 1.002
#> intrusion-flash   0.341 0.002 0.064 1035.398 1.001
#> intrusion-upset   0.189 0.082 0.008   98.786 1.000
#> intrusion-physior 0.186 0.079 0.008   97.763 1.013
#> dreams-flash      0.497 0.002 0.062 1339.509 1.005
#> dreams-upset      0.228 0.002 0.056  645.291 1.002
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean    sd  mcse n0->0 n0->1 n1->0 n1->1  n_eff
#> intrusion-dreams  1.000 0.000           0     0     0  1999       
#> intrusion-flash   1.000 0.000           0     0     0  1999       
#> intrusion-upset   0.909 0.288 0.038   173     9     9  1808 55.925
#> intrusion-physior 0.905 0.293 0.038   180    10    10  1799 59.901
#> dreams-flash      1.000 0.000           0     0     0  1999       
#> dreams-upset      1.000 0.000           0     0     0  1999       
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.001
#> intrusion-physior 1.087
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
#> intrusion  0.4863416 -1.875767 -4.791863  -9.417115
#> dreams    -0.5986160 -3.786460 -7.110609 -11.538625
#> flash     -0.1129889 -2.590796 -5.406574  -9.738771
#> upset      0.4163331 -1.307319 -3.369964  -7.022134
#> physior   -0.5866610 -3.114764 -6.128465 -10.417837
#> 
#> $pairwise
#>           intrusion      dreams      flash      upset     physior
#> intrusion 0.0000000 0.630266467 0.34134255 0.18935058 0.186344985
#> dreams    0.6302665 0.000000000 0.49680046 0.22817892 0.007084831
#> flash     0.3413426 0.496800456 0.00000000 0.01389899 0.306203789
#> upset     0.1893506 0.228178918 0.01389899 0.00000000 0.705432395
#> physior   0.1863450 0.007084831 0.30620379 0.70543240 0.000000000
#> 
#> $indicator
#>           intrusion dreams  flash  upset physior
#> intrusion     0.000 1.0000 1.0000 0.9090  0.9050
#> dreams        1.000 0.0000 1.0000 1.0000  0.0715
#> flash         1.000 1.0000 0.0000 0.1205  1.0000
#> upset         0.909 1.0000 0.1205 0.0000  1.0000
#> physior       0.905 0.0715 1.0000 1.0000  0.0000
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
