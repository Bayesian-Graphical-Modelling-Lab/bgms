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
#> intrusion (1)  0.481 0.007 0.237 1164.310 1.001
#> intrusion (2) -1.895 0.013 0.342  642.933 1.003
#> intrusion (3) -4.830 0.024 0.565  565.276 1.004
#> intrusion (4) -9.496 0.039 0.918  563.758 1.005
#> dreams (1)    -0.583 0.007 0.198  756.158 1.026
#> dreams (2)    -3.759 0.015 0.369  626.240 1.039
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.316 0.001 0.034  989.114            1.000
#> intrusion-flash   0.167 0.001 0.032 1266.799            1.000
#> intrusion-upset   0.098 0.003 0.037  198.019     157.21 1.016
#> intrusion-physior 0.098 0.002 0.032  464.427    415.751 1.003
#> dreams-flash      0.251 0.001 0.031 1149.700            1.003
#> dreams-upset      0.104 0.014 0.042  577.799      8.936 1.265
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  1999           
#> intrusion-flash   1.000       0.000     0     0     0  1999           
#> intrusion-upset   0.946 0.024 0.226    99     9     9  1882     92.152
#> intrusion-physior 0.982 0.011 0.133    31     5     5  1958    152.199
#> dreams-flash      1.000       0.000     0     0     0  1999           
#> dreams-upset      0.905 0.121 0.293   189     1     1  1808           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset       1
#> intrusion-physior 1.196
#> dreams-flash           
#> dreams-upset      1.428
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
#> intrusion  0.4808374 -1.895005 -4.830296  -9.495580
#> dreams    -0.5828464 -3.758748 -7.065042 -11.465775
#> flash     -0.1003941 -2.559979 -5.368881  -9.674206
#> upset      0.4310768 -1.281712 -3.328146  -6.967497
#> physior   -0.6189765 -3.179397 -6.241250 -10.597811
#> 
#> $pairwise
#>            intrusion      dreams       flash       upset     physior
#> intrusion 0.00000000 0.316397896 0.166553964 0.097903603 0.098247571
#> dreams    0.31639790 0.000000000 0.251074831 0.104326848 0.006629126
#> flash     0.16655396 0.251074831 0.000000000 0.006808563 0.152099811
#> upset     0.09790360 0.104326848 0.006808563 0.000000000 0.355232997
#> physior   0.09824757 0.006629126 0.152099811 0.355232997 0.000000000
#> 
#> $indicator
#>           intrusion dreams flash upset physior
#> intrusion     0.000 1.0000 1.000 0.946  0.9820
#> dreams        1.000 0.0000 1.000 0.905  0.1125
#> flash         1.000 1.0000 0.000 0.119  1.0000
#> upset         0.946 0.9050 0.119 0.000  1.0000
#> physior       0.982 0.1125 1.000 1.000  0.0000
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
