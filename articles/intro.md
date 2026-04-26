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
#> intrusion (1)  0.491 0.009 0.230  613.795 1.000
#> intrusion (2) -1.872 0.018 0.354  401.556 1.001
#> intrusion (3) -4.782 0.030 0.583  366.925 1.002
#> intrusion (4) -9.416 0.046 0.926  405.375 1.002
#> dreams (1)    -0.598 0.005 0.191 1458.976 1.000
#> dreams (2)    -3.801 0.010 0.330 1004.464 1.002
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.315 0.001 0.032 1404.718            1.000
#> intrusion-flash   0.170 0.001 0.031 1263.303            1.004
#> intrusion-upset   0.090 0.005 0.042  106.701     79.281 1.022
#> intrusion-physior 0.100 0.003 0.034  275.329    183.842 1.004
#> dreams-flash      0.250 0.001 0.029 1932.660            1.000
#> dreams-upset      0.116 0.001 0.028  527.426    435.499 1.003
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  1999           
#> intrusion-flash   1.000       0.000     0     0     0  1999           
#> intrusion-upset   0.884 0.045 0.320   221    10    10  1758     50.174
#> intrusion-physior 0.962 0.023 0.191    71     5     5  1918     70.811
#> dreams-flash      1.000       0.000     0     0     0  1999           
#> dreams-upset      0.997 0.004 0.055     5     1     1  1992           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.052
#> intrusion-physior     1
#> dreams-flash           
#> dreams-upset      1.293
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
#> intrusion  0.4905223 -1.872491 -4.781571  -9.415589
#> dreams    -0.5976886 -3.801467 -7.137578 -11.583421
#> flash     -0.1070827 -2.578464 -5.394010  -9.719818
#> upset      0.4245807 -1.289478 -3.344967  -6.987112
#> physior   -0.6120216 -3.165189 -6.211773 -10.556802
#> 
#> $pairwise
#>            intrusion     dreams       flash       upset    physior
#> intrusion 0.00000000 0.31485357 0.169868976 0.090325114 0.10001514
#> dreams    0.31485357 0.00000000 0.249759272 0.115708387 0.00302638
#> flash     0.16986898 0.24975927 0.000000000 0.006345009 0.15244124
#> upset     0.09032511 0.11570839 0.006345009 0.000000000 0.35442935
#> physior   0.10001514 0.00302638 0.152441240 0.354429352 0.00000000
#> 
#> $indicator
#>           intrusion dreams  flash  upset physior
#> intrusion    0.0000  1.000 1.0000 0.8845   0.962
#> dreams       1.0000  0.000 1.0000 0.9970   0.067
#> flash        1.0000  1.000 0.0000 0.1245   1.000
#> upset        0.8845  0.997 0.1245 0.0000   1.000
#> physior      0.9620  0.067 1.0000 1.0000   0.000
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
