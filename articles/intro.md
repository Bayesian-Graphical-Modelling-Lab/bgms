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
#> intrusion (1)  0.489 0.005 0.230 1945.784 1.000
#> intrusion (2) -1.876 0.012 0.340  821.458 1.001
#> intrusion (3) -4.787 0.022 0.555  625.757 1.000
#> intrusion (4) -9.421 0.032 0.888  759.944 1.000
#> dreams (1)    -0.600 0.003 0.191 3506.159 1.000
#> dreams (2)    -3.802 0.007 0.335 2173.957 1.001
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.316 0.001 0.032 3384.873            1.000
#> intrusion-flash   0.169 0.000 0.031 4167.302            1.001
#> intrusion-upset   0.094 0.003 0.040  186.359    242.303 1.013
#> intrusion-physior 0.096 0.003 0.038  282.208    209.412 1.027
#> dreams-flash      0.249 0.000 0.030 4971.520            1.000
#> dreams-upset      0.114 0.001 0.029  882.114    566.428 1.004
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  3999           
#> intrusion-flash   1.000       0.000     0     0     0  3999           
#> intrusion-upset   0.909 0.024 0.287   339    23    23  3614    144.777
#> intrusion-physior 0.927 0.024 0.259   275    15    15  3694    114.736
#> dreams-flash      1.000       0.000     0     0     0  3999           
#> dreams-upset      0.992 0.008 0.089    30     2     2  3965           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.021
#> intrusion-physior 1.077
#> dreams-flash           
#> dreams-upset      1.251
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
#> intrusion  0.4894592 -1.875703 -4.786752  -9.420814
#> dreams    -0.5996093 -3.801909 -7.133717 -11.578781
#> flash     -0.1029148 -2.567592 -5.372861  -9.683905
#> upset      0.4177050 -1.301367 -3.363430  -7.012984
#> physior   -0.6060231 -3.157062 -6.203702 -10.538650
#> 
#> $pairwise
#>            intrusion      dreams       flash       upset     physior
#> intrusion 0.00000000 0.315580608 0.169339357 0.093581361 0.096196683
#> dreams    0.31558061 0.000000000 0.249172935 0.113988604 0.005056459
#> flash     0.16933936 0.249172935 0.000000000 0.005299238 0.152998822
#> upset     0.09358136 0.113988604 0.005299238 0.000000000 0.355013101
#> physior   0.09619668 0.005056459 0.152998822 0.355013101 0.000000000
#> 
#> $indicator
#>           intrusion  dreams flash  upset physior
#> intrusion    0.0000 1.00000   1.0 0.9095 0.92750
#> dreams       1.0000 0.00000   1.0 0.9920 0.10175
#> flash        1.0000 1.00000   0.0 0.1000 1.00000
#> upset        0.9095 0.99200   0.1 0.0000 1.00000
#> physior      0.9275 0.10175   1.0 1.0000 0.00000
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
