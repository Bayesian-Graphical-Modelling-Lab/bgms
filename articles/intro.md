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
#> intrusion (1)  0.482 0.008 0.247 938.346 1.002
#> intrusion (2) -1.890 0.015 0.357 579.689 1.003
#> intrusion (3) -4.819 0.026 0.573 481.272 1.013
#> intrusion (4) -9.481 0.041 0.915 503.153 1.013
#> dreams (1)    -0.600 0.006 0.192 899.725 1.007
#> dreams (2)    -3.811 0.013 0.354 740.361 1.009
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.316 0.001 0.034  999.191            1.002
#> intrusion-flash   0.168 0.001 0.032 1053.382            1.001
#> intrusion-upset   0.094 0.004 0.040  111.831     93.457 1.016
#> intrusion-physior 0.099 0.003 0.036  559.379    152.522 1.054
#> dreams-flash      0.251 0.001 0.032 1294.255            1.000
#> dreams-upset      0.116 0.001 0.028  800.878            1.002
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  1999           
#> intrusion-flash   1.000       0.000     0     0     0  1999           
#> intrusion-upset   0.915 0.038 0.279   162     8     8  1821      52.79
#> intrusion-physior 0.956 0.026 0.205    83     5     4  1907     60.699
#> dreams-flash      1.000       0.000     0     0     0  1999           
#> dreams-upset      1.000       0.000     0     0     0  1999           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.077
#> intrusion-physior 1.318
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
#> intrusion  0.4816437 -1.890120 -4.818863  -9.480881
#> dreams    -0.6002222 -3.810807 -7.165686 -11.629661
#> flash     -0.1029018 -2.558923 -5.358136  -9.659810
#> upset      0.4295512 -1.291327 -3.349722  -7.002167
#> physior   -0.6147936 -3.172239 -6.217856 -10.560732
#> 
#> $pairwise
#>            intrusion      dreams       flash       upset     physior
#> intrusion 0.00000000 0.316103177 0.168396751 0.094425256 0.099384945
#> dreams    0.31610318 0.000000000 0.251241439 0.115638791 0.003917364
#> flash     0.16839675 0.251241439 0.000000000 0.003606312 0.152476833
#> upset     0.09442526 0.115638791 0.003606312 0.000000000 0.354029027
#> physior   0.09938494 0.003917364 0.152476833 0.354029027 0.000000000
#> 
#> $indicator
#>           intrusion dreams flash upset physior
#> intrusion     0.000  1.000 1.000 0.915   0.956
#> dreams        1.000  0.000 1.000 1.000   0.082
#> flash         1.000  1.000 0.000 0.077   1.000
#> upset         0.915  1.000 0.077 0.000   1.000
#> physior       0.956  0.082 1.000 1.000   0.000
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
