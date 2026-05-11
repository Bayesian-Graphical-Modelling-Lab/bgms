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
#>                 mean  mcse    sd    n_eff Rhat
#> intrusion (1)  0.466 0.004 0.225 4115.848    1
#> intrusion (2) -1.925 0.006 0.310 2791.903    1
#> intrusion (3) -4.887 0.009 0.498 2909.743    1
#> intrusion (4) -9.579 0.016 0.819 2691.660    1
#> dreams (1)    -0.583 0.003 0.192 3503.796    1
#> dreams (2)    -3.773 0.006 0.330 3183.012    1
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean  mcse    sd    n_eff n_eff_mixt  Rhat
#> intrusion-dreams  0.315 0.000 0.032 4608.610            1.003
#> intrusion-flash   0.168 0.000 0.031 4080.983            1.000
#> intrusion-upset   0.101 0.001 0.030 2306.829    781.608 1.007
#> intrusion-physior 0.100 0.001 0.029 2356.718    619.274 1.015
#> dreams-flash      0.250 0.000 0.029 4875.430            1.000
#> dreams-upset      0.114 0.000 0.027 3267.706            1.003
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean  mcse    sd n0->0 n0->1 n1->0 n1->1 n_eff_mixt
#> intrusion-dreams  1.000       0.000     0     0     0  3999           
#> intrusion-flash   1.000       0.000     0     0     0  3999           
#> intrusion-upset   0.991 0.009 0.097    36     2     2  3959           
#> intrusion-physior 0.989 0.011 0.104    42     2     2  3953           
#> dreams-flash      1.000       0.000     0     0     0  3999           
#> dreams-upset      1.000       0.000     0     0     0  3999           
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset     1.3
#> intrusion-physior 1.301
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
#> intrusion  0.4659369 -1.924887 -4.887214  -9.578937
#> dreams    -0.5826513 -3.772680 -7.087175 -11.504794
#> flash     -0.1000993 -2.559686 -5.345971  -9.645600
#> upset      0.4072480 -1.326603 -3.398759  -7.077606
#> physior   -0.6170481 -3.166350 -6.217356 -10.557282
#> 
#> $pairwise
#>           intrusion     dreams       flash       upset    physior
#> intrusion 0.0000000 0.31459734 0.167979157 0.101166857 0.10040077
#> dreams    0.3145973 0.00000000 0.249691205 0.113681735 0.00142724
#> flash     0.1679792 0.24969121 0.000000000 0.003241342 0.15338880
#> upset     0.1011669 0.11368174 0.003241342 0.000000000 0.35421411
#> physior   0.1004008 0.00142724 0.153388799 0.354214108 0.00000000
#> 
#> $indicator
#>           intrusion  dreams   flash   upset physior
#> intrusion    0.0000 1.00000 1.00000 0.99050 0.98900
#> dreams       1.0000 0.00000 1.00000 1.00000 0.03475
#> flash        1.0000 1.00000 0.00000 0.06925 1.00000
#> upset        0.9905 1.00000 0.06925 0.00000 1.00000
#> physior      0.9890 0.03475 1.00000 1.00000 0.00000
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
