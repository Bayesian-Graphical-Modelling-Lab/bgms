# Getting Started with bgms

## Introduction

The **bgms** package implements Bayesian methods for analyzing graphical
models of binary and ordinal variables. It estimates main effects
(category thresholds) and pairwise interactions in an ordinal Markov
random field (MRF), with optional Bayesian edge selection via
spike–and–slab priors. The package provides two main entry points:

- [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  for one-sample designs (single network),
- [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  for independent-sample designs (group comparisons).

This vignette walks through the basic workflow: fitting a model,
summarizing posterior output, and visualizing results.

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
#> intrusion (1)  0.475 0.007 0.223 1157.788 1.001
#> intrusion (2) -1.906 0.011 0.315  880.403 1.005
#> intrusion (3) -4.846 0.017 0.514  875.204 1.007
#> intrusion (4) -9.519 0.028 0.845  904.952 1.008
#> dreams (1)    -0.588 0.006 0.192 1159.841 1.001
#> dreams (2)    -3.782 0.012 0.353  871.629 1.000
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean    sd  mcse    n_eff  Rhat
#> intrusion-dreams  0.629 0.002 0.067 1382.782 1.005
#> intrusion-flash   0.339 0.002 0.060 1444.285 1.002
#> intrusion-upset   0.202 0.064 0.005  187.988 1.003
#> intrusion-physior 0.189 0.071 0.006  137.040 1.000
#> dreams-flash      0.500 0.002 0.060 1124.700 1.002
#> dreams-upset      0.226 0.053 0.002 1128.334 1.000
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean    sd  mcse n0->0 n0->1 n1->0 n1->1   n_eff
#> intrusion-dreams  1.000 0.000           0     0     0  1999        
#> intrusion-flash   1.000 0.000           0     0     0  1999        
#> intrusion-upset   0.979 0.143  0.02    40     2     2  1955  49.853
#> intrusion-physior 0.941 0.237 0.029   112     7     7  1873  64.566
#> dreams-flash      1.000 0.000           0     0     0  1999        
#> dreams-upset      0.999 0.039 0.002     2     1     1  1995 400.722
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.311
#> intrusion-physior 1.016
#> dreams-flash           
#> dreams-upset      1.292
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
#> intrusion  0.4752393 -1.905757 -4.846385  -9.519480
#> dreams    -0.5884216 -3.781525 -7.109964 -11.530285
#> flash     -0.1051154 -2.566799 -5.365407  -9.674749
#> upset      0.4216365 -1.308488 -3.370565  -7.037020
#> physior   -0.6037679 -3.140033 -6.167547 -10.489248
#> 
#> $pairwise
#>           intrusion      dreams       flash       upset     physior
#> intrusion 0.0000000 0.629338089 0.338657791 0.202376428 0.188939686
#> dreams    0.6293381 0.000000000 0.499816233 0.225595383 0.007497679
#> flash     0.3386578 0.499816233 0.000000000 0.004566734 0.308420955
#> upset     0.2023764 0.225595383 0.004566734 0.000000000 0.707339037
#> physior   0.1889397 0.007497679 0.308420955 0.707339037 0.000000000
#> 
#> $indicator
#>           intrusion dreams  flash  upset physior
#> intrusion    0.0000 1.0000 1.0000 0.9790  0.9405
#> dreams       1.0000 0.0000 1.0000 0.9985  0.0820
#> flash        1.0000 1.0000 0.0000 0.0545  1.0000
#> upset        0.9790 0.9985 0.0545 0.0000  1.0000
#> physior      0.9405 0.0820 1.0000 1.0000  0.0000
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

## Next steps

- For comparing groups, see
  [`?bgmCompare`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  or the *Model Comparison* vignette.
- For diagnostics and convergence checks, see the *Diagnostics*
  vignette.
