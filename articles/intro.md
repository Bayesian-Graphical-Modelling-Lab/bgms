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
#> intrusion (1)  0.471 0.005 0.232 2192.768 1.002
#> intrusion (2) -1.918 0.009 0.332 1497.613 1.004
#> intrusion (3) -4.868 0.015 0.538 1328.778 1.004
#> intrusion (4) -9.550 0.023 0.875 1461.059 1.003
#> dreams (1)    -0.592 0.004 0.198 1989.068 1.002
#> dreams (2)    -3.787 0.010 0.361 1416.890 1.001
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean    sd  mcse    n_eff  Rhat
#> intrusion-dreams  0.631 0.001 0.068 2381.502 1.002
#> intrusion-flash   0.339 0.001 0.062 2799.433 1.001
#> intrusion-upset   0.202 0.066 0.003  486.370 1.001
#> intrusion-physior 0.190 0.073 0.005  211.371 1.002
#> dreams-flash      0.499 0.001 0.061 2234.576 1.001
#> dreams-upset      0.226 0.056 0.002  974.685 1.001
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean    sd  mcse n0->0 n0->1 n1->0 n1->1   n_eff
#> intrusion-dreams  1.000 0.000           0     0     0  3999        
#> intrusion-flash   1.000 0.000           0     0     0  3999        
#> intrusion-upset   0.973 0.161 0.013    99     8     8  3884 159.781
#> intrusion-physior 0.933 0.250 0.024   254    13    13  3719  107.14
#> dreams-flash      1.000 0.000           0     0     0  3999        
#> dreams-upset      0.994 0.076 0.006    21     2     2  3974 182.918
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.067
#> intrusion-physior 1.021
#> dreams-flash           
#> dreams-upset       1.24
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
#> intrusion  0.4713730 -1.917789 -4.867864  -9.550006
#> dreams    -0.5922825 -3.786948 -7.116203 -11.546537
#> flash     -0.1005127 -2.557873 -5.353125  -9.655656
#> upset      0.4192987 -1.312709 -3.381854  -7.055391
#> physior   -0.6030730 -3.145634 -6.174537 -10.495866
#> 
#> $pairwise
#>           intrusion      dreams       flash       upset     physior
#> intrusion 0.0000000 0.631353688 0.338976477 0.202322026 0.189896271
#> dreams    0.6313537 0.000000000 0.498735904 0.226326133 0.006758395
#> flash     0.3389765 0.498735904 0.000000000 0.004674305 0.307737254
#> upset     0.2023220 0.226326133 0.004674305 0.000000000 0.708467755
#> physior   0.1898963 0.006758395 0.307737254 0.708467755 0.000000000
#> 
#> $indicator
#>           intrusion  dreams   flash   upset physior
#> intrusion   0.00000 1.00000 1.00000 0.97325 0.93325
#> dreams      1.00000 0.00000 1.00000 0.99425 0.07175
#> flash       1.00000 1.00000 0.00000 0.05525 1.00000
#> upset       0.97325 0.99425 0.05525 0.00000 1.00000
#> physior     0.93325 0.07175 1.00000 1.00000 0.00000
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
