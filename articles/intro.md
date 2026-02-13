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
#> intrusion (1)  0.476 0.005 0.226 2288.548 1.000
#> intrusion (2) -1.908 0.008 0.326 1519.036 1.002
#> intrusion (3) -4.852 0.014 0.537 1427.050 1.003
#> intrusion (4) -9.530 0.023 0.875 1417.951 1.003
#> dreams (1)    -0.590 0.005 0.195 1783.234 1.000
#> dreams (2)    -3.779 0.010 0.359 1349.062 1.000
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean    sd  mcse    n_eff  Rhat
#> intrusion-dreams  0.630 0.001 0.067 2320.130 1.001
#> intrusion-flash   0.339 0.001 0.060 2537.984 1.000
#> intrusion-upset   0.202 0.066 0.003  461.062 1.001
#> intrusion-physior 0.190 0.070 0.004  284.454 1.000
#> dreams-flash      0.501 0.001 0.062 2209.910 1.002
#> dreams-upset      0.223 0.059 0.003  541.946 1.002
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean    sd  mcse n0->0 n0->1 n1->0 n1->1   n_eff
#> intrusion-dreams  1.000 0.000           0     0     0  3999        
#> intrusion-flash   1.000 0.000           0     0     0  3999        
#> intrusion-upset   0.975 0.156 0.013    93     7     7  3892 148.937
#> intrusion-physior 0.946 0.227  0.02   205    13    13  3768  130.25
#> dreams-flash      1.000 0.000           0     0     0  3999        
#> dreams-upset      0.986 0.119  0.01    53     4     4  3938 147.635
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.062
#> intrusion-physior 1.005
#> dreams-flash           
#> dreams-upset      1.208
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
#> intrusion  0.4762725 -1.907620 -4.852093  -9.530335
#> dreams    -0.5902608 -3.778569 -7.100904 -11.522637
#> flash     -0.1071002 -2.570685 -5.380333  -9.700127
#> upset      0.4181016 -1.307836 -3.372591  -7.036536
#> physior   -0.6030401 -3.144014 -6.173321 -10.496004
#> 
#> $pairwise
#>           intrusion      dreams       flash       upset     physior
#> intrusion 0.0000000 0.630460950 0.338828174 0.201870274 0.190222148
#> dreams    0.6304609 0.000000000 0.500702270 0.222584657 0.007375532
#> flash     0.3388282 0.500702270 0.000000000 0.007425009 0.307834698
#> upset     0.2018703 0.222584657 0.007425009 0.000000000 0.707539083
#> physior   0.1902221 0.007375532 0.307834698 0.707539083 0.000000000
#> 
#> $indicator
#>           intrusion  dreams   flash   upset physior
#> intrusion    0.0000 1.00000 1.00000 0.97500 0.94550
#> dreams       1.0000 0.00000 1.00000 0.98575 0.07575
#> flash        1.0000 1.00000 0.00000 0.07225 1.00000
#> upset        0.9750 0.98575 0.07225 0.00000 1.00000
#> physior      0.9455 0.07575 1.00000 1.00000 0.00000
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
