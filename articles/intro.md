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
#> intrusion (1)  0.484 0.005 0.236 2422.715 1.000
#> intrusion (2) -1.889 0.009 0.341 1468.477 1.001
#> intrusion (3) -4.824 0.016 0.555 1132.811 1.001
#> intrusion (4) -9.477 0.026 0.886 1188.985 1.002
#> dreams (1)    -0.588 0.004 0.190 2485.622 1.001
#> dreams (2)    -3.782 0.008 0.343 2085.935 1.001
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean    sd  mcse    n_eff  Rhat
#> intrusion-dreams  0.629 0.068 0.001 3176.979 1.001
#> intrusion-flash   0.337 0.001 0.062 2687.296 1.003
#> intrusion-upset   0.191 0.076 0.006  152.560 1.004
#> intrusion-physior 0.199 0.065 0.003  468.944 1.001
#> dreams-flash      0.500 0.001 0.059 4266.768 1.001
#> dreams-upset      0.227 0.062 0.005  183.504 1.004
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean    sd  mcse n0->0 n0->1 n1->0 n1->1   n_eff
#> intrusion-dreams  0.999 0.027 0.001     2     1     1  3995 800.721
#> intrusion-flash   1.000 0.000           0     0     0  3999        
#> intrusion-upset   0.930 0.256 0.029   271    10    10  3708  78.047
#> intrusion-physior 0.972 0.165 0.012   102    10    10  3877 192.561
#> dreams-flash      1.000 0.000           0     0     0  3999        
#> dreams-upset      0.981 0.137 0.019    75     2     2  3920  53.679
#>                    Rhat
#> intrusion-dreams  1.292
#> intrusion-flash        
#> intrusion-upset   1.067
#> intrusion-physior 1.073
#> dreams-flash           
#> dreams-upset      1.339
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
#> intrusion  0.4837334 -1.889465 -4.824416  -9.477097
#> dreams    -0.5879116 -3.781981 -7.106266 -11.537586
#> flash     -0.1049549 -2.566646 -5.362953  -9.665535
#> upset      0.4158049 -1.305286 -3.367153  -7.027617
#> physior   -0.6123752 -3.172404 -6.221835 -10.579970
#> 
#> $pairwise
#>           intrusion      dreams       flash       upset     physior
#> intrusion 0.0000000 0.628937007 0.337201523 0.191010300 0.198537214
#> dreams    0.6289370 0.000000000 0.499527156 0.226692666 0.007061984
#> flash     0.3372015 0.499527156 0.000000000 0.008840418 0.305169126
#> upset     0.1910103 0.226692666 0.008840418 0.000000000 0.711333857
#> physior   0.1985372 0.007061984 0.305169126 0.711333857 0.000000000
#> 
#> $indicator
#>           intrusion  dreams   flash   upset physior
#> intrusion   0.00000 0.99925 1.00000 0.92975 0.97200
#> dreams      0.99925 0.00000 1.00000 0.98075 0.06675
#> flash       1.00000 1.00000 0.00000 0.08125 1.00000
#> upset       0.92975 0.98075 0.08125 0.00000 1.00000
#> physior     0.97200 0.06675 1.00000 1.00000 0.00000
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
