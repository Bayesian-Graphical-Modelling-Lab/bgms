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

Note: During fitting, progress bars are shown in interactive sessions.
In this vignette, they are suppressed for clarity. Sampling can take a
while; the progress bars usually help track progress.

## Posterior summaries

``` r
summary(fit)
#> Posterior summaries from Bayesian estimation:
#> 
#> Category thresholds:
#>                 mean  mcse    sd    n_eff  Rhat
#> intrusion (1)  0.489 0.006 0.233 1596.771 1.010
#> intrusion (2) -1.888 0.010 0.340 1053.243 1.024
#> intrusion (3) -4.828 0.019 0.556  827.554 1.025
#> intrusion (4) -9.490 0.031 0.897  860.034 1.023
#> dreams (1)    -0.597 0.005 0.190 1610.037 1.002
#> dreams (2)    -3.802 0.010 0.341 1174.409 1.003
#> ... (use `summary(fit)$main` to see full output)
#> 
#> Pairwise interactions:
#>                    mean    sd  mcse    n_eff  Rhat
#> intrusion-dreams  0.632 0.002 0.064 1716.164 1.000
#> intrusion-flash   0.338 0.001 0.062 1823.683 1.000
#> intrusion-upset   0.191 0.076 0.006  165.946 1.033
#> intrusion-physior 0.198 0.066 0.004  323.018 1.005
#> dreams-flash      0.498 0.001 0.060 2264.802 1.000
#> dreams-upset      0.231 0.056 0.002  721.220 1.008
#> ... (use `summary(fit)$pairwise` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an 
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise` still contains the NA values.
#> 
#> Inclusion probabilities:
#>                    mean    sd  mcse n0->0 n0->1 n1->0 n1->1   n_eff
#> intrusion-dreams  1.000 0.000           0     0     0  3999        
#> intrusion-flash   1.000 0.000           0     0     0  3999        
#> intrusion-upset   0.925 0.263 0.028   287    12    12  3688  88.677
#> intrusion-physior 0.966 0.182 0.016   129     8     8  3854 124.701
#> dreams-flash      1.000 0.000           0     0     0  3999        
#> dreams-upset      0.996 0.063 0.004    14     2     2  3981  267.81
#>                    Rhat
#> intrusion-dreams       
#> intrusion-flash        
#> intrusion-upset   1.173
#> intrusion-physior 1.059
#> dreams-flash           
#> dreams-upset        1.3
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
#> intrusion  0.4889848 -1.887789 -4.827802  -9.490176
#> dreams    -0.5974490 -3.801634 -7.133124 -11.578501
#> flash     -0.1063310 -2.567018 -5.359779  -9.659237
#> upset      0.4144289 -1.303246 -3.365035  -7.015933
#> physior   -0.6147997 -3.170385 -6.211313 -10.550982
#> 
#> $pairwise
#>           intrusion      dreams       flash       upset     physior
#> intrusion 0.0000000 0.631924323 0.338297975 0.190695963 0.198176315
#> dreams    0.6319243 0.000000000 0.498040695 0.230776860 0.005254907
#> flash     0.3382980 0.498040695 0.000000000 0.006462176 0.307138582
#> upset     0.1906960 0.230776860 0.006462176 0.000000000 0.707867139
#> physior   0.1981763 0.005254907 0.307138582 0.707867139 0.000000000
#> 
#> $indicator
#>           intrusion  dreams   flash   upset physior
#> intrusion   0.00000 1.00000 1.00000 0.96575  0.0560
#> dreams      1.00000 0.00000 0.92525 1.00000  0.0655
#> flash       1.00000 0.92525 0.00000 0.99600  1.0000
#> upset       0.96575 1.00000 0.99600 0.00000  1.0000
#> physior     0.05600 0.06550 1.00000 1.00000  0.0000
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
