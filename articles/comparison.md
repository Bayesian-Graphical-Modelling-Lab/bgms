# Model Comparison with bgmCompare

## Introduction

The function
[`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
extends
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
to independent-sample designs. It estimates whether edge weights and
category thresholds differ across groups in an ordinal Markov random
field (MRF).

Posterior inclusion probabilities indicate how plausible it is that a
group difference exists in a given parameter. These can be converted to
Bayes factors for hypothesis testing.

## ADHD dataset

We illustrate with a subset from the `ADHD` dataset included in
**bgms**.

``` r
library(bgms)

?ADHD
data_adhd = ADHD[ADHD$group == 1, -1]
data_adhd = data_adhd[, 1:5]
data_no_adhd = ADHD[ADHD$group == 0, -1]
data_no_adhd = data_no_adhd[, 1:5]
```

## Fitting a model

``` r
fit = bgmCompare(x = data_adhd, y = data_no_adhd, seed = 1234)
```

## Posterior summaries

The summary shows both baseline effects and group differences:

``` r
summary(fit)
#> Posterior summaries from Bayesian grouped MRF estimation (bgmCompare):
#> 
#> Category thresholds:
#>      parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid (1) -2.655 0.011 0.372 1244.698 1.001
#> 2 closeatt (1) -2.253 0.010 0.375 1453.516 1.000
#> 3 distract (1) -0.465 0.012 0.336  735.474 1.003
#> 4   forget (1) -1.587 0.010 0.333 1043.249 1.001
#> 5 instruct (1) -2.424 0.015 0.397  738.351 1.007
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.957 0.016 0.462  793.466 1.002
#> 2    avoid-distract  1.690 0.010 0.347 1226.443 1.001
#> 3      avoid-forget  0.518 0.015 0.396  718.012 1.002
#> 4    avoid-instruct  0.370 0.021 0.494  574.555 1.000
#> 5 closeatt-distract -0.260 0.010 0.390 1477.546 1.004
#> 6   closeatt-forget  0.141 0.008 0.310 1502.834 1.005
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean    sd  mcse n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000 0.000           0     0     0  1999
#>  avoid-closeatt (pairwise) 0.779 0.415 0.018   294   148   148  1409
#>  avoid-distract (pairwise) 0.408 0.492 0.013   762   421   421   395
#>    avoid-forget (pairwise) 0.850 0.358 0.016   199   101   102  1597
#>  avoid-instruct (pairwise) 0.987 0.113 0.005    14    12    12  1961
#>            closeatt (main) 1.000 0.000           0     0     0  1999
#>    n_eff  Rhat
#>               
#>  547.602 1.007
#>  1545.49     1
#>  494.857 1.001
#>   610.32 1.025
#>               
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant (all 0 or all 1) across all iterations, so sd/mcse/n_eff/Rhat
#> are undefined; `summary(fit)$indicator` still contains the NA values.
#> 
#> Group differences (main effects):
#>            parameter   mean    sd mcse n_eff  Rhat
#>     avoid (diff1; 1) -2.565 0.735            1.000
#>  closeatt (diff1; 1) -3.019 0.753            1.000
#>  distract (diff1; 1) -2.588 0.686            1.002
#>    forget (diff1; 1) -2.861 0.685            1.000
#>  instruct (diff1; 1) -2.394 0.912            1.004
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean    sd  mcse   n_eff  Rhat
#>     avoid-closeatt (diff1)  1.235 0.933 0.034 737.078 1.002
#>     avoid-distract (diff1)  0.233 0.377 0.012 961.022 1.001
#>       avoid-forget (diff1)  1.315 0.839 0.032 693.493 1.000
#>     avoid-instruct (diff1) -2.840 1.094 0.053 420.571 1.001
#>  closeatt-distract (diff1) -0.178 0.364 0.012 916.006 1.002
#>    closeatt-forget (diff1)  0.166 0.334 0.013 672.762 1.001
#> ... (use `summary(fit)$pairwise_diff` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$pairwise_diff` still contains the NA values.
#> 
#> Use `summary(fit)$<component>` to access full results.
#> See the `easybgm` package for other summary and plotting tools.
```

You can extract posterior means and inclusion probabilities:

``` r
coef(fit)
#> $main_effects_raw
#>                baseline     diff1
#> avoid(c1)    -2.6550620 -2.564608
#> closeatt(c1) -2.2531722 -3.018548
#> distract(c1) -0.4646796 -2.587562
#> forget(c1)   -1.5866867 -2.860841
#> instruct(c1) -2.4239183 -2.393771
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.9574121  1.2345545
#> avoid-distract     1.6896390  0.2334367
#> avoid-forget       0.5176610  1.3153691
#> avoid-instruct     0.3703210 -2.8402576
#> closeatt-distract -0.2599956 -0.1778824
#> closeatt-forget    0.1411610  0.1659322
#> closeatt-instruct  1.5817824  0.6172343
#> distract-forget    0.3954378  0.2197663
#> distract-instruct  1.2536329  1.2806550
#> forget-instruct    1.1397938  0.8513931
#> 
#> $main_effects_groups
#>                  group1    group2
#> avoid(c1)    -1.3727581 -3.937366
#> closeatt(c1) -0.7438980 -3.762446
#> distract(c1)  0.8291014 -1.758461
#> forget(c1)   -0.1562662 -3.017107
#> instruct(c1) -1.2270329 -3.620804
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.34013484  1.5746894
#> avoid-distract     1.57292063  1.8063574
#> avoid-forget      -0.14002356  1.1753456
#> avoid-instruct     1.79044981 -1.0498078
#> closeatt-distract -0.17105434 -0.3489368
#> closeatt-forget    0.05819485  0.2241271
#> closeatt-instruct  1.27316530  1.8903996
#> distract-forget    0.28555460  0.5053209
#> distract-instruct  0.61330533  1.8939604
#> forget-instruct    0.71409731  1.5654904
#> 
#> $indicators
#>           avoid closeatt distract forget instruct
#> avoid    1.0000   0.7790   0.4085 0.8495   0.9870
#> closeatt 0.7790   1.0000   0.3760 0.3900   0.6025
#> distract 0.4085   0.3760   1.0000 0.4020   0.8470
#> forget   0.8495   0.3900   0.4020 1.0000   0.7480
#> instruct 0.9870   0.6025   0.8470 0.7480   1.0000
```

## Visualizing group networks

We can use the output to plot the network for the ADHD group:

``` r
library(qgraph)

adhd_network = matrix(0, 5, 5)
adhd_network[lower.tri(adhd_network)] = coef(fit)$pairwise_effects_groups[, 1]
adhd_network = adhd_network + t(adhd_network)
colnames(adhd_network) = colnames(data_adhd)
rownames(adhd_network) = colnames(data_adhd)

qgraph(adhd_network,
  theme = "TeamFortress",
  maximum = 1,
  fade = FALSE,
  color = c("#f0ae0e"), vsize = 10, repulsion = .9,
  label.cex = 1, label.scale = "FALSE",
  labels = colnames(data_adhd)
)
```

![](comparison_files/figure-html/unnamed-chunk-7-1.png)

## Next steps

- For a one-sample analysis, see the *Getting Started* vignette.
- For diagnostics and convergence checks, see the *Diagnostics*
  vignette.
- For additional analysis tools and more advanced plotting options,
  consider using the **easybgm** package, which integrates smoothly with
  **bgms** objects.
