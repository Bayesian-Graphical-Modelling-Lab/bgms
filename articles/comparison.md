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
#> 1    avoid (1) -2.665 0.008 0.377 2475.955 1.001
#> 2 closeatt (1) -2.256 0.006 0.374 3418.632 1.000
#> 3 distract (1) -0.472 0.009 0.335 1433.373 1.001
#> 4   forget (1) -1.592 0.007 0.329 2188.826 1.000
#> 5 instruct (1) -2.421 0.010 0.391 1520.108 1.002
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.960 0.012 0.460 1529.980 1.000
#> 2    avoid-distract  1.693 0.007 0.350 2592.085 1.000
#> 3      avoid-forget  0.524 0.009 0.384 1749.304 1.001
#> 4    avoid-instruct  0.379 0.013 0.485 1503.087 1.000
#> 5 closeatt-distract -0.257 0.007 0.392 2933.834 1.001
#> 6   closeatt-forget  0.150 0.005 0.300 3234.641 1.003
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean    sd  mcse n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000 0.000           0     0     0  3999
#>  avoid-closeatt (pairwise) 0.777 0.416 0.012   584   306   307  2802
#>  avoid-distract (pairwise) 0.406 0.491 0.009  1541   833   834   791
#>    avoid-forget (pairwise) 0.859 0.348 0.011   364   200   200  3235
#>  avoid-instruct (pairwise) 0.989 0.104 0.003    25    19    19  3936
#>            closeatt (main) 1.000 0.000           0     0     0  3999
#>     n_eff  Rhat
#>                
#>  1136.654 1.001
#>  3042.969     1
#>  1040.431 1.002
#>  1117.125 1.022
#>                
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant (all 0 or all 1) across all iterations, so sd/mcse/n_eff/Rhat
#> are undefined; `summary(fit)$indicator` still contains the NA values.
#> 
#> Group differences (main effects):
#>            parameter   mean    sd mcse n_eff  Rhat
#>     avoid (diff1; 1) -2.546 0.750            1.002
#>  closeatt (diff1; 1) -2.982 0.751            1.001
#>  distract (diff1; 1) -2.567 0.682            1.001
#>    forget (diff1; 1) -2.837 0.680            1.001
#>  instruct (diff1; 1) -2.361 0.889            1.002
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean    sd  mcse    n_eff  Rhat
#>     avoid-closeatt (diff1)  1.214 0.927 0.024 1505.070 1.001
#>     avoid-distract (diff1)  0.232 0.377 0.009 1874.807 1.000
#>       avoid-forget (diff1)  1.317 0.823 0.022 1435.206 1.000
#>     avoid-instruct (diff1) -2.824 1.063 0.033 1012.104 1.001
#>  closeatt-distract (diff1) -0.181 0.354 0.009 1669.822 1.001
#>    closeatt-forget (diff1)  0.161 0.321 0.008 1546.368 1.000
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
#>               baseline     diff1
#> avoid(c1)    -2.665212 -2.546188
#> closeatt(c1) -2.255548 -2.982387
#> distract(c1) -0.472410 -2.567425
#> forget(c1)   -1.592422 -2.836587
#> instruct(c1) -2.421141 -2.361171
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.9603631  1.2141235
#> avoid-distract     1.6925403  0.2319958
#> avoid-forget       0.5242297  1.3165539
#> avoid-instruct     0.3785191 -2.8238631
#> closeatt-distract -0.2568894 -0.1808469
#> closeatt-forget    0.1497185  0.1606943
#> closeatt-instruct  1.5764951  0.6177448
#> distract-forget    0.4011882  0.2253859
#> distract-instruct  1.2557711  1.2598830
#> forget-instruct    1.1270862  0.8230455
#> 
#> $main_effects_groups
#>                  group1    group2
#> avoid(c1)    -1.3921184 -3.938306
#> closeatt(c1) -0.7643550 -3.746742
#> distract(c1)  0.8113023 -1.756122
#> forget(c1)   -0.1741286 -3.010715
#> instruct(c1) -1.2405558 -3.601727
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.35330133  1.5674249
#> avoid-distract     1.57654242  1.8085382
#> avoid-forget      -0.13404730  1.1825066
#> avoid-instruct     1.79045061 -1.0334125
#> closeatt-distract -0.16646598 -0.3473129
#> closeatt-forget    0.06937133  0.2300656
#> closeatt-instruct  1.26762266  1.8853675
#> distract-forget    0.28849521  0.5138811
#> distract-instruct  0.62582956  1.8857126
#> forget-instruct    0.71556346  1.5386090
#> 
#> $indicators
#>            avoid closeatt distract  forget instruct
#> avoid    1.00000  0.77725  0.40625 0.85900  0.98900
#> closeatt 0.77725  1.00000  0.37800 0.37525  0.60225
#> distract 0.40625  0.37800  1.00000 0.40825  0.84525
#> forget   0.85900  0.37525  0.40825 1.00000  0.73325
#> instruct 0.98900  0.60225  0.84525 0.73325  1.00000
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
