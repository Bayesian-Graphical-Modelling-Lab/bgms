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
#> 1    avoid (1) -2.657 0.007 0.380 2656.575 1.001
#> 2 closeatt (1) -2.249 0.007 0.373 3057.829 1.000
#> 3 distract (1) -0.475 0.009 0.335 1495.267 1.001
#> 4   forget (1) -1.595 0.007 0.334 2328.192 1.001
#> 5 instruct (1) -2.435 0.010 0.402 1483.464 1.003
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.954 0.011 0.456 1675.219 1.001
#> 2    avoid-distract  1.684 0.007 0.353 2627.475 1.001
#> 3      avoid-forget  0.506 0.009 0.386 1735.184 1.001
#> 4    avoid-instruct  0.386 0.013 0.479 1400.903 1.001
#> 5 closeatt-distract -0.264 0.007 0.393 2863.346 1.001
#> 6   closeatt-forget  0.148 0.005 0.311 3219.868 1.002
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean    sd  mcse n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000 0.000           0     0     0  3999
#>  avoid-closeatt (pairwise) 0.773 0.419 0.012   590   318   318  2773
#>  avoid-distract (pairwise) 0.397 0.489 0.009  1573   840   840   746
#>    avoid-forget (pairwise) 0.834 0.373 0.012   440   226   226  3107
#>  avoid-instruct (pairwise) 0.993 0.085 0.002    14    15    15  3955
#>            closeatt (main) 1.000 0.000           0     0     0  3999
#>     n_eff  Rhat
#>                
#>  1171.632 1.002
#>  3128.528     1
#>  1022.432 1.002
#>  1409.132 1.106
#>                
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant (all 0 or all 1) across all iterations, so sd/mcse/n_eff/Rhat
#> are undefined; `summary(fit)$indicator` still contains the NA values.
#> 
#> Group differences (main effects):
#>            parameter   mean    sd mcse n_eff  Rhat
#>     avoid (diff1; 1) -2.527 0.734            1.002
#>  closeatt (diff1; 1) -2.988 0.741            1.001
#>  distract (diff1; 1) -2.569 0.679            1.001
#>    forget (diff1; 1) -2.815 0.673            1.003
#>  instruct (diff1; 1) -2.359 0.913            1.002
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean    sd  mcse    n_eff  Rhat
#>     avoid-closeatt (diff1)  1.197 0.926 0.023 1572.974 1.001
#>     avoid-distract (diff1)  0.226 0.374 0.009 1811.493 1.000
#>       avoid-forget (diff1)  1.261 0.828 0.022 1414.174 1.002
#>     avoid-instruct (diff1) -2.805 1.040 0.030 1227.945 1.001
#>  closeatt-distract (diff1) -0.180 0.356 0.009 1721.366 1.000
#>    closeatt-forget (diff1)  0.171 0.333 0.009 1432.740 1.000
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
#> avoid(c1)    -2.6574124 -2.527189
#> closeatt(c1) -2.2486194 -2.988018
#> distract(c1) -0.4746695 -2.568726
#> forget(c1)   -1.5952835 -2.815087
#> instruct(c1) -2.4354222 -2.358658
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.9544429  1.1969727
#> avoid-distract     1.6844564  0.2260742
#> avoid-forget       0.5056182  1.2605800
#> avoid-instruct     0.3859527 -2.8050810
#> closeatt-distract -0.2637434 -0.1799602
#> closeatt-forget    0.1482965  0.1708644
#> closeatt-instruct  1.5751194  0.6187862
#> distract-forget    0.3970111  0.2139821
#> distract-instruct  1.2682923  1.2608154
#> forget-instruct    1.1359777  0.8168080
#> 
#> $main_effects_groups
#>                  group1    group2
#> avoid(c1)    -1.3938177 -3.921007
#> closeatt(c1) -0.7546105 -3.742628
#> distract(c1)  0.8096933 -1.759032
#> forget(c1)   -0.1877399 -3.002827
#> instruct(c1) -1.2560933 -3.614751
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.35595657  1.5529292
#> avoid-distract     1.57141931  1.7974935
#> avoid-forget      -0.12467179  1.1359082
#> avoid-instruct     1.78849317 -1.0165878
#> closeatt-distract -0.17376335 -0.3537235
#> closeatt-forget    0.06286433  0.2337287
#> closeatt-instruct  1.26572629  1.8845125
#> distract-forget    0.29002001  0.5040021
#> distract-instruct  0.63788463  1.8987000
#> forget-instruct    0.72757370  1.5443817
#> 
#> $indicators
#>            avoid closeatt distract  forget instruct
#> avoid    1.00000  0.77300  0.39675 0.83350  0.99275
#> closeatt 0.77300  1.00000  0.37400 0.38825  0.60800
#> distract 0.39675  0.37400  1.00000 0.39425  0.84525
#> forget   0.83350  0.38825  0.39425 1.00000  0.72600
#> instruct 0.99275  0.60800  0.84525 0.72600  1.00000
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
