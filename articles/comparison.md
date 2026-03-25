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
#> 1    avoid (1) -2.590 0.011 0.378 1199.128 1.001
#> 2 closeatt (1) -2.237 0.011 0.369 1105.739 1.000
#> 3 distract (1) -0.460 0.013 0.327  677.883 1.001
#> 4   forget (1) -1.542 0.011 0.323  884.620 1.002
#> 5 instruct (1) -2.340 0.012 0.389  980.089 1.009
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.797 0.017 0.445  686.463 1.000
#> 2    avoid-distract  1.639 0.009 0.357 1422.468 1.002
#> 3      avoid-forget  0.452 0.013 0.372  851.489 1.002
#> 4    avoid-instruct  0.416 0.014 0.416  943.388 1.000
#> 5 closeatt-distract -0.181 0.010 0.363 1301.033 1.001
#> 6   closeatt-forget  0.151 0.008 0.291 1222.684 1.004
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean  mcse    sd n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000       0.000     0     0     0  1999
#>  avoid-closeatt (pairwise) 0.674 0.018 0.469   434   217   217  1131
#>  avoid-distract (pairwise) 0.401 0.012 0.490   762   435   435   367
#>    avoid-forget (pairwise) 0.822 0.016 0.383   227   129   129  1514
#>  avoid-instruct (pairwise) 0.993 0.003 0.086     4    11    11  1973
#>            closeatt (main) 1.000       0.000     0     0     0  1999
#>  n_eff_mixt  Rhat
#>                  
#>     656.594     1
#>    1655.647     1
#>     565.541     1
#>    1171.778 1.019
#>                  
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;
#> `summary(fit)$indicator` still contains all computed values.
#> 
#> Group differences (main effects):
#>            parameter   mean mcse    sd   n_eff n_eff_mixt  Rhat
#>     avoid (diff1; 1) -2.541      0.738 954.171            1.000
#>  closeatt (diff1; 1) -2.883      0.720 937.369            1.002
#>  distract (diff1; 1) -2.575      0.673 592.628            1.003
#>    forget (diff1; 1) -2.849      0.658 717.155            1.005
#>  instruct (diff1; 1) -2.377      0.891 491.647            1.009
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean  mcse    sd    n_eff n_eff_mixt
#>     avoid-closeatt (diff1)  0.894 0.033 0.839  483.038    656.029
#>     avoid-distract (diff1)  0.216 0.012 0.362 1101.519    987.765
#>       avoid-forget (diff1)  1.190 0.030 0.803  511.288    715.328
#>     avoid-instruct (diff1) -2.665 0.034 0.970  769.289    832.794
#>  closeatt-distract (diff1) -0.163 0.012 0.324  961.690    763.733
#>    closeatt-forget (diff1)  0.166 0.012 0.324 1059.494    729.607
#>   Rhat
#>  1.003
#>  1.003
#>  1.000
#>  1.001
#>  1.001
#>  1.002
#> ... (use `summary(fit)$pairwise_diff` to see full output)
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
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
#> avoid(c1)    -2.5897605 -2.540850
#> closeatt(c1) -2.2371412 -2.882947
#> distract(c1) -0.4598562 -2.574924
#> forget(c1)   -1.5422258 -2.849481
#> instruct(c1) -2.3399811 -2.377159
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.7965309  0.8944052
#> avoid-distract     1.6392659  0.2163085
#> avoid-forget       0.4522752  1.1902215
#> avoid-instruct     0.4158560 -2.6651747
#> closeatt-distract -0.1805078 -0.1628199
#> closeatt-forget    0.1508933  0.1664507
#> closeatt-instruct  1.4830774  0.5053316
#> distract-forget    0.3773759  0.2330934
#> distract-instruct  1.2050536  1.2045313
#> forget-instruct    1.0588387  0.7514445
#> 
#> $main_effects_groups
#>                  group1    group2
#> avoid(c1)    -1.3193357 -3.860185
#> closeatt(c1) -0.7956675 -3.678615
#> distract(c1)  0.8276056 -1.747318
#> forget(c1)   -0.1174855 -2.966966
#> instruct(c1) -1.1514015 -3.528561
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.34932830  1.2437335
#> avoid-distract     1.53111163  1.7474202
#> avoid-forget      -0.14283557  1.0473860
#> avoid-instruct     1.74844335 -0.9167313
#> closeatt-distract -0.09909781 -0.2619177
#> closeatt-forget    0.06766800  0.2341187
#> closeatt-instruct  1.23041161  1.7357432
#> distract-forget    0.26082915  0.4939226
#> distract-instruct  0.60278797  1.8073192
#> forget-instruct    0.68311642  1.4345609
#> 
#> $indicators
#>           avoid closeatt distract forget instruct
#> avoid    1.0000   0.6745   0.4010 0.8220   0.9925
#> closeatt 0.6745   1.0000   0.3675 0.3730   0.5635
#> distract 0.4010   0.3675   1.0000 0.4035   0.8170
#> forget   0.8220   0.3730   0.4035 1.0000   0.7060
#> instruct 0.9925   0.5635   0.8170 0.7060   1.0000
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
