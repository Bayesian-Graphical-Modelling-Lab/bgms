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
#> 1    avoid (1) -2.576 0.009 0.383 2000.000 1.002
#> 2 closeatt (1) -2.240 0.008 0.348 2019.399 1.003
#> 3 distract (1) -0.442 0.009 0.319 1201.528 1.003
#> 4   forget (1) -1.530 0.007 0.325 2021.753 1.005
#> 5 instruct (1) -2.326 0.010 0.365 1415.353 1.002
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.837 0.014 0.454 1086.078 1.000
#> 2    avoid-distract  1.630 0.008 0.355 1968.617 1.004
#> 3      avoid-forget  0.457 0.010 0.368 1308.338 1.000
#> 4    avoid-instruct  0.365 0.010 0.420 1804.310 1.007
#> 5 closeatt-distract -0.176 0.007 0.352 2641.312 1.001
#> 6   closeatt-forget  0.137 0.006 0.268 2230.406 1.001
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean  mcse    sd n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000       0.000     0     0     0  1999
#>  avoid-closeatt (pairwise) 0.728 0.018 0.445   361   183   184  1271
#>  avoid-distract (pairwise) 0.419 0.012 0.493   721   441   441   396
#>    avoid-forget (pairwise) 0.827 0.017 0.378   229   117   116  1537
#>  avoid-instruct (pairwise) 0.993 0.005 0.083    10     4     4  1981
#>            closeatt (main) 1.000       0.000     0     0     0  1999
#>  n_eff_mixt  Rhat
#>                  
#>     602.231     1
#>    1657.644     1
#>     513.077 1.003
#>     336.079 1.133
#>                  
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;
#> `summary(fit)$indicator` still contains all computed values.
#> 
#> Group differences (main effects):
#>            parameter   mean mcse    sd    n_eff n_eff_mixt  Rhat
#>     avoid (diff1; 1) -2.573      0.748 1597.684            1.000
#>  closeatt (diff1; 1) -2.933      0.741 2173.052            1.000
#>  distract (diff1; 1) -2.599      0.642 1283.952            1.008
#>    forget (diff1; 1) -2.854      0.656 1334.369            1.000
#>  instruct (diff1; 1) -2.397      0.852 1115.259            1.000
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean  mcse    sd    n_eff n_eff_mixt
#>     avoid-closeatt (diff1)  1.017 0.030 0.880  892.555    879.528
#>     avoid-distract (diff1)  0.229 0.011 0.369 1459.609   1140.047
#>       avoid-forget (diff1)  1.208 0.028 0.793  716.949    795.849
#>     avoid-instruct (diff1) -2.795 0.025 0.981 1488.140   1516.295
#>  closeatt-distract (diff1) -0.164 0.010 0.335 1753.220   1026.454
#>    closeatt-forget (diff1)  0.145 0.009 0.280 1302.618    947.503
#>   Rhat
#>  1.001
#>  1.005
#>  1.003
#>  1.008
#>  1.001
#>  1.000
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
#>               baseline     diff1
#> avoid(c1)    -2.575979 -2.572790
#> closeatt(c1) -2.239661 -2.933408
#> distract(c1) -0.442476 -2.598615
#> forget(c1)   -1.529868 -2.853765
#> instruct(c1) -2.325596 -2.397421
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.8370791  1.0170685
#> avoid-distract     1.6298165  0.2294306
#> avoid-forget       0.4567097  1.2076873
#> avoid-instruct     0.3647285 -2.7951886
#> closeatt-distract -0.1764299 -0.1642858
#> closeatt-forget    0.1371961  0.1446977
#> closeatt-instruct  1.4916620  0.5402549
#> distract-forget    0.3690295  0.2202865
#> distract-instruct  1.1998664  1.2600559
#> forget-instruct    1.0508252  0.7518878
#> 
#> $main_effects_groups
#>                  group1    group2
#> avoid(c1)    -1.2895843 -3.862375
#> closeatt(c1) -0.7729574 -3.706365
#> distract(c1)  0.8568314 -1.741783
#> forget(c1)   -0.1029852 -2.956750
#> instruct(c1) -1.1268850 -3.524306
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.32854486  1.3456134
#> avoid-distract     1.51510123  1.7445318
#> avoid-forget      -0.14713390  1.0605534
#> avoid-instruct     1.76232281 -1.0328657
#> closeatt-distract -0.09428695 -0.2585728
#> closeatt-forget    0.06484720  0.2095449
#> closeatt-instruct  1.22153450  1.7617894
#> distract-forget    0.25888626  0.4791728
#> distract-instruct  0.56983849  1.8298944
#> forget-instruct    0.67488129  1.4267691
#> 
#> $indicators
#>           avoid closeatt distract forget instruct
#> avoid    1.0000   0.7275   0.4190 0.8270   0.9930
#> closeatt 0.7275   1.0000   0.3660 0.3340   0.5620
#> distract 0.4190   0.3660   1.0000 0.3905   0.8455
#> forget   0.8270   0.3340   0.3905 1.0000   0.6740
#> instruct 0.9930   0.5620   0.8455 0.6740   1.0000
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
