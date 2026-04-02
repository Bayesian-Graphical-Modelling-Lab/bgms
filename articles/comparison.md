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
#> 1    avoid (1) -2.573 0.012 0.388 1088.359 1.000
#> 2 closeatt (1) -2.232 0.012 0.374 1006.900 1.003
#> 3 distract (1) -0.441 0.014 0.331  591.792 1.001
#> 4   forget (1) -1.566 0.011 0.329  936.229 1.000
#> 5 instruct (1) -2.327 0.016 0.405  652.123 1.002
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.821 0.020 0.462  557.673 1.000
#> 2    avoid-distract  1.621 0.011 0.366 1136.131 1.000
#> 3      avoid-forget  0.464 0.013 0.359  790.028 1.001
#> 4    avoid-instruct  0.352 0.016 0.455  795.739 1.005
#> 5 closeatt-distract -0.195 0.012 0.376  978.194 1.001
#> 6   closeatt-forget  0.134 0.007 0.285 1484.529 1.003
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean  mcse    sd n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000       0.000     0     0     0  1999
#>  avoid-closeatt (pairwise) 0.686 0.019 0.464   424   204   203  1168
#>  avoid-distract (pairwise) 0.410 0.012 0.492   737   443   442   377
#>    avoid-forget (pairwise) 0.809 0.017 0.393   250   130   131  1488
#>  avoid-instruct (pairwise) 0.987 0.006 0.115    19     8     8  1964
#>            closeatt (main) 1.000       0.000     0     0     0  1999
#>  n_eff_mixt  Rhat
#>                  
#>     619.357 1.003
#>    1686.997     1
#>     536.493     1
#>      353.43   1.1
#>                  
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;
#> `summary(fit)$indicator` still contains all computed values.
#> 
#> Group differences (main effects):
#>            parameter   mean mcse    sd   n_eff n_eff_mixt  Rhat
#>     avoid (diff1; 1) -2.545      0.721 865.191            1.002
#>  closeatt (diff1; 1) -2.949      0.739 934.786            1.001
#>  distract (diff1; 1) -2.597      0.674 589.140            1.001
#>    forget (diff1; 1) -2.809      0.664 617.994            1.001
#>  instruct (diff1; 1) -2.431      0.887 542.986            1.000
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean  mcse    sd    n_eff n_eff_mixt
#>     avoid-closeatt (diff1)  0.984 0.035 0.893  435.921    636.983
#>     avoid-distract (diff1)  0.217 0.012 0.355 1072.937    810.018
#>       avoid-forget (diff1)  1.179 0.032 0.814  398.148    649.006
#>     avoid-instruct (diff1) -2.793 0.041 1.047  611.529    659.032
#>  closeatt-distract (diff1) -0.170 0.012 0.346 1268.260    817.227
#>    closeatt-forget (diff1)  0.142 0.011 0.290 1117.873    745.419
#>   Rhat
#>  1.000
#>  1.000
#>  1.000
#>  1.006
#>  1.005
#>  1.013
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
#> avoid(c1)    -2.5726946 -2.544880
#> closeatt(c1) -2.2321390 -2.948739
#> distract(c1) -0.4405119 -2.597299
#> forget(c1)   -1.5663952 -2.809333
#> instruct(c1) -2.3274717 -2.431185
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.8209105  0.9841867
#> avoid-distract     1.6214650  0.2166373
#> avoid-forget       0.4640816  1.1793304
#> avoid-instruct     0.3515582 -2.7927843
#> closeatt-distract -0.1945874 -0.1696452
#> closeatt-forget    0.1338066  0.1422537
#> closeatt-instruct  1.4990671  0.5493630
#> distract-forget    0.3923724  0.2040284
#> distract-instruct  1.1974072  1.2864408
#> forget-instruct    1.0640331  0.7475858
#> 
#> $main_effects_groups
#>                  group1    group2
#> avoid(c1)    -1.3002547 -3.845135
#> closeatt(c1) -0.7577698 -3.706508
#> distract(c1)  0.8581375 -1.739161
#> forget(c1)   -0.1617288 -2.971061
#> instruct(c1) -1.1118793 -3.543064
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.32881713  1.3130038
#> avoid-distract     1.51314638  1.7297837
#> avoid-forget      -0.12558358  1.0537469
#> avoid-instruct     1.74795041 -1.0448339
#> closeatt-distract -0.10976480 -0.2794100
#> closeatt-forget    0.06267974  0.2049334
#> closeatt-instruct  1.22438565  1.7737486
#> distract-forget    0.29035823  0.4943866
#> distract-instruct  0.55418683  1.8406277
#> forget-instruct    0.69024020  1.4378261
#> 
#> $indicators
#>           avoid closeatt distract forget instruct
#> avoid    1.0000   0.6860    0.410 0.8095   0.9865
#> closeatt 0.6860   1.0000    0.385 0.3480   0.5725
#> distract 0.4100   0.3850    1.000 0.3730   0.8510
#> forget   0.8095   0.3480    0.373 1.0000   0.6925
#> instruct 0.9865   0.5725    0.851 0.6925   1.0000
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
