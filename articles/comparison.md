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
#> 1    avoid (1) -2.578 0.009 0.373 1850.336 1.003
#> 2 closeatt (1) -2.226 0.007 0.358 2872.778 1.000
#> 3 distract (1) -0.434 0.008 0.319 1441.611 1.002
#> 4   forget (1) -1.540 0.007 0.313 1869.054 1.003
#> 5 instruct (1) -2.318 0.011 0.400 1447.159 1.006
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.825 0.014 0.456 1011.859 1.002
#> 2    avoid-distract  1.634 0.008 0.359 2221.883 1.002
#> 3      avoid-forget  0.462 0.010 0.362 1404.882 1.001
#> 4    avoid-instruct  0.353 0.010 0.420 1742.657 1.000
#> 5 closeatt-distract -0.192 0.008 0.352 2174.121 1.000
#> 6   closeatt-forget  0.139 0.006 0.289 2430.323 1.000
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean  mcse    sd n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000       0.000     0     0     0  1999
#>  avoid-closeatt (pairwise) 0.713 0.018 0.452   376   197   198  1228
#>  avoid-distract (pairwise) 0.401 0.012 0.490   760   438   438   363
#>    avoid-forget (pairwise) 0.826 0.016 0.380   226   123   123  1527
#>  avoid-instruct (pairwise) 0.995 0.004 0.071     7     3     3  1986
#>            closeatt (main) 1.000       0.000     0     0     0  1999
#>  n_eff_mixt  Rhat
#>                  
#>     636.183 1.004
#>     1677.91     1
#>     542.881     1
#>     355.031 1.295
#>                  
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;
#> `summary(fit)$indicator` still contains all computed values.
#> 
#> Group differences (main effects):
#>            parameter   mean mcse    sd    n_eff n_eff_mixt  Rhat
#>     avoid (diff1; 1) -2.549      0.706 1691.722            1.000
#>  closeatt (diff1; 1) -2.955      0.711 1668.513            1.000
#>  distract (diff1; 1) -2.621      0.647 1361.614            1.004
#>    forget (diff1; 1) -2.856      0.628 1231.545            1.001
#>  instruct (diff1; 1) -2.455      0.880  880.365            1.005
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean  mcse    sd    n_eff n_eff_mixt
#>     avoid-closeatt (diff1)  0.992 0.029 0.869  924.062    902.040
#>     avoid-distract (diff1)  0.213 0.010 0.352 1644.140   1241.465
#>       avoid-forget (diff1)  1.222 0.028 0.817  796.582    845.062
#>     avoid-instruct (diff1) -2.809 0.026 0.962 1444.597   1346.241
#>  closeatt-distract (diff1) -0.162 0.010 0.324 1718.624   1064.849
#>    closeatt-forget (diff1)  0.159 0.011 0.320 1499.768    919.930
#>   Rhat
#>  1.002
#>  1.002
#>  1.000
#>  1.002
#>  1.000
#>  1.001
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
#> avoid(c1)    -2.578117 -2.549394
#> closeatt(c1) -2.225795 -2.955245
#> distract(c1) -0.434413 -2.620691
#> forget(c1)   -1.539815 -2.856442
#> instruct(c1) -2.317686 -2.455046
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.8252890  0.9916223
#> avoid-distract     1.6344618  0.2125706
#> avoid-forget       0.4620798  1.2218174
#> avoid-instruct     0.3527761 -2.8087362
#> closeatt-distract -0.1916711 -0.1615219
#> closeatt-forget    0.1392677  0.1585576
#> closeatt-instruct  1.5042358  0.5796393
#> distract-forget    0.3774852  0.2197390
#> distract-instruct  1.1821343  1.3008091
#> forget-instruct    1.0792849  0.7833525
#> 
#> $main_effects_groups
#>                  group1    group2
#> avoid(c1)    -1.3034206 -3.852814
#> closeatt(c1) -0.7481726 -3.703417
#> distract(c1)  0.8759325 -1.744759
#> forget(c1)   -0.1115942 -2.968037
#> instruct(c1) -1.0901630 -3.545209
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.32947780  1.3211001
#> avoid-distract     1.52817646  1.7407471
#> avoid-forget      -0.14882888  1.0729885
#> avoid-instruct     1.75714415 -1.0515920
#> closeatt-distract -0.11091014 -0.2724321
#> closeatt-forget    0.05998892  0.2185465
#> closeatt-instruct  1.21441613  1.7940554
#> distract-forget    0.26761565  0.4873547
#> distract-instruct  0.53172977  1.8325389
#> forget-instruct    0.68760865  1.4709611
#> 
#> $indicators
#>           avoid closeatt distract forget instruct
#> avoid    1.0000    0.713    0.401 0.8255    0.995
#> closeatt 0.7130    1.000    0.355 0.3810    0.606
#> distract 0.4010    0.355    1.000 0.3920    0.860
#> forget   0.8255    0.381    0.392 1.0000    0.715
#> instruct 0.9950    0.606    0.860 0.7150    1.000
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
