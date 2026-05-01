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
#> 1    avoid (1) -2.581 0.006 0.385 4000.000 1.000
#> 2 closeatt (1) -2.236 0.005 0.369 4786.987 1.000
#> 3 distract (1) -0.432 0.007 0.324 2002.995 1.001
#> 4   forget (1) -1.531 0.006 0.319 3047.088 1.003
#> 5 instruct (1) -2.329 0.008 0.378 2508.437 1.000
#> 
#> Pairwise interactions:
#>           parameter   mean  mcse    sd    n_eff  Rhat
#> 1    avoid-closeatt  0.822 0.010 0.454 2200.373 1.000
#> 2    avoid-distract  1.636 0.005 0.355 4307.465 1.000
#> 3      avoid-forget  0.473 0.007 0.361 2755.046 1.000
#> 4    avoid-instruct  0.375 0.008 0.426 2936.511 1.001
#> 5 closeatt-distract -0.190 0.006 0.362 4286.014 1.001
#> 6   closeatt-forget  0.141 0.005 0.280 3724.811 1.001
#> ... (use `summary(fit)$pairwise` to see full output)
#> 
#> Inclusion probabilities:
#>                  parameter  mean  mcse    sd n0->0 n0->1 n1->0 n1->1
#>               avoid (main) 1.000       0.000     0     0     0  3999
#>  avoid-closeatt (pairwise) 0.703 0.012 0.457   769   419   419  2392
#>  avoid-distract (pairwise) 0.403 0.009 0.491  1513   874   874   738
#>    avoid-forget (pairwise) 0.849 0.011 0.359   385   220   221  3173
#>  avoid-instruct (pairwise) 0.996 0.002 0.067     9     9     9  3972
#>            closeatt (main) 1.000       0.000     0     0     0  3999
#>  n_eff_mixt  Rhat
#>                  
#>    1339.566 1.001
#>    3328.246     1
#>    1091.492 1.002
#>    1341.384 1.086
#>                  
#> ... (use `summary(fit)$indicator` to see full output)
#> Note: NA values are suppressed in the print table. They occur when an indicator
#> was constant or had fewer than 5 transitions, so n_eff_mixt is unreliable;
#> `summary(fit)$indicator` still contains all computed values.
#> 
#> Group differences (main effects):
#>            parameter   mean mcse    sd    n_eff n_eff_mixt  Rhat
#>     avoid (diff1; 1) -2.545      0.737 2901.285            1.000
#>  closeatt (diff1; 1) -2.931      0.717 3598.442            1.001
#>  distract (diff1; 1) -2.628      0.650 2067.904            1.003
#>    forget (diff1; 1) -2.878      0.644 2484.297            1.002
#>  instruct (diff1; 1) -2.408      0.869 1581.227            1.000
#> Note: NA values are suppressed in the print table. They occur here when an
#> indicator was zero across all iterations, so mcse/n_eff/n_eff_mixt/Rhat are undefined;
#> `summary(fit)$main_diff` still contains the NA values.
#> 
#> Group differences (pairwise effects):
#>                  parameter   mean  mcse    sd    n_eff n_eff_mixt Rhat
#>     avoid-closeatt (diff1)  0.972 0.020 0.866 1791.403   1806.502    1
#>     avoid-distract (diff1)  0.220 0.008 0.366 3113.310   2226.557    1
#>       avoid-forget (diff1)  1.254 0.019 0.797 1604.936   1755.915    1
#>     avoid-instruct (diff1) -2.777 0.020 0.958 2319.811   2324.516    1
#>  closeatt-distract (diff1) -0.165 0.007 0.331 3158.751   2070.827    1
#>    closeatt-forget (diff1)  0.156 0.007 0.305 2919.663   1950.611    1
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
#> avoid(c1)    -2.581259 -2.544563
#> closeatt(c1) -2.236258 -2.931231
#> distract(c1) -0.432147 -2.628363
#> forget(c1)   -1.531361 -2.878325
#> instruct(c1) -2.328738 -2.408032
#> 
#> $pairwise_effects_raw
#>                     baseline      diff1
#> avoid-closeatt     0.8218393  0.9719960
#> avoid-distract     1.6359213  0.2203003
#> avoid-forget       0.4725672  1.2537522
#> avoid-instruct     0.3753502 -2.7774496
#> closeatt-distract -0.1904266 -0.1648557
#> closeatt-forget    0.1409541  0.1561317
#> closeatt-instruct  1.4900705  0.5221248
#> distract-forget    0.3684304  0.2420162
#> distract-instruct  1.1891455  1.2797568
#> forget-instruct    1.0677543  0.7527294
#> 
#> $main_effects_groups
#>                   group1    group2
#> avoid(c1)    -1.30897769 -3.853541
#> closeatt(c1) -0.77064212 -3.701874
#> distract(c1)  0.88203470 -1.746329
#> forget(c1)   -0.09219866 -2.970523
#> instruct(c1) -1.12472172 -3.532754
#> 
#> $pairwise_effects_groups
#>                        group1     group2
#> avoid-closeatt     0.33584128  1.3078373
#> avoid-distract     1.52577117  1.7460715
#> avoid-forget      -0.15430893  1.0994432
#> avoid-instruct     1.76407502 -1.0133746
#> closeatt-distract -0.10799877 -0.2728545
#> closeatt-forget    0.06288827  0.2190199
#> closeatt-instruct  1.22900815  1.7511329
#> distract-forget    0.24742228  0.4894385
#> distract-instruct  0.54926715  1.8290239
#> forget-instruct    0.69138956  1.4441190
#> 
#> $indicators
#>           avoid closeatt distract forget instruct
#> avoid    1.0000  0.70300  0.40300 0.8485  0.99550
#> closeatt 0.7030  1.00000  0.37525 0.3655  0.56500
#> distract 0.4030  0.37525  1.00000 0.4185  0.86275
#> forget   0.8485  0.36550  0.41850 1.0000  0.70250
#> instruct 0.9955  0.56500  0.86275 0.7025  1.00000
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
