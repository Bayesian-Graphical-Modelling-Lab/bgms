# Extract Coefficients from a bgms Object

Returns the posterior mean main effects, pairwise effects, and edge
inclusion indicators from a `bgms` model fit.

## Usage

``` r
# S3 method for class 'bgms'
coef(object, ...)
```

## Arguments

- object:

  An object of class `bgms`.

- ...:

  Ignored.

## Value

A list with the following components:

- main:

  Posterior mean of the main-effect parameters. `NULL` for GGM models
  (no main effects). For OMRF models this is a numeric matrix (p x
  max_categories) of category thresholds. For mixed MRF models this is a
  list with `$discrete` (p x max_categories matrix) and `$continuous` (q
  x 1 matrix of means).

- pairwise:

  Posterior mean of the partial association matrix (zero diagonal). Use
  [`extract_precision()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_precision.md)
  for the full precision matrix.

- indicator:

  Posterior mean of the edge inclusion indicators (if available).

## See also

[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md),
[`print.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgms.md),
[`summary.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)

Other posterior-methods:
[`coef.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgmCompare.md),
[`print.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgmCompare.md),
[`print.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgms.md),
[`summary.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgmCompare.md),
[`summary.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)

## Examples

``` r
# \donttest{
fit = bgm(x = Wenchuan[, 1:3])
#> 2 rows with missing values excluded (n = 360 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 130/2000 (6.5%)
#> Chain 3 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 151/2000 (7.5%)
#> Chain 4 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 112/2000 (5.6%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 493/8000 (6.2%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2000 (20.0%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 254/2000 (12.7%)
#> Chain 3 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 158/2000 (7.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 415/2000 (20.8%)
#> Total   (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1227/8000 (15.3%)
#> Elapsed: 1s | ETA: 6s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2000 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 453/2000 (22.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 454/2000 (22.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 1024/2000 (51.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2431/8000 (30.4%)
#> Elapsed: 2s | ETA: 5s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/2000 (42.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 945/2000 (47.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 1019/2000 (50.9%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1828/2000 (91.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 4642/8000 (58.0%)
#> Elapsed: 3s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1250/2000 (62.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1215/2000 (60.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1323/2000 (66.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5788/8000 (72.4%)
#> Elapsed: 4s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/2000 (90.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1517/2000 (75.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1629/2000 (81.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6946/8000 (86.8%)
#> Elapsed: 4s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 5s | ETA: 0s
coef(fit)
#> $main
#>               cat (1)   cat (2)   cat (3)    cat (4)
#> intrusion  0.81851546 -1.163735 -3.518419  -7.412488
#> dreams    -0.40277245 -3.297814 -6.362015 -10.458782
#> flash      0.07455844 -2.124759 -4.570361  -8.404449
#> 
#> $pairwise
#>           intrusion    dreams     flash
#> intrusion 0.0000000 0.3500149 0.2123476
#> dreams    0.3500149 0.0000000 0.2806789
#> flash     0.2123476 0.2806789 0.0000000
#> 
#> $indicator
#>           intrusion dreams flash
#> intrusion         0      1     1
#> dreams            1      0     1
#> flash             1      1     0
#> 
# }
```
