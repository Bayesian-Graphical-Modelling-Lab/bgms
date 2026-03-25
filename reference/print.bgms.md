# Print method for `bgms` objects

Minimal console output for `bgms` fit objects.

## Usage

``` r
# S3 method for class 'bgms'
print(x, ...)
```

## Arguments

- x:

  An object of class `bgms`.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## See also

[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md),
[`summary.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md),
[`coef.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgms.md)

Other posterior-methods:
[`coef.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgmCompare.md),
[`coef.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgms.md),
[`print.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgmCompare.md),
[`summary.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgmCompare.md),
[`summary.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)

## Examples

``` r
# \donttest{
fit = bgm(x = Wenchuan[, 1:3])
#> 2 rows with missing values excluded (n = 360 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 97/2000 (4.9%)
#> Chain 3 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 153/2000 (7.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 91/2000 (4.5%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 441/8000 (5.5%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2000 (20.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 439/2000 (21.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 489/2000 (24.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 413/2000 (20.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1741/8000 (21.8%)
#> Elapsed: 1s | ETA: 4s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 811/2000 (40.6%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 838/2000 (41.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 805/2000 (40.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 3254/8000 (40.7%)
#> Elapsed: 1s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/2000 (55.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1118/2000 (55.9%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2000 (57.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 1153/2000 (57.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4521/8000 (56.5%)
#> Elapsed: 2s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1400/2000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1407/2000 (70.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1452/2000 (72.6%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1474/2000 (73.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5733/8000 (71.7%)
#> Elapsed: 2s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/2000 (85.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1695/2000 (84.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1735/2000 (86.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1771/2000 (88.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6901/8000 (86.3%)
#> Elapsed: 3s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 3s | ETA: 0s
print(fit)
#> Bayesian Edge Selection using a Bernoulli prior on edge inclusion 
#>  Model: OMRF (Ordinal Markov Random Field)
#>  Number of variables: 3
#>  Number of cases: 360
#>  Number of post-burnin MCMC iterations: 4000
#>  Number of MCMC chains: 4
#> Use the `summary()` function for posterior summaries and chain diagnostics.
#> See the `easybgm` package for summary and plotting tools.
# }
```
