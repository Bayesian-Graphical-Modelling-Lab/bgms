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
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 195/2000 (9.8%)
#> Chain 3 (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 201/2000 (10.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 143/2000 (7.1%)
#> Total   (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 639/8000 (8.0%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/2000 (22.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 561/2000 (28.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 575/2000 (28.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 509/2000 (25.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2094/8000 (26.2%)
#> Elapsed: 1s | ETA: 3s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 921/2000 (46.1%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 955/2000 (47.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 897/2000 (44.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3573/8000 (44.7%)
#> Elapsed: 1s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2000 (57.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1284/2000 (64.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1317/2000 (65.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1267/2000 (63.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 5018/8000 (62.7%)
#> Elapsed: 2s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1500/2000 (75.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 1651/2000 (82.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1686/2000 (84.3%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1642/2000 (82.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 6479/8000 (81.0%)
#> Elapsed: 3s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1850/2000 (92.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 7850/8000 (98.1%)
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
