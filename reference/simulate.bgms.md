# Simulate Data from a Fitted bgms Model

Generates new observations from the Markov Random Field model using the
estimated parameters from a fitted `bgms` object.

## Usage

``` r
# S3 method for class 'bgms'
simulate(
  object,
  nsim = 500,
  seed = NULL,
  method = c("posterior-mean", "posterior-sample"),
  ndraws = NULL,
  iter = 1000,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  ...
)
```

## Arguments

- object:

  An object of class `bgms`.

- nsim:

  Number of observations to simulate. Default: `500`.

- seed:

  Optional random seed for reproducibility.

- method:

  Character string specifying which parameter estimates to use:

  `"posterior-mean"`

  :   Use posterior mean parameters (faster, single simulation).

  `"posterior-sample"`

  :   Sample from posterior draws, producing one dataset per draw
      (accounts for parameter uncertainty). This method uses parallel
      processing when `cores > 1`.

- ndraws:

  Number of posterior draws to use when `method = "posterior-sample"`.
  If `NULL`, uses all available draws.

- iter:

  Number of Gibbs iterations for equilibration before collecting
  samples. Default: `1000`.

- cores:

  Number of CPU cores for parallel execution when
  `method = "posterior-sample"`. Default:
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).

- display_progress:

  Character string specifying the type of progress bar. Options:
  `"per-chain"`, `"total"`, `"none"`. Default: `"per-chain"`.

- ...:

  Additional arguments (currently ignored).

## Value

If `method = "posterior-mean"`: A matrix with `nsim` rows and `p`
columns containing simulated observations.

If `method = "posterior-sample"`: A list of matrices, one per posterior
draw, each with `nsim` rows and `p` columns.

## Details

This function uses the estimated interaction and threshold parameters to
generate new data via Gibbs sampling. When
`method = "posterior-sample"`, parameter uncertainty is propagated to
the simulated data by using different posterior draws. Parallel
processing is available for this method via the `cores` argument.

## See also

[`predict.bgms`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgms.md)
for computing conditional probabilities,
[`simulate_mrf`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate_mrf.md)
for simulation with user-specified parameters.

## Examples

``` r
# \donttest{
# Fit a model
fit <- bgm(x = Wenchuan[, 1:5])
#> 7 rows with missing values excluded (n = 355 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2000 (2.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 48/2000 (2.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 37/2000 (1.8%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 41/2000 (2.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 176/8000 (2.2%)
#> Elapsed: 1s | ETA: 44s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 86/2000 (4.3%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 78/2000 (3.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 94/2000 (4.7%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 358/8000 (4.5%)
#> Elapsed: 2s | ETA: 43s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/2000 (10.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 136/2000 (6.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 94/2000 (4.7%)
#> Chain 4 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 158/2000 (7.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 588/8000 (7.3%)
#> Elapsed: 3s | ETA: 38s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 350/2000 (17.5%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 267/2000 (13.4%)
#> Chain 3 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 167/2000 (8.3%)
#> Chain 4 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 273/2000 (13.7%)
#> Total   (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1057/8000 (13.2%)
#> Elapsed: 4s | ETA: 26s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2000 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 396/2000 (19.8%)
#> Chain 3 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 308/2000 (15.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 434/2000 (21.7%)
#> Total   (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1638/8000 (20.5%)
#> Elapsed: 4s | ETA: 16s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 650/2000 (32.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 530/2000 (26.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 449/2000 (22.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 578/2000 (28.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2207/8000 (27.6%)
#> Elapsed: 5s | ETA: 13s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 694/2000 (34.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 603/2000 (30.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 750/2000 (37.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2847/8000 (35.6%)
#> Elapsed: 6s | ETA: 11s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2000 (45.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 817/2000 (40.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 732/2000 (36.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 857/2000 (42.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3306/8000 (41.3%)
#> Elapsed: 6s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1050/2000 (52.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 929/2000 (46.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 871/2000 (43.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 986/2000 (49.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 3836/8000 (47.9%)
#> Elapsed: 7s | ETA: 8s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1200/2000 (60.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 1070/2000 (53.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1000/2000 (50.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1130/2000 (56.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4400/8000 (55.0%)
#> Elapsed: 8s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1350/2000 (67.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 1210/2000 (60.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1129/2000 (56.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1270/2000 (63.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4959/8000 (62.0%)
#> Elapsed: 8s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1500/2000 (75.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1329/2000 (66.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1239/2000 (62.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 1406/2000 (70.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 5474/8000 (68.4%)
#> Elapsed: 9s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2000 (82.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1476/2000 (73.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1355/2000 (67.8%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1540/2000 (77.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 6021/8000 (75.3%)
#> Elapsed: 10s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/2000 (90.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 1605/2000 (80.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1483/2000 (74.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1684/2000 (84.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6572/8000 (82.2%)
#> Elapsed: 10s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1950/2000 (97.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1743/2000 (87.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1628/2000 (81.4%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1843/2000 (92.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7164/8000 (89.5%)
#> Elapsed: 11s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 13s | ETA: 0s

# Simulate 100 new observations using posterior means
new_data <- simulate(fit, nsim = 100)

# Simulate with parameter uncertainty (10 datasets)
new_data_list <- simulate(fit, nsim = 100, method = "posterior-sample", ndraws = 10)
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 10/10 (100.0%)
#> Elapsed: 0s | ETA: 0s

# Use parallel processing for faster simulation
new_data_list <- simulate(fit, nsim = 100, method = "posterior-sample",
                          ndraws = 100, cores = 4)
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/100 (50.0%)
#> Elapsed: 0s | ETA: 0s
#> 
#> User interrupt detected. Exiting gracefully. It may take a few seconds before all chains are terminated.
#> All chains terminated.
# }
```
