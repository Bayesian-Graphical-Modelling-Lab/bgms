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
fit <- bgm(x = Wenchuan[, 1:5], chains = 2)
#> 7 rows with missing values excluded (n = 355 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 150/2000 (7.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 239/2000 (11.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 389/4000 (9.7%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/2000 (22.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 565/2000 (28.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1015/4000 (25.4%)
#> Elapsed: 1s | ETA: 3s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 924/2000 (46.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 1724/4000 (43.1%)
#> Elapsed: 1s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/2000 (55.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1231/2000 (61.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 2331/4000 (58.3%)
#> Elapsed: 2s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1400/2000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1543/2000 (77.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 2943/4000 (73.6%)
#> Elapsed: 2s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/2000 (85.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1847/2000 (92.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 3547/4000 (88.7%)
#> Elapsed: 3s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Elapsed: 3s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Elapsed: 3s | ETA: 0s
#> NUTS issues:
#>   - Divergences: 1 (0.050%) - check R-hat and ESS 

# Simulate 100 new observations using posterior means
new_data <- simulate(fit, nsim = 100)

# Simulate with parameter uncertainty (10 datasets)
new_data_list <- simulate(fit, nsim = 100, method = "posterior-sample", ndraws = 10)
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 10/10 (100.0%)
#> Elapsed: 0s | ETA: 0s

# Use parallel processing for faster simulation
new_data_list <- simulate(fit, nsim = 100, method = "posterior-sample",
                          ndraws = 100, cores = 2)
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/100 (50.0%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/100 (100.0%)
#> Elapsed: 1s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/100 (100.0%)
#> Elapsed: 1s | ETA: 0s
# }
```
