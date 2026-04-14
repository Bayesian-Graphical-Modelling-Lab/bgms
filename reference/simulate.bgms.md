# Simulate Data from a Fitted bgms Model

Generates new observations from the Markov Random Field model using the
estimated parameters from a fitted `bgms` object. Supports ordinal,
Blume-Capel, continuous (GGM), and mixed MRF models.

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

For mixed MRF models, discrete columns contain non-negative integers and
continuous columns contain real-valued observations, ordered as in the
original data.

## Details

This function uses the estimated interaction and threshold parameters to
generate new data via Gibbs sampling. When
`method = "posterior-sample"`, parameter uncertainty is parameter
uncertainty is propagated to the simulated data by using different
posterior draws. Parallel processing is available for this method via
the `cores` argument.

## See also

[`predict.bgms`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgms.md)
for computing conditional probabilities,
[`simulate_mrf`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate_mrf.md)
for simulation with user-specified parameters.

Other prediction:
[`predict.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgmCompare.md),
[`predict.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgms.md),
[`simulate.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgmCompare.md),
[`simulate_mrf()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate_mrf.md)

## Examples

``` r
# \donttest{
# Fit a model
fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#> 7 rows with missing values excluded (n = 355 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 79/2000 (4.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 179/4000 (4.5%)
#> Elapsed: 1s | ETA: 21s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/2000 (10.0%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 161/2000 (8.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 361/4000 (9.0%)
#> Elapsed: 2s | ETA: 20s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2000 (15.0%)
#> Chain 2 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 311/2000 (15.6%)
#> Total   (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 611/4000 (15.3%)
#> Elapsed: 3s | ETA: 17s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 600/2000 (30.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 607/2000 (30.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1207/4000 (30.2%)
#> Elapsed: 3s | ETA: 7s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2000 (45.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 918/2000 (45.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1818/4000 (45.5%)
#> Elapsed: 4s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2000 (57.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1186/2000 (59.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 2336/4000 (58.4%)
#> Elapsed: 4s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1400/2000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1442/2000 (72.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 2842/4000 (71.0%)
#> Elapsed: 5s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2000 (82.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 1704/2000 (85.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3354/4000 (83.9%)
#> Elapsed: 6s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2000 (95.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 1969/2000 (98.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3869/4000 (96.7%)
#> Elapsed: 6s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Elapsed: 6s | ETA: 0s

# Simulate 100 new observations using posterior means
new_data = simulate(fit, nsim = 100)

# Simulate with parameter uncertainty (10 datasets)
new_data_list = simulate(
  fit,
  nsim = 100,
  method = "posterior-sample", ndraws = 10
)
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 10/10 (100.0%)
#> Elapsed: 0s | ETA: 0s

# Use parallel processing for faster simulation
new_data_list = simulate(fit,
  nsim = 100, method = "posterior-sample",
  ndraws = 100, cores = 2
)
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/100 (100.0%)
#> Elapsed: 1s | ETA: 0s
# }
```
