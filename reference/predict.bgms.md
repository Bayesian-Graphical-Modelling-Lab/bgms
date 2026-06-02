# Predict Conditional Probabilities from a Fitted bgms Model

Computes conditional probability distributions for one or more variables
given the observed values of other variables in the data. Supports
ordinal, Blume-Capel, continuous (GGM), and mixed MRF models.

## Usage

``` r
# S3 method for class 'bgms'
predict(
  object,
  newdata,
  variables = NULL,
  type = c("probabilities", "response"),
  method = c("posterior-mean", "posterior-sample"),
  ndraws = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `bgms`.

- newdata:

  A matrix or data frame with `n` rows and `p` columns containing the
  observed data. Must have the same variables (columns) as the original
  data used to fit the model.

- variables:

  Which variables to predict. Can be:

  - A character vector of variable names

  - An integer vector of column indices

  - `NULL` (default) to predict all variables

- type:

  Character string specifying the type of prediction:

  `"probabilities"`

  :   Return the full conditional probability distribution for each
      variable and observation.

  `"response"`

  :   Return the predicted category (mode of the conditional
      distribution).

- method:

  Character string specifying which parameter estimates to use:

  `"posterior-mean"`

  :   Use posterior mean parameters.

  `"posterior-sample"`

  :   Average predictions over posterior draws.

- ndraws:

  Number of posterior draws to use when `method = "posterior-sample"`.
  If `NULL`, uses all available draws.

- seed:

  Optional random seed for reproducibility when
  `method = "posterior-sample"`.

- ...:

  Additional arguments (currently ignored).

## Value

**Ordinal models:**

For `type = "probabilities"`: A named list with one element per
predicted variable. Each element is a matrix with `n` rows and
`num_categories + 1` columns containing \\P(X_j = c \| X\_{-j})\\ for
each observation and category.

For `type = "response"`: A matrix with `n` rows and `length(variables)`
columns containing predicted categories.

When `method = "posterior-sample"`, probabilities are averaged over
posterior draws, and an attribute `"sd"` is included containing the
standard deviation across draws.

**GGM (continuous) models:**

For `type = "probabilities"`: A named list with one element per
predicted variable. Each element is a matrix with `n` rows and 2 columns
(`"mean"` and `"sd"`) containing the conditional Gaussian parameters
\\E(X_j \| X\_{-j})\\ and \\\text{SD}(X_j \| X\_{-j})\\.

For `type = "response"`: A matrix with `n` rows and `length(variables)`
columns containing conditional means.

When `method = "posterior-sample"`, conditional parameters are averaged
over posterior draws, and an attribute `"sd"` is included.

**Mixed MRF models:**

For mixed models, the return list contains elements for both discrete
and continuous predicted variables. Discrete variables return
probability matrices (as in ordinal models); continuous variables return
conditional mean and SD matrices (as in GGM models).

## Details

For each observation, the function computes the conditional distribution
of the target variable(s) given the observed values of all other
variables. This is the same conditional distribution used internally by
the Gibbs sampler.

For GGM (continuous) models, the conditional distribution of \\X_j \|
X\_{-j}\\ is Gaussian with mean \\-\omega\_{jj}^{-1} \sum\_{k \neq j}
\omega\_{jk} x_k\\ and variance \\\omega\_{jj}^{-1}\\, where \\\Omega\\
is the precision matrix.

## See also

[`simulate.bgms`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgms.md)
for generating new data from the model.

Other prediction:
[`predict.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/predict.bgmCompare.md),
[`simulate.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgmCompare.md),
[`simulate.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgms.md),
[`simulate_mrf()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate_mrf.md)

## Examples

``` r
# \donttest{
# Fit a model
fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#> 7 rows with missing values excluded (n = 355 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/4000 (1.2%)
#> Chain 2 (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 42/4000 (1.1%)
#> Total   (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 92/8000 (1.1%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/4000 (2.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/4000 (2.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/8000 (2.5%)
#> Elapsed: 2s | ETA: 1m 18s
#> Chain 1 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/4000 (6.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 155/4000 (3.9%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 405/8000 (5.1%)
#> Elapsed: 2s | ETA: 38s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/4000 (7.5%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 322/4000 (8.1%)
#> Total   (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 622/8000 (7.8%)
#> Elapsed: 3s | ETA: 36s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/4000 (12.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 453/4000 (11.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 953/8000 (11.9%)
#> Elapsed: 4s | ETA: 30s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/4000 (17.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 462/4000 (11.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1162/8000 (14.5%)
#> Elapsed: 4s | ETA: 24s
#> Chain 1 (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 950/4000 (23.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 661/4000 (16.5%)
#> Total   (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1611/8000 (20.1%)
#> Elapsed: 5s | ETA: 20s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1200/4000 (30.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 852/4000 (21.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2052/8000 (25.7%)
#> Elapsed: 6s | ETA: 17s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1450/4000 (36.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 861/4000 (21.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2311/8000 (28.9%)
#> Elapsed: 6s | ETA: 15s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/4000 (42.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1032/4000 (25.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2732/8000 (34.2%)
#> Elapsed: 7s | ETA: 13s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/4000 (47.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1279/4000 (32.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3179/8000 (39.7%)
#> Elapsed: 7s | ETA: 11s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 2050/4000 (51.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1491/4000 (37.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3541/8000 (44.3%)
#> Elapsed: 8s | ETA: 10s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 2250/4000 (56.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1654/4000 (41.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3904/8000 (48.8%)
#> Elapsed: 9s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 2450/4000 (61.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1664/4000 (41.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4114/8000 (51.4%)
#> Elapsed: 9s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 2650/4000 (66.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1807/4000 (45.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 4457/8000 (55.7%)
#> Elapsed: 10s | ETA: 8s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 2850/4000 (71.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1955/4000 (48.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 4805/8000 (60.1%)
#> Elapsed: 10s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 3050/4000 (76.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2161/4000 (54.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 5211/8000 (65.1%)
#> Elapsed: 11s | ETA: 6s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 3250/4000 (81.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2385/4000 (59.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 5635/8000 (70.4%)
#> Elapsed: 12s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 3450/4000 (86.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2587/4000 (64.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 6037/8000 (75.5%)
#> Elapsed: 12s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3600/4000 (90.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2768/4000 (69.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6368/8000 (79.6%)
#> Elapsed: 13s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3800/4000 (95.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2992/4000 (74.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6792/8000 (84.9%)
#> Elapsed: 13s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 3228/4000 (80.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 7228/8000 (90.3%)
#> Elapsed: 14s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 16s | ETA: 0s

# Compute conditional probabilities for all variables
probs = predict(fit, newdata = Wenchuan[1:10, 1:5])

# Predict the first variable only
probs_v1 = predict(fit, newdata = Wenchuan[1:10, 1:5], variables = 1)

# Get predicted categories
pred_class = predict(fit, newdata = Wenchuan[1:10, 1:5], type = "response")
# }
```
