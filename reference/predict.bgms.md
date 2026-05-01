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
#> Chain 2 (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 43/4000 (1.1%)
#> Total   (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 93/8000 (1.2%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/4000 (2.5%)
#> Chain 2 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 102/4000 (2.5%)
#> Total   (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 202/8000 (2.5%)
#> Elapsed: 2s | ETA: 1m 17s
#> Chain 1 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/4000 (6.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 156/4000 (3.9%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 406/8000 (5.1%)
#> Elapsed: 3s | ETA: 56s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/4000 (7.5%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 350/4000 (8.8%)
#> Total   (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 650/8000 (8.1%)
#> Elapsed: 4s | ETA: 45s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/4000 (12.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 454/4000 (11.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 954/8000 (11.9%)
#> Elapsed: 4s | ETA: 30s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/4000 (17.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 465/4000 (11.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1165/8000 (14.6%)
#> Elapsed: 5s | ETA: 29s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/4000 (22.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 656/4000 (16.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1556/8000 (19.4%)
#> Elapsed: 5s | ETA: 21s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/4000 (27.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/4000 (21.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1950/8000 (24.4%)
#> Elapsed: 6s | ETA: 19s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1350/4000 (33.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 860/4000 (21.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2210/8000 (27.6%)
#> Elapsed: 6s | ETA: 16s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1600/4000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 965/4000 (24.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2565/8000 (32.1%)
#> Elapsed: 7s | ETA: 15s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/4000 (45.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1211/4000 (30.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3011/8000 (37.6%)
#> Elapsed: 8s | ETA: 13s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/4000 (50.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1459/4000 (36.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3459/8000 (43.2%)
#> Elapsed: 8s | ETA: 11s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 2150/4000 (53.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/4000 (41.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3800/8000 (47.5%)
#> Elapsed: 9s | ETA: 10s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━⦘ 2350/4000 (58.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1661/4000 (41.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 4011/8000 (50.1%)
#> Elapsed: 9s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2500/4000 (62.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 1742/4000 (43.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 4242/8000 (53.0%)
#> Elapsed: 10s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 2650/4000 (66.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1852/4000 (46.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4502/8000 (56.3%)
#> Elapsed: 10s | ETA: 8s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2800/4000 (70.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1997/4000 (49.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4797/8000 (60.0%)
#> Elapsed: 11s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 2950/4000 (73.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2180/4000 (54.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5130/8000 (64.1%)
#> Elapsed: 12s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 3150/4000 (78.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━⦘ 2402/4000 (60.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5552/8000 (69.4%)
#> Elapsed: 12s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3300/4000 (82.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2583/4000 (64.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 5883/8000 (73.5%)
#> Elapsed: 13s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 3450/4000 (86.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2758/4000 (69.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 6208/8000 (77.6%)
#> Elapsed: 13s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3600/4000 (90.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 2940/4000 (73.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6540/8000 (81.8%)
#> Elapsed: 14s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 3750/4000 (93.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 3121/4000 (78.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 6871/8000 (85.9%)
#> Elapsed: 14s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3900/4000 (97.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 3301/4000 (82.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 7201/8000 (90.0%)
#> Elapsed: 15s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 17s | ETA: 0s

# Compute conditional probabilities for all variables
probs = predict(fit, newdata = Wenchuan[1:10, 1:5])

# Predict the first variable only
probs_v1 = predict(fit, newdata = Wenchuan[1:10, 1:5], variables = 1)

# Get predicted categories
pred_class = predict(fit, newdata = Wenchuan[1:10, 1:5], type = "response")
# }
```
