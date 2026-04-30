# Bayesian Estimation or Edge Selection for Markov Random Fields

The `bgm` function estimates the pseudoposterior distribution of the
parameters of a Markov Random Field (MRF) for binary, ordinal,
continuous, or mixed (discrete and continuous) variables. Depending on
the variable types, the model is an ordinal MRF, a Gaussian graphical
model (GGM), or a mixed MRF. Optionally, it performs Bayesian edge
selection using spike-and-slab priors to infer the network structure.

## Usage

``` r
bgm(
  x,
  variable_type = "ordinal",
  baseline_category,
  iter = 1000,
  warmup = 1000,
  interaction_prior = cauchy_prior(scale = 1),
  threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
  means_prior = normal_prior(scale = 1),
  precision_scale_prior = gamma_prior(shape = 1, rate = 1),
  edge_selection = TRUE,
  edge_prior = bernoulli_prior(0.5),
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis"),
  target_accept,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  standardize = FALSE,
  verbose = getOption("bgms.verbose", TRUE),
  progress_callback = NULL,
  pairwise_scale,
  main_alpha,
  main_beta,
  inclusion_probability,
  beta_bernoulli_alpha,
  beta_bernoulli_beta,
  beta_bernoulli_alpha_between,
  beta_bernoulli_beta_between,
  dirichlet_alpha,
  lambda,
  interaction_scale,
  burnin,
  save,
  threshold_alpha,
  threshold_beta
)
```

## Arguments

- x:

  A data frame or matrix with `n` rows and `p` columns. Columns may
  contain binary, ordinal, or continuous variables (see
  `variable_type`). Discrete variables are automatically recoded to
  non-negative integers (`0, 1, ..., m`); for regular ordinal variables,
  unobserved categories are collapsed, while Blume–Capel variables
  retain all categories. Continuous variables are column-centered
  internally so that the GGM likelihood is formulated with a zero-mean
  assumption.

- variable_type:

  Character or character vector. Specifies the type of each variable in
  `x`. Allowed values: `"ordinal"`, `"blume-capel"`, or `"continuous"`.
  A single string applies to all variables. A per-variable vector that
  mixes discrete (`"ordinal"` / `"blume-capel"`) and `"continuous"`
  types fits a mixed MRF. Binary variables are automatically treated as
  `"ordinal"`. Default: `"ordinal"`.

- baseline_category:

  Integer or vector. Baseline category used in Blume–Capel variables.
  Can be a single integer (applied to all) or a vector of length `p`.
  Required if at least one variable is of type `"blume-capel"`.

- iter:

  Integer. Number of post–burn-in iterations (per chain). Default:
  `1e3`.

- warmup:

  Integer. Number of warmup iterations before collecting samples. A
  minimum of 1000 iterations is enforced, with a warning if a smaller
  value is requested. Default: `1e3`.

- interaction_prior:

  A prior specification object for pairwise interaction parameters,
  created by one of the prior constructor functions:

  - [`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md):
    Cauchy(0, scale) prior (default).

  - [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md):
    Normal(0, scale) prior.

  - [`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md):
    Beta-prime prior.

  Default: `cauchy_prior(scale = 1)`.

- threshold_prior:

  A prior specification object for threshold (main effect) parameters,
  created by one of the prior constructor functions:

  - [`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md):
    Beta-prime prior (default).

  - [`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md):
    Cauchy(0, scale) prior.

  - [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md):
    Normal(0, scale) prior.

  Default: `beta_prime_prior(alpha = 0.5, beta = 0.5)`.

- means_prior:

  A prior specification object for continuous variable means (mixed MRF
  models only), created by one of the prior constructor functions:

  - [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md):
    Normal(0, scale) prior (default).

  - [`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md):
    Cauchy(0, scale) prior.

  - [`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md):
    Beta-prime prior.

  Only used when the model includes continuous variables. Ignored for
  pure ordinal or pure continuous (GGM) models. Default:
  `normal_prior(scale = 1)`.

- precision_scale_prior:

  A prior specification object for the diagonal elements of the
  precision matrix, created by one of:

  - [`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md):
    Gamma(shape, rate) prior (default).

  - [`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md):
    Exponential(rate) prior.

  Only used for models with continuous variables (GGM and mixed MRF).
  Ignored for pure ordinal models. Default:
  `gamma_prior(shape = 1, rate = 1)`.

- edge_selection:

  Logical. Whether to perform Bayesian edge selection. If `FALSE`, the
  model estimates all edges. Default: `TRUE`.

- edge_prior:

  An edge prior specification object, or a character string
  (deprecated). Specifies the prior for edge inclusion. Preferred: pass
  an object from one of:

  - [`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md):
    Fixed inclusion probability (default).

  - [`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md):
    Beta-distributed inclusion.

  - [`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md):
    Stochastic Block Model.

  Legacy character strings `"Bernoulli"`, `"Beta-Bernoulli"`,
  `"Stochastic-Block"` are still accepted but deprecated. Default:
  `bernoulli_prior(0.5)`.

- na_action:

  Character. Specifies missing data handling. Either `"listwise"` (drop
  rows with missing values) or `"impute"` (perform single imputation
  during sampling). Default: `"listwise"`.

- update_method:

  Character. Specifies how the MCMC sampler updates the model
  parameters:

  "adaptive-metropolis"

  :   Componentwise adaptive Metropolis–Hastings with Robbins–Monro
      proposal adaptation.

  "nuts"

  :   The No-U-Turn Sampler with RATTLE constrained integration for
      Gaussian models with edge selection.

  Default: `"nuts"`.

- target_accept:

  Numeric between 0 and 1. Target acceptance rate for the sampler.
  Defaults are set automatically if not supplied: `0.44` for adaptive
  Metropolis and `0.80` for NUTS.

- nuts_max_depth:

  Integer. Maximum tree depth in NUTS. Must be positive. Default: `10`.

- learn_mass_matrix:

  Logical. If `TRUE`, adapt a diagonal mass matrix during warmup (NUTS
  only). If `FALSE`, use the identity matrix. Default: `TRUE`.

- chains:

  Integer. Number of parallel chains to run. Default: `4`.

- cores:

  Integer. Number of CPU cores for parallel execution. Default:
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).

- display_progress:

  Character. Controls progress reporting during sampling. Options:
  `"per-chain"` (separate bar per chain), `"total"` (single combined
  bar), or `"none"` (no progress). Default: `"per-chain"`.

- seed:

  Optional integer. Random seed for reproducibility. Must be a single
  non-negative integer.

- standardize:

  Logical. If `TRUE`, the prior scale for each pairwise interaction is
  adjusted based on the range of response scores. Variables with more
  response categories have larger score products \\x_i \cdot x_j\\,
  which typically correspond to smaller interaction effects
  \\\sigma\_{ij}\\. Without standardization, a fixed prior scale is
  relatively wide for these smaller effects, resulting in less shrinkage
  for high-category pairs and more shrinkage for low-category pairs.
  Standardization scales the prior proportionally to the maximum score
  product, ensuring equivalent relative shrinkage across all pairs.
  After internal recoding, regular ordinal variables have scores \\0, 1,
  \ldots, m\\. The adjusted scale for the interaction between variables
  \\i\\ and \\j\\ is `pairwise_scale * m_i * m_j`, so that
  `pairwise_scale` itself applies to the unit interval case (binary
  variables where \\m_i = m_j = 1\\). For Blume-Capel variables with
  reference category \\b\\, scores are centered as \\-b, \ldots, m-b\\,
  and the adjustment uses the maximum absolute product of the score
  endpoints. For mixed pairs, ordinal variables use raw score endpoints
  \\(0, m)\\ and Blume-Capel variables use centered score endpoints
  \\(-b, m-b)\\. Default: `FALSE`.

- verbose:

  Logical. If `TRUE`, prints informational messages during data
  processing (e.g., missing data handling, variable recoding). Defaults
  to `getOption("bgms.verbose", TRUE)`. Set
  `options(bgms.verbose = FALSE)` to suppress messages globally.

- progress_callback:

  An optional R function with signature `function(completed, total)`
  that is called at regular intervals during sampling, where `completed`
  is the number of iterations completed across all chains and `total` is
  the total number of iterations. Useful for external front-ends (e.g.,
  JASP) that supply their own progress reporting. When `NULL` (the
  default), no callback is invoked.

- pairwise_scale:

  **\[deprecated\]** Double. Scale of the Cauchy prior for pairwise
  interaction parameters. Use `interaction_prior` instead. Default: `1`.

- main_alpha, main_beta:

  **\[deprecated\]** Double. Shape parameters of the beta-prime prior
  for threshold parameters. Use `threshold_prior` instead. Defaults:
  `main_alpha = 0.5` and `main_beta = 0.5`.

- inclusion_probability:

  **\[deprecated\]** Numeric scalar. Use
  `edge_prior = bernoulli_prior(inclusion_probability)` instead.
  Default: `0.5`.

- beta_bernoulli_alpha, beta_bernoulli_beta:

  **\[deprecated\]** Double. Use
  `edge_prior = beta_bernoulli_prior(alpha, beta)` instead. Defaults:
  `1`.

- beta_bernoulli_alpha_between, beta_bernoulli_beta_between:

  **\[deprecated\]** Double. Use
  `edge_prior = sbm_prior(alpha_between, beta_between)` instead.
  Defaults: `1`.

- dirichlet_alpha:

  **\[deprecated\]** Double. Use
  `edge_prior = sbm_prior(dirichlet_alpha = ...)` instead. Default: `1`.

- lambda:

  **\[deprecated\]** Double. Use `edge_prior = sbm_prior(lambda = ...)`
  instead. Default: `1`.

- interaction_scale, burnin, save, threshold_alpha, threshold_beta:

  **\[deprecated\]** Deprecated arguments as of **bgms 0.1.6.0**. Use
  `pairwise_scale`, `warmup`, `main_alpha`, and `main_beta` instead.

## Value

A list of class `"bgms"` with posterior summaries, posterior mean
matrices, and access to raw MCMC draws. The object can be passed to
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`coef()`](https://rdrr.io/r/stats/coef.html).

Main components include:

- `posterior_summary_main`: Data frame with posterior summaries (mean,
  sd, MCSE, ESS, Rhat) for main-effect parameters. For OMRF models these
  are category thresholds; for mixed MRF models these are discrete
  thresholds and continuous means. `NULL` for GGM models (no main
  effects).

- `posterior_summary_quadratic`: Data frame with posterior summaries for
  the residual variance parameters (GGM and mixed MRF). `NULL` for OMRF
  models.

- `posterior_summary_pairwise`: Data frame with posterior summaries for
  partial association parameters.

- `posterior_summary_indicator`: Data frame with posterior summaries for
  edge inclusion indicators (if `edge_selection = TRUE`).

- `posterior_mean_main`: Posterior mean of main-effect parameters.
  `NULL` for GGM models. For OMRF: a matrix (p x max_categories) of
  category thresholds. For mixed MRF: a list with `$discrete` (threshold
  matrix) and `$continuous` (q x 1 matrix of means).

- `posterior_mean_pairwise`: Symmetric matrix of posterior mean partial
  associations (zero diagonal). For continuous variables these are
  unstandardized partial correlations; for discrete variables these are
  half the log adjacent-category odds ratio. Use
  [`extract_precision()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_precision.md),
  [`extract_partial_correlations()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_partial_correlations.md),
  or
  [`extract_log_odds()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/extract_log_odds.md)
  to convert to interpretable scales.

- `posterior_mean_residual_variance`: Named numeric vector of posterior
  mean residual variances \\1/\Theta\_{ii}\\. Present for GGM and mixed
  MRF models; `NULL` for OMRF.

- `posterior_mean_indicator`: Symmetric matrix of posterior mean
  inclusion probabilities (if edge selection was enabled).

- Additional summaries returned when `edge_prior = "Stochastic-Block"`.
  For more details about this prior see Sekulovski et al. (2025) .

  - `posterior_summary_pairwise_allocations`: Data frame with posterior
    summaries (mean, sd, MCSE, ESS, Rhat) for the pairwise cluster
    co-occurrence of the nodes. This serves to indicate whether the
    estimated posterior allocations,co-clustering matrix and posterior
    cluster probabilities (see blow) have converged.

  - `posterior_coclustering_matrix`: a symmetric matrix of pairwise
    proportions of occurrence of every variable. This matrix can be
    plotted to visually inspect the estimated number of clusters and
    visually inspect nodes that tend to switch clusters.

  - `posterior_mean_allocations`: A vector with the posterior mean of
    the cluster allocations of the nodes. This is calculated using the
    method proposed in Dahl (2009) .

  - `posterior_mode_allocations`: A vector with the posterior mode of
    the cluster allocations of the nodes.

  - `posterior_num_blocks`: A data frame with the estimated posterior
    inclusion probabilities for all the possible number of clusters.

- `raw_samples`: A list of raw MCMC draws per chain:

  - `main`:

    List of main effect samples.

  - `pairwise`:

    List of pairwise effect samples.

  - `indicator`:

    List of indicator samples (if edge selection enabled).

  - `allocations`:

    List of cluster allocations (if SBM prior used).

  - `nchains`:

    Number of chains.

  - `niter`:

    Number of post–warmup iterations per chain.

  - `parameter_names`:

    Named lists of parameter labels.

- `arguments`: A list of function call arguments and metadata (e.g.,
  number of variables, warmup, sampler settings, package version).

The [`summary()`](https://rdrr.io/r/base/summary.html) method prints
formatted posterior summaries, and
[`coef()`](https://rdrr.io/r/stats/coef.html) extracts posterior mean
matrices.

NUTS diagnostics (tree depth, divergences, energy, E-BFMI) are included
in `fit$nuts_diag` if `update_method = "nuts"`.

## Details

Depending on the variable types, the model is an ordinal MRF, a Gaussian
graphical model (GGM), or a mixed MRF. Both regular ordinal variables
and Blume–Capel ordinal variables (with a baseline category) are
supported.

Edge selection uses spike-and-slab priors with Bernoulli,
Beta-Bernoulli, or Stochastic-Block inclusion priors. Parameters are
sampled with NUTS (default) or adaptive Metropolis–Hastings, with a
multi-stage warmup schedule. Missing data can be handled via listwise
deletion or Gibbs imputation.

For full details on model specification, prior choices, warmup, and
output interpretation, see the package website at
<https://bayesian-graphical-modelling-lab.github.io/bgms-docs/>.

## References

Dahl DB (2009). “Modal clustering in a class of product partition
models.” *Bayesian Analysis*, **4**(2), 243–264.
[doi:10.1214/09-BA409](https://doi.org/10.1214/09-BA409) .  
  
Sekulovski N, Arena G, Haslbeck JMB, Huth KBS, Friel N, Marsman M
(2025). “A Stochastic Block Prior for Clustering in Graphical Models.”
*Retrieved from <https://osf.io/preprints/psyarxiv/29p3m_v1>*. OSF
preprint.

## See also

[`vignette("intro", package = "bgms")`](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/intro.md)
for a worked example.

Other model-fitting:
[`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)

## Examples

``` r
# \donttest{
# Run bgm on subset of the Wenchuan dataset
fit = bgm(x = Wenchuan[, 1:5], chains = 2)
#> 7 rows with missing values excluded (n = 355 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2000 (2.5%)
#> Chain 2 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 51/2000 (2.5%)
#> Total   (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 101/4000 (2.5%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 144/2000 (7.2%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 244/4000 (6.1%)
#> Elapsed: 2s | ETA: 31s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2000 (15.0%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 255/2000 (12.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 555/4000 (13.9%)
#> Elapsed: 2s | ETA: 12s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 550/2000 (27.5%)
#> Chain 2 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 301/2000 (15.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 851/4000 (21.3%)
#> Elapsed: 3s | ETA: 11s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 546/2000 (27.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1346/4000 (33.7%)
#> Elapsed: 3s | ETA: 6s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1000/2000 (50.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 800/2000 (40.0%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1800/4000 (45.0%)
#> Elapsed: 4s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1250/2000 (62.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1049/2000 (52.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2299/4000 (57.5%)
#> Elapsed: 4s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1450/2000 (72.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1267/2000 (63.3%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 2717/4000 (67.9%)
#> Elapsed: 5s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/2000 (85.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1538/2000 (76.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━⦘ 3238/4000 (81.0%)
#> Elapsed: 5s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2000 (95.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1758/2000 (87.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3658/4000 (91.5%)
#> Elapsed: 6s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Elapsed: 6s | ETA: 0s

# Posterior inclusion probabilities
summary(fit)$indicator
#>                     mean       mcse         sd n0->0 n0->1 n1->0 n1->1
#> intrusion-dreams  1.0000         NA 0.00000000     0     0     0  1999
#> intrusion-flash   1.0000         NA 0.00000000     0     0     0  1999
#> intrusion-upset   0.9785 0.01686819 0.14504396    40     3     3  1953
#> intrusion-physior 0.9310 0.03161494 0.25345414   130     8     8  1853
#> dreams-flash      1.0000         NA 0.00000000     0     0     0  1999
#> dreams-upset      0.9905 0.01312940 0.09700387    18     1     1  1979
#> dreams-physior    0.1070 0.01967130 0.30911325  1743    42    42   172
#> flash-upset       0.0705 0.01150602 0.25598779  1806    52    52    89
#> flash-physior     1.0000         NA 0.00000000     0     0     0  1999
#> upset-physior     1.0000         NA 0.00000000     0     0     0  1999
#>                   n_eff_mixt     Rhat
#> intrusion-dreams          NA       NA
#> intrusion-flash           NA       NA
#> intrusion-upset     73.93709 1.158511
#> intrusion-physior   64.27084 1.000519
#> dreams-flash              NA       NA
#> dreams-upset        54.58691 1.299312
#> dreams-physior     246.92731 1.023105
#> flash-upset        494.98105 1.000179
#> flash-physior             NA       NA
#> upset-physior             NA       NA

# Posterior pairwise effects
summary(fit)$pairwise
#>                          mean         mcse         sd     n_eff
#> intrusion-dreams  0.314503857 0.0007346667 0.03277968 1990.8048
#> intrusion-flash   0.168737949 0.0007595850 0.03218127 1794.9545
#> intrusion-upset   0.100814888 0.0020004680 0.03220720  512.8043
#> intrusion-physior 0.094454205 0.0032675270 0.03653265  156.3676
#> dreams-flash      0.249154680 0.0005951231 0.03007591 2554.0180
#> dreams-upset      0.111944591 0.0016570651 0.02891952  610.7461
#> dreams-physior    0.005161725 0.0009817649 0.01519639  208.5053
#> flash-upset       0.002813570 0.0004886639 0.01040716  475.4309
#> flash-physior     0.153835523 0.0005817084 0.02690679 2139.5038
#> upset-physior     0.354728946 0.0007122993 0.02942063 1706.0001
#>                   n_eff_mixt      Rhat
#> intrusion-dreams          NA 1.0016355
#> intrusion-flash           NA 0.9995401
#> intrusion-upset     259.2047 1.0090286
#> intrusion-physior   125.0040 1.0003540
#> dreams-flash              NA 1.0004376
#> dreams-upset        304.5812 1.0043932
#> dreams-physior      239.5884 1.0151082
#> flash-upset         453.5695 0.9997131
#> flash-physior             NA 0.9997872
#> upset-physior             NA 0.9996597
# }
```
