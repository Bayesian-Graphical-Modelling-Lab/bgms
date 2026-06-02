# Bayesian Estimation or Edge Selection for Markov Random Fields

The `bgm` function estimates the pseudoposterior distribution of the
parameters of a Markov Random Field (MRF) for binary, ordinal,
continuous, or mixed (discrete and continuous) variables. Optionally, it
performs Bayesian edge selection using discrete spike-and-slab priors to
infer the network structure.

## Usage

``` r
bgm(
  x,
  variable_type = "ordinal",
  baseline_category,
  iter = 2000,
  warmup = 2000,
  interaction_prior = cauchy_prior(scale = 1),
  threshold_prior = beta_prime_prior(alpha = 0.5, beta = 0.5),
  means_prior = normal_prior(scale = 1),
  precision_scale_prior = gamma_prior(shape = 1, rate = 1),
  delta = NULL,
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
  `2e3`.

- warmup:

  Integer. Number of warmup iterations before collecting samples. Short
  warmups trigger progressive warnings (NUTS only); see
  `validate_sampler()` for the thresholds. Default: `2e3`.

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

- delta:

  Non-negative numeric, or `NULL` for the dimension- adaptive default.
  Determinant-tilt exponent on the continuous-block precision matrix
  \\K\\ (GGM) or \\K\_{yy}\\ (mixed MRF): multiplies the prior by
  \\\|K\|^{\delta}\\, softly repelling the chain from the
  positive-definite cone boundary. `delta = NULL` (default)
  auto-resolves to \\0.5 \log(p)\\ where \\p\\ is the dimension of the
  continuous precision matrix (the number of variables for GGM, the
  number of continuous variables for mixed MRF). The rule is the simple
  form of the dimension-adaptive scaling \\\delta(p) = c \log p\\ with
  \\c \in (0.3, 0.6)\\ discussed in the companion paper on
  determinant-tilted spike-and-slab priors (Marsman et al., in
  preparation). Pass `delta = 0` for the untilted prior (the
  companion-paper baseline) or a non-negative numeric to override. Both
  NUTS and adaptive-Metropolis update paths apply the tilt. Not allowed
  for pure ordinal models (no precision matrix to tilt).

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
Beta-Bernoulli, or Stochastic-Block priors on the edge inclusion
indicators. Parameters are sampled with NUTS (default) or adaptive
Metropolis–Hastings, with a multi-stage warmup schedule. Missing data
can be handled via listwise deletion or Gibbs imputation.

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
#> Chain 1 (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/4000 (1.2%)
#> Chain 2 (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 49/4000 (1.2%)
#> Total   (Warmup): ⦗╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 99/8000 (1.2%)
#> Elapsed: 1s | ETA: 1m 19s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/4000 (2.5%)
#> Chain 2 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 131/4000 (3.3%)
#> Total   (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 231/8000 (2.9%)
#> Elapsed: 2s | ETA: 1m 7s
#> Chain 1 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 250/4000 (6.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 251/4000 (6.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 501/8000 (6.3%)
#> Elapsed: 3s | ETA: 45s
#> Chain 1 (Warmup): ⦗━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 450/4000 (11.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 260/4000 (6.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 710/8000 (8.9%)
#> Elapsed: 3s | ETA: 31s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/4000 (12.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 579/4000 (14.5%)
#> Total   (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1079/8000 (13.5%)
#> Elapsed: 5s | ETA: 32s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/4000 (17.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 797/4000 (19.9%)
#> Total   (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1497/8000 (18.7%)
#> Elapsed: 5s | ETA: 22s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/4000 (22.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 972/4000 (24.3%)
#> Total   (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1872/8000 (23.4%)
#> Elapsed: 6s | ETA: 20s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/4000 (27.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1194/4000 (29.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2294/8000 (28.7%)
#> Elapsed: 6s | ETA: 15s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1300/4000 (32.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1386/4000 (34.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2685/8000 (33.6%)
#> Elapsed: 7s | ETA: 14s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1500/4000 (37.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1590/4000 (39.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3090/8000 (38.6%)
#> Elapsed: 7s | ETA: 11s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1700/4000 (42.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 1735/4000 (43.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3435/8000 (42.9%)
#> Elapsed: 9s | ETA: 12s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1850/4000 (46.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1871/4000 (46.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3721/8000 (46.5%)
#> Elapsed: 10s | ETA: 11s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/4000 (50.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 2016/4000 (50.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 4016/8000 (50.2%)
#> Elapsed: 10s | ETA: 10s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/4000 (55.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 2201/4000 (55.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 4401/8000 (55.0%)
#> Elapsed: 11s | ETA: 9s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2400/4000 (60.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2396/4000 (59.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4796/8000 (60.0%)
#> Elapsed: 12s | ETA: 8s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2600/4000 (65.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2589/4000 (64.7%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5189/8000 (64.9%)
#> Elapsed: 12s | ETA: 7s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2800/4000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2783/4000 (69.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5583/8000 (69.8%)
#> Elapsed: 13s | ETA: 6s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3000/4000 (75.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2984/4000 (74.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5984/8000 (74.8%)
#> Elapsed: 14s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3200/4000 (80.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3183/4000 (79.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6383/8000 (79.8%)
#> Elapsed: 14s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 3350/4000 (83.8%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 3345/4000 (83.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 6695/8000 (83.7%)
#> Elapsed: 15s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3500/4000 (87.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3488/4000 (87.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 6988/8000 (87.4%)
#> Elapsed: 15s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3700/4000 (92.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3684/4000 (92.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7384/8000 (92.3%)
#> Elapsed: 16s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3900/4000 (97.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3870/4000 (96.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7770/8000 (97.1%)
#> Elapsed: 17s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 17s | ETA: 0s

# Posterior inclusion probabilities
summary(fit)$indicator
#>                      mean        mcse        sd n0->0 n0->1 n1->0
#> intrusion-dreams  1.00000          NA 0.0000000     0     0     0
#> intrusion-flash   1.00000          NA 0.0000000     0     0     0
#> intrusion-upset   0.95225 0.018254152 0.2132368   179    12    12
#> intrusion-physior 0.93500 0.020477792 0.2465259   243    17    17
#> dreams-flash      1.00000          NA 0.0000000     0     0     0
#> dreams-upset      0.98825 0.008032146 0.1077587    43     4     4
#> dreams-physior    0.10075 0.011693132 0.3009974  3493   103   103
#> flash-upset       0.10175 0.013458813 0.3023193  3511    81    82
#> flash-physior     1.00000          NA 0.0000000     0     0     0
#> upset-physior     1.00000          NA 0.0000000     0     0     0
#>                   n1->1 n_eff_mixt     Rhat
#> intrusion-dreams   3999         NA       NA
#> intrusion-flash    3999         NA       NA
#> intrusion-upset    3796   136.4586 1.011704
#> intrusion-physior  3722   144.9301 0.999751
#> dreams-flash       3999         NA       NA
#> dreams-upset       3948   179.9872 1.103209
#> dreams-physior      300   662.6189 1.015119
#> flash-upset         325   504.5662 1.000324
#> flash-physior      3999         NA       NA
#> upset-physior      3999         NA       NA

# Posterior pairwise effects
summary(fit)$pairwise
#>                          mean         mcse         sd     n_eff
#> intrusion-dreams  0.315060550 0.0005236290 0.03277394 3917.5080
#> intrusion-flash   0.169199168 0.0004891054 0.03108406 4038.9684
#> intrusion-upset   0.099373511 0.0021064173 0.03640799  254.2589
#> intrusion-physior 0.095257062 0.0022101178 0.03644230  203.2061
#> dreams-flash      0.249247938 0.0004284155 0.02961637 4778.9573
#> dreams-upset      0.112501497 0.0010925338 0.02969141  944.9056
#> dreams-physior    0.005097723 0.0007212473 0.01554461  464.7216
#> flash-upset       0.005328535 0.0007188026 0.01608590  446.0611
#> flash-physior     0.153312042 0.0005254989 0.02665881 2573.5825
#> upset-physior     0.355441688 0.0005332073 0.02968287 3098.9853
#>                   n_eff_mixt      Rhat
#> intrusion-dreams          NA 0.9999356
#> intrusion-flash           NA 0.9997602
#> intrusion-upset     298.7477 1.0036149
#> intrusion-physior   271.8822 1.0003558
#> dreams-flash              NA 1.0007538
#> dreams-upset        738.5701 1.0061278
#> dreams-physior      464.5057 1.0564844
#> flash-upset         500.8081 0.9999379
#> flash-physior             NA 0.9998513
#> upset-physior             NA 1.0000148
# }
```
