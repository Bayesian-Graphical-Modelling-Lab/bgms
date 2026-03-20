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
  pairwise_scale = 2.5,
  main_alpha = 0.5,
  main_beta = 0.5,
  edge_selection = TRUE,
  include_edge = NULL,
  edge_prior = c("Bernoulli", "Beta-Bernoulli", "Stochastic-Block"),
  inclusion_probability = 0.5,
  beta_bernoulli_alpha = 1,
  beta_bernoulli_beta = 1,
  beta_bernoulli_alpha_between = 1,
  beta_bernoulli_beta_between = 1,
  dirichlet_alpha = 1,
  lambda = 1,
  na_action = c("listwise", "impute"),
  update_method = c("nuts", "adaptive-metropolis", "hamiltonian-mc"),
  target_accept,
  hmc_num_leapfrogs = 100,
  nuts_max_depth = 10,
  learn_mass_matrix = TRUE,
  chains = 4,
  cores = parallel::detectCores(),
  display_progress = c("per-chain", "total", "none"),
  seed = NULL,
  standardize = FALSE,
  pseudolikelihood = c("conditional", "marginal"),
  verbose = getOption("bgms.verbose", TRUE),
  interaction_scale,
  burnin,
  save,
  threshold_alpha,
  threshold_beta
)
```

## Arguments

- x:

  A data frame or matrix with `n` rows and `p` columns containing binary
  and ordinal responses. Variables are automatically recoded to
  non-negative integers (`0, 1, ..., m`). For regular ordinal variables,
  unobserved categories are collapsed; for Blume–Capel variables, all
  categories are retained.

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

- pairwise_scale:

  Double. Scale of the Cauchy prior for pairwise interaction parameters.
  Default: `2.5`.

- main_alpha, main_beta:

  Double. Shape parameters of the beta-prime prior for threshold
  parameters. Must be positive. If equal, the prior is symmetric.
  Defaults: `main_alpha = 0.5` and `main_beta = 0.5`.

- edge_selection:

  Logical. Whether to perform Bayesian edge selection. If `FALSE`, the
  model estimates all edges (or only those specified by `include_edge`).
  Default: `TRUE`.

- include_edge:

  Optional binary matrix (`p x p`). Specifies which edges to include in
  the model. Must be symmetric with zero diagonal and contain only 0s
  and 1s. When `NULL` (the default), all edges are included (complete
  graph). When `edge_selection = FALSE`, only the edges marked 1 in this
  matrix are estimated; when `edge_selection = TRUE`, edge selection is
  restricted to the edges marked 1. Default: `NULL`.

- edge_prior:

  Character. Specifies the prior for edge inclusion. Options:
  `"Bernoulli"`, `"Beta-Bernoulli"`, or `"Stochastic-Block"`. Default:
  `"Bernoulli"`.

- inclusion_probability:

  Numeric scalar. Prior inclusion probability of each edge (used with
  the Bernoulli prior). Default: `0.5`.

- beta_bernoulli_alpha, beta_bernoulli_beta:

  Double. Shape parameters for the beta distribution in the
  Beta–Bernoulli and the Stochastic-Block priors. Must be positive. For
  the Stochastic-Block prior these are the shape parameters for the
  within-cluster edge inclusion probabilities. Defaults:
  `beta_bernoulli_alpha = 1` and `beta_bernoulli_beta = 1`.

- beta_bernoulli_alpha_between, beta_bernoulli_beta_between:

  Double. Shape parameters for the between-cluster edge inclusion
  probabilities in the Stochastic-Block prior. Must be positive.
  Default: `beta_bernoulli_alpha_between = 1` and
  `beta_bernoulli_beta_between = 1`

- dirichlet_alpha:

  Double. Concentration parameter of the Dirichlet prior on block
  assignments (used with the Stochastic Block model). Default: `1`.

- lambda:

  Double. Rate of the zero-truncated Poisson prior on the number of
  clusters in the Stochastic Block Model. Default: `1`.

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

  "hamiltonian-mc"

  :   Hamiltonian Monte Carlo with fixed path length (number of leapfrog
      steps set by `hmc_num_leapfrogs`).

  "nuts"

  :   The No-U-Turn Sampler, an adaptive form of HMC with dynamically
      chosen trajectory lengths.

  Default: `"nuts"`.

- target_accept:

  Numeric between 0 and 1. Target acceptance rate for the sampler.
  Defaults are set automatically if not supplied: `0.44` for adaptive
  Metropolis, `0.65` for HMC, and `0.80` for NUTS.

- hmc_num_leapfrogs:

  Integer. Number of leapfrog steps for Hamiltonian Monte Carlo. Must be
  positive. Default: `100`.

- nuts_max_depth:

  Integer. Maximum tree depth in NUTS. Must be positive. Default: `10`.

- learn_mass_matrix:

  Logical. If `TRUE`, adapt a diagonal mass matrix during warmup
  (HMC/NUTS only). If `FALSE`, use the identity matrix. Default: `TRUE`.

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

  Logical. If `TRUE`, the Cauchy prior scale for each pairwise
  interaction is adjusted based on the range of response scores.
  Variables with more response categories have larger score products
  \\x_i \cdot x_j\\, which typically correspond to smaller interaction
  effects \\\sigma\_{ij}\\. Without standardization, a fixed prior scale
  is relatively wide for these smaller effects, resulting in less
  shrinkage for high-category pairs and more shrinkage for low-category
  pairs. Standardization scales the prior proportionally to the maximum
  score product, ensuring equivalent relative shrinkage across all
  pairs. After internal recoding, regular ordinal variables have scores
  \\0, 1, \ldots, m\\. The adjusted scale for the interaction between
  variables \\i\\ and \\j\\ is `pairwise_scale * m_i * m_j`, so that
  `pairwise_scale` itself applies to the unit interval case (binary
  variables where \\m_i = m_j = 1\\). For Blume-Capel variables with
  reference category \\b\\, scores are centered as \\-b, \ldots, m-b\\,
  and the adjustment uses the maximum absolute product of the score
  endpoints. For mixed pairs, ordinal variables use raw score endpoints
  \\(0, m)\\ and Blume-Capel variables use centered score endpoints
  \\(-b, m-b)\\. Default: `FALSE`.

- pseudolikelihood:

  Character. Specifies the pseudo-likelihood approximation used for
  mixed MRF models (ignored for pure ordinal or pure continuous data).
  Options:

  `"conditional"`

  :   Conditions on the observed continuous variables when computing the
      discrete full conditionals. Faster because the discrete
      pseudo-likelihood does not depend on the continuous precision
      matrix.

  `"marginal"`

  :   Integrates out the continuous variables, giving discrete full
      conditionals that account for induced interactions through the
      continuous block. More expensive per iteration.

  Default: `"conditional"`.

- verbose:

  Logical. If `TRUE`, prints informational messages during data
  processing (e.g., missing data handling, variable recoding). Defaults
  to `getOption("bgms.verbose", TRUE)`. Set
  `options(bgms.verbose = FALSE)` to suppress messages globally.

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

- `posterior_mean_associations`: Symmetric matrix of posterior mean
  partial associations (zero diagonal). For continuous variables these
  are unstandardized partial correlations; for discrete variables these
  are half the log adjacent-category odds ratio. Use
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

This function models the joint distribution of binary, ordinal,
continuous, or mixed variables using a Markov Random Field, with support
for edge selection through Bayesian variable selection. The statistical
foundation of the model is described in Marsman et al. (2025) , where
the ordinal MRF model and its Bayesian estimation procedure were first
introduced. While the implementation in bgms has since been extended and
updated (e.g., alternative priors, parallel chains, HMC/NUTS warmup), it
builds on that original framework.

Key components of the model are described in the sections below.

## Ordinal Variables

The function supports two types of ordinal variables:

**Regular ordinal variables**: Assigns a category threshold parameter to
each response category except the lowest. The model imposes no
additional constraints on the distribution of category responses.

**Blume-Capel ordinal variables**: Assume a baseline category (e.g., a
“neutral” response) and score responses by distance from this baseline.
Category thresholds are modeled as:

\$\$\mu\_{c} = \alpha \cdot (c-b) + \beta \cdot (c - b)^2\$\$

where:

- \\\mu\_{c}\\: category threshold for category \\c\\

- \\\alpha\\: linear trend across categories

- \\\beta\\: preference toward or away from the baseline

  - If \\\beta \< 0\\, the model favors responses near the baseline
    category;

  - if \\\beta \> 0\\, it favors responses farther away (i.e.,
    extremes).

- \\b\\: baseline category

Accordingly, pairwise interactions between Blume-Capel variables are
modeled in terms of \\c-b\\ scores.

## Edge Selection

When `edge_selection = TRUE`, the function performs Bayesian variable
selection on the pairwise interactions (edges) in the MRF using
spike-and-slab priors.

Supported priors for edge inclusion:

- **Bernoulli**: Fixed inclusion probability across edges.

- **Beta-Bernoulli**: Inclusion probability is assigned a Beta prior
  distribution.

- **Stochastic-Block**: Cluster-based edge priors with Beta, Dirichlet,
  and Poisson hyperpriors.

All priors operate via binary indicator variables controlling the
inclusion or exclusion of each edge in the MRF.

## Prior Distributions

- **Pairwise effects**: Modeled with a Cauchy (slab) prior.

- **Main effects**: Modeled using a beta-prime distribution.

- **Edge indicators**: Use either a Bernoulli, Beta-Bernoulli, or
  Stochastic-Block prior (as above).

## Sampling Algorithms and Warmup

Parameters are updated within a Gibbs framework, but the conditional
updates can be carried out using different algorithms:

- **Adaptive Metropolis–Hastings**: Componentwise random–walk updates
  for main effects and pairwise effects. Proposal standard deviations
  are adapted during burn–in via Robbins–Monro updates toward a target
  acceptance rate.

- **Hamiltonian Monte Carlo (HMC)**: Joint updates of all parameters
  using fixed–length leapfrog trajectories. Step size is tuned during
  warmup via dual–averaging; the diagonal mass matrix can also be
  adapted if `learn_mass_matrix = TRUE`.

- **No–U–Turn Sampler (NUTS)**: An adaptive extension of HMC that
  dynamically chooses trajectory lengths. Warmup uses a staged
  adaptation schedule (fast–slow–fast) to stabilize step size and, if
  enabled, the mass matrix.

When `edge_selection = TRUE`, updates of edge–inclusion indicators are
carried out with Metropolis–Hastings steps. These are switched on after
the core warmup phase, ensuring that graph updates occur only once the
samplers’ tuning parameters (step size, mass matrix, proposal SDs) have
stabilized.

After warmup, adaptation is disabled. Step size and mass matrix are
fixed at their learned values, and proposal SDs remain constant.

## Warmup and Adaptation

The warmup procedure in `bgm` uses a multi-stage adaptation schedule
(Stan Development Team 2023) . Warmup iterations are split into several
phases:

- **Stage 1 (fast adaptation)**: A short initial interval where only
  step size (for HMC/NUTS) is adapted, allowing the chain to move
  quickly toward the typical set.

- **Stage 2 (slow windows)**: A sequence of expanding, memoryless
  windows where both step size and, if `learn_mass_matrix = TRUE`, the
  diagonal mass matrix are adapted. Each window ends with a reset of the
  dual–averaging scheme for improved stability.

- **Stage 3a (final fast interval)**: A short interval at the end of the
  core warmup where the step size is adapted one final time.

- **Stage 3b (proposal–SD tuning)**: Only active when
  `edge_selection = TRUE` under HMC/NUTS. In this phase, Robbins–Monro
  adaptation of proposal standard deviations is performed for the
  Metropolis steps used in edge–selection moves.

- **Stage 3c (graph selection warmup)**: Also only relevant when
  `edge_selection = TRUE`. At the start of this phase, a random graph
  structure is initialized, and Metropolis–Hastings updates for edge
  inclusion indicators are switched on.

When `edge_selection = FALSE`, the total number of warmup iterations
equals the user–specified `burnin`. When `edge_selection = TRUE` and
`update_method` is `"nuts"` or `"hamiltonian-mc"`, the schedule
automatically appends additional Stage-3b and Stage-3c intervals, so the
total warmup is strictly greater than the requested `burnin`.

After all warmup phases, the sampler transitions to the sampling phase
with adaptation disabled. Step size and mass matrix (for HMC/NUTS) are
fixed at their learned values, and proposal SDs remain constant.

This staged design improves stability of proposals and ensures that both
local parameters (step size) and global parameters (mass matrix,
proposal SDs) are tuned before collecting posterior samples.

For adaptive Metropolis–Hastings runs, step size and mass matrix
adaptation are not relevant. Proposal SDs are tuned continuously during
burn–in using Robbins–Monro updates, without staged fast/slow intervals.

## Missing Data

If `na_action = "listwise"`, rows with missing values are removed. If
`na_action = "impute"`, missing values are imputed within the MCMC loop
via Gibbs sampling.

## References

Dahl DB (2009). “Modal clustering in a class of product partition
models.” *Bayesian Analysis*, **4**(2), 243–264.
[doi:10.1214/09-BA409](https://doi.org/10.1214/09-BA409) .  
  
Marsman M, van den Bergh D, Haslbeck JMB (2025). “Bayesian analysis of
the ordinal Markov random field.” *Psychometrika*, **90**(1), 146–182.
[doi:10.1017/psy.2024.4](https://doi.org/10.1017/psy.2024.4) .  
  
Sekulovski N, Arena G, Haslbeck JMB, Huth KBS, Friel N, Marsman M
(2025). “A Stochastic Block Prior for Clustering in Graphical Models.”
*Retrieved from <https://osf.io/preprints/psyarxiv/29p3m_v1>*. OSF
preprint.  
  
Stan Development Team (2023). *Stan Modeling Language Users Guide and
Reference Manual*. Version 2.33, <https://mc-stan.org/docs/>.

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
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 157/2000 (7.8%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 257/4000 (6.4%)
#> Elapsed: 1s | ETA: 15s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 350/2000 (17.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 411/2000 (20.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 761/4000 (19.0%)
#> Elapsed: 1s | ETA: 4s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 650/2000 (32.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 725/2000 (36.2%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1375/4000 (34.4%)
#> Elapsed: 2s | ETA: 4s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2000 (45.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 968/2000 (48.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1868/4000 (46.7%)
#> Elapsed: 2s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2000 (57.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1238/2000 (61.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2388/4000 (59.7%)
#> Elapsed: 3s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1400/2000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1508/2000 (75.4%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 2908/4000 (72.7%)
#> Elapsed: 3s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1650/2000 (82.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1773/2000 (88.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 3423/4000 (85.6%)
#> Elapsed: 4s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2000 (95.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3900/4000 (97.5%)
#> Elapsed: 4s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4000/4000 (100.0%)
#> Elapsed: 5s | ETA: 0s

# Posterior inclusion probabilities
summary(fit)$indicator
#>                     mean        sd       mcse n0->0 n0->1 n1->0 n1->1
#> intrusion-dreams  1.0000 0.0000000         NA     0     0     0  1999
#> intrusion-flash   1.0000 0.0000000         NA     0     0     0  1999
#> intrusion-upset   0.9050 0.2932149 0.04548424   183     7     7  1802
#> intrusion-physior 0.9765 0.1514851 0.03227556    46     1     1  1951
#> dreams-flash      1.0000 0.0000000         NA     0     0     0  1999
#> dreams-upset      1.0000 0.0000000         NA     0     0     0  1999
#> dreams-physior    0.3295 0.4700316 0.06964994  1321    19    20   639
#> flash-upset       0.4700 0.4990992 0.09887164  1047    12    13   927
#> flash-physior     1.0000 0.0000000         NA     0     0     0  1999
#> upset-physior     0.5000 0.5000000 0.49987498   999     1     0   999
#>                      n_eff     Rhat
#> intrusion-dreams        NA       NA
#> intrusion-flash         NA       NA
#> intrusion-upset   41.55759 1.427530
#> intrusion-physior 22.02888 1.314054
#> dreams-flash            NA       NA
#> dreams-upset            NA       NA
#> dreams-physior    45.54207 2.233717
#> flash-upset       25.48181 5.139062
#> flash-physior           NA       NA
#> upset-physior      1.00050      Inf

# Posterior pairwise effects
summary(fit)$pairwise
#>                         mean           sd        mcse       n_eff
#> intrusion-dreams  0.31251788 0.0011086516 0.034593533  973.642758
#> intrusion-flash   0.16045461 0.0009040525 0.032658399 1304.975966
#> intrusion-upset   0.10480672 0.0447146517 0.006197096   52.062294
#> intrusion-physior 0.12041328 0.0386713620 0.006253218   38.244741
#> dreams-flash      0.24233305 0.0008939708 0.030951893 1198.747652
#> dreams-upset      0.12390582 0.0017188152 0.029821963  301.032974
#> dreams-physior    0.02583121 0.0379260349 0.005514738   47.296084
#> flash-upset       0.04395110 0.0483165055 0.009265344   27.193674
#> flash-physior     0.15623361 0.0011448480 0.028153156  604.725848
#> upset-physior     0.17887233 0.1795900823 0.178829050    1.008529
#>                       Rhat
#> intrusion-dreams  1.015295
#> intrusion-flash   1.208862
#> intrusion-upset   1.298104
#> intrusion-physior 1.266191
#> dreams-flash      1.066401
#> dreams-upset      1.008753
#> dreams-physior    1.499734
#> flash-upset       2.542126
#> flash-physior     1.010603
#> upset-physior           NA
# }
```
