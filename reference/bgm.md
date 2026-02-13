# Bayesian Estimation or Edge Selection for Markov Random Fields

The `bgm` function estimates the pseudoposterior distribution of
category thresholds (main effects) and pairwise interaction parameters
of a Markov Random Field (MRF) model for binary and/or ordinal
variables. Optionally, it performs Bayesian edge selection using
spike-and-slab priors to infer the network structure.

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
  `x`. Allowed values: `"ordinal"` or `"blume-capel"`. Binary variables
  are automatically treated as `"ordinal"`. Default: `"ordinal"`.

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
  model estimates all edges. Default: `TRUE`.

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

- verbose:

  Logical. If `TRUE`, prints informational messages during data
  processing (e.g., missing data handling, variable recoding). Defaults
  to `getOption("bgms.verbose", TRUE)`. Set
  `options(bgms.verbose = FALSE)` to suppress messages globally.

- interaction_scale, burnin, save, threshold_alpha, threshold_beta:

  \`r lifecycle::badge("deprecated")\` Deprecated arguments as of
  \*\*bgms 0.1.6.0\*\*. Use \`pairwise_scale\`, \`warmup\`,
  \`main_alpha\`, and \`main_beta\` instead.

## Value

A list of class `"bgms"` with posterior summaries, posterior mean
matrices, and access to raw MCMC draws. The object can be passed to
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`coef()`](https://rdrr.io/r/stats/coef.html).

Main components include:

- `posterior_summary_main`: Data frame with posterior summaries (mean,
  sd, MCSE, ESS, Rhat) for category threshold parameters.

- `posterior_summary_pairwise`: Data frame with posterior summaries for
  pairwise interaction parameters.

- `posterior_summary_indicator`: Data frame with posterior summaries for
  edge inclusion indicators (if `edge_selection = TRUE`).

- `posterior_mean_main`: Matrix of posterior mean thresholds (rows =
  variables, cols = categories or parameters).

- `posterior_mean_pairwise`: Symmetric matrix of posterior mean pairwise
  interaction strengths.

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

This function models the joint distribution of binary and ordinal
variables using a Markov Random Field, with support for edge selection
through Bayesian variable selection. The statistical foundation of the
model is described in Marsman et al. (2025) , where the ordinal MRF
model and its Bayesian estimation procedure were first introduced. While
the implementation in bgms has since been extended and updated (e.g.,
alternative priors, parallel chains, HMC/NUTS warmup), it builds on that
original framework.

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

If `na_action = "listwise"`, observations with missing values are
removed. If `na_action = "impute"`, missing values are imputed during
Gibbs sampling.

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

## Examples

``` r
# \donttest{
# Run bgm on subset of the Wenchuan dataset
fit = bgm(x = Wenchuan[, 1:5])
#> 7 rows with missing values excluded (n = 355 remaining).
#> To impute missing values instead, use na_action = "impute".
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2000 (5.0%)
#> Chain 2 (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 113/2000 (5.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 99/2000 (5.0%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 93/2000 (4.7%)
#> Total   (Warmup): ⦗━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 405/8000 (5.1%)
#> Elapsed: 0s | ETA: 0s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2000 (15.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 296/2000 (14.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 288/2000 (14.4%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 288/2000 (14.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1172/8000 (14.6%)
#> Elapsed: 1s | ETA: 6s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 500/2000 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 493/2000 (24.6%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 491/2000 (24.6%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 481/2000 (24.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1965/8000 (24.6%)
#> Elapsed: 1s | ETA: 3s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 750/2000 (37.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 732/2000 (36.6%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 742/2000 (37.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 723/2000 (36.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2947/8000 (36.8%)
#> Elapsed: 2s | ETA: 3s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 950/2000 (47.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 927/2000 (46.4%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 943/2000 (47.1%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 922/2000 (46.1%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3742/8000 (46.8%)
#> Elapsed: 3s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1100/2000 (55.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━⦘ 1075/2000 (53.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1103/2000 (55.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1076/2000 (53.8%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4354/8000 (54.4%)
#> Elapsed: 3s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1250/2000 (62.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1236/2000 (61.8%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━⦘ 1264/2000 (63.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1244/2000 (62.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4994/8000 (62.4%)
#> Elapsed: 4s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1400/2000 (70.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1392/2000 (69.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1427/2000 (71.4%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1399/2000 (70.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━⦘ 5618/8000 (70.2%)
#> Elapsed: 4s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1550/2000 (77.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1553/2000 (77.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1589/2000 (79.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 1562/2000 (78.1%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━⦘ 6254/8000 (78.2%)
#> Elapsed: 5s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/2000 (87.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1760/2000 (88.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━⦘ 1809/2000 (90.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1765/2000 (88.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 7084/8000 (88.5%)
#> Elapsed: 5s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2000 (95.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 1913/2000 (95.7%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 1962/2000 (98.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 1904/2000 (95.2%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 7679/8000 (96.0%)
#> Elapsed: 6s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2000/2000 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8000/8000 (100.0%)
#> Elapsed: 6s | ETA: 0s
#> NUTS issues:
#>   - Divergences: 2 (0.050%) - check R-hat and ESS 

# Posterior inclusion probabilities
summary(fit)$indicator
#>                      mean        sd       mcse n0->0 n0->1 n1->0 n1->1
#> intrusion-dreams  1.00000 0.0000000         NA     0     0     0  3999
#> intrusion-flash   1.00000 0.0000000         NA     0     0     0  3999
#> intrusion-upset   0.91400 0.2803640 0.02743455   328    16    16  3639
#> intrusion-physior 0.91775 0.2747452 0.02717975   313    15    15  3656
#> dreams-flash      1.00000 0.0000000         NA     0     0     0  3999
#> dreams-upset      0.98225 0.1320414 0.01412438    68     3     2  3926
#> dreams-physior    0.09100 0.2876091 0.01190917  3552    84    84   279
#> flash-upset       0.08400 0.2773878 0.01088939  3578    85    86   250
#> flash-physior     1.00000 0.0000000         NA     0     0     0  3999
#> upset-physior     1.00000 0.0000000         NA     0     0     0  3999
#>                       n_eff     Rhat
#> intrusion-dreams         NA       NA
#> intrusion-flash          NA       NA
#> intrusion-upset   104.43571 1.096348
#> intrusion-physior 102.18079 1.079932
#> dreams-flash             NA       NA
#> dreams-upset       87.39397 1.277365
#> dreams-physior    583.23307 1.031354
#> flash-upset       648.88539 1.014082
#> flash-physior            NA       NA
#> upset-physior            NA       NA

# Posterior pairwise effects
summary(fit)$pairwise
#>                          mean          sd        mcse     n_eff
#> intrusion-dreams  0.631503597 0.001455279 0.067232862 2134.3741
#> intrusion-flash   0.339602596 0.001257684 0.063250258 2529.1904
#> intrusion-upset   0.190240027 0.081675158 0.006114818  178.4074
#> intrusion-physior 0.190193773 0.078164164 0.005868097  177.4275
#> dreams-flash      0.499328447 0.001227521 0.061412693 2502.9877
#> dreams-upset      0.225564812 0.064091870 0.003796965  284.9264
#> dreams-physior    0.010428208 0.033461570 0.001560120  460.0197
#> flash-upset       0.008412741 0.028248758 0.001160336  592.6949
#> flash-physior     0.305800951 0.001143566 0.053622000 2198.6893
#> upset-physior     0.712651045 0.001360732 0.062581321 2115.1665
#>                       Rhat
#> intrusion-dreams  1.006107
#> intrusion-flash   1.000530
#> intrusion-upset   1.008680
#> intrusion-physior 1.010389
#> dreams-flash      1.001722
#> dreams-upset      1.004833
#> dreams-physior    1.012984
#> flash-upset       1.004246
#> flash-physior     1.004997
#> upset-physior     1.002660
# }
```
