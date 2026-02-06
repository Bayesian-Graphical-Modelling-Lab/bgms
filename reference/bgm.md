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

  Logical. Whether to show a progress bar during sampling. Default:
  `TRUE`.

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

The warmup procedure in `bgm` is based on the multi–stage adaptation
schedule used in Stan (Stan Development Team 2023) . Warmup iterations
are split into several phases:

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
#> Warning: There were 7 rows with missing observations in the input matrix x.
#> Since na_action = listwise these rows were excluded from the analysis.
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2200 (2.3%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 50/2200 (2.3%)
#> Chain 3 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 60/2200 (2.7%)
#> Chain 4 (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 61/2200 (2.8%)
#> Total   (Warmup): ⦗━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 221/8800 (2.5%)
#> Elapsed: 2s | ETA: 1m 17s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 100/2200 (4.5%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 154/2200 (7.0%)
#> Chain 3 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 168/2200 (7.6%)
#> Chain 4 (Warmup): ⦗━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 189/2200 (8.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 611/8800 (6.9%)
#> Elapsed: 3s | ETA: 40s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 200/2200 (9.1%)
#> Chain 2 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 276/2200 (12.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 305/2200 (13.9%)
#> Chain 4 (Warmup): ⦗━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 332/2200 (15.1%)
#> Total   (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1113/8800 (12.6%)
#> Elapsed: 4s | ETA: 28s
#> Chain 1 (Warmup): ⦗━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 300/2200 (13.6%)
#> Chain 2 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 386/2200 (17.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 445/2200 (20.2%)
#> Chain 4 (Warmup): ⦗━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 449/2200 (20.4%)
#> Total   (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1580/8800 (18.0%)
#> Elapsed: 5s | ETA: 23s
#> Chain 1 (Warmup): ⦗━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 400/2200 (18.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 494/2200 (22.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 543/2200 (24.7%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 570/2200 (25.9%)
#> Total   (Warmup): ⦗━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2007/8800 (22.8%)
#> Elapsed: 5s | ETA: 17s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 550/2200 (25.0%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 671/2200 (30.5%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 738/2200 (33.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 745/2200 (33.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2704/8800 (30.7%)
#> Elapsed: 6s | ETA: 14s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 700/2200 (31.8%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 853/2200 (38.8%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 900/2200 (40.9%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━⦘ 891/2200 (40.5%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 3344/8800 (38.0%)
#> Elapsed: 7s | ETA: 11s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━━━⦘ 850/2200 (38.6%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 965/2200 (43.9%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1017/2200 (46.2%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━⦘ 1004/2200 (45.6%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 3836/8800 (43.6%)
#> Elapsed: 7s | ETA: 9s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━━━⦘ 950/2200 (43.2%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1072/2200 (48.7%)
#> Chain 3 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1133/2200 (51.5%)
#> Chain 4 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━⦘ 1107/2200 (50.3%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 4262/8800 (48.4%)
#> Elapsed: 8s | ETA: 9s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━━━━⦘ 1050/2200 (47.7%)
#> Chain 2 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1194/2200 (54.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1259/2200 (57.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━━━━━⦘ 1229/2200 (55.9%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 4732/8800 (53.8%)
#> Elapsed: 8s | ETA: 7s
#> Chain 1 (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1150/2200 (52.3%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1319/2200 (60.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1368/2200 (62.2%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1350/2200 (61.4%)
#> Total   (Warmup): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 5187/8800 (58.9%)
#> Elapsed: 9s | ETA: 6s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1300/2200 (59.1%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1472/2200 (66.9%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━⦘ 1497/2200 (68.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1515/2200 (68.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 5784/8800 (65.7%)
#> Elapsed: 10s | ETA: 5s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━━━━⦘ 1450/2200 (65.9%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1634/2200 (74.3%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1638/2200 (74.5%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━⦘ 1660/2200 (75.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 6382/8800 (72.5%)
#> Elapsed: 10s | ETA: 4s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━━━━━⦘ 1600/2200 (72.7%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1796/2200 (81.6%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1798/2200 (81.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━━⦘ 1818/2200 (82.6%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 7012/8800 (79.7%)
#> Elapsed: 11s | ETA: 3s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1750/2200 (79.5%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1940/2200 (88.2%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━⦘ 1929/2200 (87.7%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1956/2200 (88.9%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━━━━⦘ 7575/8800 (86.1%)
#> Elapsed: 11s | ETA: 2s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 1900/2200 (86.4%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2100/2200 (95.5%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2071/2200 (94.1%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━⦘ 2100/2200 (95.5%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 8171/8800 (92.9%)
#> Elapsed: 12s | ETA: 1s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺━━⦘ 2050/2200 (93.2%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╺⦘ 8650/8800 (98.3%)
#> Elapsed: 13s | ETA: 0s
#> Chain 1 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 2 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 3 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Chain 4 (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 2200/2200 (100.0%)
#> Total   (Sampling): ⦗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━⦘ 8800/8800 (100.0%)
#> Elapsed: 13s | ETA: 0s
#> NUTS Diagnostics Summary:
#>   Total divergences:         0 
#>   Max tree depth hits:       0 
#>   Min E-BFMI across chains:  0.972 

# Posterior inclusion probabilities
summary(fit)$indicator
#>                      mean         sd        mcse n0->0 n0->1 n1->0
#> intrusion-dreams  1.00000 0.00000000          NA     0     0     0
#> intrusion-flash   1.00000 0.00000000          NA     0     0     0
#> intrusion-upset   0.92375 0.26539770 0.024546504   289    16    16
#> intrusion-physior 0.96775 0.17666335 0.015352867   121     8     8
#> dreams-flash      1.00000 0.00000000          NA     0     0     0
#> dreams-upset      0.99025 0.09825954 0.007728604    36     3     3
#> dreams-physior    0.06075 0.23887117 0.008656393  3683    73    73
#> flash-upset       0.09800 0.29731465 0.015164232  3545    62    62
#> flash-physior     1.00000 0.00000000          NA     0     0     0
#> upset-physior     1.00000 0.00000000          NA     0     0     0
#>                   n1->1    n_eff     Rhat
#> intrusion-dreams   3999       NA       NA
#> intrusion-flash    3999       NA       NA
#> intrusion-upset    3678 116.9001 1.035502
#> intrusion-physior  3862 132.4079 1.203798
#> dreams-flash       3999       NA       NA
#> dreams-upset       3957 161.6394 1.314275
#> dreams-physior      170 761.4714 1.038887
#> flash-upset         330 384.4074 1.048852
#> flash-physior      3999       NA       NA
#> upset-physior      3999       NA       NA

# Posterior pairwise effects
summary(fit)$pairwise
#>                         mean          sd         mcse     n_eff
#> intrusion-dreams  0.63322577 0.001769036 0.0657926511 1383.1892
#> intrusion-flash   0.33915123 0.001445381 0.0602123871 1735.4303
#> intrusion-upset   0.19045836 0.075104543 0.0053418243  197.6757
#> intrusion-physior 0.19920445 0.065900971 0.0037517518  308.5428
#> dreams-flash      0.50069916 0.001472836 0.0596007813 1637.5513
#> dreams-upset      0.22802296 0.058240240 0.0025973772  502.7780
#> dreams-physior    0.00561329 0.022329097 0.0008753891  650.6390
#> flash-upset       0.01060057 0.032587956 0.0017027343  366.2861
#> flash-physior     0.30819861 0.001425912 0.0523173871 1346.1906
#> upset-physior     0.71093147 0.001639468 0.0590201885 1295.9711
#>                       Rhat
#> intrusion-dreams  1.001787
#> intrusion-flash   1.002640
#> intrusion-upset   1.007915
#> intrusion-physior 1.010656
#> dreams-flash      1.004181
#> dreams-upset      1.005801
#> dreams-physior    1.009234
#> flash-upset       1.015692
#> flash-physior     1.003802
#> upset-physior     1.002067
# }
```
