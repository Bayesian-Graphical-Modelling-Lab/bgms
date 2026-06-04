# Diagnostics and Spike-and-Slab Summaries

## Introduction

This vignette illustrates how to inspect convergence diagnostics and how
to interpret spike-and-slab summaries in **bgms** models. For some of
the model variables spike-and-slab priors introduce binary indicator
variables that govern whether the effect is included or not. Their
posterior distributions can be summarized with inclusion probabilities
and Bayes factors.

## Example fit

We use a subset of the Wenchuan dataset:

``` r

library(bgms)
data = Wenchuan[, 1:5]
fit = bgm(data, seed = 1234)
```

## Convergence diagnostics

The quality of the Markov chain can be assessed with common MCMC
diagnostics:

``` r

summary(fit)$pairwise
#>                          mean         mcse          sd     n_eff
#> intrusion-dreams  0.314801118 0.0004860977 0.032025542 4340.5611
#> intrusion-flash   0.169977835 0.0004759129 0.030619258 4139.3707
#> intrusion-upset   0.095469103 0.0023457146 0.036762631  276.6104
#> intrusion-physior 0.100587283 0.0012805995 0.031273334  619.1365
#> dreams-flash      0.249816272 0.0004452531 0.030527181 4700.6665
#> dreams-upset      0.114758607 0.0006497733 0.027979791 2135.1585
#> dreams-physior    0.002218409 0.0002963578 0.009321959 1387.4272
#> flash-upset       0.004178344 0.0005480479 0.013722190  659.9317
#> flash-physior     0.152727316 0.0004325839 0.026775582 3831.2198
#> upset-physior     0.355014255 0.0004894762 0.029733884 3690.1176
#>                   n_eff_mixt      Rhat
#> intrusion-dreams          NA 1.0020851
#> intrusion-flash           NA 1.0013489
#> intrusion-upset     245.6195 1.0012452
#> intrusion-physior   596.3781 0.9998937
#> dreams-flash              NA 0.9997563
#> dreams-upset       1854.2370 1.0003545
#> dreams-physior      989.4222 1.0035589
#> flash-upset         626.9167 0.9997884
#> flash-physior             NA 1.0001972
#> upset-physior             NA 1.0010886
```

- R-hat values close to 1 (typically below 1.01) suggest convergence
  ([Vehtari et al., 2021](#ref-VehtariEtAl_2021)).
- The effective sample size (ESS) reflects the number of independent
  samples that would provide equivalent precision. Larger ESS values
  indicate more reliable estimates.
- The Monte Carlo standard error (MCSE) measures the additional
  variability introduced by using a finite number of MCMC draws. A small
  MCSE relative to the posterior standard deviation indicates stable
  estimates, whereas a large MCSE suggests that more samples are needed.

### Two ESS measures for edge-selected parameters

With edge or difference selection active, the effect parameters are
governed by spike-and-slab priors. The corresponding parameter is set to
exactly zero when the effect is excluded, rather than being removed from
the model. Because the parameter has a well-defined value at every
iteration, the full chain — including zeros — is a valid sequence for
computing ESS.

- **n_eff** is the unconditional ESS, computed from the full effect
  chain. It measures how precisely the overall posterior mean is
  estimated.
- **n_eff_mixt** is the mixture ESS. It measures how precisely the
  posterior mean of the effect is estimated while accounting for the
  additional uncertainty introduced by the spike-and-slab selection.
  When the indicator rarely switches between inclusion and exclusion
  (fewer than 5 transitions), `n_eff_mixt` is suppressed in the printed
  output.

### Traceplots

Users can inspect traceplots by extracting raw samples directly. Here is
an example for the pairwise effect parameter.

``` r

param_index = 1
chains = fit$raw_samples$pairwise
nchains = length(chains)
cols = c("firebrick", "steelblue", "darkgreen", "goldenrod")

plot(chains[[1]][, param_index],
  type = "l", col = cols[1],
  xlab = "Iteration", ylab = "Value",
  main = "Traceplot of pairwise[1]",
  ylim = range(sapply(chains, function(ch) range(ch[, param_index])))
)
if(nchains > 1) {
  for(c in 2:nchains) {
    lines(chains[[c]][, param_index], col = cols[c])
  }
}
```

![](diagnostics_files/figure-html/unnamed-chunk-5-1.png)

## Spike-and-slab summaries

The spike-and-slab prior yields posterior inclusion probabilities for
edges:

``` r

coef(fit)$indicator
#>           intrusion  dreams   flash   upset physior
#> intrusion   0.00000 1.00000 1.00000 0.93550 0.97825
#> dreams      1.00000 0.00000 1.00000 0.99900 0.05525
#> flash       1.00000 1.00000 0.00000 0.08725 1.00000
#> upset       0.93550 0.99900 0.08725 0.00000 1.00000
#> physior     0.97825 0.05525 1.00000 1.00000 0.00000
```

- Values near 1.0: strong evidence the edge is present.
- Values near 0.0: strong evidence the edge is absent.
- Values near 0.5: inconclusive (absence of evidence).

## Bayes factors

When the prior inclusion probability for an edge is equal to 0.5 (e.g.,
using `edge_prior = bernoulli_prior(0.5)`, the default, or a symmetric
Beta prior `edge_prior = beta_bernoulli_prior(alpha, beta)` with
`alpha == beta`), we can directly transform inclusion probabilities into
Bayes factors for edge presence vs absence:

``` r

# Example for one edge
p = coef(fit)$indicator[1, 5]
BF_10 = p / (1 - p)
BF_10
#> [1] 44.97701
```

Here the Bayes factor in favor of inclusion (H1) is small, meaning that
there is little evidence for inclusion. Since the Bayes factor is
transitive, we can use it to express the evidence in favor of exclusion
(H0) as

``` r

1 / BF_10
#> [1] 0.02223358
```

This Bayes factor shows that there is strong evidence for the absence of
a network relation between the variables `intrusion` and `physior`.

## NUTS diagnostics

When using `update_method = "nuts"` (the default), additional
diagnostics are available to assess the quality of the Hamiltonian Monte
Carlo sampling. These can be accessed via `fit$nuts_diag`:

``` r

fit$nuts_diag$summary
#> $total_divergences
#> [1] 0
#> 
#> $total_non_reversible
#> [1] 0
#> 
#> $max_tree_depth_hits
#> [1] 0
#> 
#> $min_ebfmi
#> [1] 0.9441257
#> 
#> $mean_accept_prob
#> [1] 0.8453192
#> 
#> $warmup_incomplete
#> [1] TRUE
```

### E-BFMI

E-BFMI (Energy Bayesian Fraction of Missing Information) measures how
efficiently the sampler explores the posterior. It compares the typical
size of energy changes between successive samples to the overall spread
of energies. Values close to 1 indicate that the sampler moves freely
across the energy landscape; values below 0.3 suggest the sampler may be
getting stuck or that the chain has not yet settled into its stationary
distribution.

A low E-BFMI does not necessarily mean your results are wrong, but it
does warrant further investigation. In models with edge selection, the
most common cause is that the warmup period was too short for the
discrete graph structure to equilibrate. Increasing `warmup` often
resolves this.

### Divergent transitions

Divergent transitions occur when the numerical integrator encounters
regions of the posterior where the curvature changes too rapidly for the
current step size. A small number of divergences (say, fewer than 0.1%
of samples) is generally acceptable. However, many divergences indicate
that the sampler may be missing important parts of the posterior.

If you see a large number of divergences, consider increasing
`target_accept` (which makes the sampler use a smaller step size) and,
if this does not fix it, switching to
`update_method = "adaptive-metropolis"`.

### Tree depth

NUTS builds trajectories by repeatedly doubling their length until a
“U-turn” criterion is satisfied. If the trajectory frequently reaches
the maximum allowed depth (`nuts_max_depth`, default 10), it suggests
the sampler may benefit from longer trajectories to explore the
posterior efficiently. Hitting the maximum depth occasionally is normal;
hitting it on most iterations may indicate challenging posterior
geometry. If this happens, consider increasing `nuts_max_depth`.

### Non-reversible steps

For MRFs with continuous variables, the leapfrog integrator enforces
equality constraints through a projection step. After each forward step,
the integrator checks whether reversing the step returns to the starting
point. When the round-trip error exceeds a tolerance scaled by the
square of the step size, the step is flagged as non-reversible.

A small number of non-reversible steps is not a concern. A large number
indicates that the step size is too large for the constraint geometry.
Because the step size is tuned during warmup, the most effective remedy
is to increase `warmup` so the adapter has more time to find an
appropriate step size. If non-reversible steps persist after increasing
warmup, switch to `update_method = "adaptive-metropolis"`.

### Warmup and equilibration

Standard HMC/NUTS warmup is designed to tune the step size and mass
matrix for the continuous parameters. In models with edge selection, the
discrete graph structure may take longer to reach its stationary
distribution than the continuous parameters. As a result, even after
warmup completes, the first portion of the sampling phase may still show
transient behavior (i.e., non-stationarity).

The `warmup_check` component provides simple diagnostics that compare
the first and second halves of the post-warmup samples:

``` r

fit$nuts_diag$warmup_check
#> $warmup_incomplete
#> [1] FALSE  TRUE
#> 
#> $energy_slope
#>      time_idx      time_idx 
#>  0.0005399232 -0.0008251581 
#> 
#> $slope_significant
#> time_idx time_idx 
#>    FALSE     TRUE 
#> 
#> $ebfmi_first_half
#> [1] 0.9885984 0.9539132
#> 
#> $ebfmi_second_half
#> [1] 0.9438844 0.9502778
#> 
#> $var_ratio
#> [1] 0.8858937 1.0602087
```

The returned list contains the following fields (one value per chain):

- **warmup_incomplete**: A logical flag that is `TRUE` when any of the
  indicators below suggest the chain may not have reached stationarity.
- **energy_slope**: The slope of a linear regression of energy against
  iteration number. A slope near zero indicates stable energy; a
  significant negative slope suggests the chain is still drifting toward
  higher-probability regions.
- **slope_significant**: `TRUE` if the energy slope is statistically
  significant (p \< 0.01).
- **ebfmi_first_half** and **ebfmi_second_half**: E-BFMI computed
  separately for the first and second halves of the post-warmup samples.
  If the first-half value is much lower (for example, below 0.3) while
  the second-half value is healthy, the early samples were likely still
  settling.
- **var_ratio**: The ratio of energy variance in the first half to that
  in the second half. A ratio much greater than 1 (for example, above 2)
  indicates higher variability early on, consistent with transient
  behavior.

If these diagnostics suggest the chain was still settling, increase
`warmup` and re-run the model. If diagnostics remain problematic after a
substantial increase (for example, doubling or tripling `warmup`),
consider re-fitting with `update_method = "adaptive-metropolis"` and
comparing the posterior summaries. If the two samplers produce similar
results, the estimates are likely trustworthy despite the warnings; if
they differ substantially, that warrants further investigation of the
model or data.

## Next steps

- See *Getting Started* for a simple one-sample workflow.
- See *Model Comparison* for group differences.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved
$`\hat{R}`$ for assessing convergence of MCMC. *Bayesian Analysis*,
*16*(2), 667–718. <https://doi.org/10.1214/20-BA1221>
