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
#> intrusion-dreams  0.315073338 0.0007021890 0.031597763 2024.9052
#> intrusion-flash   0.168341600 0.0007134522 0.030973093 1884.6878
#> intrusion-upset   0.097807766 0.0025927525 0.035583797  249.4741
#> intrusion-physior 0.095269549 0.0035382828 0.038442828  251.6125
#> dreams-flash      0.248862725 0.0007313147 0.029855210 1666.6011
#> dreams-upset      0.113270460 0.0006943014 0.026555486 1462.8924
#> dreams-physior    0.002065929 0.0003823586 0.008723002  666.8232
#> flash-upset       0.006060051 0.0010190478 0.016802191  250.5115
#> flash-physior     0.154306636 0.0007548889 0.026856317 1265.6886
#> upset-physior     0.355636928 0.0007899810 0.029262067 1372.0713
#>                   n_eff_mixt      Rhat
#> intrusion-dreams          NA 1.0077927
#> intrusion-flash           NA 1.0006146
#> intrusion-upset     188.3573 1.0009777
#> intrusion-physior   118.0444 1.0500111
#> dreams-flash              NA 0.9999733
#> dreams-upset              NA 1.0014495
#> dreams-physior      520.4634 1.0109284
#> flash-upset         271.8584 1.0792918
#> flash-physior             NA 1.0004244
#> upset-physior             NA 0.9995018
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
#>           intrusion dreams  flash  upset physior
#> intrusion     0.000 1.0000 1.0000 0.9530  0.9220
#> dreams        1.000 0.0000 1.0000 1.0000  0.0545
#> flash         1.000 1.0000 0.0000 0.1195  1.0000
#> upset         0.953 1.0000 0.1195 0.0000  1.0000
#> physior       0.922 0.0545 1.0000 1.0000  0.0000
```

- Values near 1.0: strong evidence the edge is present.
- Values near 0.0: strong evidence the edge is absent.
- Values near 0.5: inconclusive (absence of evidence).

## Bayes factors

When the prior inclusion probability for an edge is equal to 0.5 (e.g.,
using a Bernoulli prior with `inclusion_probability = 0.5` or a
symmetric Beta prior, `main_alpha = main_beta`), we can directly
transform inclusion probabilities into Bayes factors for edge presence
vs absence:

``` r
# Example for one edge
p = coef(fit)$indicator[1, 5]
BF_10 = p / (1 - p)
BF_10
#> [1] 11.82051
```

Here the Bayes factor in favor of inclusion (H1) is small, meaning that
there is little evidence for inclusion. Since the Bayes factor is
transitive, we can use it to express the evidence in favor of exclusion
(H0) as

``` r
1 / BF_10
#> [1] 0.0845987
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
#> [1] 0.9681427
#> 
#> $mean_accept_prob
#> [1] 0.9024846
#> 
#> $warmup_incomplete
#> [1] FALSE
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
#> [1] FALSE FALSE
#> 
#> $energy_slope
#>      time_idx      time_idx 
#> -0.0000385324 -0.0005017369 
#> 
#> $slope_significant
#> time_idx time_idx 
#>    FALSE    FALSE 
#> 
#> $ebfmi_first_half
#> [1] 0.9220304 1.0638666
#> 
#> $ebfmi_second_half
#> [1] 1.027856 1.008601
#> 
#> $var_ratio
#> [1] 1.109731 0.968445
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
$\widehat{R}$ for assessing convergence of MCMC. *Bayesian Analysis*,
*16*(2), 667–718. <https://doi.org/10.1214/20-BA1221>
