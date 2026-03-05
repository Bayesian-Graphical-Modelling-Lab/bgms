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
#>                          mean          sd         mcse      n_eff
#> intrusion-dreams  0.631658334 0.001898891 0.0670343797 1246.22120
#> intrusion-flash   0.338183096 0.002012232 0.0669211465 1106.03956
#> intrusion-upset   0.179825453 0.090531049 0.0105310730   73.90092
#> intrusion-physior 0.195283873 0.076611368 0.0069415792  121.80633
#> dreams-flash      0.500506864 0.001723344 0.0597502784 1202.08670
#> dreams-upset      0.232857983 0.057805006 0.0025528180  512.73289
#> dreams-physior    0.003390355 0.016707941 0.0007761832  463.35876
#> flash-upset       0.009111470 0.030407574 0.0018550568  268.68853
#> flash-physior     0.309192966 0.001816074 0.0560624503  952.96569
#> upset-physior     0.709557416 0.002184435 0.0642879809  866.12682
#>                       Rhat
#> intrusion-dreams  1.003441
#> intrusion-flash   0.999552
#> intrusion-upset   1.043464
#> intrusion-physior 1.052697
#> dreams-flash      1.000466
#> dreams-upset      1.018838
#> dreams-physior    1.000111
#> flash-upset       1.001300
#> flash-physior     1.018683
#> upset-physior     1.001125
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

Advanced users can inspect traceplots by extracting raw samples and
using external packages such as `coda` or `bayesplot`. Here is an
example using the `coda` package to create a traceplot for a pairwise
effect parameter.

``` r
library(coda)

param_index = 1
chains = lapply(fit$raw_samples$pairwise, function(mat) mat[, param_index])
mcmc_obj = mcmc.list(lapply(chains, mcmc))

traceplot(mcmc_obj,
  col = c("firebrick", "steelblue", "darkgreen", "goldenrod"),
  main = "Traceplot of pairwise[1]"
)
```

![](diagnostics_files/figure-html/unnamed-chunk-5-1.png)

## Spike-and-slab summaries

The spike-and-slab prior yields posterior inclusion probabilities for
edges:

``` r
coef(fit)$indicator
#>           intrusion dreams  flash  upset physior
#> intrusion    0.0000 1.0000 1.0000 0.8575  0.9315
#> dreams       1.0000 0.0000 1.0000 0.9985  0.0405
#> flash        1.0000 1.0000 0.0000 0.0845  1.0000
#> upset        0.8575 0.9985 0.0845 0.0000  1.0000
#> physior      0.9315 0.0405 1.0000 1.0000  0.0000
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
#> [1] 13.59854
```

Here the Bayes factor in favor of inclusion (H1) is small, meaning that
there is little evidence for inclusion. Since the Bayes factor is
transitive, we can use it to express the evidence in favor of exclusion
(H0) as

``` r
1 / BF_10
#> [1] 0.07353731
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
#> $max_tree_depth_hits
#> [1] 0
#> 
#> $min_ebfmi
#> [1] 0.8760355
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
#> -1.064418e-03  8.434377e-05 
#> 
#> $slope_significant
#> time_idx time_idx 
#>    FALSE    FALSE 
#> 
#> $ebfmi_first_half
#> [1] 0.8225484 0.8692865
#> 
#> $ebfmi_second_half
#> [1] 0.9258201 1.0231384
#> 
#> $var_ratio
#> [1] 1.0301991 0.9474113
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
