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

Note: During fitting, progress bars are shown in interactive sessions.
In this vignette, they are suppressed for clarity. Sampling can take a
while; the progress bars usually help track progress.

## Convergence diagnostics

The quality of the Markov chain can be assessed with common MCMC
diagnostics:

``` r
summary(fit)$pairwise
#>                          mean          sd        mcse     n_eff
#> intrusion-dreams  0.631924323 0.001551325 0.064266150 1716.1643
#> intrusion-flash   0.338297975 0.001448421 0.061854243 1823.6826
#> intrusion-upset   0.190695963 0.076279341 0.005921381  165.9464
#> intrusion-physior 0.198176315 0.065555747 0.003647515  323.0185
#> dreams-flash      0.498040695 0.001270565 0.060466106 2264.8025
#> dreams-upset      0.230776860 0.056205053 0.002092867  721.2197
#> dreams-physior    0.005254907 0.021861327 0.000843593  671.5636
#> flash-upset       0.006462176 0.024675968 0.001036232  567.0671
#> flash-physior     0.307138582 0.001205035 0.053371694 1961.6559
#> upset-physior     0.707867139 0.001478310 0.059720575 1631.9869
#>                        Rhat
#> intrusion-dreams  0.9998582
#> intrusion-flash   0.9999985
#> intrusion-upset   1.0325834
#> intrusion-physior 1.0053426
#> dreams-flash      0.9999931
#> dreams-upset      1.0076821
#> dreams-physior    1.0051407
#> flash-upset       1.0115874
#> flash-physior     1.0013808
#> upset-physior     0.9999992
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
#>           intrusion  dreams   flash   upset physior
#> intrusion   0.00000 1.00000 1.00000 0.96575  0.0560
#> dreams      1.00000 0.00000 0.92525 1.00000  0.0655
#> flash       1.00000 0.92525 0.00000 0.99600  1.0000
#> upset       0.96575 1.00000 0.99600 0.00000  1.0000
#> physior     0.05600 0.06550 1.00000 1.00000  0.0000
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
#> [1] 0.05932203
```

Here the Bayes factor in favor of inclusion (H1) is small, meaning that
there is little evidence for inclusion. Since the Bayes factor is
transitive, we can use it to express the evidence in favor of exclusion
(H0) as

``` r
1 / BF_10
#> [1] 16.85714
```

This Bayes factor shows that there is strong evidence for the absence of
a network relation between the variables `intrusion` and `physior`.

## Notes on runtime

- Sampling with spike-and-slab priors can take longer.
- In interactive sessions, progress bars are displayed. In this
  vignette, they are suppressed for readability.

## Next steps

- See *Getting Started* for a simple one-sample workflow.
- See *Model Comparison* for group differences.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved
$\widehat{R}$ for assessing convergence of MCMC. *Bayesian Analysis*,
*16*(2), 667–718. <https://doi.org/10.1214/20-BA1221>
