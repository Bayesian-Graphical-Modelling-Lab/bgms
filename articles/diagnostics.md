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
#>                          mean          sd        mcse     n_eff
#> intrusion-dreams  0.629338089 0.001804519 0.067102463 1382.7822
#> intrusion-flash   0.338657791 0.001589875 0.060421205 1444.2846
#> intrusion-upset   0.202376428 0.064014992 0.004668924  187.9879
#> intrusion-physior 0.188939686 0.071198984 0.006082051  137.0400
#> dreams-flash      0.499816233 0.001803795 0.060493065 1124.7000
#> dreams-upset      0.225595383 0.053378801 0.001589095 1128.3344
#> dreams-physior    0.007497679 0.025566651 0.001247486  420.0263
#> flash-upset       0.004566734 0.019291950 0.001055790  333.8854
#> flash-physior     0.308420955 0.001426490 0.053789058 1421.8396
#> upset-physior     0.707339037 0.001706434 0.059781270 1227.3022
#>                        Rhat
#> intrusion-dreams  1.0047113
#> intrusion-flash   1.0024625
#> intrusion-upset   1.0026218
#> intrusion-physior 1.0002991
#> dreams-flash      1.0018317
#> dreams-upset      0.9999766
#> dreams-physior    1.0020293
#> flash-upset       0.9999943
#> flash-physior     0.9998950
#> upset-physior     1.0047092
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
#> intrusion    0.0000 1.0000 1.0000 0.9790  0.9405
#> dreams       1.0000 0.0000 1.0000 0.9985  0.0820
#> flash        1.0000 1.0000 0.0000 0.0545  1.0000
#> upset        0.9790 0.9985 0.0545 0.0000  1.0000
#> physior      0.9405 0.0820 1.0000 1.0000  0.0000
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
#> [1] 15.80672
```

Here the Bayes factor in favor of inclusion (H1) is small, meaning that
there is little evidence for inclusion. Since the Bayes factor is
transitive, we can use it to express the evidence in favor of exclusion
(H0) as

``` r
1 / BF_10
#> [1] 0.06326422
```

This Bayes factor shows that there is strong evidence for the absence of
a network relation between the variables `intrusion` and `physior`.

## Next steps

- See *Getting Started* for a simple one-sample workflow.
- See *Model Comparison* for group differences.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved
$\widehat{R}$ for assessing convergence of MCMC. *Bayesian Analysis*,
*16*(2), 667–718. <https://doi.org/10.1214/20-BA1221>
