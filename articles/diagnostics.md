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
#>                          mean          sd         mcse     n_eff
#> intrusion-dreams  0.631353688 0.001391632 0.0679125228 2381.5024
#> intrusion-flash   0.338976477 0.001163740 0.0615731169 2799.4335
#> intrusion-upset   0.202322026 0.066266758 0.0030047786  486.3697
#> intrusion-physior 0.189896271 0.072940538 0.0050170247  211.3710
#> dreams-flash      0.498735904 0.001296742 0.0612986788 2234.5755
#> dreams-upset      0.226326133 0.056330851 0.0018043223  974.6853
#> dreams-physior    0.006758395 0.024683811 0.0008585757  826.5465
#> flash-upset       0.004674305 0.019575299 0.0007052012  770.5320
#> flash-physior     0.307737254 0.001070298 0.0539656532 2542.2940
#> upset-physior     0.708467755 0.001283216 0.0610154116 2260.8915
#>                        Rhat
#> intrusion-dreams  1.0015037
#> intrusion-flash   1.0007634
#> intrusion-upset   1.0007649
#> intrusion-physior 1.0024785
#> dreams-flash      1.0010618
#> dreams-upset      1.0011885
#> dreams-physior    1.0018154
#> flash-upset       1.0001949
#> flash-physior     0.9998367
#> upset-physior     1.0017084
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
#> intrusion   0.00000 1.00000 1.00000 0.97325 0.93325
#> dreams      1.00000 0.00000 1.00000 0.99425 0.07175
#> flash       1.00000 1.00000 0.00000 0.05525 1.00000
#> upset       0.97325 0.99425 0.05525 0.00000 1.00000
#> physior     0.93325 0.07175 1.00000 1.00000 0.00000
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
#> [1] 13.98127
```

Here the Bayes factor in favor of inclusion (H1) is small, meaning that
there is little evidence for inclusion. Since the Bayes factor is
transitive, we can use it to express the evidence in favor of exclusion
(H0) as

``` r
1 / BF_10
#> [1] 0.07152424
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
