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
#> intrusion-dreams  0.632293951 0.001661179 0.0657939558 1568.6981
#> intrusion-flash   0.335941200 0.001418260 0.0598391890 1780.1616
#> intrusion-upset   0.189497746 0.075619906 0.0054547304  192.1877
#> intrusion-physior 0.194802150 0.071131862 0.0049449384  206.9220
#> dreams-flash      0.502768076 0.001633580 0.0620233738 1441.5502
#> dreams-upset      0.228964481 0.060644348 0.0030175433  403.8998
#> dreams-physior    0.003975021 0.017652134 0.0006208483  808.3946
#> flash-upset       0.009480741 0.031081369 0.0013771486  509.3762
#> flash-physior     0.307821055 0.001453666 0.0560608285 1487.2689
#> upset-physior     0.710876259 0.001532448 0.0597194081 1518.6560
#>                       Rhat
#> intrusion-dreams  1.003909
#> intrusion-flash   1.000981
#> intrusion-upset   1.005080
#> intrusion-physior 1.016768
#> dreams-flash      1.000174
#> dreams-upset      1.007224
#> dreams-physior    1.000268
#> flash-upset       1.009813
#> flash-physior     1.005348
#> upset-physior     1.000656
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
#> intrusion   0.00000 1.00000 1.00000 0.92700 0.94875
#> dreams      1.00000 0.00000 1.00000 0.98475 0.04975
#> flash       1.00000 1.00000 0.00000 0.08775 1.00000
#> upset       0.92700 0.98475 0.08775 0.00000 1.00000
#> physior     0.94875 0.04975 1.00000 1.00000 0.00000
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
#> [1] 18.5122
```

Here the Bayes factor in favor of inclusion (H1) is small, meaning that
there is little evidence for inclusion. Since the Bayes factor is
transitive, we can use it to express the evidence in favor of exclusion
(H0) as

``` r
1 / BF_10
#> [1] 0.05401845
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
