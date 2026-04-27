# Bernoulli Prior for Inclusion Indicators

Specifies a Bernoulli prior for inclusion indicators with a fixed
inclusion probability. Used for edge selection in
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
and difference selection in
[`bgmCompare`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md).

## Usage

``` r
bernoulli_prior(inclusion_probability = 0.5)
```

## Arguments

- inclusion_probability:

  Numeric scalar or symmetric matrix. Prior probability of each edge
  being included. A scalar applies to all edges; a matrix allows
  edge-specific probabilities. Must be in (0, 1). Default: `0.5`.

## Value

An object of class `"bgms_indicator_prior"` with `family = "Bernoulli"`.

## See also

[`beta_bernoulli_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`sbm_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)

## Examples

``` r
bernoulli_prior()
#> Edge prior: Bernoulli(0.5)
bernoulli_prior(inclusion_probability = 0.25)
#> Edge prior: Bernoulli(0.25)
```
