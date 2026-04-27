# Gamma Prior for Scale Parameters

Specifies a Gamma(shape, rate) prior for positive scale parameters such
as the diagonal elements of the precision matrix. The default
`gamma_prior(1, 1)` corresponds to an Exponential(1) distribution.

## Usage

``` r
gamma_prior(shape = 1, rate = 1)
```

## Arguments

- shape:

  Positive numeric. Shape parameter of the Gamma distribution. Default:
  `1`.

- rate:

  Positive numeric. Rate parameter of the Gamma distribution. Default:
  `1`.

## Value

An object of class `"bgms_scale_prior"` with `family = "gamma"`.

## See also

[`exponential_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)

## Examples

``` r
gamma_prior()
#> Scale prior: Gamma(shape = 1, rate = 1)
gamma_prior(shape = 2, rate = 0.5)
#> Scale prior: Gamma(shape = 2, rate = 0.5)
```
