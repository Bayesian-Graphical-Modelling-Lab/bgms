# Exponential Prior for Scale Parameters

Specifies an Exponential(rate) prior for positive scale parameters. This
is a convenience function equivalent to
`gamma_prior(shape = 1, rate = rate)`.

## Usage

``` r
exponential_prior(rate = 1)
```

## Arguments

- rate:

  Positive numeric. Rate parameter of the Exponential distribution.
  Default: `1`.

## Value

An object of class `"bgms_scale_prior"` with `family = "exponential"`.

## See also

[`gamma_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)

## Examples

``` r
exponential_prior()
#> Scale prior: Exponential(rate = 1)
exponential_prior(rate = 2)
#> Scale prior: Exponential(rate = 2)
```
