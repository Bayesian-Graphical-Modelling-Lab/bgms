# Cauchy Prior for Model Parameters

Specifies a Cauchy(0, scale) prior on model parameters. This is the
default prior for pairwise interactions in
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
and produces heavy-tailed shrinkage toward zero.

## Usage

``` r
cauchy_prior(scale = 1)
```

## Arguments

- scale:

  Positive numeric. Scale (half-width at half-maximum) of the Cauchy
  distribution. Default: `1`.

## Value

An object of class `"bgms_parameter_prior"` with `family = "cauchy"`.

## See also

[`normal_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`beta_prime_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)

## Examples

``` r
cauchy_prior()
#> Parameter prior: Cauchy(0, 1)
cauchy_prior(scale = 2.5)
#> Parameter prior: Cauchy(0, 2.5)
```
