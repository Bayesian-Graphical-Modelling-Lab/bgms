# Normal Prior for Model Parameters

Specifies a Normal(0, scale) prior on model parameters. Produces
lighter-tailed shrinkage than the Cauchy prior and is better suited for
simulation-based calibration (SBC) studies. Can be used for
interactions, thresholds, or continuous means.

## Usage

``` r
normal_prior(scale = 1)
```

## Arguments

- scale:

  Positive numeric. Standard deviation of the normal distribution.
  Default: `1`.

## Value

An object of class `"bgms_parameter_prior"` with `family = "normal"`.

## See also

[`cauchy_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`beta_prime_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)

## Examples

``` r
normal_prior()
#> Parameter prior: Normal(0, 1)
normal_prior(scale = 0.5)
#> Parameter prior: Normal(0, 0.5)
```
