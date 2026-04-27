# Beta-Prime Prior for Model Parameters

Specifies a beta-prime prior on model parameters. The parameterization
follows the logistic transformation: \\\sigma(\mu) \sim
\textrm{Beta}(\alpha, \beta)\\, so \\\mu = \textrm{logit}(Y)\\ where \\Y
\sim \textrm{Beta}(\alpha, \beta)\\.

## Usage

``` r
beta_prime_prior(alpha = 0.5, beta = 0.5)
```

## Arguments

- alpha:

  Positive numeric. First shape parameter. Default: `0.5`.

- beta:

  Positive numeric. Second shape parameter. Default: `0.5`.

## Value

An object of class `"bgms_parameter_prior"` with
`family = "beta-prime"`.

## See also

[`cauchy_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`normal_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)

## Examples

``` r
beta_prime_prior()
#> Parameter prior: Beta-prime(alpha = 0.5, beta = 0.5)
beta_prime_prior(alpha = 1, beta = 1)
#> Parameter prior: Beta-prime(alpha = 1, beta = 1)
```
