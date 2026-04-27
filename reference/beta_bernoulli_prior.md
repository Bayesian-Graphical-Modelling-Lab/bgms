# Beta-Bernoulli Prior for Inclusion Indicators

Specifies a Beta-Bernoulli prior for inclusion indicators. The inclusion
probability is drawn from a \\\textrm{Beta}(\alpha, \beta)\\
distribution and shared across all edges.

## Usage

``` r
beta_bernoulli_prior(alpha = 1, beta = 1)
```

## Arguments

- alpha:

  Positive numeric. First shape parameter of the Beta distribution.
  Default: `1`.

- beta:

  Positive numeric. Second shape parameter of the Beta distribution.
  Default: `1`.

## Value

An object of class `"bgms_indicator_prior"` with
`family = "Beta-Bernoulli"`.

## See also

[`bernoulli_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`sbm_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`sbm_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/sbm_prior.md)

## Examples

``` r
beta_bernoulli_prior()
#> Edge prior: Beta-Bernoulli(alpha = 1, beta = 1)
beta_bernoulli_prior(alpha = 2, beta = 5)
#> Edge prior: Beta-Bernoulli(alpha = 2, beta = 5)
```
