# Stochastic Block Model Prior for Inclusion Indicators

Specifies a Stochastic Block Model (SBM) prior for inclusion indicators.
Variables are assigned to latent clusters, with separate Beta priors on
within-cluster and between-cluster inclusion probabilities.

## Usage

``` r
sbm_prior(
  alpha = 1,
  beta = 1,
  alpha_between = 1,
  beta_between = 1,
  dirichlet_alpha = 1,
  lambda = 1
)
```

## Arguments

- alpha:

  Positive numeric. First shape parameter of the Beta distribution for
  within-cluster edges. Default: `1`.

- beta:

  Positive numeric. Second shape parameter of the Beta distribution for
  within-cluster edges. Default: `1`.

- alpha_between:

  Positive numeric. First shape parameter of the Beta distribution for
  between-cluster edges. Default: `1`.

- beta_between:

  Positive numeric. Second shape parameter of the Beta distribution for
  between-cluster edges. Default: `1`.

- dirichlet_alpha:

  Positive numeric. Concentration parameter of the Dirichlet prior on
  cluster assignments. Default: `1`.

- lambda:

  Positive numeric. Rate parameter of the zero-truncated Poisson prior
  on the number of clusters. Default: `1`.

## Value

An object of class `"bgms_indicator_prior"` with
`family = "Stochastic-Block"`.

## See also

[`bernoulli_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_bernoulli_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

Other prior-constructors:
[`bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bernoulli_prior.md),
[`beta_bernoulli_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_bernoulli_prior.md),
[`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md),
[`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md)

## Examples

``` r
sbm_prior()
#> Edge prior: Stochastic-Block
#>   Within:    Beta(1, 1)
#>   Between:   Beta(1, 1)
#>   Dirichlet: 1, Lambda: 1
sbm_prior(alpha = 2, beta = 1, alpha_between = 1, beta_between = 5)
#> Edge prior: Stochastic-Block
#>   Within:    Beta(2, 1)
#>   Between:   Beta(1, 5)
#>   Dirichlet: 1, Lambda: 1
```
