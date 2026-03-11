# Extract Coefficients from a bgmCompare Object

Returns posterior means for raw parameters (baseline + differences) and
group-specific effects from a `bgmCompare` fit, as well as inclusion
indicators.

## Usage

``` r
# S3 method for class 'bgmCompare'
coef(object, ...)
```

## Arguments

- object:

  An object of class `bgmCompare`.

- ...:

  Ignored.

## Value

A list with components:

- main_effects_raw:

  Posterior means of the raw main-effect parameters (variables x
  (baseline + differences)).

- pairwise_effects_raw:

  Posterior means of the raw pairwise-effect parameters (pairs x
  (baseline + differences)).

- main_effects_groups:

  Posterior means of group-specific main effects (variables x groups),
  computed as baseline plus projected differences.

- pairwise_effects_groups:

  Posterior means of group-specific pairwise effects (pairs x groups),
  computed as baseline plus projected differences.

- indicators:

  Posterior mean inclusion probabilities as a symmetric matrix, with
  diagonals corresponding to main effects and off-diagonals to pairwise
  effects.

## See also

[`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md),
[`print.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgmCompare.md),
[`summary.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgmCompare.md)

Other posterior-methods:
[`coef.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/coef.bgms.md),
[`print.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgmCompare.md),
[`print.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/print.bgms.md),
[`summary.bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgmCompare.md),
[`summary.bgms()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/summary.bgms.md)

## Examples

``` r
# \donttest{
# See ?bgmCompare for a full example
# }
```
