# Unpack a parameter prior into the flat parameters used by bgm_spec

Unpack a parameter prior into the flat parameters used by bgm_spec

## Usage

``` r
unpack_parameter_prior(prior)
```

## Arguments

- prior:

  A `bgms_parameter_prior` object.

## Value

A list with `prior_type` (character), `scale` (numeric), `alpha`
(numeric), and `beta` (numeric). Unused hyperparameters are set to
`NA_real_`.
