# Unpack a threshold prior into the flat parameters used by bgm_spec

Unpack a threshold prior into the flat parameters used by bgm_spec

## Usage

``` r
unpack_threshold_prior(prior)
```

## Arguments

- prior:

  A `bgms_parameter_prior` object.

## Value

A list with `threshold_prior_type` (character), `main_alpha`,
`main_beta`, and `threshold_scale`.
