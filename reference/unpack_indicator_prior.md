# Unpack an edge prior into the flat parameters used by bgm_spec

Unpack an edge prior into the flat parameters used by bgm_spec

## Usage

``` r
unpack_indicator_prior(prior, num_variables)
```

## Arguments

- prior:

  A `bgms_indicator_prior` object.

- num_variables:

  Integer. Number of variables (for inclusion matrix).

## Value

A list matching the fields expected by `validate_edge_prior` output.
