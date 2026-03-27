# Access elements of a bgmCompare object

Provides `$` access to S7 properties. Lazy `posterior_summary_*`
properties trigger computation on first access via S7 property getters.
Also supports legacy S3 list-based fit objects.

## Usage

``` r
# S3 method for class 'bgmCompare'
x$name

# S3 method for class 'bgmCompare'
x[[name, ...]]
```

## Arguments

- x:

  A `bgmCompare` object.

- name:

  Name of the element to access.

- ...:

  Ignored.

## Value

The requested element.
