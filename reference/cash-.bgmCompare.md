# Access elements of a bgmCompare object

Intercepts access to `posterior_summary_*` fields and triggers lazy
computation from cache when needed. All other fields pass through using
standard list extraction.

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
