# Access elements of a bgms object

Intercepts access to `posterior_summary_*` fields and triggers lazy
computation from cache when needed. All other fields pass through using
standard list extraction.

## Usage

``` r
# S3 method for class 'bgms'
x$name

# S3 method for class 'bgms'
x[[name, ...]]
```

## Arguments

- x:

  A `bgms` object.

- name:

  Name of the element to access.

- ...:

  Ignored.

## Value

The requested element.
