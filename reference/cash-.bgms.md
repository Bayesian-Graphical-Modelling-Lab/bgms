# Access elements of a bgms object

Provides `$` access to S7 properties. Lazy `posterior_summary_*`
properties trigger computation on first access via S7 property getters.
Also supports legacy S3 list-based fit objects.

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
