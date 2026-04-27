# NA

## Description

Provide a concise summary of the changes.

### Problem / Motivation

If this is a bug fix, describe the specific problem or edge case being
addressed. If this is a methodological update, explain the theoretical
gap or need for the change.

If helpful, add a minimal code example that demonstrates the problem.

``` r
# Problem demonstration
```

### Proposed Changes / New Functionality

Describe the new functionality or logic introduced. Note whether this
affects prior specifications, sampling algorithms in `src/`, summary
outputs, or user facing R functions.

Fixes \#

If helpful, add a minimal code example that demonstrates the new
behaviour.

``` r
# New behaviour
```

## Type of Change

Bug fix

New feature

Breaking change

Documentation update

Performance optimization

Statistical / methodological update

## Documentation and Release Notes

Updated function documentation for any signature or behaviour changes

Regenerated `man/*.Rd` files if roxygen comments changed

Added or updated `NEWS.md` entry if the change affects users

## Testing and Validation

Describe the checks you ran and the main results.

Added or updated tests in `tests/testthat/` when needed

Ran the project styler

``` r
source("inst/styler/bgms_style.R")
styler::style_pkg(style = bgms_style)
```

Ran `lintr::lint_package()`

Ran
[`roxygen2::roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html)
if roxygen comments changed

Ran `rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))` for
non-trivial changes

Verified numerical behaviour where relevant

## Additional Notes

List any reviewer context that is important for evaluation: open design
questions, follow-up work, expected numerical differences, or known
limits.
