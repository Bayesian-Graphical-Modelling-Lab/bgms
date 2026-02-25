# Contributing to bgms

Thank you for your interest in contributing to bgms. This document
covers the practical steps for building, testing, and submitting
changes.

## Getting started

```bash
git clone https://github.com/Bayesian-Graphical-Modelling-Lab/bgms.git
cd bgms
```

Build and check the package:

```r
devtools::document()
devtools::check()
```

## Code style

bgms uses `=` for assignment (not `<-`) and omits the space between
`if`/`for`/`while` and `(`. These rules are enforced by the custom
styler in `inst/styler/bgms_style.R`. Run before committing:

```r
source("inst/styler/bgms_style.R")
styler::style_pkg(style = bgms_style)
```

## Documentation

- **Exported R functions** get full roxygen2 blocks (`@title`,
  `@description`, `@param`, `@return`, `@examples`, `@seealso`,
  `@export`).
- **Internal R functions** use plain `#` comment blocks, not `#'`
  roxygen.
- **C++ headers** use Doxygen `/** */` blocks on classes and methods.
- **C++ implementations** use inline `//` comments only.

When modifying a function signature, update its documentation in the
same commit. When adding a new exported function, add it to
`_pkgdown.yml`.

## Tests

Run the test suite with:

```r
devtools::test()
```

Test files live in `tests/testthat/`. Each `test_that()` description
should be an imperative sentence.

## CI

Pull requests are checked automatically by GitHub Actions:

- `R-CMD-check.yaml` — builds and checks on multiple platforms
- `test-coverage.yaml` — reports test coverage
- `pkgdown.yaml` — builds the documentation site

## Commit messages

Use the format: `type: short description`. Common types: `feat`,
`fix`, `docs`, `refactor`, `test`.

## Contributors

All contributors are acknowledged in `inst/CONTRIBUTORS.md`.
