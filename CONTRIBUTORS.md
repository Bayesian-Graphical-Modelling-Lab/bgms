# Contributors to bgms

We are grateful to everyone who has contributed to the development of
**bgms**. This document acknowledges all contributions.

## Core Development Team

### Maarten Marsman (Lead Developer & Maintainer)

- Package design and architecture
- Pseudolikelihood inference framework
- Core MCMC sampling methods
- Model implementation and extensive testing
- Primary contact: <m.marsman@uva.nl>

### Don van den Bergh (Author)

- **Code Infrastructure & Refactoring**: Project-wide code refactoring
  including C++ best practices, namespace management, and architecture
  improvements for maintainability and reliability
- **Cross-Platform Build System**: Windows compilation debugging and
  fixes, resolving platform-specific issues with include ordering, type
  handling, and compiler toolchain integration; R CMD CHECK compliance
  and Makevars configuration
- **Numerical Stability & Performance**: Custom exp/log implementations
  for numerical robustness, memory efficiency optimizations, and
  algorithm-level performance improvements
- **Parallel Computing Infrastructure**: Design and implementation of
  parallel RNG infrastructure, progress bar mechanisms, and user
  interrupt handling for multi-chain sampling
- **Package Infrastructure & Testing**: Setup of pkgdown documentation
  pipeline, testing framework infrastructure, GitHub Actions workflows
  for automated testing and documentation, and code styling integration

## Contributors

### Nikola Sekulovski

- Refinements and improvements to the stochastic block model (SBM) edge
  prior
- SBM output functions and extraction methods
- SBM-related bug fixes and enhancements

### Giuseppe Arena

- Stochastic block model (SBM) prior code optimization and refactoring
- Bug fixes in SBM implementation
- Corrections to SBM output functions

### Laura Groot

- Preparation and curation of ADHD and Boredom datasets
- Dataset documentation and Roxygen help files
- Bibliography and citation management for datasets

### Gali Geller

- Robustness and tolerance testing for bgms and bgmCompare outputs
- Test infrastructure for range validation and symmetry checks

## Acknowledgments

We acknowledge Karoline Huth for contributions to the early development
of the bgms package.

## Attribution Standards

This package follows R community standards for author attribution using
the `Authors@R` field and role codes in the DESCRIPTION file:

- **`"aut"` (Author):** Substantial contributions; appears in package
  citation
- **`"cre"` (Creator/Maintainer):** Package maintainer
- **`"ctb"` (Contributor):** Smaller contributions (patches, bug fixes);
  does not appear in package citation

See the
[DESCRIPTION](https://bayesian-graphical-modelling-lab.github.io/bgms/DESCRIPTION)
file for complete role assignments. Standards follow: - Hornik, K.,
Murdoch, D., & Zeileis, A. (2012). Who Did What? The Roles of R Package
Authors and How to Refer to Them. *The R Journal*, 4(1), 64-69.
<https://doi.org/10.32614/RJ-2012-009> - R Core Team. *Writing R
Extensions*.
<https://cran.r-project.org/doc/manuals/r-release/R-exts.html>
