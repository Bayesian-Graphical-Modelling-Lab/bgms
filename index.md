![bgms](reference/figures/bgms-banner.svg)

**Bayesian analysis of graphical models**

The **bgms** package provides Bayesian estimation and edge selection for
Markov random field models of mixed binary, ordinal, and continuous
variables. The variable types in the data determine the model: an
**ordinal MRF** for ordinal data, a **Gaussian graphical model** for
continuous data, or a **mixed MRF** combining both. Posterior inference
uses Markov chain Monte Carlo, combining a Metropolis approach for
between-model moves (i.e., edge selection) with the No-U-Turn sampler
for within-model parameter updates. The package supports both
single-threaded and parallel chains, and uses a C++ backend for
computational efficiency.

## Main functions

- [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  — estimate a graphical model in a one-sample design.
- [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  — compare graphical models between groups.

Both functions support **edge selection** via spike-and-slab priors,
yielding posterior inclusion probabilities for each edge.
[`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
can additionally model **community structure**, and
[`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
can test for **group differences** in individual parameters.

## Installation

Install from CRAN:

``` r
install.packages("bgms")
```

Or install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("Bayesian-Graphical-Modelling-Lab/bgms")
```

## Citation

If you use bgms in your work, please cite:

1.  Marsman, M., van den Bergh, D., & Haslbeck, J. M. B. (2025).
    Bayesian analysis of the ordinal Markov random field.
    *Psychometrika*, *90*(1), 146–182.
    [![DOI:10.1017/psy.2024.4](https://zenodo.org/badge/DOI/10.1017/psy.2024.4.svg)](https://doi.org/10.1017/psy.2024.4)

Related methodological papers:

2.  Sekulovski, N., Keetelaar, S., Huth, K. B. S., Wagenmakers, E.-J.,
    van Bork, R., van den Bergh, D., & Marsman, M. (2024). Testing
    conditional independence in psychometric networks: An analysis of
    three Bayesian methods. *Multivariate Behavioral Research*, *59*,
    913–933.
    [![DOI:10.1080/00273171.2024.2345915](https://zenodo.org/badge/DOI/10.1080/00273171.2024.2345915.svg)](https://doi.org/10.1080/00273171.2024.2345915)

3.  Marsman, M., Waldorp, L. J., Sekulovski, N., & Haslbeck, J. M. B.
    (2025). Bayes factor tests for group differences in ordinal and
    binary graphical models. *Psychometrika*, *90*(5), 1809–1842.
    [![DOI:10.1017/psy.2025.10060](https://zenodo.org/badge/DOI/10.1017/psy.2025.10060.svg)](https://doi.org/10.1017/psy.2025.10060)

4.  Sekulovski, N., Arena, G., Haslbeck, J. M. B., Huth, K. B. S.,
    Friel, N., & Marsman, M. (2025). A stochastic block prior for
    clustering in graphical models.
    [![PsyArXiv:29p3m](https://img.shields.io/badge/PsyArXiv-29p3m-blue.svg)](https://osf.io/preprints/psyarxiv/29p3m_v3)

You can also retrieve the citation from R:

``` r
citation("bgms")
```

## Contributing

Contributions are welcome. See
[CONTRIBUTING.md](https://bayesian-graphical-modelling-lab.github.io/bgms/CONTRIBUTING.md)
for how to get started.

## Code of Conduct

This project follows the [Contributor Covenant Code of
Conduct](https://bayesian-graphical-modelling-lab.github.io/bgms/CODE_OF_CONDUCT.md).
