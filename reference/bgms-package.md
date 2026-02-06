# bgms: Bayesian Analysis of Networks of Binary and/or Ordinal Variables

The `R` package **bgms** provides tools for Bayesian analysis of the
ordinal Markov random field (MRF), a graphical model describing networks
of binary and/or ordinal variables (Marsman et al. 2025) . The
likelihood is approximated via a pseudolikelihood, and Markov chain
Monte Carlo (MCMC) methods are used to sample from the corresponding
pseudoposterior distribution of model parameters.

The main entry points are:

- **bgm**: estimation in a one-sample design.

- **bgmCompare**: estimation and group comparison in an
  independent-sample design.

Both functions support Bayesian effect selection with spike-and-slab
priors.

- In one-sample designs, `bgm` models the presence or absence of edges
  between variables. Posterior inclusion probabilities quantify the
  plausibility of each edge and can be converted into Bayes factors for
  conditional independence tests.

- `bgm` can also model communities (clusters) of variables. The
  posterior distribution of the number of clusters provides evidence for
  or against clustering (Sekulovski et al. 2025) .

- In independent-sample designs, `bgmCompare` estimates group
  differences in edge weights and category thresholds. Posterior
  inclusion probabilities quantify the evidence for differences and can
  be converted into Bayes factors for parameter equivalence tests
  (Marsman et al. 2025) .

## Tools

The package also provides:

1.  Simulation of response data from MRFs with a Gibbs sampler
    ([`mrfSampler`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/mrfSampler.md)).

2.  Posterior estimation and edge selection in one-sample designs
    ([`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)).

3.  Posterior estimation and group-difference selection in
    independent-sample designs
    ([`bgmCompare`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)).

## Vignettes

For tutorials and worked examples, see:

- [`vignette("intro", package = "bgms")`](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/intro.md)
  — Getting started.

- [`vignette("comparison", package = "bgms")`](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/comparison.md)
  — Model comparison.

- [`vignette("diagnostics", package = "bgms")`](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/diagnostics.md)
  — Diagnostics and spike-and-slab summaries.

## References

Marsman M, van den Bergh D, Haslbeck JMB (2025). “Bayesian analysis of
the ordinal Markov random field.” *Psychometrika*, **90**(1), 146–182.
[doi:10.1017/psy.2024.4](https://doi.org/10.1017/psy.2024.4) .  
  
Marsman M, Waldorp LJ, Sekulovski N, Haslbeck JMB (2025). “Bayes factor
tests for group differences in ordinal and binary graphical models.”
*Psychometrika*, **90**(5), 1809–1842.
[doi:10.1017/psy.2025.10060](https://doi.org/10.1017/psy.2025.10060) .  
  
Sekulovski N, Arena G, Haslbeck JMB, Huth KBS, Friel N, Marsman M
(2025). “A Stochastic Block Prior for Clustering in Graphical Models.”
*Retrieved from <https://osf.io/preprints/psyarxiv/29p3m_v1>*. OSF
preprint.

## See also

Useful links:

- <https://Bayesian-Graphical-Modelling-Lab.github.io/bgms/>

- <https://github.com/Bayesian-Graphical-Modelling-Lab/bgms>

- Report bugs at
  <https://github.com/Bayesian-Graphical-Modelling-Lab/bgms/issues>

## Author

**Maintainer**: Maarten Marsman <m.marsman@uva.nl>
([ORCID](https://orcid.org/0000-0001-5309-7502))

Authors:

- Don van den Bergh ([ORCID](https://orcid.org/0000-0002-9838-7308))

Other contributors:

- Nikola Sekulovski ([ORCID](https://orcid.org/0000-0001-7032-1684))
  \[contributor\]

- Giuseppe Arena ([ORCID](https://orcid.org/0000-0001-5204-3326))
  \[contributor\]

- Laura Groot \[contributor\]

- Gali Geller \[contributor\]
