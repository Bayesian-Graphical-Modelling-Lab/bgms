# bgms

[![bgms
website](reference/figures/bgms_sticker.svg)](https://bayesiangraphicalmodeling.com)

**Bayesian analysis of graphical models with binary and ordinal
variables**

The **bgms** package implements Bayesian estimation and model comparison
for **ordinal Markov random fields (MRFs)**, graphical models that
represent networks of binary and/or ordinal variables (Marsman et al.,
2025). The likelihood is approximated with a pseudolikelihood, and
Markov chain Monte Carlo (MCMC) methods are used to sample from the
corresponding pseudoposterior distribution of the model parameters.

## Main functions

The package has two main entry points:

- [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  – estimates a single network in a one-sample design.
- [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  – compares networks between groups in an independent-sample design.

## Effect selection

Both functions support **effect selection** with spike-and-slab priors:

- **Edges in one-sample designs**:
  [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  models the presence or absence of edges between variables. Posterior
  inclusion probabilities indicate the plausibility of each edge and can
  be converted into Bayes factors for conditional independence tests
  (see Marsman et al., 2025; Sekulovski et al., 2024).

- **Communities/clusters in one-sample designs**:
  [`bgm()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  can also model community structure. Posterior probabilities for the
  number of clusters quantify the plausibility of clustering solutions
  and can be converted into Bayes factors (see Sekulovski et al., 2025).

- **Group differences in independent-sample designs**:
  [`bgmCompare()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgmCompare.md)
  models differences in edge weights and category thresholds between
  groups. Posterior inclusion probabilities indicate the plausibility of
  parameter differences and can be converted into Bayes factors for
  tests of parameter equivalence (see Marsman et al., 2024).

## Learn more

For worked examples and tutorials, see the package vignettes:

- [Getting
  Started](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/intro.html)
- [Model
  Comparison](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/comparison.html)
- [Diagnostics and Spike-and-Slab
  Summaries](https://bayesian-graphical-modelling-lab.github.io/bgms/articles/diagnostics.html)

You can also access these directly from R with:

``` r
browseVignettes("bgms")
```

## Why use Markov Random Fields?

Graphical models or networks have become central in recent psychological
and psychometric research (Contreras et al., 2019; Marsman & Rhemtulla,
2022; Robinaugh et al., 2020). Most are **Markov random field (MRF)**
models, where the graph structure reflects partial associations between
variables (Kindermann & Snell, 1980).

In an MRF, a missing edge between two variables implies **conditional
independence** given the rest of the network (Lauritzen, 2004). In other
words, the remaining variables fully explain away any potential
association between the unconnected pair.

## Why use a Bayesian approach?

When analyzing an MRF, we often want to compare competing hypotheses:

- **Edge presence vs. edge absence** (conditional dependence
  vs. independence) in one-sample designs.
- **Parameter difference vs. parameter equivalence** in
  independent-sample designs.

Frequentist approaches are limited in such comparisons: they can reject
a null hypothesis, but they cannot provide evidence *for* it. As a
result, when an edge or difference is excluded, it remains unclear
whether this reflects true absence or simply insufficient power.

Bayesian inference avoids this problem. Using **inclusion Bayes
factors** (Huth et al., 2023; Sekulovski et al., 2024), we can quantify
evidence in both directions:

- **Evidence of edge presence** vs. **evidence of edge absence**, or
- **Evidence of parameter difference** vs. **evidence of parameter
  equivalence**.

This makes it possible not only to detect structure and group
differences, but also to conclude when there is an *absence of
evidence*.

## Installation

The current developmental version can be installed with

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("Bayesian-Graphical-Modelling-Lab/bgms")
```

## References

Contreras, A., Nieto, I., Valiente, C., Espinosa, R., & Vazquez, C.
(2019). The study of psychopathology from the network analysis
perspective: A systematic review. *Psychotherapy and Psychosomatics*,
*88*, 71–83. <https://doi.org/10.1159/000497425>

Huth, K., de Ron, J., Goudriaan, A. E., Luigjes, K., Mohammadi, R., van
Holst, R. J., Wagenmakers, E.-J., & Marsman, M. (2023). Bayesian
analysis of cross-sectional networks: A tutorial in R and JASP.
*Advances in Methods and Practices in Psychological Science*, *6*, 1–18.
<https://doi.org/10.1177/25152459231193334>

Kindermann, R., & Snell, J. L. (1980). *Markov random fields and their
applications* (Vol. 1). American Mathematical Society.

Lauritzen, S. L. (2004). *Graphical models*. Oxford University Press.

Marsman, M., & Rhemtulla, M. (2022). Guest editors’ introduction to the
special issue “network psychometrics in action”: Methodological
innovations inspired by empirical problems. *Psychometrika*, *87*, 1–11.
<https://doi.org/10.1007/s11336-022-09861-x>

Marsman, M., van den Bergh, D., & Haslbeck, J. M. B. (2025). Bayesian
analysis of the ordinal Markov random field. *Psychometrika*, *90*,
146–182.

Marsman, M., Waldorp, L. J., Sekulovski, N., & Haslbeck, J. M. B.
(2024). Bayes factor tests for group differences in ordinal and binary
graphical models. *Retrieved from <Https://Osf.io/Preprints/Osf/F4pk9>*.

Robinaugh, D. J., Hoekstra, R. H. A., Toner, E. R., & Borsboom, D.
(2020). The network approach to psychopathology: A review of the
literature 2008–2018 and an agenda for future research. *Psychological
Medicine*, *50*, 353–366. <https://doi.org/10.1017/S0033291719003404>

Sekulovski, N., Arena, G., Haslbeck, J. M. B., Huth, K. B. S., Friel,
N., & Marsman, M. (2025). A stochastic block prior for clustering in
graphical models. *Retrieved from
[Https://Osf.io/Preprints/Psyarxiv/29p3m_v1](https://osf.io/preprints/psyarxiv/29p3m_v1)*.

Sekulovski, N., Keetelaar, S., Huth, K. B. S., Wagenmakers, E.-J., van
Bork, R., van den Bergh, D., & Marsman, M. (2024). Testing conditional
independence in psychometric networks: An analysis of three Bayesian
methods. *Multivariate Behavioral Research*, *59*, 913–933.
<https://doi.org/10.1080/00273171.2024.2345915>
