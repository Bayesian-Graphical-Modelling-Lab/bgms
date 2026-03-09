---
title: "bgms: Bayesian Analysis of Graphical Models for Ordinal, Continuous, and Mixed Data"
tags:
  - R
  - Bayesian inference
  - graphical models
  - Markov random fields
  - Gaussian graphical models
  - psychometrics
authors:
  - name: Maarten Marsman
    orcid: 0000-0001-5309-7502
    affiliation: "1"
  - name:
      given-names: Don
      non-dropping-particle: van den
      surname: Bergh
    orcid: 0000-0002-9838-7308
    affiliation: "1"
affiliations:
  - index: 1
    name: University of Amsterdam, The Netherlands
date: 9 March 2026
bibliography: paper.bib
---

# Summary

`bgms` is an R package for Bayesian analysis of undirected graphical models
with binary, ordinal, continuous, and mixed variables. The package supports
single-group estimation with `bgm()` and multi-group comparison with
`bgmCompare()`. Core use cases include estimating network structure,
quantifying evidence for conditional dependence or independence, and comparing
parameter equality across groups.

The package combines pseudolikelihood-based model fitting with Markov chain
Monte Carlo (MCMC) sampling to obtain posterior summaries and Bayes factors for
model comparison tasks [@MarsmanVandenBerghHaslbeck_2025;
@MarsmanWaldorpSekulovskiHaslbeck_2024]. `bgms` is distributed on CRAN and
includes vignettes, tests, and C++ backends for computationally intensive
routines.

# Statement of need

Network models are widely used in psychological and psychometric research to
study conditional dependence structures among observed variables
[@ContrerasEtAl_2019; @RobinaughEtAl_2020]. In these settings, analysts often
need to answer inferential questions that go beyond point estimation: whether a
specific edge is present, whether a missing edge is plausible, and whether
network parameters differ across populations.

Frequentist workflows for graphical models are useful for estimation and
selection, but they are less direct when researchers need graded evidence for
both null and alternative hypotheses in the same analysis. `bgms` addresses this
need by offering a Bayesian workflow where posterior inclusion probabilities and
Bayes factors can be used to evaluate edge presence/absence and group
differences/equivalence within one framework [@HuthEtAl_2023_intro;
@SekulovskiEtAl_2024].

The target audience is applied researchers and methodologists who work with
cross-sectional network models in psychology, psychometrics, and related
fields, and who need a practical Bayesian toolchain in R for estimation,
comparison, prediction, simulation, and diagnostics.

# State of the field

Graphical-model software for applied research is broad, but many workflows are
either specialized to a single data type, focused on point estimation rather
than evidence quantification, or centered on one modeling task at a time
(estimation without comparison, or comparison without integrated diagnostics).

`bgms` contributes by integrating several tasks in one package: ordinal MRF
estimation, Gaussian graphical model estimation, mixed-variable MRF estimation,
and multi-group comparison with a shared inferential language based on posterior
inclusion and Bayes factors [@MarsmanVandenBerghHaslbeck_2025;
@MarsmanWaldorpSekulovskiHaslbeck_2024]. This design reduces friction between
analysis stages that are often split across multiple tools.

The package was developed as a standalone implementation rather than as a thin
wrapper because its scope required coordinated support for multiple model
families, custom prior specifications, shared diagnostics, and a compiled
backend tied to the package API. The result is a coherent research software
artifact that can be installed, tested, and reused as a single unit.

# Software design

`bgms` exposes high-level R interfaces while delegating compute-heavy operations
to C++ implementations. This architecture separates user-facing modeling
workflows from low-level sampling and likelihood code, allowing the package to
maintain a consistent R API while keeping core routines efficient.

Key design choices include:

1. A unified front-end (`bgm()`, `bgmCompare()`) across model families.
2. Shared extraction, summary, simulation, and prediction interfaces for fitted
   objects.
3. Bayesian variable/effect selection via spike-and-slab style workflows for
   edges, group differences, and (optionally) clustering components.
4. Pseudolikelihood-based modeling to make estimation practical in settings
   where full likelihood approaches are expensive.

These choices matter for research applications because users can move from model
estimation to hypothesis evaluation and reporting without re-implementing
separate pipelines for each data type.

# Research impact statement

`bgms` has been used in methodological work on Bayesian ordinal network
analysis, Bayesian group-difference testing in graphical models, and related
applications [@MarsmanVandenBerghHaslbeck_2025;
@MarsmanWaldorpSekulovskiHaslbeck_2024; @SekulovskiEtAl_2024]. The package also
supports reproducible workflows through scripted estimation, diagnostics,
simulation, and prediction, and includes public documentation and vignettes.

The software is publicly developed on GitHub with issue tracking, pull requests,
continuous integration, and automated tests. A stable CRAN release and package
website are available, providing accessible installation and user guidance. This
combination of methodological grounding, public development, and practical
interfaces supports continued adoption and citation in network-analytic
research.

# AI usage disclosure

Generative AI tools were used to assist drafting parts of this manuscript and to
support code/documentation editing during package development. Human authors
reviewed, edited, and validated all AI-assisted outputs and made all substantive
design and scientific decisions.

<!--
If needed, replace with a more specific statement before submission, for example:
- tools/models used and versions
- where used (code, docs, paper)
- verification process used by authors
-->

# Acknowledgements

We acknowledge contributors to the `bgms` codebase and documentation, and thank
reviewers and users who provided feedback through issues and pull requests.

<!-- Add funding and grant information here if applicable. -->

# References
