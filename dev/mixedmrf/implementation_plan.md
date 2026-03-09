# Mixed MRF — Implementation Plan

**Date:** 2026-02-25 (updated 2026-03-04; review amendments 2026-03-04;
Blume-Capel + PL consistency 2026-03-05; commit strategy 2026-03-05;
consistency review 2026-03-05; phases J-L added 2026-03-08)
**Branch:** `ggm_mixed` (PR #78)
**Goal:** Build a monolithic `MixedMRFModel` in C++ that supports both
conditional and marginal pseudo-likelihood, with and without edge selection.

**Reference code:** `MaartenMarsman/mixedGM` repository (R + Rcpp prototype)
**Theory:** `dev/plans/mixedMRF/A_Mixed_Graphical_Model_for_Continuous_and_Ordinal_Variables (1).pdf`

---

## Prototype Status (mixedGM)

The `mixedGM` package (`/Users/maartenmarsman/Documents/GitHub/mixedGM`) provides
a complete R + Rcpp prototype with:

### Implemented & Tested Components

| Component | Location | Status |
|-----------|----------|--------|
| Conditional OMRF likelihood | `src/log_likelihoods.cpp` | ✅ C++ (Rcpp) |
| Marginal OMRF likelihood | `src/log_likelihoods.cpp` | ✅ C++ (Rcpp) |
| Conditional GGM likelihood | `src/log_likelihoods.cpp` | ✅ C++ (Rcpp) |
| Θ computation | `src/log_likelihoods.cpp` | ✅ C++ (Rcpp) |
| Gibbs data generator | `src/mixed_gibbs.cpp` | ✅ C++ (Rcpp) |
| Cholesky update/downdate | `src/cholupdate.cpp` | ✅ C++ (Rcpp) |
| MH parameter updates | `R/cond_*_mh_update_functions.R` | R only |
| Edge selection | `R/cond_*_mh_update_functions.R` | R only |
| Main sampler loop | `R/mixed_sampler.R` | R only |
| Stan exact-likelihood model | `inst/stan/mixed_mrf_exact.stan` | ✅ Stan |
| Unit tests (likelihood) | `tests/testthat/test-likelihood-correctness.R` | ✅ |
| Parameter recovery tests | `tests/testthat/test-parameter-recovery.R` | ✅ |
| Edge selection tests | `tests/testthat/test-edge-selection.R` | ✅ |

### Bug Fixes Already Applied

The following issues identified in the plan reviews have been fixed in mixedGM:

| Issue | Status | Location |
|-------|--------|----------|
| Missing factor of 2 in Θ | ✅ Fixed | `rcpp_compute_Theta()`: `Kxx + 2.0 * Kxy * Sigma_yy * Kxy.t()` |
| `dnorm` without `log=TRUE` | ✅ Fixed | All edge selection functions use `log = TRUE` |
| Cache invalidation | ✅ Fixed | `update_Kyy_cache()` called after every Kyy change |
| Marginal PL Kyy acceptance | ✅ Fixed | Includes all p marginal OMRF terms |

### Validation Infrastructure

- **Stan exact model** (`inst/stan/mixed_mrf_exact.stan`): Enumerates all ordinal
  configurations for p ≤ 5 as gold-standard posterior. Serves as non-circular
  validation against the pseudolikelihood.
- **Simulation study** (`dev/simulation_study_plan.md`): 54-cell design covering
  p ∈ {5,10,15}, q ∈ {5,10,15}, n ∈ {250,500,1000}.

### Pseudolikelihood consistency (verified)

The discrete pseudolikelihood formulas in mixedGM and the bgms OMRF have
been compared and **match exactly** for ordinal variables.  The mixedGM code
was originally ported from bgms (`compute_denom_ordinal`); the same
category-threshold parameterization, rest-score definition, and log-sum-exp
stabilization are used.  The marginal OMRF correctly handles the
$\Theta_{ss}$ self-interaction term.  No reconciliation is needed.

The mixedGM prototype does **not** implement Blume-Capel variables.  The
Blume-Capel pseudolikelihood paths must therefore be adapted from the bgms
OMRF code directly rather than ported from mixedGM (see §1 and Phase B/C).

---

## Table of contents

1. [Model overview](#1-model-overview)
2. [Two pseudo-likelihood approaches](#2-two-pseudo-likelihood-approaches)
3. [Parameter groups and update schedule](#3-parameter-groups-and-update-schedule)
4. [File layout](#4-file-layout)
5. [Implementation phases](#5-implementation-phases)
6. [Phase A — Skeleton and data structures](#phase-a--skeleton-and-data-structures)
7. [Phase B — Conditional pseudo-likelihood](#phase-b--conditional-pseudo-likelihood)
8. [Phase B+ — Rank-1 Cholesky optimization](#phase-b--rank-1-cholesky-optimization)
9. [Phase C — Marginal pseudo-likelihood](#phase-c--marginal-pseudo-likelihood)
10. [Phase D — Edge selection](#phase-d--edge-selection)
11. [Phase E — R interface and integration](#phase-e--r-interface-and-integration)
12. [Phase F — Warmup, adaptation, and diagnostics](#phase-f--warmup-adaptation-and-diagnostics)
13. [Phase G — Simulation and prediction](#phase-g--simulation-and-prediction)
14. [Testing strategy](#testing-strategy)
15. [Reuse inventory](#reuse-inventory)
16. [Risk register](#risk-register)
17. [PR and commit strategy](#pr-and-commit-strategy)

---

## 1. Model overview

The mixed MRF models the joint distribution of $p$ discrete variables
$x$ ($x_s \in \{0, 1, \ldots, C_s\}$) and $q$ continuous variables $y$:

$$\log f(x, y) \propto \sum_s \mu_{x,s}(x_s) + x^\top K_{xx}\, x
  - \tfrac{1}{2}(y - \mu_y)^\top K_{yy}\,(y - \mu_y)
  + 2\, x^\top K_{xy}\, y$$

### Discrete variable types

Each discrete variable $x_s$ is either **ordinal** or **Blume-Capel**,
matching the existing `OMRFModel` design:

| Type | Main-effect parameterization | # free params | Reference level |
|------|------------------------------|:---:|---|
| Ordinal | $\mu_{x,s}(c) = \mu_{s,c}$, free thresholds ($c = 1,\ldots,C_s$; category 0 fixed at 0) | $C_s$ | Category 0 (fixed) |
| Blume-Capel | $\mu_{x,s}(c) = \alpha_s(c - \text{ref}_s) + \beta_s(c - \text{ref}_s)^2$ | 2 | User-specified $\text{ref}_s$ |

For Blume-Capel variables the observations are centered at
$\text{ref}_s$ in the constructor (i.e., `discrete_observations_.col(s) -= baseline_category_(s)`),
so all downstream code — rest-scores, sufficient statistics, likelihoods —
operates in a shifted coordinate system where the reference category
corresponds to zero, exactly as in `OMRFModel`.

The $x^\top K_{xx}\, x$ term uses the (centered) observation values
for both variable types; only the main-effect structure $\mu_{x,s}(\cdot)$
and the log-partition denominator differ.

Parameters:

| Symbol | Storage | Dimension | Role |
|--------|---------|-----------|------|
| $\mu_x$ | `mux_` | $p \times \max(C_s)$ | Ordinal thresholds **or** Blume-Capel ($\alpha, \beta$) |
| $K_{xx}$ | `Kxx_` | $p \times p$ symmetric, zero diag | Discrete pairwise interactions |
| $K_{yy}$ | `Kyy_` | $q \times q$ SPD | Continuous precision matrix |
| $K_{xy}$ | `Kxy_` | $p \times q$ | Cross-type interactions |
| $\mu_y$ | `muy_` | $q$-vector | Continuous means |
| $G$ | `edge_indicators_` | $(p+q) \times (p+q)$ | Edge inclusion indicators |

The factor 2 on $x^\top K_{xy}\, y$ reflects a symmetric parameterization:
the bilinear coupling between ordinal variable $x_s$ and continuous variable
$y_j$ contributes once as $(x_s, y_j)$ and once as $(y_j, x_s)$ in the
sufficient statistics sum over node pairs.  Absorbing both contributions into
a single $K_{xy}$ matrix introduces the factor 2.  This is a pure convention
(not a free parameter) and must be applied consistently in all likelihoods,
conditional means, and $\Theta$ computations.

### Conventions carried into C++

- $K_{yy}$ stores the **positive-definite precision** $\Sigma_{yy}^{-1}$. The
  LaTeX note writes $K_{yy} = -\tfrac{1}{2}\Sigma_{yy}^{-1}$; we absorb the
  $-\tfrac{1}{2}$ into the log-density and always work with SPD matrices.
- The joint density is written as $-\tfrac{1}{2}(y - \mu_y)^\top K_{yy}
  (y - \mu_y)$ so that $\mu_y$ is the literal continuous mean. This matches
  the R prototype and keeps the conditional mean expression short.
- With this convention, all marginal rest-scores that arise from integrating
  out $y$ must include $K_{yy}^{-1}$ explicitly. Whenever the LaTeX note shows
  $\Sigma_{yy}$, substitute $K_{yy}^{-1}$.

---

## 2. Two pseudo-likelihood approaches

Both approaches approximate the intractable joint $f(x, y)$ using
pseudo-likelihoods. They share the same GGM part but differ in how the
discrete pseudo-likelihood handles the coupling to $y$.

### 2.1 Conditional pseudo-likelihood

$$\text{PL}_{\text{cond}}(x, y) =
  \underbrace{f(y \mid x)}_{\text{conditional GGM}} \cdot
  \prod_{s=1}^{p} \underbrace{f(x_s \mid x_{-s}, y)}_{\text{conditional OMRF}}$$

**Conditional OMRF** — full conditional of $x_s$ given $x_{-s}$ **and** $y$:

$$r_s = x_{-s}^\top K_{xx,-s,s} + 2\,y^\top K_{yx,s}$$
$$\log f(x_s = c \mid x_{-s}, y) = \mu_{x,s,c} + c \cdot r_s
  - \log\!\Bigl(1 + \sum_{c'=1}^{C_s} \exp(\mu_{x,s,c'} + c' \cdot r_s)\Bigr)$$

The rest-score $r_s$ depends on $K_{xx}$ and $K_{xy}$ but **not** on $K_{yy}$.

The formula above shows the ordinal form.  For Blume-Capel variables the
main-effect term $\mu_{x,s,c}$ is replaced by
$\alpha_s(c - \text{ref}_s) + \beta_s(c - \text{ref}_s)^2$ and the
denominator uses `compute_denom_blume_capel()` (see §1 and Phase B.1).

**Conditional GGM** — $y \mid x$ is multivariate Gaussian:

$$y \mid x \sim N\bigl(\mu_y + 2\,x\,K_{xy}\,K_{yy}^{-1},\; K_{yy}^{-1}\bigr)$$
$$\log f(y \mid x) = \frac{n}{2}\log|K_{yy}|
  - \frac{1}{2}\sum_{v=1}^{n}(y_v - M_v)^\top K_{yy}\,(y_v - M_v)$$

where $M = \mathbf{1}\mu_y^\top + 2\,x\,K_{xy}\,K_{yy}^{-1}$ is the
$n \times q$ conditional mean matrix.

### 2.2 Marginal pseudo-likelihood

$$\text{PL}_{\text{marg}}(x, y) =
  \underbrace{f(y \mid x)}_{\text{conditional GGM}} \cdot
  \prod_{s=1}^{p} \underbrace{f(x_s \mid x_{-s})}_{\text{marginal OMRF}}$$

**Marginal OMRF** — full conditional of $x_s$ after integrating out $y$:

$$\Theta = K_{xx} + 2\,K_{xy}\,K_{yy}^{-1}\,K_{yx}$$
$$r_s = \bigl(x^\top \Theta_{\cdot,s}\bigr) - x_s\,\Theta_{ss} + 2\,(K_{xy})_{s\cdot}\,K_{yy}^{-1}\,\mu_y$$

The self-interaction $x_s\,\Theta_{ss}$ must be subtracted because the rest-score
conditions only on $x_{-s}$.  In practice:
- Compute $r^{\text{row}} = x^\top \Theta_{\cdot,s}$ (matrix-vector product over all observations)
- Subtract `x_dbl.col(s) * Theta_(s,s)` from the result
- Add the scalar bias `2.0 * arma::dot(Kxy_.row(s), Kyy_inv_ * muy_)` (same for every observation)

This matches the mixedGM implementation in `rcpp_log_pl_marginal_omrf` exactly.

Same categorical form as conditional PL, but the effective interaction
matrix $\Theta$ absorbs the continuous variables. This means:
- Changing $K_{yy}$ or $K_{xy}$ requires recomputing $\Theta$
- Changing $\mu_y$ changes all rest-scores through the
  $K_{yy}^{-1}\,\mu_y$ term

### 2.3 Which parameters affect which likelihoods

| Parameter | Conditional OMRF | Conditional GGM | Marginal OMRF |
|-----------|:---:|:---:|:---:|
| $\mu_x$ | ✓ | | ✓ |
| $K_{xx}$ | ✓ | | ✓ (via $\Theta$) |
| $K_{yy}$ | | ✓ | ✓ (via $\Theta$) |
| $K_{xy}$ | ✓ | ✓ | ✓ (via $\Theta$) |
| $\mu_y$ | | ✓ | ✓ |

---

## 3. Parameter groups and update schedule

Each MCMC iteration updates 5 parameter groups in sequence, matching
the R prototype. The table shows which log-likelihood components enter
the MH acceptance ratio for each group.

### 3.1 Conditional PL mode

| Step | Parameter | Components in acceptance ratio |
|------|-----------|-------------------------------|
| 1 | $\mu_{x,s}$ (one main effect) | `cond_omrf(s)` + prior |
| 2 | $\mu_{y,j}$ (one mean) | `cond_ggm()` + prior |
| 3 | $K_{xx,ij}$ (one pair) | `cond_omrf(i) + cond_omrf(j)` + prior |
| 4 | $K_{yy}$ (one element) | `cond_ggm()` + prior (Cholesky proposal) |
| 5 | $K_{xy,ij}$ (one element) | `cond_omrf(i) + cond_ggm()` + prior |

### 3.2 Marginal PL mode

| Step | Parameter | Components in acceptance ratio |
|------|-----------|-------------------------------|
| 1 | $\mu_{x,s}$ | `marg_omrf(s)` + prior |
| 2 | $\mu_{y,j}$ | `cond_ggm()` + $\sum_s$ `marg_omrf(s)` + prior |
| 3 | $K_{xx,ij}$ | `marg_omrf(i) + marg_omrf(j)` + prior |
| 4 | $K_{yy}$ | `cond_ggm()` + $\sum_s$ `marg_omrf(s)` + prior (Cholesky proposal) |
| 5 | $K_{xy,ij}$ | $\sum_s$ `marg_omrf(s)` + `cond_ggm()` + prior |

**Key difference:** In marginal mode, updating $\mu_y$ requires evaluating
ALL $p$ marginal OMRF terms because $\mu_y$ enters every rest-score. The
same is true for $K_{yy}$ and every $K_{xy,ij}$ because $\Theta$ changes
globally with those parameters. Marginal PL is therefore much more
expensive per iteration.

### 3.3 Edge selection (post-warmup)

Three independent edge-selection sweeps per iteration.  The acceptance
components differ by PL mode because $\Theta$ depends on $K_{yy}$:

| Edge type | Indicator | RJ proposal | Conditional PL | Marginal PL |
|-----------|-----------|-------------|----------------|-------------|
| Discrete-discrete | $G_{xx}$ | Toggle + spike-and-slab | `cond_omrf(i) + cond_omrf(j)` | `marg_omrf(i) + marg_omrf(j)` |
| Continuous-continuous | $G_{yy}$ | Cholesky-based toggle | `cond_ggm()` | `cond_ggm()` + $\sum_s$ `marg_omrf(s)` |
| Cross | $G_{xy}$ | Toggle + spike-and-slab | `cond_omrf(i) + cond_ggm()` | `cond_ggm()` + $\sum_s$ `marg_omrf(s)` |

---

## 4. File layout

```
src/
  models/
    mixed/
      mixed_mrf_model.h            # Class declaration (follows GGMModel pattern)
      mixed_mrf_model.cpp          # Constructor, clone, vectorization
      mixed_mrf_likelihoods.cpp    # Port from mixedGM/src/log_likelihoods.cpp
      mixed_mrf_metropolis.cpp     # Port MH updates from mixedGM R code
      mixed_mrf_edge_selection.cpp # Port edge selection from mixedGM R code
      mixed_mrf_cholesky.cpp       # Cholesky permute/R() (port from mixedGM R)
  sample_mixed.cpp                 # Rcpp interface (copy sample_ggm.cpp pattern)
  mrf_simulation.cpp               # Extend with mixed_gibbs_generate()
  mrf_prediction.cpp             # Add mixed MRF prediction
R/
  bgm.R                          # Extend bgm() to dispatch mixed data
  validate_data.R                # Add mixed data validation
  validate_model.R               # Add mixed model validation
tests/
  testthat/
    test-mixed-mrf-likelihood.R   # Likelihood correctness vs R prototype
    test-mixed-mrf-sampling.R     # Recovery tests
    test-mixed-mrf-edge-sel.R     # Edge selection tests
```

`configure` + `inst/generate_makevars_sources.R` already glob every `.cpp`
under `src/`, but Phase A must still run the script (or re-run `configure`)
so the new translation units show up in `sources.mk`/`Makevars`. Likewise,
add the new Rcpp export to `src/RcppExports.cpp`, `R/RcppExports.R`, and
`NAMESPACE` when `sample_mixed.cpp` lands.

---

## 5. Implementation phases

| Phase | What | Depends on | Deliverable |
|-------|------|------------|-------------|
| **A** | Skeleton: class, data structures, `BaseModel` overrides | — | Compiles, no sampling | ✅ |
| **B** | Conditional PL: all 5 MH updates, no edge selection | A | Recovery test passes (cond PL, estimation only) | ✅ |
| **B+** | Rank-1 Cholesky optimization for Kyy updates | B | Same correctness, $O(q^2)$ per Kyy move instead of $O(q^3)$ | ✅ |
| **C** | Marginal PL: $\Theta$ caching, marginal OMRF, $\mu_y$ full sweep | B+ | Recovery test passes (marg PL, estimation only) | ✅ |
| **D** | Edge selection: 3 RJ sweeps | B+ | Structure recovery test passes | ✅ |
| **E** | R interface: `bgm()` dispatch, output formatting | B+ | End-to-end `bgm(mixed_data)` works | ✅ |
| **F** | Warmup schedule, adaptation, diagnostics | E | Full warmup pipeline | ✅ |
| **G** | Simulation and prediction | E | `simulate.bgms` and `predict.bgms` for mixed | ✅ |
| **H** | Stochastic block model edge prior | G | SBM prior wired into mixed MRF | ✅ |
| **I** | Missing data imputation | G | `na_action = "impute"` for mixed data | ✅ |
| **J** | Performance profiling | G | Benchmark report, hotspot fixes | partial |
| **K** | Code deduplication audit | G | Shared helpers extracted, dead code removed | ✅ |
| **L** | Model output R code inspection | G | Output format verified, bugs fixed |
| **M** | Test suite cleanup and completion | G | 4 new test files, promoted dev tests, full T1–T30 coverage |
| **N** | Alternative proposals for GGM edge updates | G | New proposal mechanisms for GGM edges |

---

## Phase A — Skeleton and data structures ✅

### A.1 Create `mixed_mrf_model.h`

```cpp
class MixedMRFModel : public BaseModel {
public:
    // Construction
    MixedMRFModel(
        const arma::imat& discrete_observations,  // n × p discrete (0-based categories)
        const arma::mat& continuous_observations, // n × q continuous
        const arma::ivec& num_categories, // p-vector
        const arma::uvec& is_ordinal_variable, // 1 = ordinal, 0 = Blume-Capel
        const arma::ivec& baseline_category,   // reference category per variable
        bool edge_selection,
        const std::string& pseudolikelihood, // "conditional" or "marginal"
        int seed
    );

    MixedMRFModel(const MixedMRFModel& other);

    // BaseModel overrides (all 13 pure virtuals)
    void do_one_metropolis_step(int iteration = -1) override;
    void update_edge_indicators() override;
    size_t parameter_dimension() const override;
    arma::vec get_vectorized_parameters() const override;
    void set_vectorized_parameters(const arma::vec& params) override;
    arma::ivec get_vectorized_indicator_parameters() override;
    size_t full_parameter_dimension() const override;
    arma::vec get_full_vectorized_parameters() const override;
    void set_seed(int seed) override;
    std::unique_ptr<BaseModel> clone() const override;
    SafeRNG& get_rng() override;
    const arma::imat& get_edge_indicators() const override;
    arma::mat& get_inclusion_probability() override;
    int get_num_variables() const override;
    int get_num_pairwise() const override;
    void prepare_iteration() override;
    void set_edge_selection_active(bool active) override;
    void initialize_graph() override;
    void init_metropolis_adaptation(const WarmupSchedule& schedule) override;
    void tune_proposal_sd(int iteration, const WarmupSchedule& schedule) override;
    bool has_missing_data() const override;
    void impute_missing() override;

    // Capability queries
    bool has_edge_selection() const override;
    bool has_adaptive_metropolis() const override;

private:
    // --- Data ---
    arma::imat discrete_observations_;        // n × p discrete observations
                                              //   Blume-Capel columns centered at baseline_category_
                                              //   in the constructor (same convention as OMRFModel).
    arma::mat discrete_observations_dbl_;     // n × p  double version (post-centering)
    arma::mat continuous_observations_;       // n × q continuous observations
    int n_, p_, q_;
    arma::ivec num_categories_;       // p-vector
    int max_cats_;                    // max(num_categories)
    arma::uvec is_ordinal_variable_;  // 1 = ordinal, 0 = Blume-Capel
    arma::ivec baseline_category_;    // reference category per discrete variable

    // --- Sufficient statistics ---
    arma::imat counts_per_category_;  // (max_cats+1) × p  category counts (ordinal vars only)
    arma::imat blume_capel_stats_;    // 2 × p  linear and quadratic sums (BC vars only)

    // --- Parameters ---
    arma::mat mux_;                   // p × max_cats thresholds / Blume-Capel coefficients
                                     //   Ordinal: mux_(s, c) = threshold for category c+1;
                                     //     category 0 is the reference (fixed at 0).
                                     //   Blume-Capel: mux_(s, 0) = linear α_s,
                                     //     mux_(s, 1) = quadratic β_s.
    arma::vec muy_;                   // q-vector continuous means
    arma::mat Kxx_;                   // p × p discrete interactions (diagonal always zero;
                                     //   enforced by construction — not a free parameter)
    arma::mat Kyy_;                   // q × q SPD precision
    arma::mat Kxy_;                   // p × q cross interactions

    // --- Edge indicators ---
    // Single combined (p+q)×(p+q) matrix (Decision: Option A).
    //   Gxx block : rows [0,p),    cols [0,p)    — symmetric, zero diag
    //   Gyy block : rows [p,p+q),  cols [p,p+q) — symmetric, zero diag
    //   Gxy block : rows [0,p),    cols [p,p+q) — full p×q rectangle
    //              (lower-left mirror [p,p+q)×[0,p) unused; Gxy is not symmetric)
    //
    // Accessor helpers (prefer over raw index arithmetic throughout):
    //   int& gxx(int i, int j)  { return edge_indicators_(i, j); }
    //   int& gyy(int i, int j)  { return edge_indicators_(p_+i, p_+j); }
    //   int& gxy(int i, int j)  { return edge_indicators_(i, p_+j); }
    //
    // Serialization order for get_vectorized_indicator_parameters():
    //   1. upper-tri(Gxx) row-major — length p(p-1)/2
    //   2. upper-tri(Gyy) row-major — length q(q-1)/2
    //   3. full Gxy      row-major — length p*q
    //   Total length: p(p-1)/2 + q(q-1)/2 + p*q
    arma::imat edge_indicators_;      // (p+q) × (p+q)
    arma::mat inclusion_probability_; // (p+q) × (p+q)
    bool edge_selection_;
    bool edge_selection_active_;

    // --- Proposal SDs (Robbins-Monro) ---
    arma::mat prop_sd_Kxx_;           // p × p
    arma::mat prop_sd_Kyy_;           // q × q
    arma::mat prop_sd_Kxy_;           // p × q
    arma::mat prop_sd_mux_;           // p × max_cats
    arma::vec prop_sd_muy_;           // q-vector

    // --- Cached quantities ---
    arma::mat Kyy_inv_;              // q × q  inverse of Kyy (always maintained)
    arma::mat Kyy_chol_;             // q × q  upper Cholesky of Kyy
    double Kyy_log_det_;             // log|Kyy|
    arma::mat Theta_;                // p × p  Kxx + 2 * Kxy * Kyy_inv * Kyx
                                     //   (marginal PL only)
    arma::mat conditional_mean_;     // n × q  muy + 2 * x * Kxy * Kyy_inv

    // --- Configuration ---
    bool use_marginal_pl_;           // true = marginal, false = conditional
    SafeRNG rng_;

    // --- Edge update order ---
    arma::uvec edge_order_;          // shuffled pair indices
};
```

### A.2 Implement trivial overrides

Implement `parameter_dimension`, `get/set_vectorized_parameters`,
`get_full_vectorized_parameters`, `clone`, `set_seed`, `get_rng`,
`get_edge_indicators`, `get_inclusion_probability`, `get_num_variables`,
`get_num_pairwise`, `has_edge_selection`, `has_adaptive_metropolis`.

**Constructor must center Blume-Capel observations:**
```cpp
for (int s = 0; s < p_; ++s) {
    if (!is_ordinal_variable_(s)) {
        discrete_observations_.col(s) -= baseline_category_(s);
    }
}
discrete_observations_dbl_ = arma::conv_to<arma::mat>::from(discrete_observations_);
```
This mirrors `OMRFModel`'s constructor exactly and ensures that rest-scores,
sufficient statistics, and the conditional mean all use the same coordinate
system regardless of variable type.

**Sufficient statistics** (computed in the constructor, identical to `OMRFModel`):
- `counts_per_category_`: for each ordinal variable, count observations per
  category $c \in \{0, \ldots, C_s\}$.
- `blume_capel_stats_`: for each Blume-Capel variable, store
  $\sum_i x_{is}$ (row 0, linear) and $\sum_i x_{is}^2$ (row 1, quadratic),
  where $x$ is already centered.

**Parameter vectorization order (free parameters only):**
1. `mux_`: For each variable $s = 0,\ldots,p-1$:
   - If ordinal: copy columns $0 \ldots C_s-1$ ($C_s$ thresholds).
   - If Blume-Capel: copy columns 0 and 1 ($\alpha_s, \beta_s$; always 2 entries).
   Total count: $\sum_{s \in \text{ord}} C_s + 2 \cdot |\text{BC}|$.
   Variables produce runs of different length — **do not use a fixed stride**.
   ```
   idx = 0
   for s in 0..p-1:
       if is_ordinal_variable_[s]:
           for c in 0..num_categories_[s]-1:
               out[idx++] = mux_(s, c)
       else:  // Blume-Capel
           out[idx++] = mux_(s, 0)  // linear α
           out[idx++] = mux_(s, 1)  // quadratic β
   ```
2. `Kxx_`: strictly upper-triangular entries (symmetry supplies the lower
   half). Count = $p(p-1)/2$. Diagonal entries are always zero and excluded.
3. `muy_`: all $q$ means.
4. `Kyy_`: upper-triangle **including** the diagonal. Count = $q(q+1)/2$.
5. `Kxy_`: all $p\times q$ cross entries, row-major.

**`prepare_iteration()`** should be included in Phase A.2's list of trivial
overrides: it shuffles `edge_order_` (same pattern as `OMRFModel`) so that
RJ sweeps have no order bias.

`parameter_dimension()` returns the number of currently **active**
parameters under edge selection (mirrors `OMRFModel`).
`full_parameter_dimension()` always returns the total count above so that
sample buffers have a fixed width even while RJ toggles edges on/off.

`get_vectorized_indicator_parameters()` serializes in three contiguous
blocks (matching the layout documented in the class header comment):
1. Upper triangle of $G_{xx}$ (rows/cols $[0,p)$) — length $p(p-1)/2$
2. Upper triangle of $G_{yy}$ (rows/cols $[p,p+q)$) — length $q(q-1)/2$
3. Full $G_{xy}$ block (rows $[0,p)$, cols $[p,p+q)$) row-major — length $pq$

Total length: $p(p-1)/2 + q(q-1)/2 + pq$.
The diagonal of `edge_indicators_` is always zero and excluded throughout.

### A.3 Testing checkpoint

- Class compiles and links
- `parameter_dimension()` returns correct count
- Round-trip: `set_vectorized_parameters(get_vectorized_parameters())` is identity

---

## Phase B — Conditional pseudo-likelihood ✅

### B.1 Implement log-likelihood functions ✅

Implemented in `src/models/mixed/mixed_mrf_likelihoods.cpp`.
Declared as private methods in `mixed_mrf_model.h` with a `friend`
declaration for the test helper.

Tests: `tests/testthat/test-mixed-mrf-likelihoods.R` (14 assertions)
validated against pure-R reference implementations covering ordinal,
Blume-Capel, GGM, zero-params, nonzero-params, minimal (p=1,q=1),
and all-BC configurations.

In `mixed_mrf_likelihoods.cpp`:

#### `log_conditional_omrf(int s)` — per-variable discrete PL

Computes $\log f(x_s \mid x_{-s}, y)$ summed over all $n$ observations.
The function **branches on `is_ordinal_variable_(s)`** to select the
appropriate main-effect structure and denominator.

**Rest-score (same for both variable types):**
```
rest_score = discrete_observations_dbl_ * Kxx_.col(s) - discrete_observations_dbl_.col(s) * Kxx_(s,s)
             + 2 * continuous_observations_ * Kxy_.row(s).t()
           = n-vector
```
(`Kxx_(s,s) = 0` by construction, so the self-interaction subtraction is
a no-op for $K_{xx}$; it becomes relevant in the marginal case for $\Theta$.)

**Ordinal path:**
```
For each observation v:
  For c = 1..C_s:
    eta[v,c] = mux_(s, c-1) + c * rest_score[v]
  log_Z[v] = log(1 + sum_c exp(eta[v,c]))   (log-sum-exp stabilized)
  ll += mux_(s, x[v,s]-1) + x[v,s] * rest_score[v] - log_Z[v]
        (only if x[v,s] > 0 for the threshold part)
```
Use `compute_denom_ordinal()` from `variable_helpers.h`.

**Blume-Capel path:**
```
alpha = mux_(s, 0);  beta = mux_(s, 1);  ref = baseline_category_(s)
For each observation v:
  For c = 0..C_s:
    theta_c = alpha * (c - ref) + beta * (c - ref)^2
    eta[v,c] = theta_c + (c - ref) * rest_score[v]
  log_Z[v] = log_sum_exp over c
  ll += alpha * x[v,s] + beta * x[v,s]^2 + x[v,s] * rest_score[v] - log_Z[v]
        (x already centered at ref)
```
Use `compute_denom_blume_capel()` from `variable_helpers.h`.

**Reuse opportunity:** Both denominator functions already exist in
`src/utils/variable_helpers.{h,cpp}` and are battle-tested in the OMRF.
The mixed MRF calls them directly—no duplication needed.

#### `log_conditional_ggm()` — conditional GGM log-likelihood

Computes $\log f(y \mid x)$ using cached `Kyy_inv_`, `Kyy_log_det_`,
and `conditional_mean_`.

```
conditional_mean_ = 1*muy_^T + 2 * discrete_observations_ * Kxy_ * Kyy_inv_
D = continuous_observations_ - conditional_mean_
quad_sum = sum((D * Kyy_) .* D)       // trace of Kyy * D^T * D
ll = n/2 * (-q * log(2*pi) + Kyy_log_det_) - quad_sum / 2
```

**Reuse opportunity:** This is structurally identical to `GGMModel::log_density_impl`
but with a non-zero, observation-dependent conditional mean.

### B.2 Implement cache maintenance ✅

All three functions implemented in Phase A (`mixed_mrf_model.cpp`).
The proposed-state bookkeeping rules below are the design contract
consumed by B.3's MH updates.

#### `recompute_conditional_mean()`
Recompute `conditional_mean_ = 1*muy^T + 2 * discrete_observations_ * Kxy_ * Kyy_inv_`.
Called after any change to `muy_`, `Kxy_`, or `Kyy_`.

#### `recompute_Kyy_decomposition()`
Recompute `Kyy_chol_`, `Kyy_inv_`, `Kyy_log_det_` from `Kyy_`.
Called after any change to `Kyy_`.

#### `recompute_Theta()` (marginal PL only)
Recompute `Theta_ = Kxx_ + 2 * Kxy_ * Kyy_inv_ * Kxy_^T`.
Called after any change to `Kxx_`, `Kxy_`, or `Kyy_`.

**Cache update order — proposed-state bookkeeping:**

The guiding rule is: **all caches must be consistent with the proposed
parameter value before evaluating the proposed likelihood**. On rejection,
temporaries are discarded cheaply (Armadillo move semantics).

1. **Kyy proposal:** Build proposed `Kyy_chol_prop`, `Kyy_inv_prop`,
   `Kyy_log_det_prop`, and `conditional_mean_prop` (and `Theta_prop` if
   marginal PL) into temporaries first. Evaluate proposed likelihood with
   these temporaries. On acceptance, swap into cached members in the order:
   `Kyy_` → `Kyy_chol_`/`Kyy_inv_`/`Kyy_log_det_` → `conditional_mean_`
   → `Theta_` (if needed). On rejection, discard temporaries.
2. **Kxx proposal (marginal PL mode):** In marginal PL mode, `Theta_`
   changes when `Kxx_` changes. Build `Theta_prop` from the proposed Kxx
   into a temporary before evaluating `log_marginal_omrf`. On acceptance,
   swap `Kxx_` and `Theta_`. On rejection, discard `Theta_prop`. In
   conditional PL mode, no cache depends on `Kxx_`, so no temporary is
   needed.
3. **Kxy proposal:** Build `conditional_mean_prop` (both modes) and
   `Theta_prop` (marginal mode) into temporaries before evaluating the
   proposed likelihood. On acceptance, swap `Kxy_`, `conditional_mean_`,
   and (if marginal) `Theta_`. On rejection, discard temporaries.
4. **muy proposal (marginal PL mode):** `conditional_mean_` depends on
   `muy_` directly. Build `conditional_mean_prop` (and re-evaluate all
   rest-score bias terms `2 * dot(Kxy_.row(s), Kyy_inv_ * muy_prop)`) before
   evaluating the proposed likelihood. Swap on acceptance.
5. **Every accepted RJ edge toggle (Phase D)** must trigger the same
   cache-refresh logic as an MH acceptance for that parameter type.

### B.3 Implement MH updates (conditional PL) ✅

Each update function follows the pattern:
1. Save current parameter value
2. Propose new value from $N(\text{current}, \text{prop\_sd})$
3. Compute log acceptance ratio (log-likelihood change + log-prior change)
4. Accept/reject
5. Robbins-Monro update of proposal SD

#### `update_main_effect(int s, int c)` — one main-effect parameter

- **Ordinal** ($c \in \{0,\ldots,C_s-1\}$): updates threshold $\mu_{s,c}$.
  - Propose: $\mu'_{s,c} \sim N(\mu_{s,c}, \sigma)$
  - Acceptance: `log_conditional_omrf(s)` at proposed vs current + prior
  - Prior: Beta-type (same as `OMRFModel`: $\alpha\,x - (\alpha+\beta)\log(1+e^x)$)
- **Blume-Capel** ($c \in \{0, 1\}$: linear, quadratic):
  updates $\alpha_s$ or $\beta_s$.
  - Propose: same Normal random-walk
  - Acceptance: `log_conditional_omrf(s)` at proposed vs current + prior
  - Prior: same Beta-type prior on each coefficient

The function signature is the same for both types; `is_ordinal_variable_(s)`
selects which main-effect slot to perturb.  Matches `OMRFModel`'s
`do_one_metropolis_step` loop which iterates over categories for ordinal
variables and over `{0, 1}` for Blume-Capel variables.

#### `update_continuous_mean(int j)` — one mean $\mu_{y,j}$
- Propose: $\mu'_{y,j} \sim N(\mu_{y,j}, \sigma_{\mu_y,j})$
- Acceptance: `log_conditional_ggm()` at proposed vs current + prior ratio
- Must update `conditional_mean_` before evaluating proposed likelihood
- Prior: Normal(0, $\sigma^2$)

#### `update_Kxx(int i, int j)` — one discrete interaction
- Propose: $K'_{xx,ij} \sim N(K_{xx,ij}, \sigma_{K_{xx},ij})$
- Set both $K_{xx,ij}$ and $K_{xx,ji}$ (symmetric)
- Acceptance: `log_conditional_omrf(i) + log_conditional_omrf(j)` + prior
- Prior: Cauchy(0, scale)
- Only if $G_{xx,ij} = 1$

#### `update_Kyy_element(int i, int j)` — one precision element (Cholesky)
- Uses the same Cholesky-based proposal as `GGMModel`:
  1. Permute rows/cols so (i,j) maps to (q-1, q) position
  2. Compute Cholesky, extract constants
  3. Propose on the Cholesky scale: $\phi' \sim N(\phi, \sigma_{K_{yy},ij})$
  4. Rebuild $K'_{yy}$ maintaining positive definiteness
  5. Unpermute
- Acceptance: `log_conditional_ggm()` at proposed vs current + priors
- Priors: Cauchy(0, scale) on off-diagonal; Gamma(shape, rate) on diagonal
- Must update all Kyy-dependent caches after acceptance

Unlike `GGMModel`, no rank-2 determinant lemma shortcut is available: the
conditional mean depends on $K_{yy}^{-1}$, so each proposal evaluates the
full `log_conditional_ggm()` with freshly computed `Kyy_inv_` and
`conditional_mean_`.

**Status (Phase B):** The current implementation follows the mixedGM
permute-then-Cholesky approach ($O(q^3)$ per Kyy move). Phase B+ replaces
this with the GGM's rank-1 infrastructure to bring it down to $O(q^2)$.
See [Phase B+](#phase-b--rank-1-cholesky-optimization) for the detailed plan.

**Jacobian for diagonal Kyy proposals.** The diagonal element is proposed
on the log scale to guarantee positivity:
```
// Diagonal element i:
theta_curr = log(L(i, i))              // L is the Cholesky factor
theta_prop = rnorm(theta_curr, prop_sd)
L_prop(i, i) = exp(theta_prop)
// Include log-Jacobian in acceptance:
ln_alpha += theta_prop - theta_curr
```
Omitting this Jacobian biases the diagonal distribution. Off-diagonal
elements are proposed linearly and need no Jacobian.

**Cholesky update strategy.** `GGMModel` uses rank-1 Cholesky
update/downdate (`cholupdate.h`) applied to a single off-diagonal entry.
The mixed MRF follows the `mixedGM` approach instead: **permute** the
target row/column pair to the last two positions, perform the 2×2 block
update, then unpermute. These are distinct algorithms. Place the permute
helpers in `mixed_mrf_cholesky.cpp`, ported from
`mixedGM/R/continuous_variable_helper.R`. Do **not** attempt to reuse
`GGMModel`'s rank-1 routines for this purpose.

#### `update_Kxy(int i, int j)` — one cross interaction
- Propose: $K'_{xy,ij} \sim N(K_{xy,ij}, \sigma_{K_{xy},ij})$
- Acceptance: `log_conditional_omrf(i) + log_conditional_ggm()` + prior
- Must update `conditional_mean_` after modifying `Kxy_`
- Prior: Cauchy(0, scale)
- Only if $G_{xy,ij} = 1$

### B.4 Implement `do_one_metropolis_step(int iteration)` ✅

```cpp
void MixedMRFModel::do_one_metropolis_step(int iteration) {
    // Step 1: Update all main effects (ordinal thresholds or BC α/β)
    for (int s = 0; s < p_; ++s) {
        if (is_ordinal_variable_(s)) {
            for (int c = 0; c < num_categories_[s]; ++c)
                update_main_effect(s, c, iteration);
        } else {
            update_main_effect(s, 0, iteration);  // linear α
            update_main_effect(s, 1, iteration);  // quadratic β
        }
    }

    // Step 2: Update all continuous means
    for (int j = 0; j < q_; ++j)
        update_continuous_mean(j, iteration);

    // Step 3: Update Kxx (upper triangle, edge-gated)
    for (int i = 0; i < p_ - 1; ++i)
        for (int j = i + 1; j < p_; ++j)
            if (!edge_selection_active_ || gxx(i, j) == 1)
                update_Kxx(i, j, iteration);

    // Step 4: Update Kyy (off-diag + diagonal, edge-gated)
    for (int i = 0; i < q_ - 1; ++i)
        for (int j = i + 1; j < q_; ++j)
            if (!edge_selection_active_ || gyy(i, j) == 1)
                update_Kyy_offdiag(i, j, iteration);
    for (int i = 0; i < q_; ++i)
        update_Kyy_diag(i, iteration);  // diagonal always active

    // Step 5: Update Kxy (edge-gated)
    for (int i = 0; i < p_; ++i)
        for (int j = 0; j < q_; ++j)
            if (!edge_selection_active_ || gxy(i, j) == 1)
                update_Kxy(i, j, iteration);
}
```

### B.5 Testing checkpoint — conditional PL estimation

**Test 1 (T1, T2, T10): Likelihood agreement + cross-validation**
- Generate data from `mixed_gibbs_generate()` in R
- Compute `log_conditional_omrf(s)` and `log_conditional_ggm()` in both
  R (prototype) and C++ at the same parameter values
- Assert agreement to machine precision
- Cross-validate against `mixedGM::rcpp_log_pl_conditional_omrf` and
  `mixedGM::rcpp_log_pl_conditional_ggm` for a three-way check
  (bgms C++, R reference, mixedGM C++)

**Test 2 (T13): Recovery (estimation only, no edge selection)**
- Generate data from known parameters (p=3 ordinal, q=2 continuous)
- Run mixed sampler with conditional PL, no edge selection
- Check posterior means recover true parameters (correlation > 0.9)
- Use `dev/plans/mixedMRF/mixedGM/dev/conditional_vs_marginal_pl.R`
  as template

---

## Phase B+ — Rank-1 Cholesky optimization ✅

The Phase B implementation ports the mixedGM permute-then-Cholesky approach
for $K_{yy}$ updates: each proposal calls `arma::chol` ($O(q^3)$), evaluates
the full `log_conditional_ggm()` twice (proposed and current), and recomputes
all caches from scratch on acceptance. The GGM model avoids this entirely
with rank-1 Cholesky updates, log-likelihood ratios via the matrix
determinant lemma, and cached covariance. Phase B+ brings the same
infrastructure to the mixed model.

**Design decision: Cholesky over Woodbury.**  An alternative is to
maintain only $(K_{yy}, \Sigma_{yy}, \log|K_{yy}|)$ and update
$\Sigma_{yy}$ directly via the Woodbury identity, dropping the Cholesky
factor entirely (fewer members, simpler code).  The asymptotic cost is
identical — both eliminate the $O(q^3)$ bottleneck — but the Cholesky
factor provides a structural positive-definiteness guarantee between
periodic full recomputes, while a Woodbury-only covariance can silently
lose symmetry or PD from floating-point drift.  We therefore keep the
Cholesky factor.  The Woodbury-only variant remains a future option to
investigate if profiling shows the Cholesky maintenance constant matters.

### B+.0 Why this matters now

Every subsequent phase (C, D, E, F) exercises $K_{yy}$ updates. Marginal PL
(Phase C) adds $\sum_s$ `marg_omrf(s)` to every $K_{yy}$ acceptance ratio,
and edge selection (Phase D) adds Cholesky-based birth/death moves for
$G_{yy}$. Optimizing $K_{yy}$ after those phases are built means
retrofitting the rank-1 machinery into more code paths. Doing it now, while
only `update_Kyy_offdiag` and `update_Kyy_diag` exist, is cleaner.

### B+.1 Add GGM-style member variables

New members in `MixedMRFModel` (mirroring `GGMModel`):

```cpp
// --- Cholesky workspace (Phase B+) ---
arma::mat covariance_yy_;            // q × q  Kyy^{-1}, replaces Kyy_inv_
arma::mat inv_cholesky_yy_;          // q × q  R^{-1} where Kyy = R'R
std::array<double, 6> kyy_constants_; // reparametrization constants
arma::mat precision_yy_proposal_;    // q × q  scratch for proposed Kyy
arma::vec v1_, v2_;                  // 2-vectors for rank-2 decomposition
arma::vec vf1_, vf2_;               // q-vectors (full-size update vectors)
arma::vec u1_, u2_;                  // q-vectors (decomposed rank-1 vectors)
```

`Kyy_inv_` is renamed to `covariance_yy_` for clarity (it is the covariance
$\Sigma_{yy} = K_{yy}^{-1}$).  `Kyy_chol_` is retained as the upper
Cholesky factor.

### B+.2 Port `get_constants` and `constrained_diagonal` (no-permutation form)

Replace the current permutation-based `get_constants(const arma::mat& L,
int q)` and `constrained_diagonal(double omega, const CholConstants& C)`
with the GGM's direct-indexing approach that reads constants from
`covariance_yy_` without permuting:

```cpp
void MixedMRFModel::get_kyy_constants(int i, int j) {
    double logdet = get_log_det(Kyy_chol_);
    double log_adj_ii = logdet + std::log(std::abs(covariance_yy_(i, i)));
    double log_adj_ij = logdet + std::log(std::abs(covariance_yy_(i, j)));
    double log_adj_jj = logdet + std::log(std::abs(covariance_yy_(j, j)));

    double inv_sub = compute_inv_submatrix_i(covariance_yy_, i, j, j);
    double log_abs_inv_sub = log_adj_ii + std::log(std::abs(inv_sub));

    double Phi_q1q  = (2 * std::signbit(covariance_yy_(i, j)) - 1)
                    * std::exp(log_adj_ij - (log_adj_jj + log_abs_inv_sub) / 2);
    double Phi_q1q1 = std::exp((log_adj_jj - log_abs_inv_sub) / 2);

    kyy_constants_[0] = Phi_q1q;
    kyy_constants_[1] = Phi_q1q1;
    kyy_constants_[2] = Kyy_(i, j) - Phi_q1q * Phi_q1q1;   // c1
    kyy_constants_[3] = Phi_q1q1;                           // c2
    kyy_constants_[4] = Kyy_(j, j) - Phi_q1q * Phi_q1q;    // c3
    kyy_constants_[5] = kyy_constants_[4]
                      + kyy_constants_[2] * kyy_constants_[2]
                      / (kyy_constants_[3] * kyy_constants_[3]); // c4
}
```

This is `GGMModel::get_constants` operating on the $q \times q$ precision
and covariance pair instead of the GGM's $p \times p$ pair.  The helper
`compute_inv_submatrix_i` and `get_log_det` are reused directly from
`GGMModel` (they are pure algebra with no model-specific state).

### B+.3 Implement `log_ggm_ratio_edge(int i, int j)` — rank-2 log-LR

The GGM's `log_density_impl_edge` computes the log-likelihood **ratio**
for a 2-element change $\Delta\omega_{ij}, \Delta\omega_{jj}$ using the
matrix determinant lemma. The formula involves:

1. **Log-determinant ratio** from the $2 \times 2$ Schur complement.
2. **Trace ratio** from the sufficient statistic $S = X^\top X$.

For the mixed model's conditional GGM, the sufficient statistic is
$S^* = D^\top D$ where $D = Y - M$.  The key question is whether $S^*$
remains fixed when only $K_{yy}$ changes.

**It does not.**  The conditional mean $M = \mathbf{1}\mu_y^\top + 2\,X\,K_{xy}\,K_{yy}^{-1}$ depends on $K_{yy}^{-1}$.  However, the
Cholesky-based proposal changes $K_{yy}$ by a structured rank-2
perturbation (off-diagonal + constrained diagonal), and the
corresponding rank-2 change in $K_{yy}^{-1}$ (via the Woodbury identity)
can be computed in $O(q)$ from the update vectors.

The `log_ggm_ratio_edge` function therefore:

1. Computes the log-det ratio using the same $2 \times 2$ determinant
   lemma as the GGM: $O(q)$ from covariance entries.
2. Computes the proposed covariance $\Sigma'_{yy}$ via rank-2 Woodbury
   update on `covariance_yy_`: $O(q^2)$.
3. Computes the proposed conditional mean $M' = \mathbf{1}\mu_y^\top +
   2\,X\,K_{xy}\,\Sigma'_{yy}$: $O(npq)$ (matrix multiply, but this is
   the same cost as `recompute_conditional_mean` so no regression).
4. Computes the proposed quadratic form $\text{tr}(K'_{yy} (Y-M')^\top (Y-M'))$
   using the proposed precision and residuals.
5. Returns the full log-likelihood ratio including both log-det and
   quadratic-form changes.

The initial implementation recomputes the full conditional mean
$M'$, making step 3 the bottleneck at $O(npq)$.  Step B+.10 eliminates
this by exploiting the rank-2 structure of $\Sigma'_{yy} - \Sigma_{yy}$.

### B+.4 Implement `log_ggm_ratio_diag(int i)` — rank-1 log-LR

Same approach as B+.3 but for a diagonal-only change.  The rank-1 case
is simpler (1 × 1 Schur complement, single Woodbury vector).  Mirrors
`GGMModel::log_density_impl_diag`.

### B+.5 Implement rank-1 Cholesky update after acceptance

Port `GGMModel::cholesky_update_after_edge` and
`cholesky_update_after_diag` to operate on `Kyy_chol_`,
`inv_cholesky_yy_`, and `covariance_yy_`:

```cpp
void MixedMRFModel::cholesky_update_after_kyy_edge(
    double omega_ij_old, double omega_jj_old, int i, int j)
{
    // Decompose rank-2 Omega change into 2 rank-1 ops (same as GGM)
    v2_[0] = omega_ij_old - precision_yy_proposal_(i, j);
    v2_[1] = (omega_jj_old - precision_yy_proposal_(j, j)) / 2;
    vf1_[i] = v1_[0];  vf1_[j] = v1_[1];
    vf2_[i] = v2_[0];  vf2_[j] = v2_[1];
    u1_ = (vf1_ + vf2_) / std::sqrt(2.0);
    u2_ = (vf1_ - vf2_) / std::sqrt(2.0);

    cholesky_update(Kyy_chol_, u1_);
    cholesky_downdate(Kyy_chol_, u2_);

    arma::inv(inv_cholesky_yy_, arma::trimatu(Kyy_chol_));
    covariance_yy_ = inv_cholesky_yy_ * inv_cholesky_yy_.t();
    Kyy_log_det_ = get_log_det(Kyy_chol_);

    vf1_[i] = 0; vf1_[j] = 0;
    vf2_[i] = 0; vf2_[j] = 0;
}
```

The Cholesky update/downdate calls are $O(q^2)$ each.  Inverting the
triangular factor and forming `covariance_yy_` is also $O(q^2)$ (back-
substitution + symmetric product).  Total: $O(q^2)$ vs. the current
$O(q^3)$.

After the Cholesky update, recompute `conditional_mean_` (and `Theta_`
if marginal PL) as before — those are independent of the Cholesky
strategy.

### B+.6 Rewrite `update_Kyy_offdiag` and `update_Kyy_diag`

Replace the permute-then-Cholesky implementation with:

```
update_Kyy_offdiag(i, j):
    get_kyy_constants(i, j)
    propose phi' ~ N(Phi_q1q, sigma)
    map to omega'_ij, omega'_jj via constrained_diagonal
    fill precision_yy_proposal_
    ln_alpha = log_ggm_ratio_edge(i, j)  // rank-2 log-LR
    + prior ratio (Cauchy off-diag + Gamma diag)
    if accept:
        old_ij = Kyy_(i,j); old_jj = Kyy_(j,j)
        Kyy_ = precision_yy_proposal_ (only 3 entries changed)
        cholesky_update_after_kyy_edge(old_ij, old_jj, i, j)
        recompute_conditional_mean()
        if marginal: recompute_Theta()
```

The current save-evaluate-restore-evaluate-compare pattern is replaced
by a single log-ratio call.  No full `log_conditional_ggm()` evaluation
is needed.

### B+.7 Remove permutation helpers

Delete `make_perm_offdiag`, `make_perm_diag`, `permute_matrix`,
`get_constants(const arma::mat& L, int q)`, `constrained_diagonal(double,
CholConstants)`, and the `CholConstants` struct from
`mixed_mrf_metropolis.cpp`.  They are fully superseded by the GGM-style
functions.

### B+.8 Functions reused from GGM

| Function | Source | Reuse form |
|----------|--------|------------|
| `cholesky_update` / `cholesky_downdate` | `cholupdate.h` | Direct call |
| `get_log_det` | `GGMModel` | Copy as free function or method |
| `compute_inv_submatrix_i` | `GGMModel` | Copy as free function or method |
| `constrained_diagonal` (GGM form) | `GGMModel` | Adapted as `kyy_constrained_diagonal` |

`get_log_det` and `compute_inv_submatrix_i` are currently private methods
of `GGMModel`.  They are pure algebra with no model state.  Options:

- **(a)** Extract to a shared utility header (`src/math/cholesky_helpers.h`)
  and include from both models.
- **(b)** Duplicate as private methods on `MixedMRFModel`.

Option (a) is preferred (avoids code duplication and makes the shared
algebra explicit), but option (b) is acceptable if the extraction touches
too many files for this phase.

### B+.9 Testing checkpoint

**Correctness (non-negotiable):**
- Existing `test-mixed-mrf-likelihoods.R` must still pass unchanged
  (likelihood functions are not modified).
- Existing `test-mixed-mrf-sampling.R` must still pass (recovery results
  must be identical at the same seed).
- New test: for a known $K_{yy}$ proposal, verify that `log_ggm_ratio_edge`
  equals `log_conditional_ggm(proposed) - log_conditional_ggm(current)`
  to machine precision (T28).
- New test: verify rank-1 Cholesky update produces the same
  `Kyy_chol_`, `covariance_yy_`, `Kyy_log_det_` as a full recompute (T29).

**Performance (informational):**
- Benchmark `update_Kyy_offdiag` before and after on a $q = 20$ problem.
  Expected speedup: order of magnitude for the Cholesky portion.

### B+.10 Rank-2 quadratic-form shortcut

After B+.3–B+.6 are working, the remaining bottleneck in
`log_ggm_ratio_edge` is the $O(npq)$ conditional-mean recompute (step 3)
followed by the $O(nq^2)$ quadratic form (step 4).  Both can be reduced
by exploiting the rank-2 structure of the covariance change.

The Woodbury identity gives:

$$\Sigma'_{yy} - \Sigma_{yy} = \Delta\Sigma$$

where $\Delta\Sigma$ is a rank-2 matrix (it comes from the rank-2
perturbation of $K_{yy}$).  The conditional-mean change is:

$$\Delta M = M' - M = 2\,X\,K_{xy}\,\Delta\Sigma$$

Since $\Delta\Sigma$ is rank 2, $\Delta M$ is an $n \times q$ matrix of
rank at most 2.  Write $\Delta\Sigma = a\,b^\top + b\,a^\top$ (or its
actual Woodbury form) and define the two $n$-vectors
$g_1 = 2\,X\,K_{xy}\,a$ and $g_2 = 2\,X\,K_{xy}\,b$.  Then:

$$\Delta M = g_1\,b^\top + g_2\,a^\top$$

The proposed residual is $D' = D - \Delta M$.  The quadratic-form
difference becomes:

$$\text{tr}(K'_{yy}\,D'^\top D') - \text{tr}(K_{yy}\,D^\top D)
  = \text{tr}(\Delta K_{yy}\,D^\top D)
  - 2\,\text{tr}(K'_{yy}\,D^\top \Delta M)
  + \text{tr}(K'_{yy}\,\Delta M^\top \Delta M)$$

Each term can be evaluated without forming the full $n \times q$
matrices $M'$ or $D'$:

- $\text{tr}(\Delta K_{yy}\,D^\top D)$: $\Delta K_{yy}$ is nonzero at
  3 entries, and $D^\top D$ is the cached sufficient statistic $S^*$.
  Cost: $O(1)$.
- $\text{tr}(K'_{yy}\,D^\top \Delta M)$: Compute $D^\top g_1$ and
  $D^\top g_2$ ($O(nq)$ each), then contract with $K'_{yy} b$ and
  $K'_{yy} a$ ($O(q^2)$ each).  Total: $O(nq)$.
- $\text{tr}(K'_{yy}\,\Delta M^\top \Delta M)$: Inner products of
  $g_1, g_2$ ($O(n)$ each) contracted with a $2 \times 2$ block from
  $a^\top K'_{yy}\,b$ etc.  Total: $O(n + q^2)$.

Overall cost: $O(nq)$ instead of $O(npq + nq^2)$.  This also eliminates
the need to store the $n \times q$ proposed conditional mean.

**Pre-condition:** $D^\top D$ (the current residual cross-product) must be
cached.  Add a member `suf_stat_ggm_` ($q \times q$) updated after every
accepted move that changes $K_{yy}$, $K_{xy}$, or $\mu_y$.

**Testing:**
- T30: For a known proposal, verify the shortcut quadratic-form difference
  matches the brute-force computation to machine precision.
- Existing T28 still covers end-to-end log-ratio agreement.

---

## Phase C — Marginal pseudo-likelihood ✅

### C.1 Implement `log_marginal_omrf(int s)`

Same structure as `log_conditional_omrf(s)` but uses $\Theta$ instead
of $K_{xx}$ and adds the $\mu_y$ bias term.  The self-interaction must be
excluded explicitly:

```
Theta_ = Kxx_ + 2 * Kxy_ * Kyy_inv_ * Kxy_^T     (cached; see C.3)

// n-vector of rest-scores:
rest_score = discrete_observations_dbl_ * Theta_.col(s)             // include all n obs, all p vars
           - discrete_observations_dbl_.col(s) * Theta_(s, s)        // subtract self-interaction
           + 2.0 * arma::dot(Kxy_.row(s),
                             Kyy_inv_ * muy_)   // scalar bias, same for all obs
```

This matches `rcpp_log_pl_marginal_omrf` in `mixedGM/src/log_likelihoods.cpp`
exactly for ordinal variables.

**Blume-Capel variables in marginal PL:** The rest-score is the same as
above. The difference is in the numerator and denominator:
- The numerator adds $\alpha_s x_{vs} + \beta_s x_{vs}^2 +
  x_{vs} \cdot r_{vs}$ (using centered observations) plus
  $\Theta_{ss} \cdot x_{vs}^2$ for the self-interaction.
- The denominator uses `compute_denom_blume_capel()` with `col_offset(c)`
  set to $\theta_c + (c - \text{ref})^2 \cdot \Theta_{ss}$ instead of
  $\mu_{s,c} + c^2 \cdot \Theta_{ss}$.

This parallels the ordinal-vs-Blume-Capel branching already present in
`OMRFModel::log_pseudoposterior_with_state()`.

### C.2 Modify `update_continuous_mean(int j)` for marginal mode

In marginal PL mode, changing $\mu_y$ affects all $p$ marginal OMRF terms.
The acceptance ratio becomes:

```
sum_{s=1}^{p} [log_marginal_omrf(s, proposed) - log_marginal_omrf(s, current)]
+ [log_conditional_ggm(proposed) - log_conditional_ggm(current)]
+ log_prior_ratio
```

This is more expensive; patch the R prototype (factor-of-two fix +
`dnorm(log = TRUE)`) before generating fixtures so C++ and R target the
same expression.

### C.3 Cache invalidation for marginal PL

**Per-proposal proposed-Theta rule** (see Phase B.2): When proposing a
$K_{xx}$, $K_{xy}$, or $K_{yy}$ move in marginal PL mode, the proposed
marginal OMRF likelihood must be evaluated with a proposed $\Theta$
computed from the proposed parameter values.  This means each of those
update functions must build a local `Theta_prop` before calling
`log_marginal_omrf(s)` with the proposed value.  See Phase B.2 for the
exact temporary-variable protocol.

**Θ recompute granularity.** Recomputing $\Theta$ after every individual
element change is $O(p^2 q + pq^2)$ per move, which is expensive for
large $p, q$. The mitigation already noted in the risk register covers
the accepted-Theta path:  after **accepting** a $K_{xx}$/$K_{xy}$/$K_{yy}$
move, update the cached `Theta_` once from the new parameter state.  Do
not recompute `Theta_` more than once per accepted move.

For initial implementation, full recompute
is simpler and correct. Rank-1 shortcuts for single $K_{xy,ij}$ changes
are a future optimization.

### C.4 Configuration dispatch

`do_one_metropolis_step` dispatches between conditional and marginal
via `use_marginal_pl_`. The dispatch happens inside each update function,
not at the loop level — most updates have different acceptance ratios
in the two modes.

### C.5 Testing checkpoint — marginal PL estimation

**Test 3 (T3, T11): Marginal likelihood agreement**
- Same data as Test 1, compute `log_marginal_omrf(s)` in R and C++
- Assert agreement

**Test 4 (T14): Recovery (marginal PL)**
- Same setup as Test 2 but `pseudolikelihood = "marginal"`
- Check posterior means recover true parameters

**Test 5 (T16): Conditional vs marginal agreement**
- Run both modes on the same data, compare posterior means
- They should be similar (not identical — different approximations)

---

## Phase D — Edge selection ✅

### D.1 Discrete edge selection (`update_edge_indicator_Kxx`) ✅

For each pair $(i, j)$ with $i < j$:
1. Propose $G'_{xx,ij} = 1 - G_{xx,ij}$
2. If **birth** ($G_{xx,ij}=0 \to 1$): propose $K'_{xx,ij} \sim N(K_{xx,ij}, \sigma)$;
   `k_curr = Kxx_(i,j)`, `k_prop = rnorm(k_curr, sigma)` (k_curr = 0 on a true birth).
3. If **death** ($G_{xx,ij}=1 \to 0$): set $K'_{xx,ij} = 0$.
4. Compute log acceptance ratio:
   - Likelihood: `omrf(i) + omrf(j)` at proposed vs current parameters
   - Slab prior: on birth add `log_slab_prior(k_prop)`, on death subtract `log_slab_prior(k_curr)`
   - **Hastings ratio** (proposal asymmetry):
     - On birth: subtract `dnorm(k_prop, k_curr, sigma, log=true)` (cost of generating the proposed value)
     - On death: add `dnorm(k_curr, k_prop, sigma, log=true)` = `dnorm(k_curr, 0, sigma, log=true)`
       (cost of the reverse birth, which would propose k_curr from a Normal centred on 0)
   - Inclusion prior: $\log(\pi) - \log(1 - \pi)$ ratio (or reverse on death)
5. After an **accepted** move, refresh `Theta_` (marginal mode).

All log-density calls must use `log = true` to keep all terms on the log scale.
Follows `cond_omrf_update_association_indicator_pair` in the R prototype.
Transplant the Hastings terms verbatim from that R code; do not re-derive.

### D.2 Continuous edge selection (`update_edge_indicator_Kyy`) ✅

Cholesky-based birth/death as in `GGMModel`:
1. Permute, Cholesky, extract constants
2. If birth: propose $\phi' \sim N(0, \sigma)$, rebuild
3. If death: set off-diagonal to 0, rebuild
4. Acceptance: `cond_ggm()` + $\sum_s$ `marg_omrf(s)` (marginal mode only)
  + inclusion prior + slab/proposal density
5. On acceptance, run the same cache pipeline as an MH precision update
  (`Kyy` → decompositions → `conditional_mean_` → `Theta_`).

Follows `cond_ggm_update_precision_indicator_pair` in the R prototype.

**Reuse opportunity:** Directly reuse `GGMModel::update_edge_indicator_parameter_pair`
logic, adapted for the mixed model's conditional GGM likelihood.

### D.3 Cross edge selection (`update_edge_indicator_Kxy`) ✅

Cross-edges $G_{xy,ij}$ share the same Bernoulli inclusion prior $\pi$ as
$G_{xx}$ and $G_{yy}$ (decided).  A single $\pi$ keeps the SBM prior
uniform across all edge types; if different sparsity assumptions are needed
for cross-type edges in the future, a separate hyperparameter can be added
then.

For each pair $(i, j)$ where $i \in \{1..p\}$, $j \in \{1..q\}$:
1. Propose $G'_{xy,ij} = 1 - G_{xy,ij}$
2. If birth: propose $K'_{xy,ij} \sim N(K_{xy,ij}, \sigma)$
3. If death: set $K'_{xy,ij} = 0$
4. Acceptance: `omrf(i) + cond_ggm()` in conditional mode;
  `cond_ggm() + \sum_s marg_omrf(s)` in marginal mode; spike-and-slab priors
  plus log-dense Hastings terms.
5. On acceptance, update `conditional_mean_` (both modes) and `Theta_`
  (marginal mode) before the next likelihood evaluation.

Follows `cond_omrf_update_cross_association_indicator_pair` in the R prototype.

### D.4 Implement `update_edge_indicators()` ✅

```cpp
void MixedMRFModel::update_edge_indicators() {
    if (!edge_selection_active_) return;

    // Shuffle edge order
    // ...

    // Discrete-discrete edges
    for (int i = 0; i < p_ - 1; ++i)
        for (int j = i + 1; j < p_; ++j)
            update_edge_indicator_Kxx(i, j);

    // Continuous-continuous edges
    for (int i = 0; i < q_ - 1; ++i)
        for (int j = i + 1; j < q_; ++j)
            update_edge_indicator_Kyy(i, j);

    // Cross edges
    for (int i = 0; i < p_; ++i)
        for (int j = 0; j < q_; ++j)
            update_edge_indicator_Kxy(i, j);
}
```

### D.5 Blume-Capel and edge selection

Edge selection for $K_{xx}$ pairs is **type-agnostic**: the acceptance
ratio evaluates `log_conditional_omrf(i) + log_conditional_omrf(j)`
(or marginal equivalents), and those functions already branch on
`is_ordinal_variable_`.  No additional edge-selection logic is needed
for Blume-Capel variables—this mirrors the `OMRFModel` design where
`update_edge_indicator_parameter_pair()` calls
`compute_log_likelihood_ratio_for_variable()`, which dispatches to the
correct denominator internally.

Cross-edges $G_{xy}$ involving a Blume-Capel discrete variable also
require no special handling: the rest-score uses the centered
observations, and the likelihood function dispatches to the Blume-Capel
denominator as needed.

### D.6 Testing checkpoint — edge selection ✅

**Test 6 (T15): Structure recovery**
- Generate data from a sparse mixed graph (some edges zero)
- Run with edge selection, check posterior inclusion probabilities
  recover the true structure (true edges have high PIP, false edges low)

### D.7 Pairwise prior standardization

The pure OMRF scales the Cauchy slab prior per edge pair using
`pairwise_scaling_factors_`, computed in `compute_scaling_factors()`.
For two ordinal variables with $M_i$ and $M_j$ categories respectively,
the factor is $M_i \times M_j$; for Blume-Capel pairs it is the maximum
absolute product of the endpoint ranges; for mixed ordinal/Blume-Capel
pairs it is the corresponding cross product.  This makes the prior
comparable across variable pairs with different numbers of categories.

The mixed MRF currently uses a single `pairwise_scale_` for all Kxx
edges and does not apply per-pair scaling.  This must be fixed:

1. **Add `pairwise_scaling_factors_xx_`** (p × p matrix) to
   `MixedMRFModel`.  Kxy and Kyy edges are not scaled (continuous
   variables have no category-count dependence; cross-type edges
   couple a score in {0..M} with a continuous value, and the OMRF
   does not scale cross-type edges either).
2. **Compute the scaling factors** in the R interface layer (Phase E)
   using the existing `compute_scaling_factors()` applied to the
   discrete block, and pass them to `MixedMRFModel`.
3. **Apply** `pairwise_scale_ * pairwise_scaling_factors_xx_(i, j)` in
   `update_Kxx()`, `update_edge_indicator_Kxx()`, and anywhere
   else the Cauchy slab prior appears for Kxx edges.
4. **Test** that a model with variables of differing category counts
   produces sensible edge estimates (no systematic bias toward
   high-category pairs).

---

## Phase E — R interface and integration ✅

### E.0 Remove `src/test_mixed_mrf.cpp` ✅

Delete `src/test_mixed_mrf.cpp` once `sample_mixed.cpp` provides the
production Rcpp interface (E.1).  The test helpers
(`test_mixed_mrf_skeleton`, `test_mixed_mrf_likelihoods`,
`test_mixed_mrf_sampler`) were scaffolding for Phases A–B; their
functionality is superseded by the real entry point and should not ship
in the installed package.  Update `src/RcppExports.cpp` and
`R/RcppExports.R` accordingly (re-run `Rcpp::compileAttributes()`).

### E.1 Create `sample_mixed.cpp` ✅

Rcpp interface function `sample_mixed_mrf_cpp(...)`:
- Takes R data (integer matrix `x`, numeric matrix `y`)
- Creates `MixedMRFModel`
- Creates edge prior
- Calls `run_mcmc_sampler()`
- Returns results as `Rcpp::List`

Follow the pattern of `sample_ggm.cpp` and `sample_omrf.cpp`.

### E.2 Extend `bgm()` in R ✅

The user interface uses **Option A** (decided): a single data frame plus a
`variable_type` argument:

```r
bgm(data, variable_type = c("ordinal", "blume-capel", "continuous", ...))
```

- `data`: an $n \times (p+q)$ data frame or matrix with all variables.
- `variable_type`: character vector, length $p+q$; values `"ordinal"`,
  `"blume-capel"`, or `"continuous"`. Column order must match
  `variable_type`. Abbreviated forms (`"o"`, `"b"`, `"c"`) may be
  supported via `match.arg` for convenience.
- `baseline_category`: integer vector, length $p+q$ (entries for
  continuous and ordinal columns are ignored). Specifies the reference
  category for each Blume-Capel variable. Passed through to
  `sample_mixed_mrf_cpp()`.

`bgm()` splits `data` into the integer matrix `x` (ordinal + Blume-Capel
columns) and numeric matrix `y` (continuous columns), constructs
`is_ordinal_variable` and `baseline_category` vectors for the discrete
subset, then dispatches to `sample_mixed_mrf_cpp()`. The split indices
and variable-type metadata must be stored in the output object so that
`coef`, `predict`, and `simulate` methods can reconstruct the original
column order.

**Missing data** is an explicit non-goal for this PR.
`has_missing_data()` returns `false` and `impute_missing()` is an empty
override so the shared sampler pipeline compiles. Future imputation support
for mixed data should be a separate phase (Phase H).

### E.3 Extend `build_output.R` ✅

The output structure needs to accommodate:
- Separate interaction matrices: `Kxx`, `Kyy`, `Kxy` (or a combined one)
- Threshold samples: `mux` array
- Mean samples: `muy` matrix
- Edge indicators decomposed by type

### E.4 Testing checkpoint — end-to-end ✅

**Test 7 (T19): `bgm()` with mixed data**
- Call `bgm()` with `variable_type` containing both ordinal and continuous
- Check output structure matches expected format
- Verify S3 methods work (`print`, `summary`, `coef`, `predict`)

---

## Phase F — Warmup, adaptation, and diagnostics ✅

### F.1 Warmup schedule ✅

**Pre-condition check (before Phase F begins):** Verify that `WarmupSchedule`
and `ChainRunner` support a Metropolis-only model without NUTS mass-matrix
or step-size hooks.  If `WarmupSchedule::stage_2_windows()` triggers
NUTS-specific adaptation that `MixedMRFModel` cannot satisfy, the class
must either no-op those hooks or a simplified schedule must be added.
Check `GGMModel`'s implementation of `init_metropolis_adaptation` and
`tune_proposal_sd` as the reference pattern.

The mixed MRF is Metropolis-only (no NUTS/HMC), matching the GGM model.
Use the same warmup schedule:
- Stage 1: Initial fast adaptation (75 iterations)
- Stage 2: Doubling windows for covariance adaptation
- Stage 3a: Terminal fast adaptation
- Stage 3b: Proposal SD tuning with edge selection (if enabled)
- Stage 3c: Step-size re-adaptation with edge selection active

### F.2 Robbins-Monro adaptation ✅

Per-parameter proposal SD adaptation, matching the R prototype:
```
sigma_new = sigma_old + (acceptance - target) * weight
weight = (1/iter)^0.6
target = 0.44
```

Clamp to [0.001, 2.0].

### F.3 Init metropolis adaptation ✅

Override `init_metropolis_adaptation(const WarmupSchedule&)` to store
the schedule for use in `tune_proposal_sd()`.

Override `tune_proposal_sd(int iteration, const WarmupSchedule&)` — the
Robbins-Monro is already embedded in each update function.

---

## Phase G — Simulation and prediction ✅

### G.1 Mixed MRF simulation

Extend `mrf_simulation.cpp` to support mixed data generation via
block Gibbs sampling, matching `mixed_gibbs_generate()`:

1. Sample $x_s \mid x_{-s}, y$: categorical from log-sum-exp stabilized
   probabilities.
   - For ordinal variables use `compute_probs_ordinal()`.
   - For Blume-Capel variables use `compute_probs_blume_capel()`
     (both already available in `variable_helpers.h`).
2. Sample $y \mid x$: multivariate Normal with conditional mean
   $\mu_y + 2\,x\,K_{xy}\,K_{yy}^{-1}$ and covariance $K_{yy}^{-1}$

Note: the mixedGM Gibbs generator handles ordinal variables only.
The Blume-Capel Gibbs sampling step must be added using the existing
bgms probability helpers.

#### Implementation

| File | Component |
|------|-----------|
| `src/mrf_simulation.cpp` | `simulate_mixed_mrf()` (core Gibbs), `sample_mixed_mrf_gibbs()` (single-draw Rcpp export), `MixedSimulationWorker` (RcppParallel worker with SafeRNG), `run_mixed_simulation_parallel()` (parallel Rcpp export) |
| `R/simulate_predict.R` | `simulate_bgms_mixed()` dispatched from `simulate.bgms()`, `reconstruct_mixed_params_from_means()`, `split_mixed_raw_samples()`, `combine_mixed_result()` |

### G.2 Mixed MRF prediction

Extend `mrf_prediction.cpp` for posterior predictive checks on mixed data.

#### Implementation

| File | Component |
|------|-----------|
| `src/mrf_prediction.cpp` | `compute_conditional_mixed()` (Rcpp export): discrete probs via rest-score + ordinal/BC helpers, continuous conditional mean and SD from precision matrix |
| `R/simulate_predict.R` | `predict_bgms_mixed()` dispatched from `predict.bgms()`, `reconstruct_mixed_params_from_row()`, `format_mixed_predictions()`, `format_mixed_response()` |

### G.3 Cycle-test validation

Simulation-recovery cycle tests (fit → simulate → refit → correlate)
confirm that estimation and simulation code are consistent.

| PL mode | Pairwise r | Main r | Max |diff| (main) |
|---------|-----------|--------|---------------------|
| Conditional | 1.000 | 1.000 | 0.50 |
| Marginal | 0.994 | 0.994 | 1.02 |

The marginal PL main-effect offsets are a property of the parameterization,
not a code bug. The marginal PL estimates $\mu_x$ under $\Theta = K_{xx} -
K_{xy} K_{yy}^{-1} K_{xy}^\top$, but the Gibbs generator uses conditional
parameters ($K_{xx}$ directly). Both bgms and the mixedGM prototype show the
same pattern (mixedGM max |diff| = 2.4 on the same data). The pairwise
interactions are unaffected because $K_{xx}$, $K_{yy}$, and $K_{xy}$ are
shared across both parameterizations.

Test files:

| File | Content |
|------|---------|
| `tests/testthat/test-mixed-mrf-simulate-predict.R` | 77 tests: Gibbs sanity, parallel simulation, conditional prediction, S3 methods, edge cases (p=1, q=1) |
| `tests/testthat/test-simulate-predict-regression.R` | Mixed MRF fixtures added to existing regression test infrastructure |
| `dev/tests/test-simulation-recovery.R` | `run_simrec_test_mixed()` + Boredom cycle test |
| `dev/tests/plot_cycle_scatterplots.R` | Scatterplots for conditional and marginal PL cycle tests |
| `dev/tests/compare_bgms_mixedgm_cycle.R` | Cross-package comparison (bgms vs mixedGM) |

---

## Phase H — Stochastic block model edge prior ✅

The SBM edge prior is already implemented for GGM and OMRF models
(`src/priors/sbm_edge_prior.h`, `sbm_edge_prior_interface.cpp`).
The unified `BaseEdgePrior` / `ChainRunner` architecture means
the mixed MRF required **no additional code** — all five scope
items were already satisfied:

1. ✅ `bgm_spec()` and `build_spec_mixed_mrf()` forward
   `edge_prior` / SBM hyperparameters (6 params).
2. ✅ `run_sampler_mixed_mrf()` passes them to `sample_mixed_mrf()`.
3. ✅ `sample_mixed.cpp` creates the edge prior via
   `create_edge_prior()` and passes it to `run_mcmc_sampler()`.
4. ✅ `ChainRunner` calls `edge_prior.update()` after each
   edge-indicator sweep, modifying `inclusion_probability_` in
   place via `MixedMRFModel::get_inclusion_probability()`.
   Allocation samples are stored automatically.
5. ✅ `build_output_mixed_mrf()` extracts allocations,
   co-clustering matrix, posterior mean/mode allocations,
   and block-count summaries.

Verified with a smoke test: `bgm(mixed_data, edge_prior = "Stochastic-Block")`
produces correct allocation, co-clustering, and block-count outputs
across 4 chains.

---

## Phase I — Missing data imputation ✅

Missing data imputation wires the existing `BaseModel` virtual interface
(`has_missing_data()`, `impute_missing()`) into the mixed MRF model,
mirroring the GGM and OMRF implementations.

### Implementation

Four layers were connected:

1. **R spec** (`bgm_spec.R`): Removed the `stop()` guards that rejected
   `NA` entries for mixed data. Added listwise-deletion and imputation
   handling that processes discrete and continuous sub-matrices
   separately via `handle_impute()`. The spec stores
   `missing_index_discrete` and `missing_index_continuous` (0-based
   row/col index matrices).
2. **R sampler** (`run_sampler.R`): Passes `na_impute`,
   `missing_index_discrete`, and `missing_index_continuous` to the C++
   entry point.
3. **C++ entry** (`sample_mixed.cpp`): Accepts the three new
   parameters, calls `model.set_missing_data()`, and sets
   `config.na_impute`.
4. **C++ model** (`mixed_mrf_model.h/.cpp`):
   - `has_missing_data()` returns the stored flag (was hardcoded
     `false`).
   - `set_missing_data()` stores the index matrices and sets the flag.
   - `impute_missing()` performs Gibbs imputation in three phases:
     - Phase 1 — Discrete: samples each missing entry from its full
       conditional using rest scores and category probabilities
       (ordinal or Blume-Capel). Updates sufficient statistics
       (`counts_per_category_` or `blume_capel_stats_`).
     - Phase 2 — Refreshes `conditional_mean_` when both discrete
       and continuous entries are missing.
     - Phase 3 — Continuous: samples each missing entry from the
       conditional Normal $N(\mu^*, 1/K_{yy,jj})$. Invalidates the
       gradient cache.

`ChainRunner` already calls `impute_missing()` each iteration when
`config.na_impute && model.has_missing_data()`.

The `RcppExports` were regenerated, which also fixed a pre-existing
bug where the R wrapper for `sample_mixed_mrf` had only 17 parameters
while the C++ function expected 21.

### Files changed

| File | Change |
|------|--------|
| `R/bgm_spec.R` | Removed `stop()` guards; added listwise/impute handling for mixed data |
| `R/run_sampler.R` | Passes missing-data indices and `na_impute` flag to C++ |
| `R/RcppExports.R` | Regenerated (also fixed stale 17-param wrapper) |
| `src/sample_mixed.cpp` | Three new parameters; wires missing data to model |
| `src/models/mixed/mixed_mrf_model.h` | New data members, `set_missing_data()`, `impute_missing()` declaration |
| `src/models/mixed/mixed_mrf_model.cpp` | `impute_missing()` implementation, copy constructor updated |
| `src/RcppExports.cpp` | Regenerated |
---

## Phase J — Performance profiling (partial)

Profile the mixed MRF sampler to identify bottlenecks and verify
that the rank-1 Cholesky optimizations (Phase B+) deliver the
expected speedup over naive recomputation.

### Completed

- `MY_LOG`/`MY_EXP` macros applied across GGM and OMRF models
  (commit `cad3426`).
- Benchmark baseline saved (`dev/fixtures/benchmark_baseline.rds`);
  rounds 1–3 saved. Optimization rounds 1–8 stashed (`stash@{0}`).

### Remaining

- Formal scaling-curve plots and profiling report not yet produced.
- Marginal vs conditional PL cost comparison not yet documented.

### Scope

1. **Wall-clock benchmarks.** Time full `bgm()` runs at several
   network sizes (p = 4/8/16, q = 3/6/12) with fixed iteration
   counts. Record total time, per-iteration time, and breakdown
   by update type (Kxx, Kyy, Kxy, Kxy, main effects, edge
   selection when active). Compare MH vs NUTS.
2. **Cache-hit verification.** Under a profiler (or instrumented
   C++ timers), confirm that the Cholesky rank-1 update path
   dominates over any fallback full-recompute path. Quantify
   the fraction of iterations that hit the fast path.
3. **Scaling curves.** Plot per-iteration cost against q (fixing
   p) and against p (fixing q) to confirm empirical $O(q^2)$
   and $O(p)$ scaling respectively.
4. **Marginal vs conditional PL cost.** Conditional PL evaluates
   one OMRF term per discrete variable per move; marginal PL
   evaluates all $p$ terms for mean/precision moves. Benchmark
   both and document the crossover point.
5. **Comparison with OMRF and GGM.** Run the same-sized pure
   ordinal and pure continuous models and compare iteration cost
   to the mixed MRF.

### Deliverables

| File | Content |
|------|---------|
| `dev/tests/validation/group9_profiling.R` | Benchmark script |
| `dev/tests/validation/output/group9_profiling.pdf` | Scaling curves and timing tables |
| `dev/ggm_speed.md` or `dev/ggm_performance.md` | Updated with mixed MRF results |

---

## Phase K — Code deduplication audit ✅

Review the mixed MRF implementation and the rest of the package for
repeated code patterns that could be shared. The goal is to reduce
maintenance burden without over-abstracting.

### Completed

- **Part 1** (`676b49f`): Renamed all MixedMRF member variables to
  match GGM/OMRF naming conventions (`mux_` → `main_effects_discrete_`,
  `Kxx_` → `pairwise_effects_discrete_`, etc.). Replaced 6 `std::log`
  calls with `MY_LOG` in `mixed_mrf_metropolis.cpp`.
- **Part 2** (`d58848b`): Moved `cholupdate.h`/`.cpp` from
  `src/models/ggm/` to `src/math/`. Deduplicated R output-builder
  code in `build_output.R` and `run_sampler.R`.
- **Part 3** (`b944cd4`): Documented overflow-guard bitwise break,
  regenerated HMC fixtures, restored progress interval.
- Full audit in `dev/mixedmrf/deduplication_audit.md`.

### Known candidates

1. **Matrix assembly pattern in `build_output.R`.**
   `build_output_mixed_mrf()` contains four near-identical nested
   loops that fill an (p+q) × (p+q) matrix from offset-indexed
   flat vectors (discrete main effects, continuous main effects,
   pairwise matrix, indicator matrix). Extract a shared helper.
2. **SBM summarization computed twice.**
   `summarize_alloc_pairs()` is called identically for the summary
   table and the co-clustering matrix. Cache the result and reuse.
3. **Raw-samples list assembly.**
   `build_output_bgm()` and `build_output_mixed_mrf()` construct
   `results$raw_samples` with identical structure. Extract to a
   shared helper.
4. **Parameter index computation.**
   The 50-line block in `build_output_mixed_mrf()` that computes
   mux/Kxx/muy/Kxy/Kyy offsets is inline and hard to audit.
   Extract to a named function so it can be unit-tested.
5. **Cholesky header location.**
   `cholupdate.h` / `cholupdate.cpp` live in `src/models/ggm/`
   but are included by the mixed MRF model. Evaluate moving to
   `src/utils/` for discoverability.
6. **Metropolis step skeleton.**
   GGM, OMRF, and Mixed MRF each implement Robbins-Monro
   adaptation with target acceptance 0.234 / 0.44. The adaptation
   arithmetic is short and varies per model, so duplication is
   low, but verify no copy-paste drift has occurred.

### Scope

- Audit, not refactor. This phase produces a ranked list of
  deduplication opportunities with effort estimates. Actual
  refactoring happens in follow-up commits.
- For each candidate, determine: lines saved, risk of introducing
  bugs, and whether a shared helper already exists in the codebase.

### Deliverables

| File | Content |
|------|---------|
| `dev/mixedmrf/deduplication_audit.md` | Ranked list of candidates with code pointers |

---

## Phase L — Model output R code inspection

Systematic review of the R code that assembles, summarises, and
exposes the mixed MRF model output. The mixed MRF builder was
written alongside the C++ port and has accumulated inline logic
that may be fragile or inconsistent with the GGM/OMRF paths.

### Inspection points

1. **`build_output_mixed_mrf()` correctness.**
   - Verify that the flat-to-matrix index mapping for mux, Kxx,
     muy, Kxy, Kyy matches the C++ parameter ordering exactly.
   - Trace the NUTS diagnostic column naming (`__` suffix
     convention) and confirm it matches `build_output_bgm()`.
   - Check that posterior_mean_main for BC variables correctly
     stores alpha/beta in columns 1:2 with the rest NA.
   - Confirm that `posterior_mean_pairwise` symmetry and sign
     conventions match GGM/OMRF output.
2. **Posterior mean type consistency.**
   Mixed MRF returns `posterior_mean_main` as a list with
   `$discrete` and `$continuous` sub-elements, while GGM and
   OMRF return a matrix. Determine whether this inconsistency
   causes problems downstream (in `coef.bgms()`, print methods,
   user code) and propose a resolution.
3. **`coef.bgms()` with mixed MRF.**
   Verify that `coef()` on a mixed fit returns a structure that
   users can work with. Check edge cases: all-ordinal mixed MRF
   (q=0 should not occur but the code path should be guarded),
   all-BC discrete variables, single variable of each type.
4. **`simulate.bgms()` with mixed MRF.**
   Confirm that simulation from posterior means and from posterior
   samples both produce correctly ordered data frames (original
   column order, not internal discrete-first order).
5. **`predict.bgms()` with mixed MRF.**
   Confirm that conditional prediction handles mixed types and
   that the reordering from internal to user column order is
   correct.
6. **NUTS diagnostic wiring.**
   The `update_method` string in the mixed MRF path is
   `"hybrid-nuts"` while GGM uses `"nuts"`. Verify that
   downstream code (Rhat computation, ESS, trace extraction)
   handles both strings correctly.
7. **`summary.bgms()` and `print.bgms()`.**
   Confirm these methods produce sensible output for a mixed fit.

### Deliverables

| File | Content |
|------|---------|
| `dev/mixedmrf/output_inspection.md` | Findings, bugs found, and fix list |

---

## Phase M — Test suite cleanup and completion

The planned test structure (§4 file layout, commit strategy) envisioned
four dedicated test files in `tests/testthat/`:

```
test-mixed-mrf-skeleton.R     # T4, T5 (dimension, vectorization round-trip)
test-mixed-mrf-likelihood.R   # T1–T3, T7–T12 (likelihood unit tests)
test-mixed-mrf-sampling.R     # T13–T14, T16–T18 (recovery, reproducibility)
test-mixed-mrf-edge-sel.R     # T15, T26–T27 (edge selection, BC combos)
```

What actually shipped is different: the scaffolding test helpers were
removed in commit `04b9562`, the likelihood tests were removed with
them, and the remaining coverage landed in two files plus scattered
additions to existing files. This phase reconciles the plan with
reality, fills coverage gaps, and consolidates test organization.

### M.1 Current state audit

**`tests/testthat/` (run by `R CMD check`):**

| File | Mixed-MRF `test_that` blocks | Covers |
|------|:---:|---|
| `test-mixed-mrf-simulate-predict.R` | 15 | Gibbs generator, simulate, predict, p=1/q=1 edge cases |
| `test-bgm.R` (mixed section) | 8 | End-to-end `bgm()`, output dimensions, symmetry, marginal PL |
| `test-validate-variable-types.R` | 13 | Type validation (not model-specific) |

**`dev/tests/` (manual only, not run by `R CMD check`):**

| File | Blocks | Covers |
|------|:---:|---|
| `test-numerical-gradient.R` | 4 | Finite-diff vs analytical gradient (cond + marg, ord + BC) |
| `test-simulation-recovery.R` | 1 mixed | Mixed MRF sim-recovery cycle |
| `compare_bgms_mixedgm_cycle.R` | 0 | Script, no assertions — bgms vs mixedGM comparison |
| `validation/group1–7` | 0 | PDF-generating scripts, no programmatic assertions |

### M.2 Gap analysis against plan's test matrix

#### Covered (have a `test_that` block in `tests/testthat/`)

| Test | Description | Covered by |
|------|-------------|------------|
| T17 | Reproducibility (same seed) | `test-bgm.R` "bgm mixed MRF is reproducible" |
| T19 | End-to-end `bgm()` | `test-bgm.R` (8 blocks) |
| T22 | Degenerate p=1 | `test-mixed-mrf-simulate-predict.R` "single discrete variable" |
| T23 | Degenerate q=1 | `test-mixed-mrf-simulate-predict.R` "single continuous variable" |
| T25 | Gibbs generator sanity | `test-mixed-mrf-simulate-predict.R` (dimension, range, SD, seed) |
| R1 | GGM regression | Existing `test-bgm.R` GGM section |
| R2 | OMRF regression | Existing `test-bgm.R` OMRF sections |
| R3 | bgmCompare regression | `test-bgmCompare.R` |

#### Covered in `dev/tests/` only (not run by `R CMD check`)

| Test | Description | Covered by |
|------|-------------|------------|
| T1 | `log_conditional_omrf(s)` | Was in scaffolding tests (removed); partially in `validation/group4` |
| T2 | `log_conditional_ggm()` | Same |
| T3 | `log_marginal_omrf(s)` | Same |
| T13 | Conditional PL recovery | `validation/group1`, `test-simulation-recovery.R` |
| T14 | Marginal PL recovery | `validation/group1` |
| T15 | Edge selection recovery | `validation/group5` |
| T16 | Cond vs marginal agreement | `validation/group3` |

#### Covered in `dev/tests/` only (gradient-specific)

| Test | Description | Covered by |
|------|-------------|------------|
| NUTS-7.1 | Numerical gradient check | `test-numerical-gradient.R` (4 blocks) |
| NUTS-7.2 | Posterior equivalence MH vs NUTS | `validation/group2` |
| NUTS-7.3 | ESS comparison | `validation/group6` |
| NUTS-7.4 | Edge selection with hybrid | `validation/group5` |

#### Not covered anywhere

| Test | Description | Status |
|------|-------------|--------|
| T4 | `parameter_dimension()` | Was in scaffolding tests, removed |
| T5 | Vectorization round-trip | Was in scaffolding tests, removed |
| T6 | Cholesky permutation involution | Superseded by B+ (permutation removed) |
| T7 | Analytic Gaussian ($K_{xy}=0$) | Was in scaffolding tests, removed |
| T8 | Fixture replay (bit-for-bit) | Never implemented |
| T9 | Cache freshness after parameter tweak | Never implemented |
| T10 | BC conditional OMRF vs standalone OMRF | Never implemented |
| T11 | BC marginal OMRF vs standalone OMRF | Never implemented |
| T12 | BC observation centering consistency | Never implemented |
| T18 | Multi-chain (4 chains, R-hat < 1.1) | Not in `tests/testthat/` |
| T20 | Kyy PD invariant under stress | Never implemented |
| T21 | Cache consistency under RJ | Never implemented |
| T24 | Binary-only ordinal matches logistic | Never implemented |
| T26 | Mixed ordinal + BC recovery | Never implemented |
| T27 | BC-only discrete recovery + edge selection | Never implemented |
| T28 | `log_ggm_ratio_edge` agreement | Was in scaffolding tests, removed |
| T29 | Cholesky update fidelity | Was in scaffolding tests, removed |
| T30 | Rank-2 quadratic shortcut | Never implemented (B+.10 not implemented) |

#### Additional gaps not in original matrix

| Gap | Description |
|-----|-------------|
| G1 | No `bgm_spec()` test for `mixed_mrf` model type |
| G2 | No `test-spec-wiring.R` coverage for mixed pipeline |
| G3 | No `test-methods.R` coverage (`coef`, `print`, `summary` for mixed fixture) |
| G4 | No `predict.bgms` posterior-draw-path test (only posterior-mean tested) |
| G5 | No multi-chain test for mixed MRF (chain-combining logic) |
| G6 | No `na_action = "listwise"` test for mixed data |
| G7 | No Blume-Capel + continuous combo end-to-end test |
| G8 | `compare_bgms_mixedgm_cycle.R` is unstructured script (no `test_that`) |
| G9 | Validation scripts (group1–7) have no programmatic pass/fail |
| G10 | `test-numerical-gradient.R` and `test-simulation-recovery.R` use `<-` assignment (style violation) |

### M.3 Reorganization plan

#### Target file structure

Consolidate mixed-MRF tests into four files matching the original plan,
plus keep the simulate/predict file which has grown organically:

```
tests/testthat/
  test-mixed-mrf-skeleton.R         # T4, T5 (dimension, round-trip)
  test-mixed-mrf-likelihood.R       # T1–T3, T7, T10–T12 (likelihood unit)
  test-mixed-mrf-sampling.R         # T13, T14, T16–T18, T26, T27 (recovery)
  test-mixed-mrf-edge-selection.R   # T15, T20, T21 (edge selection, PD, cache)
  test-mixed-mrf-simulate-predict.R # (existing — keep as-is)
  test-mixed-mrf-gradient.R         # Numerical gradient (promote from dev/)
```

Mixed-MRF blocks currently in `test-bgm.R` stay there — they test the
`bgm()` interface, not the model internals. Add a mixed-MRF fixture to
`test-methods.R` and `test-bgm-spec.R`.

#### What to promote from `dev/tests/` to `tests/testthat/`

| Source | Target | Notes |
|--------|--------|-------|
| `test-numerical-gradient.R` (4 blocks) | `test-mixed-mrf-gradient.R` | Wrap in `skip_on_cran()`. Fix `<-` to `=`. |
| `test-simulation-recovery.R` (mixed block) | `test-mixed-mrf-sampling.R` | Merge into recovery section. Fix `<-` to `=`. |

The validation suite (`dev/tests/validation/`) stays in `dev/`. These
are long-running visual diagnostics, not `R CMD check` tests. No
changes needed except documenting their role.

`compare_bgms_mixedgm_cycle.R` stays in `dev/` as a manual diagnostic
script. Add a header comment clarifying it is not a testthat test.

#### What to add (new tests)

**Tier 1 — essential for correctness (add in this phase):**

| File | Tests | What |
|------|-------|------|
| `test-mixed-mrf-skeleton.R` | T4, T5 | Parameter count for known (p,q, mixed ord+BC); vectorization round-trip |
| `test-mixed-mrf-likelihood.R` | T1, T2, T3, T7 | Likelihood values vs R reference at known params; analytic Gaussian check |
| `test-mixed-mrf-likelihood.R` | T10, T11, T12 | BC likelihood vs standalone OMRF; BC observation centering |
| `test-mixed-mrf-sampling.R` | T18 | Multi-chain (2 chains, R-hat < 1.2; `skip_on_cran()`) |
| `test-mixed-mrf-sampling.R` | T26, T27 | Mixed ord+BC recovery; BC-only recovery (`skip_on_cran()`) |
| `test-mixed-mrf-edge-selection.R` | T15 | Structure recovery on sparse graph (`skip_on_cran()`) |
| `test-bgm-spec.R` | G1 | Spec builder produces `mixed_mrf` for `c("ordinal","continuous")` input |
| `test-methods.R` | G3 | `coef`, `print`, `summary` on mixed fixture |

**Tier 2 — defensive (add if time permits):**

| File | Tests | What |
|------|-------|------|
| `test-mixed-mrf-likelihood.R` | T9 | Cache freshness: tweak one param, verify cache matches recompute |
| `test-mixed-mrf-edge-selection.R` | T20 | Kyy PD invariant under 500 iters with wide proposals |
| `test-mixed-mrf-edge-selection.R` | T21 | Cache consistency after RJ birth/death |
| `test-mixed-mrf-sampling.R` | T24 | Binary-only ordinal vs logistic rest-scores |
| `test-mixed-mrf-sampling.R` | G4 | Predict with posterior draws (not just posterior-mean) |
| `test-mixed-mrf-sampling.R` | G5 | Multi-chain output combining |
| `test-bgm.R` | G7 | `bgm()` with Blume-Capel + continuous data |

**Dropped (superseded or not applicable):**

| Test | Reason |
|------|--------|
| T6 | Cholesky permutation removed in Phase B+; rank-1 approach has no permutation |
| T8 | Fixture replay requires stable binary fixtures; deferred to future maintenance |
| T28, T29 | Rank-1 log-ratio and Cholesky fidelity: internal C++ state not easily testable from R without scaffolding exports. Cover implicitly via T13/T14 recovery. |
| T30 | B+.10 quadratic shortcut not implemented |

### M.4 Style fixes

Apply to all promoted and new test files:

1. Replace `<-` with `=` for assignment (except inside `expect_*()`
   calls where the result must be captured:
   `expect_message(result <- foo(), "pat")`).
2. Remove `#'` roxygen-style comments; use plain `#` comments.
3. Use `if(` not `if (`.
4. Run `styler::style_file()` with `bgms_style` on each file.

### M.5 Deliverables

| # | Commit message | What ships |
|:-:|----------------|------------|
| 1 | `test: add mixed MRF skeleton and likelihood tests` | `test-mixed-mrf-skeleton.R` (T4, T5), `test-mixed-mrf-likelihood.R` (T1–T3, T7, T10–T12) |
| 2 | `test: add mixed MRF sampling recovery tests` | `test-mixed-mrf-sampling.R` (T13, T14, T16–T18, T26, T27) |
| 3 | `test: add mixed MRF edge selection tests` | `test-mixed-mrf-edge-selection.R` (T15) |
| 4 | `test: promote gradient tests to testthat` | `test-mixed-mrf-gradient.R` (from `dev/tests/`); style fixes |
| 5 | `test: add mixed MRF coverage to spec and methods tests` | `test-bgm-spec.R` (G1), `test-methods.R` (G3) |
| 6 | `test: defensive tests (tier 2)` | T9, T20, T21, T24, G4, G5, G7 |

### M.6 Reference implementation for likelihood tests

Likelihood tests (T1–T3, T10–T12) need R-level reference values.
The scaffolding test helpers that originally provided these were
Rcpp-exported functions in `test_mixed_mrf.cpp` (deleted in
commit `04b9562`). Two options:

**(a) Pure-R reference functions.** Write compact R functions that
compute `log_conditional_omrf(s)`, `log_conditional_ggm()`, and
`log_marginal_omrf(s)` using base-R matrix operations. These are
10–20 lines each and serve as the gold standard. Place in
`tests/testthat/helper-mixed-mrf-reference.R`.

**(b) Cross-validate against mixedGM.** Call
`mixedGM::rcpp_log_pl_conditional_omrf()` etc. from the test. Requires
`mixedGM` as a `Suggests` dependency or a `skip_if_not_installed()`
guard.

Option (a) is preferred: no external dependency, self-contained, and
exercises a different code path than the C++ implementation.

---

## Phase N — Alternative proposals for GGM edge updates

The current edge birth/death proposals for the continuous part of the
mixed MRF (and the standalone GGM) use an adaptive random-walk
Metropolis-Hastings scheme: the birth move draws
$\epsilon \sim \mathcal{N}(0,\sigma_{ij})$ on a reparameterized scale,
and the death move proposes setting the off-diagonal to zero. Proposal
standard deviations $\sigma_{ij}$ are tuned via Robbins-Monro during
warmup.

Van den Bergh, Clyde, Raftery, and Marsman (in preparation) show that
the MoMS (mixtures of mutually singular distributions) framework used in
bgms is equivalent to RJMCMC in terms of the Metropolis-Hastings
acceptance probability. However, the paper also demonstrates that
RJMCMC can achieve higher effective sample size (ESS) per iteration
in practice, because it allows more flexible proposal distributions
that better match the target. In principle, the MoMS framework admits
any proposal mechanism — the current random-walk proposal is a
convenience, not a requirement.

This phase investigates alternative proposal distributions for the
conditional GGM edge updates ($K_{yy}$) and cross-type edge updates
($K_{xy}$), with the goal of improving edge-indicator mixing and ESS.
The same ideas apply to the standalone `GGMModel`.

### N.1 Current proposal anatomy

Both `GGMModel` and `MixedMRFModel` use the same proposal pattern
for $K_{yy}$ edges:

**Birth** ($\gamma_{ij}: 0 \to 1$):
1. Draw $\epsilon \sim \mathcal{N}(0, \sigma_{ij})$.
2. Set $\omega_{ij}^* = C_3 \cdot \epsilon$ (reparameterized via
   Cholesky constants).
3. Constrain $\omega_{jj}^*$ to maintain positive definiteness.
4. Accept with MH ratio combining likelihood ratio, Cauchy slab
   prior, and $\mathcal{N}(0, \sigma_{ij})$ proposal density.

**Death** ($\gamma_{ij}: 1 \to 0$):
1. Set $\omega_{ij}^* = 0$.
2. Constrain $\omega_{jj}^*$.
3. Accept with the reverse MH ratio.

The proposal is symmetric on the $\epsilon$ scale. The
Robbins-Monro adaptation targets acceptance rate 0.44.

### N.2 Candidate alternatives

#### (a) Posterior-informed birth proposals

Instead of drawing from $\mathcal{N}(0, \sigma_{ij})$, use a
proposal centered on a point estimate of $\omega_{ij}$ conditional
on the current state:

$$\epsilon \sim \mathcal{N}\!\left(\hat{\omega}_{ij}^{\text{cond}},\,\sigma_{ij}\right)$$

where $\hat{\omega}_{ij}^{\text{cond}}$ is derived from the
conditional posterior mode. For the Gaussian likelihood with
Cauchy slab, the conditional mode given the sufficient statistics
$S = X^\top X$ and the rest of $\Omega$ can be computed in
closed form as a function of $S_{ij}$, $\Sigma_{ij}^{(-ij)}$,
and $n$. This turns the birth proposal into a near-optimal
independent proposal.

The death proposal remains deterministic (set to zero), so
the MH ratio adjusts by replacing the reverse-proposal
$q(\omega_{ij} \mid 0)$ with
$\mathcal{N}(\omega_{ij} \mid \hat{\omega}_{ij}^{\text{cond}}, \sigma_{ij})$.

#### (b) Exchange-type proposals

For the within-model parameter update (when $\gamma_{ij} = 1$),
replace the random-walk proposal with a Gibbs draw from the
conditional posterior of $\omega_{ij}$ given the rest of $\Omega$.
Under the Gaussian likelihood with a Cauchy slab, this conditional
posterior is not available in closed form, but a single Newton step
or a Laplace approximation to the conditional can serve as a
high-quality proposal.

#### (c) Block proposals for $K_{xy}$

The current $K_{xy}$ birth/death proposals update one element at a
time. Because each $K_{xy}$ change affects the conditional mean
$\mu_y + 2\,K_{xy}^\top x$, nearby $K_{xy}$ elements can be
strongly correlated in the posterior. A block proposal that updates
an entire row $K_{xy,i\cdot}$ simultaneously may improve mixing.
The NUTS sampler already handles $K_{xy}$ as part of the
gradient-based block, so this mainly benefits the Metropolis path.

### N.3 Evaluation plan

1. **Baseline ESS.** Run the current sampler on three benchmark
   problems (p=4/q=3, p=8/q=6, p=4/q=6) with 50,000 iterations
   after warmup. Record per-parameter ESS, per-edge-indicator ESS,
   and bulk/tail ESS.
2. **Implement candidate (a).** Add `update_edge_indicator_Kyy_informed()`
   as an alternative to `update_edge_indicator_Kyy()`. Selection via
   a compile-time or constructor flag (`proposal_type`).
3. **Comparison.** Run the same benchmarks with the informed proposal.
   Compare ESS, ESS/second, and acceptance rates. Focus on the
   edge-indicator ESS (the primary bottleneck).
4. **Evaluate (b) and (c)** if (a) shows improvement. Otherwise,
   document findings and close the phase.
5. **Apply to standalone GGM.** If the improved proposal works for
   the mixed MRF $K_{yy}$, port it to `GGMModel` which uses the
   same Cholesky reparameterization.

### N.4 Deliverables

| File | Content |
|------|---------|
| `dev/mixedmrf/proposal_study.md` | Design notes, derivation of informed proposal |
| `dev/tests/validation/group10_proposal_comparison.R` | Benchmark script |
| `dev/tests/validation/output/group10_proposal_comparison.pdf` | ESS comparison plots |
| `src/models/mixed/mixed_mrf_metropolis.cpp` | New proposal function(s) if adopted |
| `src/models/ggm/ggm_model.cpp` | Port to standalone GGM if adopted |

---

## Testing strategy



### Existing test infrastructure (from mixedGM)

The mixedGM prototype already provides comprehensive tests that can be adapted
for the bgms C++ implementation:

| mixedGM Test File | What it covers | Adaptation for bgms |
|-------------------|----------------|---------------------|
| `test-likelihood-correctness.R` | C++ likelihoods vs hand-computed values | Use as reference for bgms unit tests |
| `test-parameter-recovery.R` | Posterior mean recovery (p=2,q=1 and larger) | Port test scenarios |
| `test-edge-selection.R` | Edge birth/death moves, PIP calibration | Port test scenarios |
| `test-edge-cases.R` | p=1, q=1, binary-only ordinal | Port edge cases |
| `test-cholesky-update.R` | Cholesky update/downdate correctness | Already covered by bgms GGM tests |
| `test-data-generator.R` | Gibbs generator sanity checks | Port |
| `test-mcmc-diagnostics.R` | ESS, R-hat convergence checks | Port |
| `test-pl-comparison.R` | Conditional vs marginal PL agreement | Port |

### Test fixture availability

The mixedGM prototype implements all C++ likelihood functions as Rcpp exports.
This means test fixtures can be generated by running mixedGM R code and
comparing against the bgms C++ port:

```r
# Generate fixture in mixedGM
library(mixedGM)
set.seed(42)   # ALWAYS use a documented seed
ll_cond_omrf = rcpp_log_pl_conditional_omrf(x, y, Kxx, Kxy, mux, num_categories, i)
# Compare against bgms MixedMRFModel::log_conditional_omrf(i)
```

All fixture generation must use `set.seed(<documented value>)` before any
RNG-dependent step.  The seed and the mixedGM package version must be
recorded in a comment at the top of each fixture file.  The mixedGM
tests use `set.seed(42)` consistently; use the same seed for bgms fixtures
unless a specific test requires otherwise.

No need to patch the prototype — all math bugs have already been fixed.

Cross-reference the existing mixedGM test files for each scenario:
- `tests/testthat/test-likelihood-correctness.R` (likelihood unit tests)
- `tests/testthat/test-parameter-recovery.R` (recovery scenarios)
- `tests/testthat/test-edge-selection.R` (structure recovery scenarios)
- `tests/testthat/test-edge-cases.R` (p=1, q=1, binary ordinal)

### Unit tests (per-function correctness)

| Test | What | How |
|------|------|-----|
| T1 | `log_conditional_omrf(s)` | Compare C++ vs R prototype at known parameters |
| T2 | `log_conditional_ggm()` | Compare C++ vs R prototype at known parameters |
| T3 | `log_marginal_omrf(s)` | Compare C++ vs R prototype at known parameters |
| T4 | `parameter_dimension()` | Check count for known p, q, num_categories (ordinal + BC mix) |
| T5 | Vectorization round-trip | `set(get(params)) == params` (include BC variables) |
| T6 | Cholesky permutation | Verify permute is involution, PD maintained |
| T7 | Analytic Gaussian check | Set $K_{xy}=0$ so `log_conditional_ggm()` reduces to a standard MVN and compare against closed form |
| T8 | Fixture replay | Load deterministic R fixtures and ensure C++ reproduces saved log-likelihood components bit-for-bit |
| T9 | Cache freshness | After each parameter tweak, recompute `conditional_mean_` and (if needed) `Theta_` from scratch and compare with cached copies |
| T10 | BC conditional OMRF | Compare `log_conditional_omrf(s)` for a Blume-Capel variable against OMRF-only model at same parameters |
| T11 | BC marginal OMRF | Same as T10 but for `log_marginal_omrf(s)` |
| T12 | BC observation centering | Verify that centered observations + `compute_denom_blume_capel` produce the same result as the standalone `OMRFModel` |

### Rank-1 Cholesky optimization tests (Phase B+)

| Test | What | How |
|------|------|-----|
| T28 | Log-ratio agreement | For a known $K_{yy}$ off-diag proposal, verify `log_ggm_ratio_edge(i,j)` equals `log_conditional_ggm(proposed) - log_conditional_ggm(current)` to machine precision |
| T29 | Cholesky update fidelity | After rank-1 update, verify `Kyy_chol_`, `covariance_yy_`, `Kyy_log_det_` match a full `arma::chol` / `arma::inv` recompute |
| T30 | Rank-2 quadratic shortcut | For a known proposal, verify the B+.10 shortcut quadratic-form difference matches brute-force computation to machine precision |

### Integration tests (sampling correctness)

| Test | What | How |
|------|------|-----|
| T13 | Conditional PL recovery | Generate → estimate → check cor > 0.9 |
| T14 | Marginal PL recovery | Generate → estimate → check cor > 0.9 |
| T15 | Edge selection recovery | Sparse graph → check PIP > 0.5 for true, < 0.5 for false |
| T16 | Cond vs marginal agreement | Both modes → posterior means close |
| T17 | Reproducibility | Same seed → identical output |
| T18 | Multi-chain | 4 chains → R-hat < 1.1 |
| T19 | End-to-end `bgm()` | Call `bgm()` with mixed data; check output structure and S3 methods |

### Edge-case & stress tests

| Test | What | How |
|------|------|-----|
| T20 | Kyy PD invariant | Run 1k iterations with aggressive proposals; Cholesky of $K_{yy}$ must succeed after every accepted move |
| T21 | Cache consistency under RJ | After edge births/deaths, recompute caches and assert equality (debug build) |
| T22 | Degenerate $p=1$ | Single ordinal variable, ensure sampler runs and recovers $K_{xy}, K_{yy}$ |
| T23 | Degenerate $q=1$ | Scalar precision, ensure Cholesky permutation code handles the trivial case |
| T24 | Binary-only ordinal | All $C_s=1$, verify conditional PL matches logistic rest-scores |
| T25 | Gibbs generator sanity | Compare empirical moments from `mixed_gibbs_generate()` against theoretical targets before using it for fixtures |
| T26 | Mixed ordinal + BC | Graph with both ordinal and Blume-Capel discrete variables; verify recovery |
| T27 | BC-only discrete | All discrete variables are Blume-Capel; verify recovery and edge selection |

### Regression tests

| Test | What |
|------|------|
| R1 | Existing GGM tests still pass |
| R2 | Existing OMRF tests still pass |
| R3 | Existing bgmCompare tests still pass |

---

## Reuse inventory

### From mixedGM prototype (port to bgms)

| Component | mixedGM Location | Target in bgms | Action |
|-----------|-----------------|----------------|--------|
| Conditional OMRF likelihood | `src/log_likelihoods.cpp` | `mixed_mrf_likelihoods.cpp` | Port (remove Rcpp exports) |
| Marginal OMRF likelihood | `src/log_likelihoods.cpp` | `mixed_mrf_likelihoods.cpp` | Port |
| Conditional GGM likelihood | `src/log_likelihoods.cpp` | `mixed_mrf_likelihoods.cpp` | Port |
| Θ computation | `src/log_likelihoods.cpp` | `mixed_mrf_likelihoods.cpp` | Port |
| Log-sum-exp helper | `src/log_likelihoods.cpp` | `mixed_mrf_likelihoods.cpp` | Port |
| Gibbs data generator | `src/mixed_gibbs.cpp` | `mrf_simulation.cpp` | Extend |
| Cholesky helpers | `src/cholupdate.cpp` | Already in bgms | N/A |
| Cholesky permute/R() | `R/continuous_variable_helper.R` | `mixed_mrf_cholesky.cpp` | Port to C++ |
| MH Kxx updates | `R/cond_omrf_mh_update_functions.R` | `mixed_mrf_metropolis.cpp` | Port to C++ |
| MH Kyy updates | `R/cond_ggm_mh_update_functions.R` | `mixed_mrf_metropolis.cpp` | Port to C++ |
| MH Kxy updates | `R/cond_omrf_mh_update_functions.R` | `mixed_mrf_metropolis.cpp` | Port to C++ |
| Edge selection (Kxx) | `R/cond_omrf_mh_update_functions.R` | `mixed_mrf_edge_selection.cpp` | Port to C++ |
| Edge selection (Kyy) | `R/cond_ggm_mh_update_functions.R` | `mixed_mrf_edge_selection.cpp` | Port to C++ |
| Edge selection (Kxy) | `R/cond_omrf_mh_update_functions.R` | `mixed_mrf_edge_selection.cpp` | Port to C++ |
| Test fixtures | `tests/testthat/test-*.R` | `tests/testthat/test-mixed-*.R` | Adapt |

### From bgms existing infrastructure (direct reuse)

| Component | Source | Reuse type |
|-----------|--------|------------|
| Cholesky update/downdate | `src/models/ggm/cholupdate.h` | Direct include |
| `get_log_det` | `GGMModel` | Extract to shared utility or copy (Phase B+) |
| `compute_inv_submatrix_i` | `GGMModel` | Extract to shared utility or copy (Phase B+) |
| Rank-2 log-LR pattern | `GGMModel::log_density_impl_edge` | Adapt for observation-dependent mean (Phase B+) |
| Rank-2 Cholesky decomposition | `GGMModel::cholesky_update_after_edge` | Port to $q \times q$ Kyy workspace (Phase B+) |
| Log-sum-exp stabilization | `OMRFModel::compute_logZ_*` | Adapt pattern |
| Ordinal denominator | `src/utils/variable_helpers.{h,cpp}` `compute_denom_ordinal` | Direct call |
| Blume-Capel denominator | `src/utils/variable_helpers.{h,cpp}` `compute_denom_blume_capel` | Direct call |
| Ordinal probabilities | `src/utils/variable_helpers.{h,cpp}` `compute_probs_ordinal` | Direct call (Gibbs generator) |
| Blume-Capel probabilities | `src/utils/variable_helpers.{h,cpp}` `compute_probs_blume_capel` | Direct call (Gibbs generator) |
| Robbins-Monro adaptation | `OMRFModel::robbins_monro_*` | Direct reuse |
| Edge prior (SBM) | `src/priors/sbm_edge_prior.h` | Direct reuse |
| WarmupSchedule | `src/mcmc/execution/warmup_schedule.h` | Direct reuse |
| ChainRunner | `src/mcmc/execution/chain_runner.h` | Direct reuse |
| SafeRNG | `src/rng/rng_utils.h` | Direct reuse |
| Progress manager | `src/utils/progress_manager.h` | Direct reuse |
| R output builder | `R/build_output.R` | Extend |
| Validation functions | `R/validate_data.R` | Extend |
| Rcpp interface pattern | `src/sample_ggm.cpp` | Copy pattern |

---

## Risk register

| Risk | Impact | Mitigation |
|------|--------|------------|
| Kyy inversion per MH step is $O(q^3)$ | Slow for large $q$ | Phase B+: rank-1 Cholesky updates + GGM-style log-ratio shortcut bring per-move cost to $O(q^2)$ |
| Marginal PL $\mu_y$ update evaluates all $p$ OMRF terms | Slow for large $p$ | Cache rest-scores; provide conditional PL as faster default |
| $\Theta$ recomputation after every $K_{yy}$/$K_{xy}$ change | $O(p^2 q + pq^2)$ | Defer $\Theta$ recompute to once per accepted move, not per-element; proposed-state uses a temporary (Phase B.2, C.3) |
| Edge selection order effects | Bias | Shuffle edge order each iteration (already done in OMRF) |
| **Sticky edge indicators** | Poor PIP mixing | Shuffle cross-edge and xx-edge update order each iteration; monitor PIP stability and autocorrelation across chains |
| Numerical instability in log-sum-exp | NaN/Inf | Use stabilized version from OMRF (subtract max) |
| Factor 2 convention mismatch | Wrong posteriors | Document consistently; unit-test against R prototype |
| PD violation during Kyy proposals | Crash | Cholesky-based proposals guarantee PD by construction |
| Large parameter space mixing | Poor ESS | Per-parameter Robbins-Monro; future: block updates or HMC |
| Cholesky maintenance constant | Moderate overhead for small $q$ | Investigate Woodbury-only covariance update (same asymptotics, fewer members, but weaker PD guarantee; see Phase B+ design decision) |
| Blume-Capel not in prototype | No mixedGM reference for BC paths | Adapt from bgms `OMRFModel`; validated via OMRF-vs-mixed comparison tests (T10-T12) |

---

## PR and commit strategy

**Branch:** `ggm_mixed` (currently at the same commit as `main`).
**PR:** A single draft PR opened early for visibility, merged only after
all phases are complete and CI passes.

### Principles

1. **Every commit compiles and passes `R CMD check`.** No broken
   intermediate states on the branch. If a phase needs multiple
   sub-steps before the package can build, stage them locally and
   squash into one commit that restores a green state.
2. **One commit per logical unit of work.** A "logical unit" is the
   smallest change that is self-contained and reviewable on its own:
   a new file, a completed function group, a test suite, a wiring
   change. Avoid mixing unrelated concerns in one commit.
3. **Tests accompany the code they exercise.** When a commit adds
   a function, the corresponding unit test belongs in the same commit
   (or the immediately following one if fixtures must be generated
   first). Do not defer all tests to the end.
4. **Commit messages use the project's conventional-commit prefixes:**
   `feat:`, `fix:`, `test:`, `refactor:`, `docs:`, `dev:`.
   Keep the first line under 72 characters; add a body paragraph only
   when the "why" is non-obvious.

### Planned commits

The table below maps implementation phases to commits. Each row is one
`git commit`. Rows may be split further if a phase turns out large, or
merged if two steps are trivially small — the guiding rule is always
principle 1 (green build) and principle 2 (one logical concern).

| # | Prefix | Message | Phase | What ships |
|:-:|--------|---------|:-----:|------------|
| 1 | `feat` | `feat: add MixedMRFModel skeleton and data structures` | A | Header, constructor (incl. BC centering), trivial overrides, vectorization, sufficient statistics. Stubs for sampling functions (`do_one_metropolis_step` etc.) so the class compiles. No sampling yet. |
| 2 | `test` | `test: add mixed MRF skeleton tests` | A | `test-mixed-mrf-skeleton.R`: compile check, `parameter_dimension`, vectorization round-trip (ordinal + BC mix). |
| 3 | `feat` | `feat: add mixed MRF likelihood functions` | B.1 | `mixed_mrf_likelihoods.cpp`: `log_conditional_omrf`, `log_conditional_ggm`, cache helpers. Ordinal + Blume-Capel paths. |
| 4 | `test` | `test: add mixed MRF likelihood unit tests` | B.1 | `test-mixed-mrf-likelihood.R`: T1–T3 (compare C++ vs R fixtures), T7 (analytic Gaussian), T10–T12 (BC-specific). |
| 5 | `feat` | `feat: add conditional PL Metropolis updates` | B.3–B.4 | `mixed_mrf_metropolis.cpp`: 6 MH update functions (main effect, muy, Kxx, Kyy off-diag, Kyy diag, Kxy) + static Cholesky helpers + `do_one_metropolis_step` 5-step sweep. |
| 6 | `test` | `test: add conditional PL recovery test` | B.5 | `test-mixed-mrf-sampling.R`: T13 (conditional PL parameter recovery, estimation only). |
| 6a | `refactor` | `refactor: rank-1 Cholesky optimization for Kyy updates` | B+ | Replace permute-based $O(q^3)$ Kyy updates with GGM-style rank-1 infrastructure. Shared utility for `get_log_det` / `compute_inv_submatrix_i`. |
| 6b | `test` | `test: add rank-1 Cholesky correctness tests` | B+ | T28 (log-ratio agreement), T29 (Cholesky update fidelity). Verify existing tests still pass at same seed. |
| 6c | `refactor` | `refactor: rank-2 quadratic-form shortcut for Kyy log-ratio` | B+.10 | Exploit rank-2 structure of $\Delta\Sigma$ to reduce log-ratio cost from $O(npq + nq^2)$ to $O(nq)$. Cache `suf_stat_ggm_`. |
| 6d | `test` | `test: add quadratic shortcut correctness test` | B+.10 | T30 (shortcut vs brute-force agreement). |
| 7 | `feat` | `feat: add marginal PL mode and Theta caching` | C | `log_marginal_omrf`, Theta cache, marginal-mode branches in MH updates. |
| 8 | `test` | `test: add marginal PL recovery and comparison tests` | C.5 | T14 (marginal PL recovery), T16 (cond vs marginal agreement). |
| 9 | `feat` | `feat: add mixed MRF edge selection` | D | `mixed_mrf_edge_selection.cpp`: RJ sweeps for Kxx, Kyy, Kxy. `update_edge_indicators`, `initialize_graph`, `prepare_iteration`. |
| 10 | `test` | `test: add mixed MRF edge selection tests` | D | T15 (structure recovery), T26–T27 (BC edge cases). |
| 11 | `feat` | `feat: add Rcpp interface and bgm() mixed dispatch` | E | `sample_mixed.cpp`, `bgm.R` changes, `validate_data.R`/`validate_model.R` extensions, `build_output.R` mixed output. |
| 12 | `test` | `test: add end-to-end bgm() mixed model tests` | E | T19 (end-to-end), S3 method checks (`print`, `summary`, `coef`). |
| 13 | `feat` | `feat: add mixed MRF warmup and adaptation` | F | `init_metropolis_adaptation`, `tune_proposal_sd`, warmup schedule integration. |
| 14 | `feat` | `feat: add mixed MRF simulation and prediction` | G | Gibbs generator in `mrf_simulation.cpp`, `simulate.bgms` + `predict.bgms` extensions. |
| 15 | `test` | `test: add simulation and prediction tests` | G | T25 (Gibbs sanity), generator tests, prediction tests. |
| 16 | `test` | `test: verify GGM, OMRF, bgmCompare regression suite` | — | Run full existing test suite, fix any regressions. R1–R3. |
| 17 | `docs` | `docs: add mixed MRF documentation and pkgdown entries` | — | Roxygen for new exported functions, `_pkgdown.yml` entries, NEWS.md entry. |

### Conventions for this branch

- **Do not rebase onto main mid-flight.** Merge main into `ggm_mixed`
  only if a conflict blocks progress, and note it in the PR
  description.
- **Fixture files** (`.rds`, `.rda`) generated for tests go in
  `dev/fixtures/mixed/` and are committed alongside the test that
  consumes them.
- **Plan updates** during implementation get their own `docs:` commit
  (e.g., `docs: update mixed MRF plan — Phase B lessons`).
- **Squash-merge** is used when merging the PR into `main`. The full
  per-phase history remains on the branch for reference, but `main`
  gets a single clean merge commit.
