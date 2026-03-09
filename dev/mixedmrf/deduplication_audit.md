# Phase K — Code Deduplication & Architectural Conformance Audit

**Date:** 2026-03-09
**Branch:** `ggm_mixed`
**Scope:** Audit only — no refactoring in this phase.

---

## Part 1: Architectural Conformance

The GGM model is the prototype for how model classes should be structured.
The OMRF follows GGM conventions in most areas. This section documents
where the MixedMRF deviates from the GGM/OMRF pattern and recommends
alignment.

### 1.1 Edge selection: embedded vs separated

| Model | Pattern | Location |
|-------|---------|----------|
| GGM | Embedded in `do_one_metropolis_step` | `ggm_model.cpp` L341–346 |
| OMRF | Separated: `update_edge_indicators()` called by ChainRunner | `omrf_model.cpp` L1296–1305 |
| MixedMRF | Separated: `update_edge_indicators()` with three sub-sweeps | `mixed_mrf_model.cpp` |

**Finding:** GGM embeds edge selection inside the Metropolis step; OMRF
and MixedMRF separate it. The OMRF/MixedMRF pattern is the better
design (matches the ChainRunner contract where `update_edge_indicators()`
is called as a distinct step). GGM is the outlier here.

**Decision:** MixedMRF correctly follows the OMRF pattern. No change.
GGM should align upward to this pattern in a future refactor.

### 1.2 Edge-order shuffling

| Model | Shuffle | Where |
|-------|---------|-------|
| GGM | No shuffle — deterministic i,j loop | `do_one_metropolis_step` |
| OMRF | Shuffles `shuffled_edge_order_` in `prepare_iteration()` | `omrf_model.cpp` |
| MixedMRF | Shuffles three vectors in `prepare_iteration()` | `mixed_mrf_model.cpp` |

**Finding:** GGM does not shuffle edge update order at all. OMRF and
MixedMRF do, which is better practice (reduces order bias in RJ sweeps).

**Decision:** OMRF/MixedMRF pattern is preferred. MixedMRF is correct.
The three separate shuffle vectors (xx, yy, xy) are necessary because
MixedMRF has three distinct edge-type blocks. GGM should add shuffling
in a future refactor.

### 1.3 Metropolis adaptation mechanism

| Model | Mechanism | Controller object? |
|-------|-----------|-------------------|
| GGM | Inline Robbins-Monro in each update function | No |
| OMRF | Delegated to `MetropolisAdaptationController` | Yes (2 instances: main, pairwise) |
| MixedMRF | Inline Robbins-Monro in each update function | No |

**Adaptation formula (GGM and MixedMRF):**
```cpp
if(iteration >= 1 && iteration < total_warmup_) {
    double rm_weight = std::pow(iteration, -0.75);
    prop_sd_(i, j) = update_proposal_sd_with_robbins_monro(
        prop_sd_(i, j), ln_alpha, rm_weight, 0.44);
}
```

**Finding:** GGM and MixedMRF use identical inline Robbins-Monro (target
0.44, exponent -0.75). OMRF uses a controller object that collects
acceptance probabilities into matrices and applies batch updates. The
OMRF approach is more structured but also more complex.

**Recommendation:** MixedMRF matches GGM exactly. No copy-paste drift
detected. The OMRF's controller pattern was designed for its specific
batch-collection workflow and is not obviously better for the mixed
model's five-block sweep. Status quo is fine.

### 1.4 Proposal SD storage

| Model | Storage |
|-------|---------|
| GGM | `proposal_sds_` — flat `arma::vec` of length dim_ |
| OMRF | `proposal_sd_main_` (p x max_cats), `proposal_sd_pairwise_` (p x p) |
| MixedMRF | `prop_sd_mux_`, `prop_sd_muy_`, `prop_sd_Kxx_`, `prop_sd_Kyy_`, `prop_sd_Kxy_` |

**Finding:** GGM flattens everything into one vector; OMRF uses two
matrices; MixedMRF uses five separate storage objects. The MixedMRF
approach is natural for a model with five parameter blocks. No issue.

**Recommendation:** No change needed. The naming abbreviation `prop_sd_`
vs OMRF's `proposal_sd_` is a minor style divergence; consider
standardizing to `proposal_sd_` for consistency with OMRF.

### 1.5 Naming conventions

The naming scheme uses a two-level hierarchy: role (`main_effects_`,
`pairwise_effects_`) qualified by variable type (`discrete_`,
`continuous_`, `cross_`). Both category thresholds and continuous
means are main effects — the variable type is the qualifier, not the
statistical role.

#### Agreed rename mapping

**Parameters:**

| Current | New |
|---------|-----|
| `mux_` | `main_effects_discrete_` |
| `muy_` | `main_effects_continuous_` |
| `Kxx_` | `pairwise_effects_discrete_` |
| `Kyy_` | `pairwise_effects_continuous_` |
| `Kxy_` | `pairwise_effects_cross_` |

**Proposal SDs:**

| Current | New |
|---------|-----|
| `prop_sd_mux_` | `proposal_sd_main_discrete_` |
| `prop_sd_muy_` | `proposal_sd_main_continuous_` |
| `prop_sd_Kxx_` | `proposal_sd_pairwise_discrete_` |
| `prop_sd_Kyy_` | `proposal_sd_pairwise_continuous_` |
| `prop_sd_Kxy_` | `proposal_sd_pairwise_cross_` |

**Derived quantities (GGM pattern for Cholesky names):**

| Current | New |
|---------|-----|
| `Kyy_chol_` | `cholesky_of_precision_` |
| `inv_cholesky_yy_` | `inv_cholesky_of_precision_` |
| `covariance_yy_` | `covariance_continuous_` |

**Cross-model consistency:**

| Concept | GGM | OMRF | MixedMRF (new) |
|---------|-----|------|----------------|
| Main effects | N/A | `main_effects_` | `main_effects_discrete_`, `main_effects_continuous_` |
| Pairwise effects | `precision_matrix_` | `pairwise_effects_` | `pairwise_effects_discrete_`, `pairwise_effects_continuous_`, `pairwise_effects_cross_` |
| Cholesky | `cholesky_of_precision_` | N/A | `cholesky_of_precision_` |
| Inv Cholesky | `inv_cholesky_of_precision_` | N/A | `inv_cholesky_of_precision_` |
| Covariance | `covariance_matrix_` | N/A | `covariance_continuous_` |
| Proposal SDs | `proposal_sds_` | `proposal_sd_*` | `proposal_sd_*` |
| Observations (int) | N/A | `observations_` | `discrete_observations_` (OK — mixed needs qualifier) |
| Observations (dbl) | `observations_` | `observations_double_` | `discrete_observations_dbl_` (OK) |

Grepping for `main_effects_` or `pairwise_effects_` finds all relevant
members across all three models. The Cholesky names parallel GGM
directly — anyone familiar with the GGM code recognizes the role.

**Decision (high priority):** Apply this rename. The change is mechanical
but touches many lines across 4 C++ files and the R output builder.

### 1.6 `MY_LOG` vs `std::log` in edge selection

| Model | Usage in edge-indicator inclusion-prior terms |
|-------|-----------------------------------------------|
| GGM | `MY_LOG(inclusion_probability_(i, j))` |
| OMRF | `MY_LOG(inclusion_probability_ij)` |
| MixedMRF | `std::log(inclusion_probability_(i, j))` |

**Finding:** GGM and OMRF use `MY_LOG`; MixedMRF uses `std::log`.
The `MY_LOG` macro was introduced for consistent fast-math behavior
across the codebase.

**Recommendation (fix):** Replace all 6 `std::log(inclusion_probability_`
occurrences in `mixed_mrf_metropolis.cpp` with `MY_LOG`. This is a
mechanical change (lines 621, 627, 702, 712, 796, 802).

### 1.7 File organization

| Model | Implementation files | Pattern |
|-------|---------------------|---------|
| GGM | 1: `ggm_model.cpp` (~470 lines) + `cholupdate.cpp` (~120 lines) | Monolithic |
| OMRF | 1: `omrf_model.cpp` (~1450 lines) | Monolithic |
| MixedMRF | 4: `mixed_mrf_model.cpp` (~850), `mixed_mrf_likelihoods.cpp` (~140), `mixed_mrf_gradient.cpp` (~500), `mixed_mrf_metropolis.cpp` (~800) | Modular split |

**Finding:** MixedMRF splits into four files while GGM and OMRF are
monolithic. The MixedMRF approach is arguably better for a 2300-line
model, but diverges from the established pattern.

**Recommendation:** No change. The split is pragmatic given the model's
complexity. If GGM or OMRF grow, they should adopt the same split
pattern rather than forcing MixedMRF back into a monolith.

### 1.8 Public API surface

**GGM-specific public methods:**
- `log_likelihood()` (2 overloads)
- `set_missing_data(const arma::imat&)`

**OMRF-specific public methods:**
- ~15 accessors/setters (`get_main_effects`, `set_pairwise_effects`, etc.)
- `set_step_size`, `get_step_size`, `set_inv_mass`, `get_inv_mass`
- `set_pairwise_scaling_factors`

**MixedMRF-specific public methods:**
- `do_kyy_metropolis_step(int)` — exposed for hybrid NUTS sampler
- `set_missing_data(const arma::imat&, const arma::imat&)` — two matrices
- `storage_dimension()`, `get_storage_vectorized_parameters()` — extra vectorization levels
- `get_active_inv_mass()` — override for NUTS mass matrix

**Finding:** MixedMRF has a relatively narrow public API compared to OMRF.
The `do_kyy_metropolis_step` public method is specific to the hybrid
NUTS+MH design and is called by `HybridNUTSSampler`. The OMRF's many
accessors exist because the NUTS sampler needs direct parameter
access. MixedMRF absorbs this into the vectorization interface instead.

**Recommendation:** The OMRF's many public accessors (`get_main_effects`,
`set_pairwise_effects`, etc.) create a wide coupling surface. The
MixedMRF's approach of using vectorization methods is cleaner. Consider
narrowing OMRF's public API in a future refactor to match.

### 1.9 `storage_dimension()` override

Only MixedMRF overrides `storage_dimension()`. BaseModel's default
returns `full_parameter_dimension()`. MixedMRF needs this because
Kyy is excluded from the NUTS gradient block but must still be stored
in the sample buffer.

**Finding:** This is a necessary architectural extension, not a
conformance issue. GGM and OMRF don't need it because they don't
have a hybrid sampler.

**Recommendation:** No change. Document in BaseModel that
`storage_dimension` exists for hybrid samplers.

### 1.10 Documentation quality

| Model | Class-level doc | Private method doc | Member field doc | Section dividers |
|-------|:---:|:---:|:---:|:---:|
| GGM | `/** */` | Minimal | `///` on most | Minimal |
| OMRF | `/** */` | Most methods | `///` on all | `//` section headers |
| MixedMRF | `/** */` extensive | All methods | `///` on all | `// ===` and `// ---` |

**Finding:** MixedMRF has the most thorough documentation. This exceeds
GGM/OMRF quality. Good.

**Recommendation:** No change for MixedMRF. Consider improving GGM's
private-method documentation to match OMRF/MixedMRF level.

---

## Part 2: Code Duplication Candidates

### 2.1 Matrix assembly in `build_output.R`

**Location:** `R/build_output.R`
- `build_output_mixed_mrf()` lines 476–569: Four nested-loop blocks
  filling (p+q)×(p+q) matrices from flat vectors (pairwise means,
  pairwise indicators, each with Kxx/Kyy/Kxy sub-blocks).
- `build_output_bgm()` lines 225–253: Analogous matrix assembly using
  `lower.tri()` + transpose (simpler because only one matrix block).

**Duplication:** The mixed MRF has 3 sub-blocks per matrix × 2 matrices
(pairwise + indicator) = 6 near-identical loop nests. Each is ~12 lines.
Total redundant code: ~60 lines.

**Extraction candidate:**
```r
fill_symmetric_matrix = function(flat_values, idx_pairs, target_matrix) {
  for(k in seq_along(flat_values)) {
    i = idx_pairs[k, 1]; j = idx_pairs[k, 2]
    target_matrix[i, j] = flat_values[k]
    target_matrix[j, i] = flat_values[k]
  }
  target_matrix
}
```

**Effort:** Low (30 min). **Risk:** Low. **Lines saved:** ~50.

### 2.2 SBM summarization computed twice

**Location:** `R/build_output.R`
- `build_output_bgm()`: `summarize_alloc_pairs()` called at lines ~176
  and ~243, producing `sbm_summary` and `co_occur_matrix` respectively.
  Same inputs both times.

**Finding confirmed:** The same function is called twice with identical
arguments. The first call extracts `$sbm_summary`, the second extracts
`$co_occur_matrix`.

**Extraction:** Call once, store result, extract both fields.

**Effort:** Trivial (5 min). **Risk:** None. **Lines saved:** ~8.

### 2.3 Raw-samples list assembly

**Location:** `R/build_output.R`
- `build_output_bgm()` lines 274–304
- `build_output_mixed_mrf()` lines 678–703

**Duplication:** Near-identical structure. Both build:
```r
list(main = ..., pairwise = ..., indicator = ..., allocations = ...,
     nchains = ..., niter = ..., parameter_names = ...)
```

**Difference:** `build_output_bgm` includes `allocations` naming;
`build_output_mixed_mrf` omits it. Otherwise identical.

**Extraction candidate:**
```r
build_raw_samples_list = function(raw, edge_selection, edge_prior,
                                   names_main, edge_names) { ... }
```

**Effort:** Low (20 min). **Risk:** Low. **Lines saved:** ~25.

### 2.4 Parameter index computation in `build_output_mixed_mrf()`

**Location:** `R/build_output.R` lines 337–429 (~90 lines)

**Finding:** This block computes mux/Kxx/muy/Kxy/Kyy offsets and
separates them into `main_idx` and `pairwise_idx`. It's inline and
hard to audit. Contains a nested loop (lines 361–375) to identify
Kyy diagonal vs off-diagonal positions.

**Extraction candidate:**
```r
compute_mixed_parameter_indices = function(num_thresholds, p, q) {
  # Returns: list(mux_range, kxx_range, muy_range, kxy_range,
  #               kyy_range, kyy_diag_abs, kyy_offdiag_abs,
  #               main_idx, pairwise_idx)
}
```

**Effort:** Medium (45 min). **Risk:** Low — can unit-test independently.
**Lines saved:** 0 (extracted, not eliminated), but testability improves.

### 2.5 `cholupdate.h` location

**Current location:** `src/models/ggm/cholupdate.h` and `cholupdate.cpp`
**Included by:** GGMModel and MixedMRFModel

**Finding:** This is pure linear algebra with no GGM-specific state.
Living under `src/models/ggm/` is misleading since MixedMRF depends on it.

**Recommendation:** Move to `src/math/cholupdate.h` (alongside
`cholesky_helpers.h` which is already in `src/math/`).

**Effort:** Low (15 min — move files, update includes, regenerate
`sources.mk`). **Risk:** Low — `#include` paths change, nothing else.

### 2.6 Robbins-Monro adaptation

**Finding:** GGM and MixedMRF use identical inline Robbins-Monro calls
with the same parameters (target=0.44, exponent=-0.75). Both call
`update_proposal_sd_with_robbins_monro()` which is a shared utility.
No copy-paste drift detected.

OMRF uses `MetropolisAdaptationController` which internally does the
same math but in a batch-collection pattern.

**Recommendation:** No extraction needed. The shared utility function
`update_proposal_sd_with_robbins_monro()` already prevents drift.

### 2.7 Beta-Bernoulli between-cluster parameter handling

**Location:** `R/run_sampler.R`
- Lines ~109–115 in `run_sampler_omrf()`
- Same pattern in `run_sampler_ggm()` and `run_sampler_mixed_mrf()`

```r
bb_alpha_between = if(is.null(p$beta_bernoulli_alpha_between)) -1.0
                   else p$beta_bernoulli_alpha_between
```

**Finding:** Identical 6-line pattern in all three sampler functions.

**Extraction candidate:** Helper function or inline in `run_sampler()`.

**Effort:** Trivial (10 min). **Risk:** None. **Lines saved:** ~12.

### 2.8 Variable-type normalization in simulate/predict

**Location:** `R/simulate_predict.R`
- Line ~1157 (simulate.bgms) and ~1627 (predict.bgms)

```r
if(length(variable_type) == 1) {
  variable_type = rep(variable_type, num_variables)
}
```

**Finding:** Same 3-line pattern duplicated in simulate and predict.

**Extraction:** Minor; could go in a shared input-normalization helper.

**Effort:** Trivial (5 min). **Risk:** None. **Lines saved:** ~3.

---

## Part 3: Ranked Summary

| # | Item | Type | Priority | Effort | Risk | Lines saved |
|:-:|------|------|:--------:|:------:|:----:|:-----------:|
| 1 | MY_LOG in edge selection (§1.6) | Bug-class fix | High | 5 min | None | 0 (correctness) |
| 2 | Naming conventions (§1.5) | Conformance | High | 2–3 hrs | Medium | 0 (rename) |
| 3 | cholupdate.h location (§2.5) | Organization | Medium | 15 min | Low | 0 (move) |
| 4 | Matrix assembly helper (§2.1) | Dedup | Medium | 30 min | Low | ~50 |
| 5 | SBM double-call (§2.2) | Dedup | Medium | 5 min | None | ~8 |
| 6 | Parameter index extraction (§2.4) | Testability | Medium | 45 min | Low | 0 (testable) |
| 7 | Raw-samples assembly (§2.3) | Dedup | Low | 20 min | Low | ~25 |
| 8 | BB between-cluster helper (§2.7) | Dedup | Low | 10 min | None | ~12 |
| 9 | Variable-type normalization (§2.8) | Dedup | Low | 5 min | None | ~3 |
| 10 | prop_sd_ → proposal_sd_ (§1.4) | Naming | Low | 30 min | Low | 0 (rename) |

### Items NOT recommended for change

- **File split (§1.7):** MixedMRF's 4-file split is better than monolithic.
- **Adaptation mechanism (§1.3):** Inline Robbins-Monro matches GGM. No drift.
- **Public API width (§1.8):** MixedMRF's narrow API is preferred.
- **storage_dimension (§1.9):** Necessary for hybrid sampler.
- **Edge selection separation (§1.1):** MixedMRF follows OMRF correctly.
- **Edge-order shuffling (§1.2):** MixedMRF follows OMRF correctly.

---

## Appendix: Structural Comparison Tables

### A. BaseModel Override Implementation

| Override | GGM | OMRF | MixedMRF |
|----------|-----|------|----------|
| `has_gradient()` | `false` | `true` | `true` |
| `has_adaptive_metropolis()` | `true` | `true` | `true` |
| `do_one_metropolis_step` | Off-diag → diag → edge sel | Pairwise → main | Main → muy → Kxx → Kyy → Kxy |
| `update_edge_indicators` | No-op (embedded) | Shuffled sweep | Three shuffled sweeps |
| `prepare_iteration` | No-op | Shuffle edges | Shuffle 3 edge vectors |
| `get_vectorized_parameters` | Upper-tri precision | Main + pairwise | mux + muy + Kxx + Kxy |
| `set_vectorized_parameters` | Not implemented | Full unpack + cache | Full unpack + cache |
| `init_metropolis_adaptation` | Store warmup count | Create 2 controllers | Store warmup count |
| `tune_proposal_sd` | No-op | Poll 2 controllers | No-op |
| `impute_missing` | Conditional Normal | Categorical probs | Categorical + Normal |
| `clone` | Copy constructor | Copy constructor | Copy constructor |

### B. Member Variable Inventory

| Category | GGM | OMRF | MixedMRF |
|----------|:---:|:----:|:--------:|
| Data members | 5 | 8 | 10 |
| Parameter members | 1 (precision) | 2 (main, pairwise) | 5 (mux, muy, Kxx, Kyy, Kxy) |
| Cache members | 3 (chol, inv_chol, cov) | 3 (residual, gradient, index) | 8 (chol, inv_chol, cov, logdet, Theta, cond_mean, gradient, index) |
| Prior members | 2 | 4 | 3 |
| Proposal SD members | 1 | 2 | 5 |
| Edge members | 3 | 3 | 3 |
| RNG + config | 2 | 4 | 3 |
| Workspace | 7 (rank-2 vectors) | 0 | 7 (rank-2 vectors, prefixed kyy_) |
| **Total** | ~23 | ~24 | ~44 |

### C. Cholesky Infrastructure Sharing

```
src/math/cholesky_helpers.h       ← Shared (get_log_det, compute_inv_submatrix_i)
src/models/ggm/cholupdate.h       ← Should be src/math/cholupdate.h
    ├── Used by GGMModel           (precision_matrix_ updates)
    └── Used by MixedMRFModel      (Kyy_ updates)

GGM reparameterization pattern (get_constants, constrained_diagonal):
    ├── GGMModel::get_constants()
    └── MixedMRFModel::get_kyy_constants()    ← Adapted copy (same math)
```
