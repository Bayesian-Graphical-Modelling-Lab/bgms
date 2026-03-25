# Future Tasks to Investigate

**Context:** Issues discovered during PR #78 refactoring work.
Not blocking the current work — collect here and address later.

## Status Snapshot (2026-03-10)

Quick triage after recent workflow and test refactors:

- Item 1 (Triple data storage in OMRFModel): **Open**
- Item 2 (Comprehensive correctness / accuracy testing): **Partially addressed**
  - GGM correctness tests are now in `tests/testthat/test-bgm.R` (active workflow).
  - References to `dev/tests/test-simulation-recovery.R` are stale and should be replaced.
- Item 3 (Direct MCMC algorithm tests): **Open**
- Item 4 (Function documentation standard): **Partially addressed**
  - Conventions exist in `.github/copilot-instructions.md`.
  - Proposed file `dev/documentation_standard.md` does not exist yet.
- Item 5 (GGM: sample the means): **Open**
- Item 6 (Detailed computational manuals): **Open**
- Item 7 (Stage 3b double-update of interactions): **Open**
- Item 8 (Model-level warmup_info()): **Open / Phase-5 candidate**
- Item 9 (Extract shared edge-selection logic from models): **Open / Phase-5 candidate**
- Item 10 (Placeholder): **N/A**
- Item 11 (Lazy MCMC diagnostics + FFT ESS/Rhat): **Open**
- Item 12 (Empirical validation of log adjacent-category odds ratios): **Open**
- Item 13 (On-demand Σ entries via triangular solves in GGM): **Open** (Sherman-Morrison rejected as unstable; use on-demand solves instead)
- Item 14 (Standardized cross-type partial associations in mixed MRF): **Open**
- Item 15 (Extractor functions for Bayes factors): **Open**
- Item 16 (Check Cholesky downdate failure sentinel): **Open**
- Item 17 (Dense / block-diagonal mass matrix for NUTS): **Open**
- Item 18 (Rao-Blackwellized ESS for edge indicators): **Open**
- Item 19 (Remove HMC): **Open**
- Item 20 (SBC with edge selection): **Open**

Notes:

- `dev/` content has been reorganized; some paths in this file refer to pre-reorg locations.
- This snapshot is informational and does not change task priority by itself.

---

## 1. Triple data storage in OMRFModel

**Date:** 2025-02-24
**Status:** Open

**Problem:** Each chain stores the observation data three times:

| Member | Type | Shape | Purpose |
|--------|------|-------|---------|
| `observations_` | `arma::imat` | n×p | Integer: category counting, BC score indexing, pseudolikelihood ratio, imputation |
| `observations_double_` | `arma::mat` | n×p | BLAS: residual updates, sufficient statistics, column references |
| `observations_double_t_` | `arma::mat` | p×n | Gradient: `observations_double_t_ * E` in `logp_and_gradient()` (4 sites) |

**Memory cost:** ~20·n·p bytes per chain (3 × n×p × 4 or 8 bytes).

**Elimination candidates:**

- **`observations_` (imat):** After the centering fix, all remaining integer-specific
  uses can work from `observations_double_` with casts. The values are always integral.
  Sites to update: category counting (~line 148), BC score indexing (~line 161),
  pseudolikelihood ratio (~lines 784-785), imputation (~lines 1378-1381).
  Saves n×p×4 bytes per chain.

- **`observations_double_t_` (mat):** Used only in 4 gradient calls
  (`observations_double_t_ * E`). Armadillo's lazy `.t()` should recognize the
  transpose and call BLAS `gemv` with the transpose flag — no temporary allocation.
  This is the NUTS/HMC gradient hot path, so benchmark before removing.
  Saves n×p×8 bytes per chain.

**Action:** Drop `observations_`, benchmark `.t()` vs pre-stored transpose.

---

## 2. Comprehensive correctness / accuracy testing

**Date:** 2025-02-24
**Status:** Partially addressed

**Problem:** The testthat suite currently tests structure, ordering, reproducibility,
input validation, and S3 methods, but posterior-accuracy coverage remains uneven
across model families.

**Partially addressed by Part D (GGM cleanup):**
Phase 3 Part D added several accuracy/correctness tests to `tests/testthat/test-bgm.R`:
- Multi-chain convergence (R̂ < 1.1)
- MLE convergence (large n, cor > 0.95 with MLE)
- Edge detection power (PIP discrimination between true/false edges)
- Conditional regression check (implied β vs OLS β)
These use `skip_on_cran()` for heavy tests. The GGM now has reasonable testthat-level
correctness coverage.

**What's still needed:**
- Design a coherent accuracy/correctness testing strategy covering both GGM and OMRF
- Decide which tests belong in testthat (with `skip_on_cran`) vs dev-only scripts
- Consider: simulation-recovery, known-parameter recovery, conditional regression
  checks, SBC, multi-chain convergence diagnostics
- Fit -> simulate -> refit cycle tests: fit on observed data, simulate from fitted
  model, refit on simulated data, check that estimates are consistent. This tests
  the full round-trip (estimation + simulation + re-estimation) and is more
  elaborate than direct parameter recovery since it exercises `simulate.bgms()` /
  `simulate_mrf()` in the loop.

**Current location (2026-03-10):**
- GGM correctness/convergence checks: `tests/testthat/test-bgm.R`
- Slow validation checks: `tests/testthat/test-validation-slow.R`
- Compliance regression checks: `tests/compliance/test_compliance.R`

**Note:** Legacy reference `dev/tests/test-simulation-recovery.R` is stale.

---

## 3. Direct MCMC algorithm tests

**Date:** 2025-02-24
**Status:** Open

**Problem:** All current tests go through the full `bgm()` pipeline (R dispatch →
C++ model → MCMC runner → sampler → output). There are no isolated tests of the
MCMC algorithms themselves (`nuts_step()`, `hmc_step()`, `metropolis_step()`,
leapfrog integrator, adaptation controllers).

**What's needed:**
- Unit tests for NUTS/HMC on simple known targets (e.g., multivariate normal)
  where we can verify step acceptance rates, trajectory lengths, energy conservation
- Leapfrog integrator accuracy: verify symplecticity on a harmonic oscillator
- Mass matrix adaptation: feed known variance, check convergence to true diagonal
- Step-size adaptation (dual averaging): verify convergence to target acceptance rate
- Robbins-Monro proposal-SD tuning: verify convergence for MH on a known target
- Warmup schedule stage transitions: verify correct handoff between stages

**Challenge:** The MCMC code is C++ only, so tests need either:
- Rcpp test exports for individual algorithm functions, or
- A lightweight C++ test harness (e.g., Catch2), or
- R-level tests that infer algorithm behavior from `bgm()` output (indirect)

---

## 4. Function documentation standard

**Date:** 2026-02-25
**Status:** Partially addressed

**Problem:** The package has no explicit standard for what function documentation
should contain, how detailed it should be, or how it should differ between
user-facing and internal code. This makes it hard for both human contributors
and AI assistants to write consistent, useful documentation.

**What's needed:** A written standard covering at minimum:

**Exported R functions** (user-facing, roxygen2):
- What each `@param` should describe (type, constraints, default behavior)
- `@return` format (structure, dimensions, names, class)
- `@examples` requirements (runnable, representative use cases)
- `@details` — when to include mathematical formulation, algorithmic notes,
  or references to papers
- `@seealso` — when to cross-reference related functions
- Error messages — should be self-explanatory to end users unfamiliar with
  the source code

**Internal R functions** (not exported, roxygen2 with `@keywords internal` or `@noRd`):
- Purpose and role within the package architecture (what calls it, why it exists)
- `@param` documentation — types, expected structure, upstream origin
- `@return` — what downstream code expects
- Side effects — options checked, warnings issued, data modified in place
- Less emphasis on examples, more on call context

**C++ functions and classes** (Doxygen `/** */`):
- Class-level doc: design intent, ownership semantics, relationship to other classes
- Public methods: `@param`, `@return`, pre/post-conditions
- Private methods: what they compute, mathematical formulas where applicable
  (see `ggm_model.h` — ~30 undocumented private methods with dense math)
- Free functions: algorithm description, `@param`, `@return`
- Reference implementations to follow: `omrf_model.h` (good), `ggm_model.h` (needs work)

**Action:** Write a concise documentation guide (e.g., `dev/documentation_standard.md`)
that can be referenced by both human developers and AI assistants working on the codebase.

**Note (2026-03-10):** Documentation conventions exist in `.github/copilot-instructions.md`,
but a dedicated `dev/documentation_standard.md` has not been created.

---

## 5. GGM: sample the means

**Date:** 2025-06-18 (moved from inline TODO in `ggm_model.h`)
**Status:** Open

**Problem:** The GGM constructor computes `suf_stat_ = X'X` directly,
which implicitly assumes the data are mean-centered. The algorithm
currently does not sample the mean vector. If users supply non-centered
data, the sufficient statistics are incorrect.

**Action:** Adjust the GGM sampling algorithm to also sample the means,
or document the centering requirement and enforce it at the R interface.

---

## 6. Detailed computational manuals

**Date:** 2026-02-25
**Status:** Open

**Problem:** The package implements several non-trivial computational methods
(NUTS/HMC for the GGM, Gibbs sampling for the OMRF, Blume-Capel pseudolikelihood,
MFM-SBM edge priors, Cholesky rank-1 updates, missing-data imputation) but the
algorithmic detail lives only in the C++ source code. Roxygen and Doxygen
documentation describe *what* functions do and *how* to call them, but not *why*
the algorithms work, what their mathematical foundations are, or how the pieces
fit together. This level of detail does not belong in API docs — it needs
dedicated long-form manuals.

**What's needed:**

- **GGM sampling manual:** Full derivation of the G-Wishart / BDgraph-style
  sampler, the NUTS/HMC integration, the precision-matrix parameterization,
  the Cholesky update strategy, and the edge selection mechanism. Cover the
  log-posterior, gradient computation, and mass matrix adaptation.

- **OMRF sampling manual:** The Gibbs/MH sampling scheme for ordinal MRFs,
  pseudolikelihood computation, the Blume-Capel model and its normalization,
  the relationship between interaction parameters and sufficient statistics.

- **Edge prior manual:** The MFM-SBM edge prior, its Gibbs sampler for block
  allocations and probabilities, the partition coefficient computation, and
  how it integrates with the model selection step.

- **Missing-data imputation manual:** The imputation step within the MCMC
  loop, how it interacts with the likelihood computation, and what assumptions
  it makes about the missingness mechanism.

- **Model comparison manual:** The bgmCompare framework, how group-level
  parameters are shared vs. separated, the difference prior, and the
  interpretation of posterior output.

**Format options:** Package vignettes (heavy, rendered), `inst/doc/` PDFs,
or `dev/manuals/` markdown files for developer reference. Vignettes are
preferable if the audience includes applied users; markdown is sufficient
if the audience is primarily developers and methodologists.

**Action:** Write one manual per topic. Prioritize the GGM and OMRF sampling
manuals since those contain the most complex algorithms.

---

## 7. Stage 3b double-update of interactions

**Date:** 2026-02-25
**Status:** Open

**Problem:** During warmup Stage 3b the OMRF sampler updates interaction
parameters twice per iteration. This was inherited from the legacy backend
and may be unnecessary, but removing it risks changing posterior behavior.
The double update is not harmful — it just does extra work.

**Action:** Investigate whether the second update changes mixing or
convergence. Benchmark single vs. double update on representative datasets.
Only remove if there is no measurable effect on posterior quality.

---

## 8. Model-level warmup_info()

**Date:** 2025-02-25
**Status:** Open (Phase-5 candidate)

**Problem:** `run_mcmc_chain()` decides whether to tune proposal SDs by
string-matching on `config.sampler_type` (`"nuts"`, `"hmc"`, etc.).
This encodes model knowledge in the runner.

**Proposal:** Add a `needs_proposal_sd_tuning()` virtual to `BaseModel`
(or a richer `warmup_info()` struct). The runner queries the model
instead of pattern-matching on sampler strings. Most valuable when
Phase 5 introduces new model types with different warmup needs.

**Origin:** Phase 4, item 4.3.

---

## 9. Extract shared edge-selection logic from models

**Date:** 2025-02-25
**Status:** Open (Phase-5 candidate)

**Problem:** Both `OMRFModel::update_edge_indicators()` and future models
follow the same pattern: loop over variable pairs, compute a likelihood
ratio, apply spike-and-slab accept/reject. Only the likelihood ratio
computation and state update differ between models. This shared logic
could be extracted into a reusable component in `mcmc/`.

**Complication:** GGM ties edge toggles to parameter updates (PD
constraint), so it cannot use the generic pattern.

**When:** Revisit when the mixed MRF (Phase 5) arrives — that will
reveal exactly which parts are truly shared.

**Origin:** Phase 4, item 4.4.

---

## 10. (Placeholder for future items)

**Status:** N/A

Add new items here as they are discovered during refactoring.

---

## 11. Lazy MCMC diagnostics with caching + FFT-based ESS/Rhat

**Date:** 2026-03-10
**Status:** Open

### Problem

`bgm()` computes full MCMC diagnostics (ESS, Rhat, MCSE) for every
parameter before returning. This post-processing dominates wall-clock
time: profiling on Wenchuan (13 variables, 100k iter, 4 chains) shows
sampling takes ~12s but `build_output_bgm()` takes **147s (76% of
total)**. The bottleneck is `coda::effectiveSize()`, which calls
`spectrum0.ar()` → `ar.yw()` → `embed()`, building an enormous lag
matrix per parameter per chain. With edge selection the code also
computes summaries multiple times for the same data:

- `summarize_fit()` calls `summarize_indicator()` and `summarize_slab()`
- `summarize_fit()` then calls `summarize_pair()`, which internally
  calls `summarize_indicator()` and `summarize_slab()` again
- `build_output_bgm()` calls `summarize_indicator()` a third time

Users who just want the fitted model (posterior means, raw samples)
wait 3+ minutes for diagnostics they may never inspect.

### Profiling evidence

| Function | Self time | Total time | % of total |
|---|---|---|---|
| `build_output_bgm` | — | 147s | 76% |
| `compute_rhat_ess` | — | 137s | 71% |
| `coda::effectiveSize` | — | 133s | 69% |
| `embed` (lag matrix) | 66s | 75s | 39% |
| `ar.yw` (AR fit) | 2s | 107s | 55% |

### Solution: two changes

**A. Lazy computation with environment-based caching.**

At fit time, compute only the cheap quantities: posterior means
(`colMeans`), posterior mean matrices, arguments, raw samples, NUTS
diagnostics. Store a cache environment in the fit object:

```r
results$cache = new.env(parent = emptyenv())
```

An internal function `ensure_summaries(fit)` checks the cache, computes
diagnostics on the first call, and stores them. Subsequent calls return
the cached result. Environments have reference semantics, so mutations
persist without the user reassigning the fit object.

Consumers that tap into the cache:
- `summary.bgms()` — reads `posterior_summary_*` fields
- `extract_ess()` — reads `$n_eff` from summaries
- `extract_rhat()` — reads `$Rhat` from summaries

Consumers that remain untouched (no diagnostics needed):
- `coef()` — only reads `posterior_mean_*`
- `print.bgms()` — only reads `$arguments`

This matches the pattern used by cmdstanr, rstan, and brms: store raw
draws at fit time, compute diagnostics on first request.

**B. Replace `coda::effectiveSize` + `coda::gelman.diag` with
FFT-based autocorrelation.**

R's built-in `fft()` computes autocovariance in O(n log n) without
building lag matrices. Benchmarks on the Wenchuan scenario (78
parameters, 100k iter, 4 chains):

| Method | Time | Speedup |
|---|---|---|
| `coda` (current) | 9.3s | 1x |
| `posterior` package | 31.4s | 0.3x |
| Direct FFT | 3.1s | **3x** |

The `posterior` package is slower and adds heavy dependencies. A direct
FFT implementation uses only base R and produces ESS/Rhat values
consistent with coda estimates.

The algorithm:
1. Per-chain autocovariance via FFT (zero-pad, FFT, multiply by
   conjugate, inverse-FFT).
2. Average autocovariances across chains.
3. Initial positive sequence estimator (Geyer 1992) to truncate
   the autocorrelation sum.
4. ESS = (nchains × niter) / tau, where tau is the integrated
   autocorrelation time.
5. Rhat = sqrt(var_hat / W), the standard Gelman-Rubin statistic.

### Expected impact

- `bgm()` returns in ~12s (sampling only) instead of ~193s.
- First `summary()` or `extract_ess()` call takes ~3s instead of ~147s.
- Subsequent diagnostic calls: instant (cached).
- `coda` dependency can potentially be dropped (check other uses first).

### Files to change

- `R/mcmc_summary.R` — replace `compute_rhat_ess()` with FFT version;
  eliminate redundant calls in `summarize_fit()` / `summarize_pair()`
- `R/build_output.R` — remove eager `summarize_fit()` / `summarize_indicator()`
  calls from `build_output_bgm()` and `build_output_mixed_mrf()`; add
  cache env; compute only posterior means eagerly
- `R/bgms-methods.R` — `summary.bgms()` calls `ensure_summaries()`
- `R/extractor_functions.R` — `extract_ess()` / `extract_rhat()` call
  `ensure_summaries()`
- `R/bgmcompare-methods.r` + compare extractors — same pattern
- `DESCRIPTION` — potentially remove `coda` from Imports

---

## 12. Empirical validation of log adjacent-category odds ratios

**Date:** 2026-03-14
**Status:** Open

**Problem:** The OMRF model implies that the log adjacent-category odds
ratio for each variable pair (i, j) is constant across all adjacent
category pairs (h, h+1) × (q, q+1) and equals 2σ_ij = 2K_ij. A
preliminary analysis on the Wenchuan data (`dev/wenchuan_logodds_check.R`)
compared empirical marginal cross-ratios from bivariate contingency
tables against the model-implied value 2K. Per-cell correlation was
near zero (r = 0.039); edge-mean correlation was moderate (r = 0.568).

The mismatch is expected: 2K is the *conditional* log cross-ratio
(controlling for all other variables), while the empirical values
computed from bivariate tables are marginal. A fair comparison requires
either (a) simulating data from the fitted model and computing marginal
cross-ratios from the simulated tables as a posterior predictive check,
or (b) computing conditional empirical cross-ratios (e.g., stratified
or residualised).

**Next steps:**

- Simulate replicate datasets from the fitted bgm posterior and compute
  marginal empirical cross-ratios from the simulated contingency tables.
  Compare these to the observed empirical cross-ratios (posterior
  predictive check).
- Explore conditional empirical cross-ratios (stratified by other
  variables or using pseudolikelihood residuals).
- Consider turning this into a vignette or diagnostic plot function.

---

## 13. Rank-2 Woodbury covariance update in GGM

**Discovered:** 2025-07 (cost analysis of element-wise MH vs HMC)

**Problem:** After each accepted Cholesky update (edge or diagonal),
the current implementation recomputes the full covariance matrix from
scratch:

```cpp
covariance_ = arma::inv(cholesky_of_precision_);
covariance_ = covariance_ * covariance_.t();
```

This is $O(p^3)$ per accepted move. Since a full sweep proposes
$O(p^2)$ moves, the per-sweep cost becomes $O(p^5)$ in the worst case
(dense graph). The design target for element-wise MH is $O(p^2)$ per
element / $O(p^4)$ per sweep.

**Solution:** Replace the full recomputation with a rank-2
Sherman-Morrison-Woodbury update to $\Sigma$ directly. A single
Cholesky column update changes at most two entries of $\Phi$
(one diagonal + one off-diagonal), which corresponds to a rank-2
perturbation of $K = \Phi^\top\Phi$. The Woodbury identity gives:

$$\Sigma_{\text{new}} = \Sigma - \Sigma U (I + V \Sigma U)^{-1} V \Sigma$$

where $U, V$ encode the rank-2 change to $K$. The inner matrix is
$2 \times 2$, so the update is $O(p^2)$.

This eliminates `inv_cholesky_of_precision_` entirely — the model
would maintain $\Sigma$ directly and update it incrementally.

**Files affected:**

- `src/bgm/ggm_model.cpp`: `cholesky_update_after_edge()`,
  `cholesky_update_after_diag()`
- `src/bgmCompare/mixed_mrf_metropolis.cpp`:
  `cholesky_update_after_precision_edge()`,
  `cholesky_update_after_precision_diag()`

**Risk:** Accumulated floating-point drift in $\Sigma$ over many
updates. May need periodic full recomputation (e.g., every $p^2$
accepted moves) as a correction step.

---

## 14. Standardized cross-type partial associations in mixed MRF

**Placeholder.** Details to be added.

---

## 13. Rank-2 Woodbury covariance update in GGM element-wise MH

**Date:** 2026-03-17
**Status:** Open

**Problem:** After each accepted off-diagonal precision update, the GGM
model recomputes the inverse Cholesky factor and full covariance matrix
from scratch:

```cpp
arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
```

The triangular inverse is `O(p^2)` but the covariance product is `O(p^3)`.
This runs after every accepted element, so the per-sweep cost is
`O(p^2 × p^3) = O(p^5)` instead of the intended `O(p^4)`.

The same issue exists in the MixedMRFModel continuous precision updates
(`cholesky_update_after_precision_edge` / `_diag`).

**Background:** The Cholesky rank-1/2 update itself (Givens/hyperbolic
rotations) is already `O(p^2)` — this was a major improvement over the
original approach in the mixedGM package, which permuted the Cholesky
matrix to move the target element to the bottom-right corner, computed
the Roverato constants from that position, updated, and permuted back.
The in-place update avoids the permutation overhead entirely. The
covariance recomputation is the remaining `O(p^3)` bottleneck.

**Rejected approach — Sherman-Morrison-Woodbury:**

After the rank-2 Cholesky update, apply two rank-1 Sherman-Morrison
updates to Sigma directly (`O(p^2)` each). However, this accumulates
floating-point drift across `p(p-1)/2` updates per sweep over thousands
of iterations. The denominator `1 + u^T Sigma u` can approach zero,
causing catastrophic cancellation. Periodic full recomputation mitigates
but does not eliminate the risk. **Not recommended.**

**Preferred approach — on-demand Sigma entries via triangular solves:**

The Roverato constants for edge (i,j) require only three entries of
Sigma: `sigma_ii`, `sigma_jj`, and `sigma_ij`. These can be extracted
from the Cholesky factor L without maintaining a full covariance matrix:

```
Solve L y = e_j  (forward substitution, O(p^2))
Solve L^T x = y  (back substitution, O(p^2))
=> x = column j of Sigma = K^{-1} e_j
```

Two column solves per edge proposal (columns i and j) give `O(p^2)` per
edge, `O(p^4)` per sweep. Each solve uses the current exact Cholesky
factor — no accumulated error.

The inverse Cholesky factor `R^{-1}` and full `covariance_matrix_` can
be dropped entirely. The `get_constants()` / `get_precision_constants()`
functions would call the on-demand solver instead of reading from
a cached covariance matrix.

**Expected impact:** Per-sweep cost drops from `O(p^5)` to `O(p^4)`,
which is the design-intended complexity. For `p = 20` this is a ~20x
reduction in the covariance maintenance cost. Numerically exact.

---

## 14. Standardized cross-type partial associations in mixed MRF

**Date:** 2026-03-18
**Status:** Open

**Problem:** The mixed MRF reports three types of partial associations
on a unified scale in `coef(fit)$pairwise`. The continuous–continuous
block has a natural standardized measure (partial correlations via
`extract_partial_correlations()`), and the discrete–discrete block has
log adjacent-category odds ratios (`extract_log_odds()`). The cross-type
block ($\boldsymbol{\Phi}_{xy}$) has no standardized counterpart —
its magnitude depends on the scale of the continuous variable.

**Proposal:** A rescaled log-odds that divides by the conditional
standard deviation of the continuous variable:

$$
\frac{2\,\phi_{ij}^{xy}}{\sqrt{\theta_{jj}}}
$$

where $1/\theta_{jj}$ is the residual variance of continuous variable
$y_j$. This measures how the log adjacent-category odds of $x_i$
change per conditional standard deviation of $y_j$.

**Properties:**

- Scale-invariant on the continuous side (rescaling $y_j$ changes both
  $\phi_{ij}^{xy}$ and $\theta_{jj}$, but the ratio is invariant).
- Interpretable: log-odds change per conditional SD of the continuous
  variable.
- Parallels partial correlations, which standardize $\theta_{ij}$ by
  $\sqrt{\theta_{ii}\,\theta_{jj}}$ — here only the continuous side
  needs standardization since the discrete side has natural category
  increments.
- Directly computable from existing MCMC output ($\phi_{ij}^{xy}$ and
  $\theta_{jj}$ are both sampled).

**Action:** Implement an extractor (e.g., `extract_standardized_cross()`)
that computes the rescaled log-odds from posterior samples. Document the
standardization in the mixed MRF guide page.

**Files to change:**
- `src/models/ggm/ggm_model.h` — remove `inv_cholesky_of_precision_`;
  add periodic recomputation counter
- `src/models/ggm/ggm_model.cpp` — replace full recomputation in
  `cholesky_update_after_edge()` and `cholesky_update_after_diag()`
  with Sherman-Morrison updates to `covariance_matrix_`
- `src/models/mixed/mixed_mrf_metropolis.cpp` — same changes in
  `cholesky_update_after_precision_edge()` and `_diag()`

## 15. Extractor functions for Bayes factors

**Date:** 2026-03-18
**Status:** Open

**Problem:** `bgm()` and `bgmCompare()` return posterior inclusion
probabilities via `coef(fit)$indicator`, but there is no dedicated
extractor that returns the Bayes factors directly. Users must compute
them manually from the inclusion probabilities and the prior odds.
For `bgmCompare()`, the same applies to group-difference Bayes factors.

**Action:** Implement extractor functions (e.g., `extract_bayes_factors()`)
that compute inclusion Bayes factors from the posterior inclusion
probabilities and the prior inclusion odds for both `bgm` and
`bgmCompare` objects. Return a matrix with the same layout as
`coef(fit)$indicator`.

---

## 16. Check Cholesky downdate failure sentinel

**Date:** 2026-03-18
**Status:** Open

**Problem:** The low-level downdate kernel sets a failure sentinel
(`R[1] = -2.0`) when a hyperbolic downdate would break
positive-definiteness, but the model-level update path does not
currently check this sentinel after calling `cholesky_downdate()`.
As a result, a failed downdate could propagate invalid factor/covariance
state silently.

Relevant locations:

- `src/math/cholupdate.cpp` (`chol_up`, downdate failure signal)
- `src/models/ggm/ggm_model.cpp` (`cholesky_update_after_edge`,
  `cholesky_update_after_diag`)
- `src/models/mixed/mixed_mrf_metropolis.cpp`
  (`cholesky_update_after_precision_edge`,
  `cholesky_update_after_precision_diag`)

**Action:** Add explicit post-downdate checks and define failure
handling semantics. Candidate strategies:

- Reject the move and restore previous state.
- Fallback to full refactorization/recompute from `precision_matrix_`.
- Raise a diagnostic warning/error counter for debugging.

Prefer deterministic recovery over silent continuation.

---

## 17. Dense / block-diagonal mass matrix for NUTS

**Date:** 2026-03-19
**Status:** Open

**Context:** The NUTS sampler uses a diagonal mass matrix. There is
no architectural barrier to a dense or block-diagonal mass matrix:
the dimension-change mechanism (extracting subvectors for active
parameters) generalizes to extracting principal submatrices, which
preserve positive-definiteness.

A diagonal mass matrix is a speed-for-accuracy trade-off. Each
leapfrog step costs $O(d)$ for velocity/momentum/kinetic-energy
with a diagonal, vs $O(d^2)$ with dense. For models where
parameters have strong cross-correlations, the improved
preconditioning from a dense mass matrix can reduce the number of
leapfrog steps enough to offset the per-step cost.

**Proposal:** Support a block-diagonal mass matrix in the NUTS
leapfrog integrator, where blocks correspond to natural parameter
groups:

- **GGM:** one dense block per column group $(f_q, \psi_q)$.
  Blocks are at most $(q{-}1) \times (q{-}1)$ — small.
- **OMRF:** one block per variable (main + pairwise effects
  involving that variable), or a single block for all pairwise
  effects. Block structure TBD based on correlation analysis.

**Cost analysis:** For the GGM, total leapfrog overhead is
$O(\sum d_q^2) \le O(p^3/3)$ per step — negligible compared to
the $O(p^3)$ gradient. Cholesky of each block: $O(d_q^3)$ per
column, only on graph change. For the OMRF, block sizes depend
on grouping strategy.

**Implementation scope:**

- `NUTSSampler` / leapfrog: accept block-diagonal structure for
  momentum sampling ($p \sim N(0, M)$), velocity ($M^{-1}p$),
  and kinetic energy ($p^T M^{-1} p$).
- `DiagMassMatrixAccumulator`: generalize to block-diagonal
  Welford accumulator (learn per-block covariance during stage 2).
- `get_active_inv_mass()`: return block structure or Cholesky
  factors of blocks.
- Dimension changes on edge toggle: extract principal submatrices
  of affected blocks, re-Cholesky.

**Motivation from GGM NUTS:** The GGM free-element Cholesky
parameterization rotates through null-space bases $N_q$ when
edges toggle. A diagonal rotation ($\sum_j (N_q)_{jk}^2
\sigma^2_{\phi_{jq}}$) captures marginal variances correctly but
discards off-diagonal correlations induced by $N_q$. The full
block rotation $N_q^T \operatorname{diag}(\sigma^2) N_q$ captures
all information available from the Welford input. The GGM NUTS v1
uses diagonal rotation; this item upgrades to block-diagonal.

---

## 18. Rao-Blackwellized ESS for edge indicators

**Date:** 2026-03-20
**Status:** Open

**Problem:** The transition-based mixture ESS (`n_eff_mixt`) breaks
down when the marginal inclusion probability is near 0 or 1. With
few transitions the estimators $\hat{a}$ and $\hat{b}$ are noisy,
producing volatile or undefined ESS. The AR-spectral ESS (`n_eff`)
handles the boundary correctly but ignores the two-state structure.
Neither measure captures the full information content of the sampler
output.

**Idea:** Construct a Rao-Blackwellized continuous chain from the
MH acceptance probabilities. At each iteration the sampler computes
an acceptance ratio $\alpha_t = \min(1, r_t)$ for the proposal to
flip the edge indicator. Define

$$\tilde{p}_t = \gamma_t(1 - \alpha_t) + (1 - \gamma_t)\,\alpha_t$$

which gives the conditional inclusion probability implied by the
current state and the proposed move. This chain has the same
expectation as $\gamma_t$ but lower variance, and standard AR-based
ESS applies without transition counting.

**Complications:**

- The current sampler discards $\alpha_t$ after the accept/reject
  step (see `update_edge_indicator_parameter_pair()` in
  `ggm_model.cpp`). Storing one float per edge per iteration adds
  $O(E \cdot T)$ memory.
- Because the sampler proposes $\omega'$ from the Cauchy slab
  (K=1 importance sample), each $\alpha_t$ is noisy: proposals
  from the heavy tails produce near-zero acceptance even when the
  marginal conditional probability is high. This noise inflates
  the variance of $\tilde{p}_t$ and deflates its ESS relative to
  what a better proposal (or K>1 importance sampling) would give.
- An importance-sampling extension drawing K>1 slab values per
  edge per iteration would reduce this noise but costs K extra
  likelihood evaluations per edge per step — expensive at
  $O(p^2)$ each.
- Adapting the importance proposal (e.g., normal centered at the
  conditional mode of $\omega$) would improve IS efficiency but
  adds mode-finding cost.

**Action:** Investigate whether storing $\alpha_t$ and computing
ESS on the Rao-Blackwellized chain $\tilde{p}_t$ gives useful
diagnostics in practice, especially near the boundary. Compare
against the current dual-ESS approach (`n_eff` + `n_eff_mixt`)
on simulated examples with known inclusion probabilities.

## 19. Remove HMC

**Date:** 2026-03-23
**Status:** Open

## 20. SBC with edge selection

**Date:** 2026-03-25
**Status:** Open

**Problem:** The current SBC tests run with `edge_selection = FALSE`.
SBC with edge selection is harder because the spike-and-slab prior
produces tied ranks at zero for excluded edges. Standard rank
uniformity tests (KS, chi-squared) break down under point masses.

**Action:** Explore whether SBC can be adapted for the edge-selection
setting. Possible approaches: (1) condition on the graph and test
only the slab components, (2) use a discrete calibration check for
inclusion indicators separately, (3) use a mixed
discrete-continuous rank statistic.
**Status:** Open