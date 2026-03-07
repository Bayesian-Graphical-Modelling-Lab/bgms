# Hybrid NUTS+MH Sampler for Mixed MRF (integrated plan)

Integrates feedback from review1.md and review2.md into the original
nuts-plan.md. Open points that need decisions are marked **OPEN POINT**.

## Overview

Split each MCMC iteration into two phases:

1. **NUTS phase** — jointly sample (μx, μy, Kxx, Kxy) using NUTS, treating Kyy as fixed
2. **Metropolis phase** — component-wise update Kyy using the existing Cholesky-based MH

This is a Gibbs-within-MCMC scheme: NUTS handles the high-dimensional
unconstrained block, Metropolis handles only the q(q+1)/2 constrained
precision matrix entries.

## Why it works

| Parameter            | Constraint | Current sampler        | Proposed              |
|----------------------|-----------|------------------------|-----------------------|
| μx (thresholds)      | none      | component MH           | **NUTS**              |
| μy (means)           | none      | component MH           | **NUTS**              |
| Kxx (discrete int.)  | none      | component MH           | **NUTS**              |
| Kxy (cross int.)     | none      | component MH           | **NUTS**              |
| Kyy off-diagonal     | SPD       | rank-2 Cholesky MH     | MH (unchanged)        |
| Kyy diagonal         | SPD       | rank-1 Cholesky MH     | MH (unchanged)        |

With Kyy fixed, the pseudo-log-posterior over the remaining parameters is
smooth and unconstrained — exactly the setting where NUTS excels.

## Gradient derivation

The conditional pseudo-log-posterior (with Kyy fixed) is:

    ℓ(θ) = Σ_s log p(x_s | x_{-s}, y)   [OMRF conditionals]
          + log p(y | x)                   [GGM conditional]
          + log π(θ)                       [priors]

### OMRF conditional gradients

For discrete variable s with rest score:

    r_{i,s} = Σ_{t≠s} Kxx_{st} x_{i,t} + 2 Σ_j Kxy_{sj} y_{i,j}

- ∂/∂μx_{s,c}: counts minus expected counts — identical to the existing
  OMRF gradient via `compute_logZ_and_probs_ordinal`.
- ∂/∂Kxx_{st}: Σ_i x_{i,t} · [(x_{i,s}+1) - E_s[c+1 | r_{i,s}]]
  (+ symmetric contribution from ℓ_t) — same pattern as OMRF pairwise gradient.
- ∂/∂Kxy_{sj}: 2 Σ_i y_{i,j} · [(x_{i,s}+1) - E_s[c+1 | r_{i,s}]]
  — same pattern but with continuous observations.

### GGM conditional gradients

With residual D = Y - M where M_i = μy + 2 Σyy Kxy' x_i:

- ∂/∂μy  = Kyy · D'1_n
- ∂/∂Kxy = 2 · X_disc' · D

The simplification Kyy Σyy = I avoids computing Σyy during the gradient.
These formulas have been verified by both reviewers.

### Prior gradients

- μx: Beta prior → sigmoid form (same as OMRF)
- μy, Kxy: Normal(0, 1) prior
- Kxx: Cauchy(0, scale) prior

Key property: `compute_logZ_and_probs_ordinal/blume_capel` already computes
both the log-normalizer and the softmax probabilities needed for the gradient
— the same shared utility used by the OMRF.

## Architecture: 7 phases

### Phase 1: NUTS parameter vector

**Files**: `mixed_mrf_model.h`, `mixed_mrf_model.cpp`

The NUTS vector excludes Kyy:

    [mux_params | active_Kxx | muy | active_Kxy]

- `get_vectorized_parameters()` / `set_vectorized_parameters()` — pack/unpack
  only NUTS params (excluding Kyy).
- `get_full_vectorized_parameters()` — unchanged, returns everything including
  Kyy (for sample storage).
- `get_active_inv_mass()` — returns diagonal mass vector sized for NUTS
  params only.
- New helper: `unvectorize_nuts_params()` to unpack into model matrices.

Edge selection changes the NUTS vector dimension (Gxx/Gxy edge toggles),
mirroring how the OMRF handles it.

#### Vectorization contract (from reviews)

The `BaseModel` interface uses `get_vectorized_parameters()` /
`set_vectorized_parameters()` as the NUTS-sampled state. For the hybrid
scheme these methods must return/accept the NUTS-only vector (excluding Kyy),
while `get_full_vectorized_parameters()` returns the full vector for sample
storage. This matches the OMRF pattern where `get_vectorized_parameters()`
returns active-edge-only params.

**RESOLVED — `full_parameter_dimension()` sizing → option 2.**
`GradientSamplerBase::do_initialize()` calls
`model.full_parameter_dimension()` to size the initial inverse mass vector.
For the hybrid scheme the mass matrix covers only the NUTS dimension, not
Kyy.

**Investigation result:** For the OMRF, `full_parameter_dimension()` returns
`num_main_ + num_pairwise_` (all params) while `parameter_dimension()`
returns only active-edge params. For the GGM, both return the same value
(`dim_`). Neither model would break from a refactoring that uses
`full_parameter_dimension()` for mass-matrix sizing: the OMRF already needs
the full (non-edge-subsetted) dimension for the mass matrix (since
`set_inv_mass` receives a full-dimension vector that the model subsets via
`get_active_inv_mass`), and the GGM has no distinction. The mixed model
would override `full_parameter_dimension()` to return the NUTS-block
dimension (num_main + p(p-1)/2 + q + pq — everything except Kyy), while
a new `storage_dimension()` returns the complete vector length including Kyy.

**Decision:** Refactor `GradientSamplerBase::do_initialize()` to use
`full_parameter_dimension()` for mass matrix sizing (this is already the
case), and make the mixed model's `full_parameter_dimension()` return the
NUTS-block size. Add `storage_dimension()` for sample storage. The OMRF and
GGM are unaffected because their `full_parameter_dimension()` already matches
the mass-matrix dimension.

**Implementation note — `set_inv_mass` / `get_active_inv_mass`.** The OMRF
pattern subsets `inv_mass_` to active-edge-only entries. The mixed-model
version must exclude Kyy entries entirely. Since `set_inv_mass` is called
from the adaptation controller with a full-dimension vector, the model needs
to know which entries correspond to Kyy and skip them. Pin down index
bookkeeping early.

### Phase 2: Gradient computation

**New file**: `mixed_mrf_gradient.cpp`

- `logp_and_gradient(params)` — returns (log-pseudo-posterior, ∇)
- Structure follows OMRF's `logp_and_gradient()` closely:
  1. Unvectorize params into temp matrices (not model state — see below)
  2. Compute rest scores (including 2 Kxy y contribution)
  3. Per variable: `compute_logZ_and_probs` gives both log Z and
     probabilities
  4. Gradient from probabilities: main effects + pairwise weights
     (using x and y)
  5. GGM conditional: residuals → ∇μy and ∇Kxy contributions
  6. Prior gradients
- `gradient(params)` — calls `logp_and_gradient` and discards the
  log-posterior.
- `has_gradient()` returns `true`.

#### Cache invalidation strategy (from reviews)

NUTS calls `logp_and_gradient()` at many points along the leapfrog
trajectory, but only calls `set_vectorized_parameters()` once with the final
accepted state. Following the OMRF pattern:

- `logp_and_gradient()` unvectorizes into **temporaries** (analogous to
  `unvectorize_to_temps` in the OMRF), computing derived quantities (rest
  scores, conditional means) locally. This avoids mutating model state during
  the trajectory.
- The final `set_vectorized_parameters()` updates model state and
  recomputes `conditional_mean_` and (for marginal PL) `Theta_`.

#### Gradient cache scope (from reviews)

Only the OMRF sufficient statistics (category counts, pairwise sums) can be
precomputed and cached. The GGM residual D = Y - M depends on Kxy and μy,
so the GGM gradient must be computed fresh each evaluation.

#### Marginal PL gradient — full plan

The gradient derivation above is for the **conditional** pseudo-likelihood.
The mixed model also supports **marginal** PL
(`pseudolikelihood = "marginal"`), where OMRF conditionals use
Theta = Kxx + 2 Kxy Σyy Kxy' instead of Kxx directly. Marginal PL is
important in practice, so the gradient implementation must support both
modes.

**Marginal PL pseudo-log-posterior (Kyy fixed):**

    ℓ_marg(θ) = Σ_s log p_marg(x_s | x_{-s})  [marginal OMRF conditionals]
              + log p(y | x)                     [GGM conditional]
              + log π(θ)                         [priors]

The marginal OMRF conditional for variable s uses Theta instead of Kxx.
The rest score and denominator structure differ from the conditional form
(see `log_marginal_omrf` in `mixed_mrf_likelihoods.cpp`):

**Rest score (marginal):**

    r_{i,s}^marg = Σ_{t≠s} Θ_{st} x_{i,t} + 2 Σ_j Kxy_{sj} μy_j

where Θ = Kxx + 2 Kxy Σyy Kxy'. Note: Σyy is fixed (Kyy is fixed).

**Denominator differences:** The marginal ordinal denominator absorbs the
quadratic self-interaction Θ_{ss} into the main-effect parameter vector:
main_param(c) = μx_{s,c} + (c+1)² Θ_{ss}. Blume-Capel uses
effective_beta = β + Θ_{ss}.

**Gradient structure (marginal PL):**

1. ∂/∂μx_{s,c}: Same structure as conditional, but using the marginal
   denominator probabilities (which include the Θ_{ss} offset). The
   `compute_logZ_and_probs` utility handles this transparently if given
   the modified main_param.

2. ∂/∂Kxx_{st}: The OMRF part uses Θ-based rest scores. Since
   ∂Θ/∂Kxx_{st} = e_s e_t' + e_t e_s' (just the identity), the gradient
   has the same form as the conditional case but evaluated at the Θ-based
   rest scores and marginal probabilities.

3. ∂/∂Kxy_{sj}: This is the complex case. Kxy appears in two places:
   - The GGM conditional (via conditional_mean_) — same gradient as
     conditional case: 2 [X' D]_{s,j}
   - The marginal OMRF conditionals (via Θ). The chain rule gives:

     ∂Θ_{ab}/∂Kxy_{sj} = 2 [Σyy]_{j,:} Kxy_{b,:}' δ_{as}
                        + 2 Kxy_{a,:} [Σyy]_{:,j} δ_{bs}

   For each OMRF variable a, `log_marginal_omrf(a)` changes through Θ.
   The gradient contribution from the marginal OMRF terms for one variable
   a involves:
   - The sensitivity of the rest score to Θ changes
   - The sensitivity of the denominator (via Θ_{aa}) to Kxy

   In matrix form, for the full Kxy gradient from OMRF marginal terms:

     ∂ℓ_marg_omrf/∂Kxy = 2 Σyy [partial through Θ rest scores]
                       + 2 [partial through Θ diagonal]

   Both partials are computable from the marginal probabilities already
   obtained from `compute_logZ_and_probs`.

4. ∂/∂μy: In the marginal form, μy enters the OMRF rest score as a
   constant bias 2 Kxy_{s,:} μy. So there is an additional μy gradient
   from the OMRF terms:

     ∂ℓ_marg_omrf/∂μy_j = Σ_s 2 Kxy_{s,j} · [obs_mean_s - expected_mean_s]

   Plus the GGM contribution Kyy · D' 1_n as before.

**Implementation plan:**

- Phase 2 implements `logp_and_gradient` with a branch on `use_marginal_pl_`.
- The conditional path is implemented and tested first (numerical gradient
  check).
- The marginal path is implemented second, reusing the same
  `compute_logZ_and_probs` calls but with modified main_param vectors
  (absorbing Θ_{ss}) and Θ-based rest scores.
- The marginal ∂Kxy gradient has 3 contributions: (a) GGM conditional,
  (b) OMRF via Θ rest scores, (c) OMRF via Θ diagonal. Each is O(n·p·q).
- A second numerical gradient check validates the marginal path.
- Testing: both PL modes should produce consistent posterior marginals on
  small test problems.

**Key insight:** With Kyy fixed, Σyy = Kyy⁻¹ is a constant matrix. So
`recompute_Theta()` during gradient evaluation is just a matrix product
Kxx + 2 Kxy Σyy Kxy' — O(p²q) per evaluation but only needed once per
`logp_and_gradient` call (not per variable). The marginal rest scores and
probabilities are then computed per variable as usual.

### Phase 3: Kyy-only Metropolis

**File**: `mixed_mrf_metropolis.cpp`

Create a dedicated `do_kyy_metropolis_step()` that only updates Kyy entries:
  - `update_Kyy_offdiag(i, j)` — existing rank-2 Cholesky MH (unchanged)
  - `update_Kyy_diag(i)` — existing rank-1 Cholesky MH (unchanged)

The hybrid sampler calls `do_kyy_metropolis_step()`.
The pure-MH sampler calls `do_one_metropolis_step()` (full sweep, unchanged).

Robbins-Monro adaptation for Kyy proposal SDs remains.

#### Edge-indicator ownership — investigation and decision

Currently edge-indicator updates happen in two places:

- `ChainRunner::run_iteration()` (line ~59) calls
  `model.update_edge_indicators()` when `schedule.selection_enabled(iter)`.
- `MixedMRFModel::do_one_metropolis_step()` (line ~477) unconditionally calls
  `update_edge_indicators()` at the end of the sweep.

**Investigation result:** The pattern differs across model types:

- **OMRF**: `do_one_metropolis_step()` does NOT call
  `update_edge_indicators()`. Edge updates are handled exclusively by
  `ChainRunner`. The OMRF has a separate `update_edge_indicators()` method
  that uses the shuffled edge order from `prepare_iteration()`.
- **GGM**: `do_one_metropolis_step()` DOES include edge-indicator updates
  inline (calls `update_edge_indicator_parameter_pair` inside the method).
  Its `update_edge_indicators()` override is a no-op. `ChainRunner` calls
  the no-op, so no double-counting.
- **Mixed**: `do_one_metropolis_step()` calls `update_edge_indicators()` at
  step 6. `ChainRunner` also calls it. **This is a double-counting bug in
  the current code** — edge-indicator proposals happen twice per iteration
  when edge selection is active.

**Decision:** Fix the mixed model to match the OMRF pattern:

1. Remove the `update_edge_indicators()` call from
   `MixedMRFModel::do_one_metropolis_step()`.
2. `ChainRunner` remains the sole owner of edge-indicator scheduling.
3. `do_kyy_metropolis_step()` does not include edge-indicator logic.
4. The pure-MH path is not broken because `ChainRunner` already calls
   `update_edge_indicators()` after the sampler step.

**This is an existing bug that should be fixed before the hybrid work.**

**Edge selection ordering**: indicator updates run **after** both NUTS and
Kyy MH complete in each iteration. Birth/death moves for Gxx/Gxy change the
NUTS vector dimension, so they cannot happen mid-trajectory.

### Phase 4: Hybrid sampler

**New file**: `hybrid_nuts_sampler.h`

```
class HybridNUTSSampler : public GradientSamplerBase {
    step(model, iter):
        1. GradientSamplerBase::step()        → NUTS for (mux, muy, Kxx, Kxy)
        2. model.do_kyy_metropolis_step(iter) → MH for Kyy
        return combined result
}
```

- Inherits NUTS adaptation (step size + diagonal mass matrix) from
  `GradientSamplerBase`.
- Kyy proposal-SD adaptation runs inside the model (Robbins-Monro, same
  as now).
- Registered in sampler factory as `"hybrid-nuts"`.

#### Adaptation controller alignment (from reviews)

`GradientSamplerBase::step()` calls `adapt_->update(full_params, ...)`
with the parameter vector from `model.get_full_vectorized_parameters()`.
The adaptation controller's mass matrix must match the NUTS dimension
(excluding Kyy), not the full storage dimension. Resolved by Open Point A:
refactor `GradientSamplerBase` to use `parameter_dimension()` for
mass-matrix sizing.

#### Kyy MH diagnostics — deferred

NUTS diagnostics (tree depth, divergences, energy) reflect only the NUTS
phase. Kyy MH acceptance rates are tracked by Robbins-Monro internally.
Kyy acceptance rates can be computed from the sampler output post-hoc.
Deferred to a later iteration.

### Phase 5: Warmup schedule integration

- Stages 1–3a: NUTS adaptation (step size + mass matrix for NUTS params) +
  Kyy proposal SD warmup.
- Stage 3b: Kyy proposal-SD tuning via Robbins-Monro.
- Stage 3c: Re-adapt NUTS step size with edge selection active.
- Edge selection: Gxx/Gxy toggles change NUTS vector dimension; Gyy
  toggles change Kyy MH sweep.

#### Warmup scheduling for `"hybrid-nuts"` — resolved

`ChainRunner` (line ~34) currently sets `learn_sd` and NUTS diagnostics
only for known sampler types (`"nuts"`, `"hmc"`, `"hamiltonian-mc"`).
Since the user-facing sampler name will be `"nuts"` (which automatically
becomes hybrid for mixed models — see Phase 6), the existing runner logic
recognizes the sampler type without modification. If the internal
implementation uses a distinct `"hybrid-nuts"` key in the sampler factory,
the runner check must be extended to include it.

#### Kyy proposal-SD adaptation timing — resolved

Kyy proposal-SD adaptation via Robbins-Monro runs throughout the entire
warmup phase (whenever `iteration < total_warmup_`), matching the current
adaptive-MH behavior. This is correct: the Kyy MH proposal SDs need
continuous adaptation because Kyy mixing depends on the evolving NUTS
parameter state. No gating by warmup stage.

### Phase 6: Sampler type selection

**Files**: R-level `bgm_spec.R`, C++ interface (`sample_mixed.cpp`)

#### End-to-end sampler wiring — resolved

Currently sampler selection is hard-wired to MH for mixed models across
multiple files:

- `R/bgm_spec.R` line 212: rejects non-`adaptive-metropolis` for `mixed_mrf`
- `R/bgm_spec.R` line 320: sets `sampler_is_continuous = is_continuous || is_mixed`
- `R/validate_sampler.R` lines 111, 119: forces `adaptive-metropolis` for `is_continuous`
- `src/sample_mixed.cpp` line 98: hard-codes `config.sampler_type = "mh"`
- `R/run_sampler.R` line 184: calls `sample_mixed_mrf()` with no sampler-type argument

**Decisions:**

- **User-facing argument:** `sampler = "nuts"` (same as existing). For
  `type = "mixed"`, `"nuts"` automatically becomes the hybrid NUTS+MH
  scheme internally. Users do not need to know about the hybrid distinction.
- **Default:** `"nuts"` is already the default `update_method`. Mixed models
  use the same default.
- **Fallback:** `sampler = "adaptive-metropolis"` remains available as
  pure-MH fallback for mixed models.

**Coordinated changes:**
  1. Remove the `mixed_mrf` rejection in `bgm_spec.R` for non-MH samplers
  2. Allow `"nuts"` for mixed models in `validate_sampler.R`
  3. Extend `sample_mixed_mrf` C++ signature with sampler type + NUTS config
  4. Pass `spec$sampler$update_method` from `run_sampler_mixed_mrf`
  5. In the C++ sampler factory, `"nuts"` + mixed model creates the
     `HybridNUTSSampler` (NUTS for unconstrained block + MH for Kyy)

### Phase 7: Testing and validation

1. **Numerical gradient check**: finite-difference vs analytical gradient.
   Implement *before* the full hybrid sampler as a standalone test calling
   `logp_and_gradient` directly. Highest-risk component.
2. **Posterior equivalence**: compare hybrid-NUTS posteriors against
   pure-MH posteriors on small problems. Consider using the Stan exact-
   likelihood model from `mixedGM` (`inst/stan/mixed_mrf_exact.stan`) as a
   non-pseudolikelihood gold standard.
3. **ESS comparison**: effective sample size per second, hybrid vs pure MH.
4. **Edge selection**: verify birth/death moves work correctly with the
   hybrid scheme.
5. **Edge cases**: q = 1 (scalar Kyy), p = 1 (single discrete variable),
   all edges off.
6. **Gradient consistency across PL types** (when marginal PL is
   implemented): verify conditional and marginal PL gradients produce the
   same posterior marginals.
7. **Vectorization roundtrip tests**: verify `get_vectorized_parameters()` →
   `set_vectorized_parameters()` → `get_vectorized_parameters()` is
   identity for the NUTS block.

## Key design decisions

1. **Both conditional and marginal pseudo-likelihood for NUTS.** The
   conditional form Σ_s p(x_s|x_{-s}, y) · p(y|x) is natural when Kyy is
   fixed: the GGM conditional is standard Gaussian and the OMRF conditionals
   just gain a continuous contribution to rest scores. The marginal form
   uses Θ = Kxx + 2 Kxy Σyy Kxy' in the OMRF conditionals; since Kyy (and
   thus Σyy) is fixed, Θ is a smooth function of Kxx and Kxy. Implement
   conditional first, then marginal, each with numerical gradient checks.

2. **Separate `do_kyy_metropolis_step()` method.** Rather than overloading
   `do_one_metropolis_step()` semantics via `has_gradient()` branching,
   create a dedicated Kyy-only method. Explicit and avoids fragile coupling.

3. **Single NUTS vector (no Kyy).** The NUTS parameter vector packs
   [μx, Kxx, μy, Kxy] with active-edge subsetting. Kyy lives entirely
   outside the NUTS dynamics.

4. **Reuse of `compute_logZ_and_probs`.** The gradient computation reuses
   the existing shared utility for joint log-normalizer + probability
   computation.

5. **`ChainRunner` owns edge-indicator updates.** Edge-indicator proposals
   are removed from `MixedMRFModel::do_one_metropolis_step()` (fixing an
   existing double-counting bug) and handled exclusively by the runner,
   matching the OMRF pattern.

6. **Mixed model's `full_parameter_dimension()` returns NUTS-block size.**
   The mixed model overrides `full_parameter_dimension()` to return the
   NUTS-block dimension (excluding Kyy). A new `storage_dimension()` method
   provides the complete dimension for sample storage. No refactoring of
   `GradientSamplerBase` needed — it already uses `full_parameter_dimension()`
   which will return the correct NUTS-block size.

## Risk and mitigation

| Risk                                        | Mitigation                                              |
|---------------------------------------------|---------------------------------------------------------|
| Gradient bugs (subtle sign errors)          | Numerical gradient check as first test (Phase 7.1)      |
| Mass matrix mismatch after edge toggling    | Follow OMRF pattern: resize + re-init heuristic         |
| Kyy MH mixing slower than NUTS expectations | Kyy is low-dimensional (q(q+1)/2); component MH is ok  |
| Marginal vs conditional PL mismatch         | Validate both decompositions give consistent posteriors  |
| Edge-indicator double-counting              | Single owner (ChainRunner) for indicator updates        |
| Cache invalidation during leapfrog          | Unvectorize to temps in gradient; update state once     |
| Hard-wired MH selection blocks hybrid       | Coordinated R + C++ wiring changes (Phase 6)            |
| Warmup schedule ignores hybrid-nuts         | Update ChainRunner feature detection                    |

## Expected impact

For a model with p discrete and q continuous variables:

- **Current**: ~O(p² + pq + q²) component-wise MH steps per iteration
- **Hybrid**: 1 NUTS trajectory (joint over p·C̄ + p(p-1)/2 + q + pq dims)
  + q(q+1)/2 MH steps
- NUTS typically produces near-independent samples with O(1) effective draws
  per trajectory, vs ~O(d) iterations for random-walk MH to decorrelate in d
  dimensions
- Biggest win when p is large relative to q

## Recommended implementation order

1. Resolve ownership contracts: indicator updates from `do_one_metropolis_step`,
   `do_kyy_metropolis_step` semantics, hybrid sampler hook shape.
2. Wire sampler selection end-to-end (R validation → R dispatch → C++ config
   → sampler factory).
3. Implement mixed-model NUTS vectorization and gradient, with numerical
   gradient checks.
4. Integrate warmup/diagnostics handling for `hybrid-nuts`.
5. Add posterior-equivalence and edge-selection regression tests.

---

## Resolution summary

| # | Topic | Decision |
|---|-------|----------|
| A | `full_parameter_dimension()` sizing | Mixed model overrides to return NUTS-block size; new `storage_dimension()` for full vector. No `GradientSamplerBase` refactoring needed. |
| B | Marginal PL gradient | Implement both conditional and marginal. Conditional first with gradient check, then marginal with its own gradient check. |
| C | Edge-indicator ownership | Fix existing double-counting bug: remove `update_edge_indicators()` from `do_one_metropolis_step()`. `ChainRunner` is sole owner. |
| D | Kyy MH diagnostics | Deferred. Computable from sampler output post-hoc. |
| E | Warmup scheduling | User-facing name is `"nuts"`, so runner recognizes it. Kyy proposal-SD adaptation runs through entire warmup (no stage gating). |
| F | End-to-end sampler wiring | `sampler = "nuts"` (default) becomes hybrid for mixed models. `"adaptive-metropolis"` as pure-MH fallback. Coordinated R + C++ changes. |

### Pre-implementation fix

The edge-indicator double-counting in `MixedMRFModel::do_one_metropolis_step()`
(Point C) should be fixed before starting the hybrid sampler work. This is a
standalone bug fix that benefits the existing pure-MH code path.
