# GGM Mixed Integration Project

**Started:** 2026-02-18  
**Goal:** Integrate postdoc's Gaussian model & class refactoring, then expand to mixed graphical model

---

## Overview

Integrating code from `ggm_mixed` branch into `main`, ensuring:
- Code style and conventions consistency
- Sound architecture/design decisions
- Proper documentation and comments
- Support for future mixed graphical model expansion

---

## Phase 1: Understand the Changes

### Status: IN PROGRESS

### Branches
- **Source branch:** `ggm_mixed` (postdoc's work)
- **Target branch:** `main` (with SBM fix in `fix/sbm-within-between-asymmetry`)

### Summary of Changes

**The branches have significantly diverged:**

| Branch | Unique Commits | Description |
|--------|----------------|-------------|
| `ggm_mixed` | 11 commits | GGM model, class refactoring, skeleton model |
| `main` | 20+ commits | Tests, NUTS optimization (~2x), bug fixes, documentation |

**Net change:** -10,067 / +3,869 lines across 110 files

### Key Additions in `ggm_mixed`
- `src/ggm_model.cpp/h` - Gaussian Graphical Model implementation
- `src/base_model.cpp/h` - Base class for model hierarchy
- `src/skeleton_model.cpp` - Skeleton model implementation  
- `src/sample_ggm.cpp` - GGM sampling interface
- `src/cholupdate.cpp/h` - Cholesky update routines
- `src/SkeletonVariables.h`, `src/adaptiveMetropolis.h` - Support files
- `src/chainResultNew.h` - New chain result structure
- `R/sampleMRF.R` - New R interface

### Removed/Missing from `ggm_mixed`
- `R/simulate_predict.R` (1254 lines) - simulate/predict functionality
- Many test files and fixtures
- Legacy testing infrastructure
- Recent NUTS optimizations (~2x speedup)
- Recent bug fixes and documentation improvements

### Architecture Notes

#### New Class Hierarchy (from `ggm_mixed`)

```
BaseModel (abstract base class)
    │
    ├── GaussianVariables    (GGM implementation)
    ├── SkeletonVariables    (template/documentation for creating new models)
    └── [Future: OrdinalVariables, MixedVariables]
```

#### BaseModel Interface (`base_model.h`)
Pure virtual methods that all models must implement:
- `logp(parameters)` → log posterior density
- `parameter_dimension()` → size of parameter vector

Optional methods (with capability queries):
- `has_gradient()` → true if model provides gradients (for HMC/NUTS)
- `has_adaptive_mh()` → true if model handles its own MH updates
- `gradient(parameters)` → gradient of log posterior
- `logp_and_gradient(parameters)` → fused computation
- `do_one_mh_step()` → one Metropolis-Hastings update
- `get_vectorized_parameters()` → flatten parameters to vector
- `get_vectorized_indicator_parameters()` → flatten edge indicators
- `clone()` → create deep copy
- `set_seed(seed)` → set RNG state

#### Key Design Decisions
1. **Self-contained models**: Each model stores its own data, parameters, priors, and sampler state
2. **Polymorphism via BaseModel**: Samplers can work with any model through the interface
3. **Adaptive MH internal**: Models handle their own proposal adaptation
4. **Clone support**: For parallel chains

#### Current Ordinal Model Architecture (`main`)
- **Not class-based**: Uses free functions with many parameters
- `run_gibbs_sampler_bgm()` - main entry point (43 parameters!)
- Functions in `bgm_logp_and_grad.cpp/h` for likelihood/gradient
- Functions in `bgm_sampler.cpp/h` for MCMC updates
- State passed around as separate matrices/vectors

#### Refactoring Needed for Ordinal Model
1. Create `OrdinalVariables` class inheriting from `BaseModel`
2. Move parameters and state into class members
3. Implement required interface methods
4. Move existing `bgm_logp_and_grad.cpp` functions into class methods
5. Adapt MCMC sampling to work through the class interface

### Critical Issue
`ggm_mixed` was branched from an older `main` and is missing ~20 commits of improvements. A naive merge will cause conflicts and potentially lose recent work.

---

## Phase 2: Create Integration Branch

### Status: NOT STARTED

- [ ] Create new branch from main (after merging SBM fix)
- [ ] Document integration strategy
- [ ] Selectively integrate changes

---

## Phase 3: Refactor & Expand

### Status: NOT STARTED

- [ ] Standardize code style
- [ ] Add/improve documentation
- [ ] Verify architecture supports mixed graphical model
- [ ] Implement mixed graphical model

---

## Phase 4: Testing

### Status: NOT STARTED

- [ ] Verify existing tests pass
- [ ] Add tests for Gaussian model
- [ ] Add tests for mixed graphical model

---

## Session Log

### 2026-02-18
- Created this tracking document
- Starting Phase 1: analyzing diff between `ggm_mixed` and `main`

---

## Questions / Decisions Needed

### Decision 1: Integration Strategy (BLOCKING)

**Options:**

**A. Cherry-pick GGM additions onto main (Recommended)**
- Start from current `main` 
- Manually port the new GGM files and class architecture
- Keeps all recent improvements, tests, and optimizations
- More work but cleaner result

**B. Merge main into ggm_mixed**
- Resolve conflicts, potentially complex
- Risk of regression in either GGM work or recent fixes
- Could work but messy

**C. Rebase ggm_mixed onto main**
- Replay GGM commits on top of current main
- Could have many conflicts to resolve per commit
- Result similar to A but with git history

**Decision:** *(pending)*

---

## References

- Postdoc's branch: `ggm_mixed`
- SBM fix branch: `fix/sbm-within-between-asymmetry`
