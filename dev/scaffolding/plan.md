# R Scaffolding Refactor — bgm_spec Architecture

**Date:** 2025-07-12 (updated 2026-02-25 after Review 1 & Review 2)  
**Branch:** `ggm_mixed` (PR #78)  
**Context:** Phase 3‑H (R-side validation for GGM) expanded into a full
R-layer restructuring to support the upcoming mixed MRF and future variable
types / design structures.

> **Review changelog (2026-02-24):** Incorporated feedback from
> `dev/scaffolding/review 1.md` and `dev/scaffolding/review 2.md`.
> Major changes: (1) added `new_bgm_spec()` constructor +
> `validate_bgm_spec()`, (2) added `collapse_categories_across_groups()`
> for bgmCompare ordinal collapsing, (3) replaced unified `build_output()`
> with three dispatched builders sharing helpers, (4) golden-snapshot
> fixtures are now a prerequisite for Phase A, (5) explicit dual-path
> exit criteria, (6) `$arguments` and simulate/predict regression tests,
> (7) effort estimate updated to ~25–30 commits / 8–10 PRs,
> (8) GGM-specific constraint placement documented.
> See §12 for open discussion points.

---

## 1. Problem Statement

### 1.1 Three structural problems

| # | Problem | Evidence |
|---|---------|----------|
| **P1** | Functions organized by entry point, not concern | `check_model()` (320 lines) and `check_compare_model()` (300 lines) share ~80% of their code. `reformat_data()` (193 lines) and `compare_reformat_data()` (235 lines) share ~70%. `prepare_output_bgm()` and `prepare_output_ggm()` share ~60%. |
| **P2** | Monolithic "check everything" functions | `check_model()` validates variable types, Blume-Capel inputs, prior scales, edge selection priors, and SBM priors — all in one 320-line function with no sub-structure. |
| **P3** | GGM data handling bypasses shared infrastructure | The GGM path in `bgm.R` handles missing-data, seeding, and data prep inline (~60 lines), completely bypassing `reformat_data()`. Adding a new variable type means touching each inline block separately. |

### 1.2 Why it matters now

The roadmap calls for:
- **Mixed MRF** (Phase 5): continuous + ordinal variables in one model
- **GGM-Compare**: extending `bgmCompare` to continuous data
- **New designs**: longitudinal, moderated networks
- **New variable types**: count, binary-specific, etc.

Under the current architecture, every new combination requires copy-pasting
validation logic into yet another monolithic function. The `bgm_spec` pattern
creates a single validated intermediate object that carries all model
information, eliminating duplication and making extensions additive rather than
multiplicative.

### 1.3 Reference patterns

| Package | Pattern | Key insight |
|---------|---------|-------------|
| **brms** | `validate_formula() → brmsterms() → validate_data() → brmsframe() → .standata()` | The `brmsframe` carries all validated information forward; downstream functions never re-validate. |
| **rstanarm** | `validate_family()`, `validate_glm_outcome_support()`, `handle_glm_prior()`, `center_x()` | Each validator is 20–60 lines, focused on one concern. Easy to test, easy to extend. |

---

## 2. Current Architecture

### 2.1 File inventory (R/ layer — 7,428 lines)

| File | Lines | Core functions | Role |
|------|-------|----------------|------|
| `simulate_predict.R` | 1,637 | `simulate_mrf`, `simulate.bgms`, `simulate.bgmCompare`, `predict.bgms`, `predict.bgmCompare` + 6 helpers | Simulation & prediction |
| `extractor_functions.R` | 913 | 20+ `extract_*()` generics + methods | Post-fit extraction |
| `mcmc_summary.R` | 883 | `summarize_fit`, `summarize_fit_compare`, `summarize_indicator`, `posterior_summary_SBM` + 12 helpers | MCMC chain summaries |
| `bgm.R` | 882 | `bgm()` | Main entry point |
| `bgmCompare.R` | 695 | `bgmCompare()` | Group comparison entry |
| `function_input_utils.R` | 671 | `check_model`, `check_compare_model`, 5 helpers | Input validation |
| `output_utils.R` | 623 | `prepare_output_bgm`, `prepare_output_ggm`, `prepare_output_bgmCompare`, `transform_*` | Output formatting |
| `data_utils.R` | 494 | `reformat_data`, `compare_reformat_data`, `data_check`, `compute_*` | Data preprocessing |
| `bgms-methods.R` | 201 | `print.bgms`, `summary.bgms`, `coef.bgms` | S3 methods |
| `bgmcompare-methods.R` | 377 | `print.bgmCompare`, `summary.bgmCompare`, `coef.bgmCompare` | S3 methods |

### 2.2 Pipeline: bgm()

```
User call
  │
  ├─ data_check(x)                          # data_utils.R
  ├─ check_model(x, variable_type, ...)     # function_input_utils.R  [320 lines]
  │    → {variable_bool, baseline_category, edge_selection,
  │       edge_prior, inclusion_probability, is_continuous}
  ├─ handle deprecated args                 # bgm.R inline  [30 lines]
  ├─ validate sampler args                  # bgm.R inline  [40 lines]
  ├─ validate na_action, progress, seed     # bgm.R inline  [20 lines]
  │
  ├─── GGM path ───────────────────────────── bgm.R inline  [60 lines]
  │  ├─ listwise missing removal (inline)
  │  ├─ sample_ggm(...)                     # → C++
  │  ├─ transform_ggm_backend_output()      # output_utils.R
  │  └─ prepare_output_ggm()               # output_utils.R  [140 lines]
  │
  └─── OMRF path ──────────────────────────
     ├─ reformat_data()                     # data_utils.R  [193 lines]
     │    → {x, num_categories, baseline_category,
     │       missing_index, na_impute}
     ├─ compute pairwise_scaling_factors    # bgm.R inline  [40 lines]
     ├─ sample_omrf(...)                    # → C++
     ├─ transform_new_backend_output()      # output_utils.R
     └─ prepare_output_bgm()               # output_utils.R  [180 lines]
```

### 2.3 Pipeline: bgmCompare()

```
User call
  │
  ├─ data_check(x), data_check(y)
  ├─ validate group_indicator / y
  ├─ check_compare_model(...)               # function_input_utils.R  [300 lines]
  │    → {x, group, variable_bool, baseline_category,
  │       difference_prior, inclusion_probability_difference}
  ├─ handle ~12 deprecated args             # bgmCompare.R inline  [80 lines]
  ├─ validate sampler args                  # bgmCompare.R inline  [30 lines]
  ├─ compare_reformat_data()                # data_utils.R  [235 lines]
  │    → {x, group, num_categories, baseline_category,
  │       missing_index, na_impute}
  ├─ compute_counts_per_category()          # data_utils.R
  ├─ compute_blume_capel_stats()            # data_utils.R
  ├─ compute_pairwise_stats()               # data_utils.R
  ├─ compute indices & scaling factors      # bgmCompare.R inline  [80 lines]
  ├─ compute projection matrix              # bgmCompare.R inline  [15 lines]
  ├─ run_bgmCompare_parallel(...)           # → C++
  └─ prepare_output_bgmCompare()            # output_utils.R  [100 lines]
```

### 2.4 Pipeline: simulate.bgms() / predict.bgms()

```
S3 method call
  │
  ├─ extract_arguments(object)              # extractor_functions.R
  │    → arguments$num_variables, $num_categories, $variable_type,
  │      $data_columnnames, $baseline_category, $is_continuous
  │
  ├─── GGM path ───
  │  ├─ simulate: simulate_bgms_ggm() → simulate_mrf() or run_ggm_simulation_parallel()
  │  └─ predict:  predict_bgms_ggm()  → reconstruct_precision() → compute_conditional_ggm()
  │
  └─── OMRF path ──
     ├─ simulate: extract posterior mean/samples → simulate_mrf() or run_simulation_parallel()
     └─ predict:  recode_data_for_prediction() → compute_conditional_probs()
```

### 2.5 Duplication map

| Concern | bgm() side | bgmCompare() side | % shared |
|---------|-----------|-------------------|----------|
| Variable type validation | `check_model` L57–140 | `check_compare_model` L395–470 | ~90% |
| Blume-Capel baseline checks | `check_model` L142–233 | `check_compare_model` L472–560 | ~95% |
| Prior scale validation | `check_model` L236–248 | `check_compare_model` L562–580 | ~90% |
| Edge/difference prior setup | `check_model` L250–370 | `check_compare_model` L582–670 | ~70% |
| Data: missing-data handling | `reformat_data` L5–72 | `compare_reformat_data` L210–290 | ~85% |
| Data: ordinal recoding | `reformat_data` L90–110 | `compare_reformat_data` L320–340 | ~95% |
| Data: BC recoding/validation | `reformat_data` L112–155 | `compare_reformat_data` L342–385 | ~95% |
| Output: MCMC summary | `prepare_output_bgm` L80–100 | `prepare_output_ggm` L220–240 | ~85% |
| Output: param naming | `prepare_output_bgm` L55–78 | `prepare_output_ggm` L200–215 | ~80% |
| Output: raw_samples assembly | `prepare_output_bgm` L150–175 | `prepare_output_ggm` L260–285 | ~90% |
| Scaling factors | `bgm.R` L695–745 | `bgmCompare.R` L455–505 | ~95% |

---

## 3. Target Architecture

### 3.1 Overview

```
User call: bgm(x, ...)  or  bgmCompare(x, y, ...)
  │
  ├─ handle_deprecated_args()               # [new] thin shim
  │
  ├─ bgm_spec() constructor                 # [new] → validated "bgm_spec" object
  │    ├─ validate_data(x)                   # [new] coerce + basic checks
  │    ├─ validate_variable_types()          # [new] parse variable_type → typed info
  │    ├─ validate_baseline_category()       # [new] BC-specific checks
  │    ├─ validate_missing_data()            # [new] shared listwise/impute, all types
  │    │    └─ GGM: enforce no impute        # [new] §4.4
  │    ├─ reformat_ordinal_data()            # [refactored] recode ordinal + BC data
  │    ├─ collapse_categories_across_groups()# [new] Compare-only group collapsing
  │    ├─ validate_edge_prior()              # [new] Bernoulli/Beta-Bernoulli/SBM
  │    ├─ validate_sampler()                 # [new] update_method, warmup, HMC/NUTS
  │    │    └─ GGM: enforce adaptive-metropolis # [new] §4.8
  │    ├─ compute_scaling_factors()           # [new] extracted from inline code
  │    ├─ new_bgm_spec()                     # [new] low-level constructor
  │    └─ validate_bgm_spec()                # [new] cross-field invariants
  │
  ├─ run_sampler(spec)                  # [new] dispatch to C++ based on spec$model_type
  │    ├─ GGM:  sample_ggm(...)
  │    ├─ OMRF: sample_omrf(...)
  │    └─ Compare: run_bgmCompare_parallel(...)
  │
  └─ build_output(spec, raw_out)        # [new] thin dispatcher
       ├─ build_output_bgm(spec, raw)     # [refactored] unified GGM + OMRF builder
       └─ build_output_compare(spec, raw) # [refactored] Compare-specific builder
       (both share helpers: build_posterior_mean_matrix(),
        build_raw_samples_list(), build_arguments(), generate_param_names())
```

> **Design rationale (build_output):** The GGM and OMRF output paths
> share ~80% of their code (MCMC summaries, pairwise/indicator matrices,
> SBM handling, raw_samples assembly). The only differences are
> **(1) parameter naming** (GGM: `"Var (precision)"` vs OMRF:
> `"Var (1)", "Var (2)", ...` per category) and **(2) the main
> posterior mean shape** (GGM: `p × 1` diagonal vs OMRF:
> `p × max_categories`). These are small enough to handle with an
> `if (is_continuous)` branch inside a single `build_output_bgm()`.
> The Compare path has genuinely different structure (group-level
> summaries, difference parameters, projection matrices) and stays
> separate. `build_output()` dispatches to one of these two builders.

### 3.2 The bgm_spec object

An S3 list of class `"bgm_spec"` with a well-defined structure:

```r
bgm_spec <- list(
  # ---- Model type ----
  model_type       = "ggm" | "omrf" | "compare",  # dispatch key

  # ---- Data ----
  data = list(
    x                = <matrix>,        # validated, recoded data
    data_columnnames = <character>,      # variable names
    num_variables    = <integer>,        # ncol(x)
    num_cases        = <integer>,        # nrow(x)
    num_categories   = <integer vector>, # max category per variable (OMRF/BC)
    # Compare-specific:
    group            = <integer vector>, # group membership (1-based)
    num_groups       = <integer>,
    group_indices    = <matrix>,         # [num_groups x 2] start/end (0-based)
    projection       = <matrix>          # eigenvector projection for differences
  ),

  # ---- Variable specification ----
  variables = list(
    variable_type     = <character vector>,  # "ordinal"/"blume-capel"/"continuous"
    is_ordinal        = <logical vector>,    # TRUE = ordinal, FALSE = BC
    is_continuous     = <logical>,           # TRUE if all continuous
    baseline_category = <integer vector>     # BC baseline (0 for ordinal/continuous)
  ),

  # ---- Missing data ----
  missing = list(
    na_action     = "listwise" | "impute",
    na_impute     = <logical>,
    missing_index = <matrix>              # [n_missing x 2] (0-based row, col)
  ),

  # ---- Priors ----
  prior = list(
    # Pairwise / main
    pairwise_scale   = <numeric>,
    main_alpha       = <numeric>,
    main_beta        = <numeric>,
    standardize      = <logical>,
    pairwise_scaling_factors = <matrix>,

    # Edge selection (bgm)
    edge_selection        = <logical>,
    edge_prior            = "Bernoulli" | "Beta-Bernoulli" | "Stochastic-Block" | "Not Applicable",
    inclusion_probability = <matrix>,
    beta_bernoulli_alpha  = <numeric>,
    beta_bernoulli_beta   = <numeric>,
    beta_bernoulli_alpha_between = <numeric>,   # SBM only
    beta_bernoulli_beta_between  = <numeric>,   # SBM only
    dirichlet_alpha       = <numeric>,           # SBM only
    lambda                = <numeric>,           # SBM only

    # Difference selection (bgmCompare)
    difference_selection  = <logical>,
    difference_prior      = <character>,
    difference_scale      = <numeric>,
    inclusion_probability_difference = <matrix>
  ),

  # ---- Sampler ----
  sampler = list(
    update_method     = "nuts" | "adaptive-metropolis" | "hamiltonian-mc",
    target_accept     = <numeric>,
    iter              = <integer>,
    warmup            = <integer>,
    chains            = <integer>,
    cores             = <integer>,
    hmc_num_leapfrogs = <integer>,
    nuts_max_depth    = <integer>,
    learn_mass_matrix = <logical>,
    seed              = <integer>,
    progress_type     = 0L | 1L | 2L
  ),

  # ---- Precomputed structures (OMRF/Compare) ----
  precomputed = list(
    counts_per_category       = <array>,     # Compare only
    blume_capel_stats         = <list>,      # Compare only
    pairwise_stats            = <list>,      # Compare only
    main_effect_indices       = <matrix>,    # Compare only
    pairwise_effect_indices   = <matrix>,    # Compare only
    interaction_index_matrix  = <matrix>,    # Compare only
    num_thresholds            = <integer>    # OMRF only
  )
)
```

The spec is a plain named list (no S3 class). The class tag can be added
later if method dispatch on the spec is ever needed.

### 3.2b The `new_bgm_spec()` constructor and `validate_bgm_spec()`

> **Added after Review 1 & 2:** The spec has ~40 fields across 6 nested
> sub-lists. Without a strict constructor, a missing field produces a
> cryptic `NULL` subscript error deep in `run_sampler()` or `build_output()`.

```r
new_bgm_spec <- function(model_type, data, variables, missing, prior,
                          sampler, precomputed = list()) {
  # Low-level constructor: asserts presence and type of every field.
  stopifnot(is.character(model_type), length(model_type) == 1L)
  stopifnot(model_type %in% c("ggm", "omrf", "compare"))
  stopifnot(is.list(data), is.matrix(data$x))
  stopifnot(is.list(variables), is.list(missing), is.list(prior))
  stopifnot(is.list(sampler))
  # ... type + length assertions for every field ...
  spec <- list(
    model_type  = model_type,
    data        = data,
    variables   = variables,
    missing     = missing,
    prior       = prior,
    sampler     = sampler,
    precomputed = precomputed
  )
  spec
}

validate_bgm_spec <- function(spec) {
  # Cross-field invariant checks:
  # - model_type == "ggm"     → variables$is_continuous == TRUE
  # - model_type == "ggm"     → sampler$update_method == "adaptive-metropolis"
  # - model_type == "ggm"     → missing$na_action != "impute" (not yet supported)
  # - model_type == "compare" → data$group is present, data$num_groups >= 2
  # - edge_selection == TRUE  → edge_prior is not "Not Applicable"
  # - All pairwise_scaling_factors dimensions match num_variables
  invisible(spec)
}
```

`bgm_spec()` is the *user-facing* constructor that calls individual
validators, then calls `new_bgm_spec()` to assemble, then
`validate_bgm_spec()` to verify cross-field invariants.

### 3.3 New file structure

| New file | Contents | Replaces |
|----------|----------|----------|
| `R/bgm_spec.R` | `bgm_spec()` constructor, `print.bgm_spec()` | — |
| `R/validate_data.R` | `validate_data()`, `validate_missing_data()`, `reformat_ordinal_data()` | `data_utils.R` (reformat_data, compare_reformat_data) |
| `R/validate_model.R` | `validate_variable_types()`, `validate_baseline_category()`, `validate_edge_prior()`, `validate_difference_prior()`, `validate_prior_scales()` | `function_input_utils.R` (check_model, check_compare_model) |
| `R/validate_sampler.R` | `validate_sampler()`, `validate_seed()`, `validate_progress()` | Inline code in bgm.R / bgmCompare.R |
| `R/compute_utils.R` | `compute_scaling_factors()`, `compute_indices()`, `compute_projection()` | Inline code in bgm.R / bgmCompare.R |
| `R/run_sampler.R` | `run_sampler()` dispatch, `run_sampler_ggm()`, `run_sampler_omrf()`, `run_sampler_compare()` | Inline C++ call blocks in bgm.R / bgmCompare.R |
| `R/build_output.R` | `build_output()` (thin dispatcher), `build_output_bgm()` (unified GGM + OMRF), `build_output_compare()`, shared helpers: `build_arguments()`, `build_posterior_mean_matrix()`, `build_raw_samples_list()`, `generate_param_names()` | `output_utils.R` (prepare_output_bgm, prepare_output_ggm, prepare_output_bgmCompare) |

Files that stay (mostly) unchanged:
- `R/bgm.R` — becomes thin: deprecated args → `bgm_spec()` → `run_sampler()` → `build_output()`
- `R/bgmCompare.R` — same thinning
- `R/simulate_predict.R` — stays as-is initially; can later accept spec
- `R/extractor_functions.R` — unchanged (works on fit objects, not specs)
- `R/mcmc_summary.R` — unchanged (internal MCMC processing)
- `R/bgms-methods.R`, `R/bgmcompare-methods.R` — unchanged
- `R/data_utils.R` — retains `data_check()`, `compute_counts_per_category()`, `compute_blume_capel_stats()`, `compute_pairwise_stats()`; big functions deleted
- `R/function_input_utils.R` — retains small helpers (`check_positive_integer`, etc.); big functions deleted
- `R/output_utils.R` — retains `transform_*` functions; `prepare_output_*` migrated to `build_output.R`

### 3.4 Validator design

Each validator is a pure function: input → validated output (or error).
Validators do NOT read from parent environments or use `hasArg()`.

```r
# Example: validate_variable_types
# - Called from bgm_spec() for both bgm and bgmCompare
# - Returns typed result; caller decides what to do with is_continuous
validate_variable_types <- function(variable_type, num_variables,
                                     allow_continuous = TRUE) {
  # Single string → replicate
  # Vector → validate each element
  # Returns: list(variable_type, is_ordinal, is_continuous)
}

# Example: validate_baseline_category
# - Only called when is_continuous is FALSE and some variables are BC
validate_baseline_category <- function(baseline_category, x, is_ordinal) {
  # Validate integer, range, replication
  # Returns: integer vector of length num_variables
}

# Example: validate_edge_prior
# - Encapsulates all Bernoulli/Beta-Bernoulli/SBM logic
validate_edge_prior <- function(edge_selection, edge_prior,
                                 inclusion_probability, num_variables,
                                 beta_bernoulli_alpha, beta_bernoulli_beta,
                                 ...) {
  # Returns: list(edge_selection, edge_prior, inclusion_probability, ...)
}
```

### 3.5 How entry points change

#### bgm() — before (882 lines) → after (~120 lines)

```r
bgm <- function(x, ...) {
  # 1. Handle deprecated arguments (~30 lines, unchanged)
  args <- handle_deprecated_bgm_args(...)

  # 2. Build spec (~1 line)
  spec <- bgm_spec(
    x = x,
    variable_type = args$variable_type,
    baseline_category = args$baseline_category,
    # ... all validated args forwarded
    model_type = "single"
  )

  # 3. Sample (~1 line)
  raw_out <- run_sampler(spec)

  # 4. Build output (~1 line)
  output <- build_output(spec, raw_out)

  # 5. NUTS diagnostics (optional)
  if (spec$sampler$update_method == "nuts") {
    output$nuts_diag <- summarize_nuts_diagnostics(raw_out, ...)
  }

  return(output)
}
```

#### bgmCompare() — before (695 lines) → after (~100 lines)

Same pattern, with `model_type = "compare"` and group handling.

---

## 4. Detailed Validator Specifications

### 4.1 validate_data(x, name = "x")

**Input:** Raw user data (matrix, data.frame, or other).  
**Output:** Numeric matrix.  
**Replaces:** `data_check()` (stays as thin alias).

```r
validate_data <- function(x, name = "x") {
  # 1. Coerce data.frame → matrix
  # 2. Check >= 2 rows, >= 2 columns
  # 3. Check all numeric
  # Returns: numeric matrix
}
```

### 4.2 validate_variable_types(variable_type, num_variables, allow_continuous)

**Input:** User's `variable_type` string or vector + column count.  
**Output:** `list(variable_type, is_ordinal, is_continuous)`.  
**Replaces:** Lines 57–140 of `check_model()` and 395–470 of `check_compare_model()`.

Key behavior:
- Single string → replicate to length `num_variables`
- Vector → validate each, check length matches
- `"continuous"` requires all-continuous (no mixing yet; future mixed MRF)
- `allow_continuous = FALSE` for bgmCompare (until GGM-Compare)
- Returns canonical character vector + derived booleans

### 4.3 validate_baseline_category(baseline_category, x, is_ordinal)

**Input:** User's `baseline_category` + data + ordinal flags.  
**Output:** Integer vector of length `num_variables`.  
**Replaces:** Lines 142–233 of `check_model()` and 472–560 of `check_compare_model()`.

Key behavior:
- Scalar → replicate
- Validate integer-ness for BC variables
- Validate within observed data range
- Ordinal variables get `baseline_category = 0L`

### 4.4 validate_missing_data(x, na_action, is_continuous)

**Input:** Data matrix + na_action string + model type.  
**Output:** `list(x, na_impute, missing_index)`.  
**Replaces:** Top section of `reformat_data()` and `compare_reformat_data()`,
plus inline GGM missing-data handling in `bgm.R`.

Key behavior:
- `"listwise"`: remove rows with NAs, warn
- `"impute"`: build missing_index, impute starting values
- GGM + impute: error (not yet supported)
- Shared across all model types

### 4.5 reformat_ordinal_data(x, is_ordinal, baseline_category)

**Input:** Validated data + variable flags.  
**Output:** `list(x, num_categories, baseline_category)`.  
**Replaces:** The per-variable recoding loop in `reformat_data()` (lines
80–195). Does **not** handle the group-conditional collapsing from
`compare_reformat_data()` — that is handled by §4.5b.

Key behavior:
- Ordinal variables: recode to 0-based contiguous categories
- BC variables: validate integer, shift to 0-start, adjust baseline
- Continuous variables: pass through (no recoding)
- Compute `num_categories` per variable

### 4.5b collapse_categories_across_groups(x, group, is_ordinal, num_categories, baseline_category)

> **Added after Review 1 & 2:** `compare_reformat_data()` (data_utils.R
> L290–430) does something fundamentally different from `reformat_data()`:
> it collapses ordinal categories that aren't observed in *all* groups
> via the `observed_scores` matrix. This group-conditional collapsing
> has no counterpart in `reformat_ordinal_data()` and is the most likely
> place for a silent behavioral difference to slip through.

**Input:** Recoded data (output of `reformat_ordinal_data()`), group
membership vector, variable flags, per-variable category counts,
baseline categories.  
**Output:** `list(x, num_categories, baseline_category)` — with
categories collapsed so that only categories observed in *all* groups
survive.  
**Replaces:** The group-aware cross-check and collapsing loop in
`compare_reformat_data()` (lines 290–430).  
**Called by:** `bgm_spec()` when `model_type == "compare"`, immediately
after `reformat_ordinal_data()`.

Key behavior:
- Build `observed_scores` matrix: which categories appear in which groups
- Collapse categories absent from any group
- Renumber remaining categories contiguously (0-based)
- Adjust `baseline_category` if the baseline was collapsed
- Warn when categories are dropped

### 4.6 validate_edge_prior(...)

**Replaces:** Lines 250–370 of `check_model()`.

### 4.7 validate_difference_prior(...)

**Replaces:** Lines 582–670 of `check_compare_model()`.

### 4.8 validate_sampler(update_method, target_accept, iter, warmup, ..., is_continuous)

**Input:** All sampler-related arguments + `is_continuous` flag.  
**Output:** `list(update_method, target_accept, iter, warmup, chains, cores, hmc_num_leapfrogs, nuts_max_depth, learn_mass_matrix, seed, progress_type)`.  
**Replaces:** Inline sampler validation in `bgm.R` and `bgmCompare.R`.

> **Added after Review 2:** GGM currently forces
> `update_method = "adaptive-metropolis"` and ignores user-supplied NUTS/HMC
> settings. This constraint must be explicitly enforced here (not left to
> inline code). When `is_continuous == TRUE`, this validator:
> - Sets `update_method = "adaptive-metropolis"` regardless of user input
> - Warns if user requested NUTS or HMC
> - Ignores `hmc_num_leapfrogs` and `nuts_max_depth`

### 4.9 compute_scaling_factors(num_variables, is_ordinal, num_categories, baseline_category, standardize)

**Replaces:** Duplicated 40-line loop in `bgm.R` (L695–745) and `bgmCompare.R` (L455–505).

---

## 5. The arguments Field: Backward Compatibility

The fit objects (`bgms`, `bgmCompare`) store an `$arguments` list used by:
- `extract_arguments()` (extractor_functions.R)
- All S3 methods (print, summary, coef)
- All simulate/predict methods
- The `easybgm` package (external dependency)

### Current arguments fields (bgms)

```r
arguments = list(
  num_variables, num_cases, na_impute, variable_type, iter, warmup,
  pairwise_scale, standardize, main_alpha, main_beta,
  edge_selection, edge_prior, inclusion_probability,
  beta_bernoulli_alpha, beta_bernoulli_beta,
  beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
  dirichlet_alpha, lambda, na_action,
  version, update_method, target_accept,
  hmc_num_leapfrogs, nuts_max_depth, learn_mass_matrix,
  num_chains, num_categories, data_columnnames,
  baseline_category, pairwise_scaling_factors,
  no_variables,                           # easybgm compat
  is_continuous                           # GGM only
)
```

### Strategy

The `build_output()` function reads directly from the `bgm_spec` to assemble
this `arguments` list. Because the spec contains a superset of all needed
fields, this is a straightforward extraction:

```r
build_arguments <- function(spec) {
  list(
    num_variables    = spec$data$num_variables,
    num_cases        = spec$data$num_cases,
    variable_type    = spec$variables$variable_type,
    # ... map every existing field from spec
    no_variables     = spec$data$num_variables  # easybgm compat
  )
}
```

The `simulate.bgms()` / `predict.bgms()` methods continue reading from
`object$arguments` as before. No downstream code needs to change.

---

## 6. Migration Strategy

### 6.1 Guiding principles

1. **Test-preserving**: All 1,788 tests pass at every checkpoint.
2. **Incremental**: One validator at a time, each in its own commit.
3. **Fixture-first**: Golden-snapshot fixtures are captured *before* any
   extraction begins (see §6.1b). Every validator must reproduce these.
4. **Dual-path**: New validators are called alongside existing code initially,
   with assertions that outputs match. Old code is deleted only after the
   dual-path exit criteria are met (see §6.1c).
5. **No public API changes**: `bgm()` and `bgmCompare()` signatures stay
   identical. Return values stay identical.
6. **No easybgm breakage**: The `$arguments` list keeps all existing fields
   including `no_variables`.
7. **simulate/predict in the regression suite**: `simulate()` and `predict()`
   calls are included in Phase C regression tests, not just `bgm()` output.

### 6.1b Golden-Snapshot Fixtures (prerequisite for Phase A)

> **Added after Review 1 & 2:** The existing 1,788 tests exercise the
> public API end-to-end and may not catch intermediate-state differences.
> Before extracting any validators, we capture the full outputs of
> `check_model()`, `check_compare_model()`, `reformat_data()`,
> `compare_reformat_data()`, and `prepare_output_*()` for representative
> inputs and store them as `.rds` files in `dev/fixtures/scaffolding/`.

**Fixture matrix (minimum 15 cases):**

| # | Model type | variable_type | edge_prior | na_action | Special |
|---|-----------|---------------|------------|-----------|--------|
| 1 | bgm/GGM | continuous | Bernoulli | listwise | — |
| 2 | bgm/GGM | continuous | Beta-Bernoulli | listwise | with NAs |
| 3 | bgm/OMRF | ordinal | Bernoulli | listwise | — |
| 4 | bgm/OMRF | ordinal | Beta-Bernoulli | impute | with NAs |
| 5 | bgm/OMRF | ordinal | SBM | listwise | — |
| 6 | bgm/OMRF | blume-capel | Bernoulli | listwise | custom baseline |
| 7 | bgm/OMRF | blume-capel | Beta-Bernoulli | impute | with NAs |
| 8 | bgm/OMRF | mixed ord+BC | Bernoulli | listwise | — |
| 9 | bgmCompare | ordinal | — | listwise | 2 groups |
| 10 | bgmCompare | ordinal | — | impute | 2 groups + NAs |
| 11 | bgmCompare | blume-capel | — | listwise | 2 groups |
| 12 | bgmCompare | ordinal | — | listwise | >2 groups |
| 13 | bgmCompare | ordinal | — | listwise | categories missing in 1 group |
| 14 | bgmCompare | mixed ord+BC | — | listwise | categories missing + BC |
| 15 | bgm/OMRF | ordinal | none | listwise | edge_selection = FALSE |

Each fixture stores:
- Input arguments
- `check_model()` / `check_compare_model()` return value
- `reformat_data()` / `compare_reformat_data()` return value
- Full `$arguments` list from the fit object

Fixture generation script: `dev/generate_scaffolding_fixtures.R`.

### 6.1c Dual-Path Exit Criteria (gate for Phase C)

> **Added after Review 1 & 2:** The plan previously said "assertions are
> removed in Phase C when the old path is deleted" without defining
> success criteria.

Phase C can only begin when **all** of the following are satisfied:

1. Phase B dual-path assertions have passed for at least **2 full CI
   cycles** (not just local runs).
2. All 15+ golden-snapshot fixtures match between old and new pipelines
   to floating-point tolerance (`testthat::expect_equal()` with default
   tolerance).
3. **`$arguments` regression test passes**: for each model type, the
   `output$arguments` list matches field-by-field between old and new
   pipelines.
4. **simulate/predict smoke tests pass**: for each model type, calling
   `simulate()` and `predict()` on the fit object produces identical
   results between old and new pipelines.
5. All 1,788+ existing tests pass unchanged.

### 6.2 Migration phases

#### Phase A-0: Create golden-snapshot fixtures (prerequisite) ✅ COMPLETE

| Step | What | Status |
|------|------|--------|
| A.0.1 | Write `dev/generate_scaffolding_fixtures.R` that runs all 15+ fixture cases | ✅ `2ca15a4` |
| A.0.2 | Store `.rds` files in `dev/fixtures/scaffolding/` | ✅ `2ca15a4` |
| A.0.3 | Write `tests/testthat/test-scaffolding-fixtures.R` that loads fixtures and verifies they match current behavior | ✅ `2ca15a4` |

**Checkpoint A-0**: ✅ Fixtures captured and verified. No code changes yet.

#### Phase A: Extract validators (no behavior change) ✅ COMPLETE

Each step extracts one focused validator from the monolithic functions.
The old function calls the new validator internally (behavioral equivalence).
Each step is a **separate PR** (not just a commit), especially A.5/A.6
which have the most subtle behavioral differences.

| Step | Extract | From | New file | Test | Status |
|------|---------|------|----------|------|--------|
| A.1 | `validate_variable_types()` | `check_model` + `check_compare_model` | `validate_model.R` | Unit tests for all variable type combos | ✅ `7108bf6` |
| A.2 | `validate_baseline_category()` | `check_model` + `check_compare_model` | `validate_model.R` | Unit tests for BC validation | ✅ `985da03` |
| A.3 | `validate_edge_prior()` | `check_model` | `validate_model.R` | Unit tests for prior setup | ✅ `fc8f0e3` |
| A.4 | `validate_difference_prior()` | `check_compare_model` | `validate_model.R` | Unit tests for diff prior | ✅ `0ab58b8` |
| A.5 | `validate_missing_data()` | `reformat_data` + inline GGM | `validate_data.R` | Unit tests for listwise/impute + GGM constraint | ✅ `d6df27b` |
| A.6 | `reformat_ordinal_data()` | `reformat_data` | `validate_data.R` | Unit tests for recoding (single-group only) | ✅ `6d98878` |
| A.6b | `collapse_categories_across_groups()` | `compare_reformat_data` | `validate_data.R` | Unit tests: multi-group collapsing, missing categories, >2 groups | ✅ `0252731` |
| A.7 | `validate_sampler()` | inline in `bgm.R` + `bgmCompare.R` | `validate_sampler.R` | Unit tests for sampler args + GGM constraints | ✅ `e3c0c13` |
| A.8 | `compute_scaling_factors()` | inline in `bgm.R` + `bgmCompare.R` | `compute_utils.R` | Unit tests for scaling | ✅ `ed33f34` |

> **Note (A.6b):** `collapse_categories_across_groups()` is extracted
> separately from `reformat_ordinal_data()` because the group-conditional
> collapsing logic in `compare_reformat_data()` (via the `observed_scores`
> matrix) has no counterpart in the single-model path. This was identified
> as the highest-risk extraction point by both reviewers.

**Checkpoint A**: ✅ All 2,285 tests pass (497 new unit tests added).
All golden-snapshot fixtures still match. Old functions now delegate to
new validators. No new public API.

New files created during Phase A:
- `R/validate_model.R` — `validate_variable_types()`, `validate_baseline_category()`, `validate_edge_prior()`, `validate_difference_prior()`
- `R/validate_data.R` — `validate_missing_data()`, `reformat_ordinal_data()`, `collapse_categories_across_groups()`
- `R/validate_sampler.R` — `validate_sampler()`
- `R/compute_utils.R` — `compute_scaling_factors()`

#### Phase B: Build bgm_spec constructor

| Step | What |
|------|------|
| B.1 | Create `new_bgm_spec()` low-level constructor with type/presence assertions |
| B.2 | Create `validate_bgm_spec()` with cross-field invariant checks |
| B.3 | Create `bgm_spec()` that calls all Phase A validators → `new_bgm_spec()` → `validate_bgm_spec()` |
| B.4 | Add `print.bgm_spec()` for debugging |
| B.5 | Write comprehensive unit tests for `bgm_spec()` construction |
| B.6 | Add `$arguments` regression test: compare `build_arguments(spec)` against current `output$arguments` field-by-field for all model types |
| B.7 | Add `bgm_spec()` call inside `bgm()` alongside existing code, assert equivalence |
| B.8 | Add `bgm_spec()` call inside `bgmCompare()` alongside existing code, assert equivalence |

**Checkpoint B**: `bgm_spec()` produces correct specs for all model types.
`build_arguments(spec)` matches `output$arguments` for all model types.
Old code still runs in parallel. All tests pass. All golden-snapshot
fixtures still match.

#### Phase C: Wire spec to backends

| Step | What |
|------|------|
| C.1 | Create `run_sampler()` dispatch that reads from spec → calls C++ |
| C.2 | Create `build_output_bgm()` (unified GGM + OMRF) and `build_output_compare()` + shared helpers |
| C.3 | Create thin `build_output()` dispatcher |
| C.4 | In `bgm()`: replace inline GGM/OMRF paths with `run_sampler(spec)` + `build_output(spec, raw)` |
| C.5 | In `bgmCompare()`: replace inline path with `run_sampler(spec)` + `build_output(spec, raw)` |
| C.6 | Add simulate/predict regression tests (call `simulate()` and `predict()` on refactored fit objects, compare to golden fixtures) |
| C.7 | **Gate check:** verify §6.1c exit criteria are met before proceeding |
| C.8 | Delete dead code from old functions |

**Checkpoint C**: Entry points are thin wrappers (~100–120 lines each).
All tests pass (including simulate/predict regression). Old monolithic
functions deleted.

#### Phase D: Clean up residual files

| Step | What |
|------|------|
| D.1 | Move remaining helpers from `function_input_utils.R` to appropriate new files |
| D.2 | Move remaining helpers from `data_utils.R` (compute_*) to `compute_utils.R` |
| D.3 | Consolidate `output_utils.R` → `build_output.R` (keep `transform_*` in place) |
| D.4 | Delete empty or near-empty old files |
| D.5 | Update NAMESPACE, roxygen |

**Checkpoint D**: Clean file structure. All tests pass.

### 6.3 Commit / PR strategy

> **Updated after Review 1 & 2:** The original estimate of ~15 commits is
> optimistic by ~2× once proper fixtures and tests are included.

- **Phase A-0**: 1 PR (fixtures).
- **Phase A**: 9 PRs (one per validator extraction, A.6 and A.6b are
  separate). Each PR gets its own review cycle, especially A.5/A.6/A.6b.
- **Phase B**: 2–3 PRs (spec constructor + integration).
- **Phase C**: 2–3 PRs (wiring + cleanup).
- **Phase D**: 1 PR (housekeeping).
- **Total: ~25–30 commits across 8–10 PRs.**

Every commit:
- Passes `R CMD check`
- Passes all testthat tests
- Passes all golden-snapshot fixture tests
- Has a descriptive message referencing this plan

---

## 7. How Future Extensions Plug In

### 7.1 Mixed MRF (Phase 5)

```r
# validate_variable_types gains:
validate_variable_types <- function(variable_type, num_variables,
                                     allow_continuous = TRUE,
                                     allow_mixed = FALSE) {
  # allow_mixed = TRUE permits continuous + ordinal in one model
}

# bgm_spec gains model_type = "mixed":
# The spec carries both ordinal and continuous variable indices
# run_sampler dispatches to sample_mixed_mrf()
```

No duplication needed — just extending the validator and adding a new
sampler dispatch.

### 7.2 GGM-Compare

```r
# bgm_spec(model_type = "compare") already supports continuous:
# validate_variable_types(allow_continuous = TRUE)  ← flip the flag
# run_sampler dispatches to the GGM-compare C++ backend
```

### 7.3 New variable types (count, binary-specific)

```r
# validate_variable_types gains new choices
# reformat_ordinal_data handles the new type's recoding
# validate_baseline_category handles any type-specific constraints
# The spec carries the new type info; run_sampler passes it to C++
```

### 7.4 New designs (longitudinal, moderated)

```r
# bgm_spec gains new fields in $data (e.g., time index, moderator)
# New validators: validate_longitudinal_design(), validate_moderator()
# run_sampler dispatches to the new C++ backend
# build_output handles the new output format
```

All extensions are **additive** — existing validators and spec fields are
unchanged. This is the core benefit of the spec pattern.

---

## 8. Testing Strategy

### 8.1 Validator unit tests

Each validator gets dedicated unit tests in `tests/testthat/test-validators.R`:

- Happy path: correct inputs produce expected outputs
- Edge cases: single variable, maximum categories, all BC, all ordinal
- Error paths: each validation error is triggered and message is checked
- Equivalence: validator output matches old function output for fixture inputs

### 8.2 bgm_spec integration tests

In `tests/testthat/test-bgm-spec.R`:

- Spec construction for each model type (GGM, OMRF, BC, Compare)
- Round-trip: `bgm_spec() → build_arguments()` reproduces existing `$arguments`
- Spec fields are complete and correctly typed

### 8.3 Regression tests

The existing 1,788 tests serve as regression tests. At every checkpoint,
they must all pass unchanged.

### 8.4 Fixture-based equivalence

During the dual-path period (Phase B), the old and new pipelines run in
parallel with assertions:

```r
# Inside bgm() during Phase B:
spec <- bgm_spec(x, ...)
# Assert spec$data$x matches reformat_data() output
# Assert spec$prior$inclusion_probability matches check_model() output
# etc.
```

These assertions are removed in Phase C when the old path is deleted.

---

## 9. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Subtle behavioral difference in refactored validators | Medium | High | Golden-snapshot fixtures + dual-path assertions (Phase B) + §6.1c exit criteria |
| Group-aware collapsing regression in bgmCompare | Medium | High | Separate `collapse_categories_across_groups()` extraction (A.6b) with dedicated tests for missing categories, >2 groups |
| easybgm breakage from changed `$arguments` | Low | High | Explicit `no_variables` compat field; `$arguments` regression test (B.6); test with easybgm |
| simulate/predict breakage from changed `$arguments` layout | Medium | High | simulate/predict regression tests in Phase C (C.6) |
| Large merge conflicts with other work | Low | Medium | Do this on `ggm_mixed` before Phase 5 starts |
| Scope creep (rewriting simulate_predict.R) | Medium | Medium | Explicitly out of scope — simulate/predict stay as-is |

---

## 10. What's Explicitly Out of Scope

- **simulate_predict.R** — stays as-is. It reads from `object$arguments`
  which is unchanged. Future work could have it accept a spec.
- **extractor_functions.R** — unchanged. Reads from fit objects.
- **mcmc_summary.R** — unchanged. Internal MCMC processing.
- **S3 methods** — unchanged. Read from fit objects.
- **C++ layer** — no changes. The spec just organizes the same arguments.
- **easybgm compatibility shim** — stays in `bgm.R` until easybgm ≥ 0.2.2.

---

## 11. Definition of Done

- [ ] Golden-snapshot fixtures captured for 15+ representative cases
- [ ] All validators extracted and unit-tested
- [ ] `collapse_categories_across_groups()` extracted with dedicated tests
- [ ] `new_bgm_spec()` constructor enforces field types/presence
- [ ] `validate_bgm_spec()` checks cross-field invariants
- [ ] `bgm_spec()` constructor builds correct specs for GGM, OMRF, BC, Compare
- [ ] `build_output_bgm()` (unified GGM + OMRF) and `build_output_compare()` each tested independently
- [ ] `build_arguments(spec)` matches `output$arguments` field-by-field for all model types
- [ ] `bgm()` and `bgmCompare()` are thin wrappers (~100–120 lines)
- [ ] All 1,788+ existing tests pass unchanged
- [ ] All golden-snapshot fixtures match between old and new pipelines
- [ ] simulate/predict regression tests pass for all model types
- [ ] New validator tests provide >95% coverage of validation logic
- [ ] Old monolithic functions deleted
- [ ] File structure matches Section 3.3
- [ ] `R CMD check` passes with 0 errors, 0 warnings, 0 notes
- [ ] Dev docs updated (roadmap.md references this plan)

---

## 12. Open Discussion Points

> **Added 2026-02-24.** These items had partial or conflicting guidance
> across the two reviews and need a decision before implementation begins.

### 12.1 S3 class for bgm_spec: now or later?

**Review 1** suggests skipping the S3 class initially ("just use a
validated named list returned by `build_spec()`; the `print.bgm_spec()`
method is yak-shaving"). **Review 2** is silent on this. The current
plan uses an S3 class.

**Options:**
- **(a)** Ship with S3 class from the start (current plan). Pro: method
  dispatch is available immediately if needed; consistent with brms
  pattern. Con: print method is low-value; adds a few lines.
- **(b)** Start as a validated named list; add `class()` later when method
  dispatch is needed. Pro: less yak-shaving; forces us to keep the spec
  as a pure data structure. Con: harder to add retroactively if code
  starts type-checking.

**Decision:** Use a plain validated list (option b). The S3 class adds
no practical value — nobody prints a spec in production, and
`validate_bgm_spec()` works fine on a plain list. Add the class tag
later only if method dispatch on the spec is needed.

### 12.2 Phase A granularity: PRs vs. commits

**Review 1** recommends each Phase A step be a separate PR, not just a
commit. **Review 2** says 8–10 PRs total (which implies grouping some
Phase A steps).

**Decision:** Decide on the go — adjust PR granularity based on review
capacity and how complex each extraction turns out to be.

### 12.3 Where do GGM constraints live?

**Review 2** noted that GGM currently forces `update_method =
"adaptive-metropolis"` and forbids `na_action = "impute"`, and the plan
didn't specify where these constraints land.

**Current decision:** Split across two validators:
- `validate_missing_data()` enforces GGM + impute → error (§4.4)
- `validate_sampler()` enforces GGM → adaptive-metropolis (§4.8)
- `validate_bgm_spec()` asserts both as cross-field invariants (§3.2b)

**Decision:** Keep double enforcement. The individual validators produce
user-facing error messages; the spec validator is a cheap safety net for
future code paths that might bypass a validator. Revisit if it becomes
a maintenance annoyance.

### 12.4 Hidden mutation in `$arguments`

**Review 1** (§6) identified that several fields mutate between the
entry point and `prepare_output_*()`: `baseline_category` gets replicated
and shifted, `inclusion_probability` starts as scalar and becomes a
matrix, `variable_type` starts as a single string and becomes a vector.

`build_arguments()` will read from the spec, which stores post-validation
(post-mutation) values. This is *probably* correct, but there's a risk
that some downstream code depends on the pre-mutation value.

**Current decision:** The `$arguments` regression test (B.6) will catch
any mismatch. No special handling is planned unless the regression test
fails. Document the mutation flow in `build_arguments()` so future
maintainers understand which spec field maps to which arguments field.
