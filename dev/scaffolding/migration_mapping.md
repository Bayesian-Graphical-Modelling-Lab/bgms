# Migration Mapping: Current Code → New Location

**Companion to:** `dev/scaffolding/plan.md`  
**Updated:** 2026-02-24 (after Review 1 & Review 2)

This file maps every block of current code to its destination in the
refactored architecture. Use it as a checklist during implementation.

> **Review changelog (2026-02-24):**
> - Added Phase A.0 (golden-snapshot fixtures, prerequisite)
> - Split A.6 (single-group recoding) from A.6b (group-aware collapsing)
> - `compare_reformat_data()` group-collapsing → `collapse_categories_across_groups()`
> - `build_output()` is now a thin dispatcher; three model-specific builders
> - Added `new_bgm_spec()` / `validate_bgm_spec()` in Phase B
> - GGM constraint placement documented in A.5 and A.7

---

## function_input_utils.R (671 lines)

| Lines | Current code | Destination | Phase |
|-------|-------------|-------------|-------|
| 1–6 | `check_positive_integer()` | keep in `function_input_utils.R` (or move to `validate_sampler.R`) | D.1 |
| 8–13 | `check_non_negative_integer()` | keep in place | D.1 |
| 16–23 | `check_logical()` | keep in place | D.1 |
| 26–34 | `check_seed()` | `validate_sampler.R` → `validate_seed()` | A.7 |
| 37–55 | `check_model()` signature + variable_type start | `validate_model.R` → `validate_variable_types()` | A.1 |
| 57–140 | `check_model()` variable type parsing | `validate_model.R` → `validate_variable_types()` | A.1 |
| 142–233 | `check_model()` baseline_category validation | `validate_model.R` → `validate_baseline_category()` | A.2 |
| 236–248 | `check_model()` prior scale checks | `validate_model.R` → `validate_prior_scales()` | A.3 |
| 250–370 | `check_model()` edge selection prior logic | `validate_model.R` → `validate_edge_prior()` | A.3 |
| 372–376 | `check_model()` return statement | absorbed into `bgm_spec()` | B.1 |
| 378–395 | `check_compare_model()` signature + group validation | `validate_model.R` → `validate_groups()` + `validate_variable_types()` | A.1 |
| 395–470 | `check_compare_model()` variable type parsing | `validate_model.R` → `validate_variable_types()` (shared) | A.1 |
| 472–560 | `check_compare_model()` baseline checks | `validate_model.R` → `validate_baseline_category()` (shared) | A.2 |
| 562–580 | `check_compare_model()` prior scale checks | `validate_model.R` → `validate_prior_scales()` (shared) | A.3 |
| 582–670 | `check_compare_model()` difference prior logic | `validate_model.R` → `validate_difference_prior()` | A.4 |
| 660–671 | `progress_type_from_display_progress()` | `validate_sampler.R` → `validate_progress()` | A.7 |

## data_utils.R (494 lines)

| Lines | Current code | Destination | Phase |
|-------|-------------|-------------|-------|
| 1–72 | `reformat_data()` missing-data handling | `validate_data.R` → `validate_missing_data()` | A.5 |
| 74–193 | `reformat_data()` per-variable recoding | `validate_data.R` → `reformat_ordinal_data()` | A.6 |
| 195–200 | `reformat_data()` return | absorbed into `bgm_spec()` | B.1 |
| 202–214 | `data_check()` | `validate_data.R` → `validate_data()` (keep alias) | A.5 |
| 216–290 | `compare_reformat_data()` missing-data | **deleted** — uses `validate_missing_data()` | A.5 |
| 290–340 | `compare_reformat_data()` per-variable recoding | **deleted** — uses `reformat_ordinal_data()` (shared, single-group) | A.6 |
| 340–430 | `compare_reformat_data()` group-conditional collapsing (`observed_scores` matrix, cross-group category check) | `validate_data.R` → `collapse_categories_across_groups()` **[new]** | A.6b |
| 432–450 | `compute_counts_per_category()` | `compute_utils.R` (move) | D.2 |
| 452–475 | `compute_blume_capel_stats()` | `compute_utils.R` (move) | D.2 |
| 477–495 | `compute_pairwise_stats()` | `compute_utils.R` (move) | D.2 |

## bgm.R (882 lines)

| Lines | Current code | Destination | Phase |
|-------|-------------|-------------|-------|
| 1–500 | roxygen docs | stays (trimmed) | C.3 |
| 501–510 | function signature | stays (unchanged) | — |
| 511–520 | `data_check(x)` | `bgm_spec()` calls `validate_data()` | C.3 |
| 521–530 | `check_model(...)` | `bgm_spec()` calls validators | C.3 |
| 531–545 | deprecated arg handling | `handle_deprecated_bgm_args()` helper or stays inline | C.3 |
| 546–590 | sampler validation (warmup, HMC, NUTS) | `bgm_spec()` calls `validate_sampler()` | C.3 |
| 591–610 | na_action, progress, seed | `bgm_spec()` calls validators | C.3 |
| 612–660 | **GGM path**: missing data + sample_ggm | `run_sampler(spec)` → `run_sampler_ggm()` | C.1 |
| 660–680 | GGM: transform + output | `build_output(spec, raw)` → dispatches to `build_output_ggm()` | C.2/C.3 |
| 682–698 | OMRF: `reformat_data()` | `bgm_spec()` calls `validate_missing_data()` + `reformat_ordinal_data()` | C.3 |
| 698–745 | OMRF: scaling factors | `bgm_spec()` calls `compute_scaling_factors()` | C.3 |
| 745–795 | OMRF: `sample_omrf()` | `run_sampler(spec)` → `run_sampler_omrf()` | C.1 |
| 795–810 | OMRF: transform + output | `build_output(spec, raw)` → dispatches to `build_output_omrf()` | C.2/C.3 |
| 810–845 | OMRF: NUTS diagnostics | stays (appended after `build_output()`) | C.3 |
| 845–882 | easybgm compat shim | stays | — |

## bgmCompare.R (695 lines)

| Lines | Current code | Destination | Phase |
|-------|-------------|-------------|-------|
| 1–195 | roxygen + signature | stays (trimmed) | C.4 |
| 246–345 | deprecated arg handling (~12 args) | `handle_deprecated_compare_args()` helper | C.4 |
| 346–375 | sampler validation | `bgm_spec()` calls `validate_sampler()` | C.4 |
| 376–400 | `check_compare_model()` | `bgm_spec()` calls validators | C.4 |
| 400–430 | `compare_reformat_data()` | `bgm_spec()` calls validators | C.4 |
| 430–510 | compute indices, scaling factors | `bgm_spec()` calls `compute_indices()`, `compute_scaling_factors()` | C.4 |
| 510–540 | compute projection, sort groups | `bgm_spec()` calls `compute_projection()` | C.4 |
| 540–580 | `run_bgmCompare_parallel()` | `run_sampler(spec)` → `run_sampler_compare()` | C.1 |
| 580–695 | `prepare_output_bgmCompare()` | `build_output(spec, raw)` → dispatches to `build_output_compare()` | C.2/C.3 |

## output_utils.R (623 lines)

| Lines | Current code | Destination | Phase |
|-------|-------------|-------------|-------|
| 1–175 | `prepare_output_bgm()` | `build_output.R` → `build_output_omrf()` | C.2 |
| 177–315 | `prepare_output_ggm()` | `build_output.R` → `build_output_ggm()` | C.2 |
| 317–370 | `transform_ggm_backend_output()` | stays in `output_utils.R` (or `build_output.R`) | C.2 |
| 372–408 | `transform_new_backend_output()` | stays | C.2 |
| 410–500 | `generate_param_names_bgmCompare()` | `build_output.R` → shared helper `generate_param_names()` | C.2 |
| 500–623 | `prepare_output_bgmCompare()` | `build_output.R` → `build_output_compare()` | C.2 |

---

## Lines deleted vs lines written (estimate)

| Category | Deleted | Written | Net |
|----------|---------|---------|-----|
| Duplicated validation logic | ~600 | ~300 | -300 |
| Duplicated data reformatting | ~250 | ~150 | -100 |
| Inline code in bgm/bgmCompare | ~300 | ~100 | -200 |
| bgm_spec constructor + new_bgm_spec + validate_bgm_spec | 0 | ~200 | +200 |
| collapse_categories_across_groups() | 0 | ~60 | +60 |
| run_sampler dispatch | 0 | ~80 | +80 |
| build_output (3 dispatched builders + shared helpers) | ~500 | ~350 | -150 |
| Golden-snapshot fixtures + generation script | 0 | ~200 | +200 |
| New unit tests (validators + regression) | 0 | ~500 | +500 |
| **Total** | **~1,650** | **~1,640** | **-10 + 700 test/fixture lines** |

Net effect: roughly neutral production code count, ~700 new lines of tests
and fixtures, and the remaining code is organized by concern rather than
by entry point. The three output builders are independent (not coupled via
branching), and the spec object is validated at construction time.
