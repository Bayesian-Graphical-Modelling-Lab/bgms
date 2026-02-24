# Review 1: R Scaffolding Refactor Plan

**Reviewer:** Senior R package developer (automated)
**Date:** 2025-02-24
**Documents reviewed:** `dev/scaffolding/plan.md`, `dev/scaffolding/migration_mapping.md`
**Source files reviewed:** `R/bgm.R`, `R/bgmCompare.R`, `R/function_input_utils.R`,
`R/data_utils.R`, `R/output_utils.R`

---

## Overall Assessment

The diagnosis is correct and well-evidenced. The duplication map in §2.5 is
honest. The `bgm_spec` pattern is the right idea — the brms/rstanarm
comparisons are apt. The phased migration with dual-path assertions (Phase B)
is the single best decision in the plan. The line-level mapping in
`migration_mapping.md` is unusually thorough for a plan document.

The main risks are not in the *direction* of the refactor but in
under-specified integration points, missing regression infrastructure, and
one design choice (`build_output()` unification) that is likely to trade
independent duplication for coupled complexity.

---

## Weakest Parts (ranked by risk)

### 1. The `bgm_spec` object is a plain list with no constructor validation

The spec has ~40 fields across 6 nested sub-lists. It's defined only in prose
(§3.2). There is no `new_bgm_spec()` constructor, no `validate_bgm_spec()`,
no assertion that required fields are present or correctly typed. Once
`run_sampler()` and `build_output()` depend on
`spec$prior$pairwise_scaling_factors` existing, a single missing field
produces a cryptic `NULL` subscript error deep in the pipeline.

brms avoids this because `brmsterms` and `brmsframe` are validated S3 classes
with strict constructors.

**Fix:** Add a `new_bgm_spec()` constructor that asserts types and lengths for
every field, and a `validate_bgm_spec()` that checks cross-field invariants
(e.g., `model_type == "ggm"` implies `is_continuous == TRUE`). Call it at the
end of `bgm_spec()`. This is 60–80 lines and catches 80% of integration bugs.

---

### 2. Group-aware category collapsing isn't mapped to `reformat_ordinal_data()`

`compare_reformat_data()` (data_utils.R L290–430) does something
fundamentally different from `reformat_data()`: it collapses ordinal
categories that aren't observed in *all* groups. The `observed_scores` matrix
and the per-group cross-check have no counterpart in the proposed
`reformat_ordinal_data()` spec (§4.5). The migration mapping says "deleted —
uses `reformat_ordinal_data()`" but doesn't account for this group-conditional
logic. This is the most likely place for a silent behavioral difference to
slip through.

**Fix:** The spec for `reformat_ordinal_data()` needs an explicit `group`
parameter and must document the group-conditional collapsing as a distinct code
path. Alternatively, factor it into a separate
`collapse_categories_across_groups()` that wraps `reformat_ordinal_data()`.

---

### 3. The dual-path assertion period has no defined gate criteria

The plan says "assertions are removed in Phase C when the old path is deleted"
but doesn't define what "passed" means. How many CI runs? Which test matrices?
Manual inspection or automated? In practice, these dual-path periods either
get removed too early (before edge cases are hit) or never get removed.

**Fix:** Define a concrete gate: "Phase B assertions stay for at least 2 full
CI cycles covering all 3 model types × 2 na_actions × edge_selection on/off.
A dedicated test fixture (stored in `dev/fixtures/`) captures the full
`bgm_spec` output for each combination. Phase C can only begin when all
fixtures match to floating-point tolerance."

---

### 4. `build_output()` unification is underspecified and likely counterproductive

The plan claims a single `build_output()` replaces three functions
(`prepare_output_bgm`, `prepare_output_ggm`, `prepare_output_bgmCompare`).
But these three differ substantially:

- **GGM:** diagonal params are "precision", naming is `"Var (precision)"`,
  main is a 1-column matrix.
- **OMRF:** main params expand by `num_categories`, naming is `"Var (cat)"`,
  main is a wide matrix.
- **Compare:** has baseline/difference split, indicator logic includes diag
  (main effects).

A "unified" function that handles all three via
`if (spec$model_type == ...)` branching is just the current three functions
stapled together, except now they share a namespace and can accidentally
interfere. This trades duplicated-but-independent code for
coupled-and-branchy code — often worse.

**Fix:** Keep three internal builders (`build_output_ggm`,
`build_output_omrf`, `build_output_compare`) that share *helper* functions
(e.g., `build_posterior_mean_matrix()`, `build_raw_samples_list()`). Let
`build_output()` be a thin dispatcher. This is what rstanarm does —
`stan_glm.fit`, `stan_glmer.fit`, etc. share helpers but have separate
top-level builders.

---

### 5. The existing test suite is treated as coverage — it isn't

The plan's safety net is "all 1,788 tests pass." But these tests almost
certainly exercise `bgm()` and `bgmCompare()` through the public API, meaning
they test the *composition* of validation + sampling + output. They won't
catch cases where a refactored validator produces the same final output but
with different intermediate state (e.g., a spec field gets the right value by
coincidence but the validator logic is wrong).

The plan acknowledges this ("each validator gets unit tests") but lists it as
Phase A-parallel, not as a prerequisite gating Phase B.

**Fix:** Write validator unit tests *before* extracting them. Characterize the
current behavior first (golden-snapshot tests against `check_model()` /
`reformat_data()` outputs for 10+ edge cases), then extract, then verify the
new validators match the snapshots.

---

### 6. `$arguments` backward-compatibility has a hidden mutation problem

`build_arguments()` maps from `bgm_spec` → `$arguments`. But the current code
doesn't build `$arguments` from a spec — it builds it by forwarding function
arguments directly. Some get mutated between the entry point and
`prepare_output_*()`:

- `baseline_category` gets replicated and shifted in `reformat_data()`
- `inclusion_probability` starts as scalar `0.5`, becomes a matrix in
  `check_model()`
- `variable_type` starts as `"ordinal"`, the arguments list stores the
  replicated vector

If `build_arguments()` reads from the spec, it gets the post-mutation values
(probably correct). But there's no test that verifies the arguments list
matches field-by-field between old and new pipelines.

**Fix:** Add an explicit regression test: run the current code, capture
`output$arguments`, run the new code, capture `output$arguments`,
`testthat::expect_equal()` on every field. Do this for all model types.

---

### 7. `simulate_predict.R` is out of scope but not out of blast radius

The plan marks `simulate_predict.R` as "out of scope" but `simulate.bgms()`
calls `extract_arguments()` which reconstructs model metadata from
`$arguments`. If the `$arguments` list changes field names, types, or adds
unexpected `NULL`s, simulation breaks. The plan doesn't include a
simulate/predict integration test as a Phase C gate.

**Fix:** Add `simulate()` and `predict()` calls to the Phase C regression
suite, not just `bgm()` output comparison.

---

## What I'd Do Differently

1. **Start with golden-snapshot fixtures, not validators.** Capture the full
   `check_model()`, `reformat_data()`, and `prepare_output_*()` outputs for
   ~15 representative inputs (all model types × edge cases). Store as `.rds`
   files in `dev/fixtures/`. Every subsequent change is validated against
   these snapshots. This costs a day but prevents a month of debugging.

2. **Don't unify `build_output()`.** Three dispatched builders sharing helpers
   is cleaner than a 300-line unified function with model-type switches. The
   duplication in output building is *incidental*, not *essential* — the three
   model types genuinely have different output structures.

3. **Skip the `bgm_spec` S3 class initially.** Just use a validated named
   list returned by `build_spec()`. The `print.bgm_spec()` method is nice but
   it's yak-shaving — nobody prints a spec in production. Add the class later
   when you need method dispatch on it.

4. **Phase A should be 8 PRs, not 8 commits.** Each validator extraction
   deserves its own review cycle, especially `validate_missing_data()` and
   `reformat_ordinal_data()` which have the most subtle behavioral
   differences between the `bgm` and `bgmCompare` paths.

5. **Add a `compare_reformat_data()` behavioral-equivalence test
   immediately.** This function has the most complex group-aware logic and is
   the most likely place for regressions. Before touching it, write 5 tests
   that exercise: all-ordinal multi-group, mixed ordinal/BC multi-group,
   categories missing in one group, >2 groups, and missing data + group
   collapsing together.

6. **The "~15 commits" estimate is optimistic by 2×.** Realistic count with
   proper testing: ~25–30 commits. Budget accordingly.

---

## Summary

| Concern | Severity | Effort to fix |
|---------|----------|---------------|
| No spec validation / constructor | High | 1 day |
| Group-aware collapsing not mapped | High | 0.5 day (design) |
| `build_output()` unification will couple code | Medium | 0 (just don't do it) |
| No snapshot fixtures | Medium | 1 day |
| Dual-path gate criteria missing | Medium | 1 hour (write it down) |
| `$arguments` regression test missing | Medium | 0.5 day |
| simulate/predict not in regression suite | Medium | 0.5 day |
| Commit-vs-PR granularity | Low | 0 (process decision) |
