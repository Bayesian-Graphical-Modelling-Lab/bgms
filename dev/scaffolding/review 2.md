# Review 2: R Scaffolding Refactor Plan

**Reviewer:** GitHub Copilot (GPT-5.1-Codex)
**Date:** 2026-02-24
**Scope:** dev/scaffolding/plan.md, dev/scaffolding/review 1.md, R/bgm.R, R/bgmCompare.R, R/function_input_utils.R, R/data_utils.R, R/output_utils.R

---

## Executive Summary

The bgm_spec direction is still solid, but the current plan underestimates how tricky the existing validation and preprocessing paths are. The riskiest areas are (1) lack of a strict spec constructor/validator, (2) losing group-aware data collapsing logic when extracting validators, and (3) unifying the three output builders into one large function. Without golden fixtures and explicit regression gates, the dual-path strategy will not deliver the promised safety net.

---

## High-Risk Gaps

1. **Spec construction is informal.** The plan defines bgm_spec as a plain list. No constructor enforces presence or type of the ~40 nested fields. Downstream code will fail with opaque `NULL` subscripts. Add `new_bgm_spec()` plus `validate_bgm_spec()` that asserts invariants such as "continuous ⇒ adaptive-metropolis" and "model_type == 'compare' ⇒ data$group is present" before anything reaches `run_sampler()`.

2. **Group-aware ordinal collapsing is unmapped.** `compare_reformat_data()` collapses categories that disappear in any group via the `observed_scores` matrix. The proposed `reformat_ordinal_data()` helper has no notion of groups, so refactoring will silently change sufficient statistics for bgmCompare. Either split this into a dedicated `collapse_categories_across_groups()` or ensure `reformat_ordinal_data()` accepts group metadata.

3. **Unified `build_output()` couples independent logic.** `prepare_output_bgm`, `prepare_output_ggm`, and `prepare_output_bgmCompare` have genuinely different naming, summaries, and SBM extras. Folding them into one dispatcher creates a 300-line branch-fest where a tweak for GGM can break bgmCompare. Keep three model-specific builders and share only naming/post-processing helpers.

4. **No fixture-based regression tests.** The plan leans on "1,788 tests" but none snapshot the outputs of `check_model()`, `reformat_data()`, or `prepare_output_*()`. Before extraction, capture `.rds` fixtures for representative inputs (continuous, ordinal, Blume–Capel, multi-group with missing categories, SBM). Every new validator must match those snapshots; otherwise subtle behavior shifts will go unnoticed.

5. **Dual-path exit criteria undefined.** Phase B says "remove assertions in Phase C" without specifying success gates. Require at least two full CI runs where bgm/bgmCompare outputs match the fixtures for all combinations (GGM/OMRF/compare × listwise/impute × selection on/off) before deleting the old path.

6. **$arguments compatibility untested.** `prepare_output_*()` currently builds `$arguments` after all recoding/scaling, and downstream code (simulate/predict, easybgm) depends on the exact fields. Introduce a regression test that fits each model type, saves `output$arguments`, and compares to the refactored pipeline field-for-field.

7. **simulate/predict blast radius ignored.** Marking `simulate_predict.R` as "out of scope" does not protect it. Those S3 methods read `$arguments`; if any layout changes slip through, simulation breaks silently. Add simulate/predict calls to the regression suite before Phase C.

8. **Continuous path guardrails unspecified.** In `bgm()` today, GGM forces `update_method = "adaptive-metropolis"` and forbids `na_action = "impute"`. The plan never states where that logic lands post-refactor. Make sure `validate_variable_types()` or `validate_missing_data()` re-imposes those constraints.

---

## Recommendations

1. **Start with golden fixtures, not code movement.** Snapshot current validator/output behavior first. This gives you mechanical oracles for dual-path comparisons.
2. **Introduce `new_bgm_spec()` early.** Fail fast when a field is missing or of the wrong type; add cross-field assertions.
3. **Factor data validators carefully.** Separate per-variable recoding from group-aware collapsing so bgmCompare keeps its semantics.
4. **Leave top-level output builders split.** Share helper utilities but keep `build_output_ggm`, `build_output_omrf`, and `build_output_compare` distinct.
5. **Track `$arguments` explicitly.** Add tests that compare the current `$arguments` list against the spec-derived version for all model types, plus simulate/predict smoke tests.
6. **Document Phase B gates.** Only remove dual-path assertions after fixtures for every model/NA combination match across multiple CI runs.

---

## Effort Adjustments

- Expect ~25–30 commits (or 8–10 reviewable PRs) once proper fixtures and tests are included, not the planned 15.
- Budget at least one full day for fixture creation and another for building the spec constructor/validator.

The refactor direction remains correct, but without these guardrails the probability of introducing silent behavior changes is high. Tighten validation, preserve model-specific output paths, and insist on fixture-backed regression before deleting the existing pipeline.
