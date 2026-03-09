# Phase L — Model Output R Code Inspection

**Date:** 2026-03-09
**Branch:** `ggm_mixed`
**Scope:** Audit of R code that assembles, summarises, and exposes
mixed MRF model output. Identifies bugs and inconsistencies.

---

## L1. `build_output_mixed_mrf()` correctness

### L1.1 Flat-to-matrix index mapping

The C++ `get_storage_vectorized_parameters()` emits parameters in this
order:

```
[main_discrete | Kxx_ut | main_continuous | Kxy_rowmajor | Kyy_ut+diag]
```

The R function `compute_mixed_parameter_indices()` splits this into:

- `main_idx`: discrete thresholds + continuous means + Kyy diagonal
- `pairwise_idx`: Kxx upper-tri + Kyy off-diagonal + Kxy

The Kyy block in C++ is stored row-major upper-triangular including
diagonal (i.e., for `j >= i`). The R code iterates `for(i in 1:q)
for(j in i:q)` and separates diagonal from off-diagonal positions.

**Finding: correct.** The R loop matches the C++ `for(i) for(j=i..q)`
ordering exactly. Diagonal and off-diagonal positions are correctly
identified and routed to `main_idx` and `pairwise_idx` respectively.

### L1.2 Indicator ordering vs pairwise ordering

C++ indicator vector: `[Gxx_ut | Gyy_ut | Gxy_rowmajor]`
R `pairwise_idx` order: `[Kxx_ut | Kyy_offdiag | Kxy]`

**Finding: correct.** The indicator vector excludes diagonal entries
(Gyy has no diagonal — self-edges are not indicators). The pairwise
vector also excludes diagonal (Kyy diagonal is routed to main). Both
use the same traversal order within each block. The edge names
generated in `build_output_mixed_mrf` follow the same block order
(discrete-discrete, continuous-continuous, cross), so
`indicator_summary` rows align with `pairwise_summary` rows.

### L1.3 NUTS diagnostic column naming

`build_output_mixed_mrf` checks `s$update_method == "hybrid-nuts"`.
`build_output_bgm` checks `s$update_method == "nuts"`.
Both call the same `summarize_nuts_diagnostics()` function.

**Finding: correct.** The `bgm_spec` pipeline remaps `"nuts"` to
`"hybrid-nuts"` for mixed MRF models (line 350–351 of `bgm_spec.R`),
so the string comparison works correctly in both builders.

### L1.4 BC main effects in posterior mean

For Blume-Capel discrete variables, `build_output_mixed_mrf` stores
alpha and beta in columns 1:2 of `pmm_disc`, with remaining columns
NA. This matches the OMRF behavior in `build_output_bgm`.

**Finding: correct.**

### L1.5 Pairwise symmetry and sign

`fill_mixed_symmetric()` places each edge value symmetrically
(`mat[oi,oj] = mat[oj,oi] = value`), mapping from internal
discrete-first indices to original user column indices via
`disc_idx` and `cont_idx`.

**Finding: correct.** Sign convention is inherited from C++: positive
values mean positive association, matching GGM/OMRF convention.

### L1.6 SBM allocation names missing in raw_samples

`build_output_mixed_mrf` calls `build_raw_samples_list(raw,
edge_selection, edge_prior, names_main, edge_names)` without the
`allocation_names` argument. When the SBM prior is active, the
allocations are stored in `$raw_samples$allocations`, but
`$raw_samples$parameter_names$allocations` is `NULL`.

In contrast, `build_output_bgm` computes `alloc_names` and passes it.

**Verdict: BUG (minor).** The allocation names are missing from the
raw samples list for mixed MRF models with SBM prior. Fix: pass
`allocation_names = all_internal_names` (the `c(disc_names,
cont_names)` vector already computed on line 586).

**Priority:** Low — only affects users who directly access
`$raw_samples$parameter_names$allocations` with an SBM prior on
mixed data.

---

## L2. Posterior mean type inconsistency

| Model | `$posterior_mean_main` type |
|-------|----------------------------|
| GGM | numeric matrix (p × 1, column = "mean") |
| OMRF | numeric matrix (p × max_cats) |
| Mixed MRF | **list** with `$discrete` (p × max_cats) and `$continuous` (q × 2) |

**Verdict: INCONSISTENCY (moderate).** This is the most significant
structural difference in the output. Any user code that does
`fit$posterior_mean_main[1, 2]` will fail on a mixed MRF fit because
the list does not support matrix indexing.

The list structure is arguably the correct design — discrete thresholds
and continuous means/precisions are different quantities that belong
in separate matrices. But it breaks the assumption that
`$posterior_mean_main` is always a matrix.

**Recommendation:** Keep the list structure for mixed MRF. Document it
clearly in `?bgm` and in `coef.bgms()`. The alternative — flattening
into one matrix with mixed row semantics — would be worse for users.

---

## L3. `coef.bgms()` with mixed MRF

```r
coef.bgms = function(object, ...) {
  out = list(main = object$posterior_mean_main, ...)
}
```

For mixed MRF, `out$main` is a list (with `$discrete` and
`$continuous`), yielding a nested list: `coef(fit)$main$discrete`.
For GGM/OMRF, `coef(fit)$main` is a matrix.

**Verdict: INCONSISTENCY (moderate).** The return type of
`coef()$main` depends on the model type. This is not documented.

**Recommendation:** Accept this inconsistency but document it:
- For GGM: matrix (p × 1)
- For OMRF: matrix (p × max_cats)
- For mixed MRF: list with `$discrete` (matrix) and `$continuous`
  (matrix)

The `@return` block of `coef.bgms` currently says "Posterior mean of
the category threshold parameters" — incorrect for mixed MRF.

---

## L4. `extract_category_thresholds()` — BROKEN for mixed MRF

```r
extract_category_thresholds.bgms = function(bgms_object) {
  vec = bgms_object$posterior_summary_main[, "mean"]
  num_vars = arguments$num_variables  # = p + q
  num_cats = arguments$num_categories # length p (discrete only)
  ...
  for(v in seq_len(num_vars)) {
    num_cats[v]  # ERROR: index out of bounds when v > p
  }
}
```

**Verdict: BUG (crash).** For mixed MRF, `num_variables = p + q`, but
`num_categories` has length `p`. When the loop reaches `v = p + 1`,
`num_cats[v]` is an out-of-bounds access. The function will error.

Additionally, even if the bounds were fixed, the `vec` from
`posterior_summary_main` for mixed MRF interleaves discrete thresholds
with continuous mean/precision entries. The positional extraction logic
(`pos = pos + k` vs `pos = pos + 2`) would misalign because it does
not know that continuous-variable rows have a different meaning.

**Fix required:** Either:
1. Guard: `if(isTRUE(arguments$is_mixed))` and extract only discrete
   thresholds using `arguments$num_discrete` and `arguments$is_ordinal`.
2. Or: extract from `object$posterior_mean_main$discrete` directly.

---

## L5. `extract_pairwise_interactions()` — WRONG COLUMN NAMES for mixed MRF

```r
extract_pairwise_interactions.bgms = function(bgms_object) {
  mats = bgms_object$raw_samples$pairwise
  mat = do.call(rbind, mats)
  edge_names = character()
  for(i in 1:(num_vars - 1)) {
    for(j in (i + 1):num_vars) {
      edge_names = c(edge_names, paste0(var_names[i], "-", var_names[j]))
    }
  }
  dimnames(mat) = list(..., edge_names)
}
```

**Verdict: BUG (silent wrong labels).** The function generates edge
names in standard upper-triangle order of the (p+q) × (p+q) matrix,
using user column order. But for mixed MRF, `raw_samples$pairwise`
columns are in block order: [Kxx_ut | Kyy_offdiag | Kxy]. These two
orderings are NOT the same when discrete and continuous variables are
interleaved in the original data.

Example: data with columns `[x1, y1, x2, y2]` (disc, cont, disc, cont).
- Standard upper-tri: (x1,y1), (x1,x2), (x1,y2), (y1,x2), (y1,y2), (x2,y2)
- Mixed block order: (x1,x2), (y1,y2), (x1,y1), (x1,y2), (x2,y1), (x2,y2)

The edge count is the same (6), but the assignment is scrambled.

Note: `raw_samples$parameter_names$pairwise` already contains the
correct edge names — the extractor should use those instead of
generating its own.

**Fix required:** Use `bgms_object$raw_samples$parameter_names$pairwise`
for mixed MRF models instead of regenerating edge names.

---

## L6. `extract_indicators()` — SAME BUG as L5

`extract_indicators.bgms` also generates its own column names via
the `parameter_names$indicator` stored in raw_samples — wait, it
actually *does* use `bgms_object$raw_samples$parameter_names$indicator`:

```r
param_names = bgms_object$raw_samples$parameter_names$indicator
colnames(indicator_samples) = param_names
```

**Finding: correct.** Unlike `extract_pairwise_interactions`, the
indicator extractor uses stored names. No bug here.

---

## L7. `simulate.bgms()` with mixed MRF

The `simulate_bgms_mixed` helper:
1. Extracts parameters from `posterior_mean_main$discrete` /
   `$continuous` and `posterior_mean_pairwise`.
2. Calls `sample_mixed_mrf_gibbs()` (C++).
3. `combine_mixed_result()` scatters x/y columns into original order
   using `disc_idx` and `cont_idx`.

**Finding: correct.** The `nsim` parameter is properly forwarded. The
`combine_mixed_result` helper uses `disc_idx` / `cont_idx` from
`arguments`, which are set by `build_arguments_mixed_mrf`. Column
reordering is correct. `drop = FALSE` protects against single-variable
edge cases.

---

## L8. `predict.bgms()` with mixed MRF

The `predict_bgms_mixed` helper:
1. Splits `newdata` by `disc_idx` / `cont_idx`.
2. Maps `predict_vars` to internal 0-based indices via `match()`.
3. Calls C++ `compute_conditional_mixed`.
4. Formats output with variable-type-aware labels.

**Finding: correct.** Both discrete and continuous conditioned
variables are supported. Column matching is positional (same as
GGM/OMRF), which is correct given users provide data in the same
column order as training data.

Minor note: no column-name validation is performed — if `newdata` has
different column names or order from training data, results will be
silently wrong. This is consistent with GGM/OMRF behavior so it is
not a mixed-MRF-specific issue.

---

## L9. NUTS diagnostic wiring

`build_output_mixed_mrf` checks for `"hybrid-nuts"`;
`build_output_bgm` checks for `"nuts"`. Both funnel to the same
`summarize_nuts_diagnostics()`.

The `bgm_spec` pipeline maps `"nuts"` → `"hybrid-nuts"` for mixed
MRF at spec-build time.

**Finding: correct.** No model-type-specific behavior is needed in
the diagnostics function itself.

---

## L10. `summary.bgms()` and `print.summary.bgms()`

`summary.bgms()` blindly returns `$posterior_summary_main` and
`$posterior_summary_pairwise` — this works for all model types.

`print.summary.bgms()` labels the main block as
**"Category thresholds:"** unconditionally.

**Verdict: MISLEADING LABEL (minor).** For mixed MRF, the main
summary includes rows like `"V4 (mean)"` and `"V4 (precision)"` under
the label "Category thresholds:", which is incorrect. For GGM, the
main effects are means (also mislabeled as thresholds).

**Recommendation:** Change the label to context-appropriate text:
- GGM: "Variable means:"
- OMRF: "Category thresholds:"
- Mixed MRF: "Main effects:" (covers both threshold and mean params)

This requires `print.summary.bgms` to read the model type from the
fit object, which it currently does not do. A lightweight fix: store
`model_type` in `$arguments` and read it in the print method.

---

## L11. `print.bgms()`

Prints "Bayesian Estimation" or "Bayesian Edge Selection" with no
indication of model type (GGM vs OMRF vs mixed MRF).

**Verdict: COSMETIC (low priority).** Adding a model-type line
(e.g., "Model: Mixed MRF (3 ordinal, 2 continuous)") would help
users confirm they're looking at the right fit.

---

## Summary of findings

### Bugs (must fix)

| # | Severity | Location | Issue |
|:-:|:--------:|----------|-------|
| B1 | **Crash** | `extract_category_thresholds.bgms` | Out-of-bounds on `num_categories` for mixed MRF (p+q > p) |
| B2 | **Wrong results** | `extract_pairwise_interactions.bgms` | Edge names in standard upper-tri order; mixed MRF samples in block order |
| B3 | **Minor** | `build_output_mixed_mrf` | SBM allocation names missing from `$raw_samples$parameter_names$allocations` |

### Inconsistencies (should fix for release)

| # | Severity | Location | Issue |
|:-:|:--------:|----------|-------|
| I1 | Moderate | `coef.bgms` | `$main` is list for mixed, matrix for others; undocumented |
| I2 | Minor | `print.summary.bgms` | "Category thresholds:" label incorrect for mixed and GGM |
| I3 | Minor | `print.bgms` | No model-type indication |
| I4 | Minor | `coef.bgms` roxygen | `@return` describes only threshold params |

### Correct (no action needed)

| # | Component | Notes |
|:-:|-----------|-------|
| OK1 | `build_output_mixed_mrf` index mapping | C++/R ordering matches exactly |
| OK2 | `fill_mixed_symmetric` | Block traversal and index scattering correct |
| OK3 | Indicator ordering | Matches pairwise block ordering |
| OK4 | `simulate.bgms` for mixed | Column reordering correct |
| OK5 | `predict.bgms` for mixed | Type-aware dispatch and formatting correct |
| OK6 | NUTS diagnostic wiring | `"hybrid-nuts"` string handled correctly |
| OK7 | `extract_indicators` | Uses stored parameter names (not regenerated) |
| OK8 | `extract_ess` | Reads from summary tables with row names; model-agnostic |
| OK9 | `summary.bgms` | Correct; label issue is in print method only |

---

## Recommended fix order

1. **B1** — `extract_category_thresholds`: guard for mixed MRF.
2. **B2** — `extract_pairwise_interactions`: use stored parameter
   names for mixed MRF.
3. **B3** — `build_output_mixed_mrf`: pass `allocation_names`.
4. **I1/I4** — Document `coef.bgms` return type for mixed MRF.
5. **I2** — Improve `print.summary.bgms` label.
6. **I3** — Add model type to `print.bgms`.
