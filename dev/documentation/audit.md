# Documentation Audit & Action Plan

Current state of bgms documentation measured against the strategy in
`documentation_strategy.md`. Each item is actionable and scoped.

---

## Tier 1: Exported R functions

### bgm() — OK ✅
Full roxygen: `@title`, `@description`, `@param` (all 37), `@return`
(detailed), `@examples`, `@seealso`, `@export`. 436-line Rd file.

**Action:** Add `@family model-fitting`.

### bgmCompare() — OK ✅
Mirrors bgm() quality. 286-line Rd. Full tags.

**Action:** Add `@family model-fitting`. Use `@inheritParams bgm` for
shared parameters (iter, warmup, chains, etc.) to reduce duplication
and prevent drift.

### simulate_mrf() — OK ✅
Three worked examples, full params, math formulas.

**Action:** Add `@family prediction`.

### mrfSampler() — Minor gaps
Missing `@examples`.

**Action:** Low priority (deprecated function). Add
`@family prediction` and a `\donttest{}` example if easy.

### predict.bgms() — OK ✅
Full tags, separate return descriptions for ordinal vs GGM.

**Action:** Add `@family prediction`, `@inheritParams bgm` for `object`.

### predict.bgmCompare() — OK ✅
**Action:** Add `@family prediction`.

### simulate.bgms() — OK ✅
**Action:** Add `@family prediction`.

### simulate.bgmCompare() — OK ✅
**Action:** Add `@family prediction`.

### coef.bgms() — Minor gaps
Missing `@examples`, `@seealso`.

**Action:** Add `@examples` (1-2 lines extracting coefficients from a
small fit), `@seealso [bgm()]`, `@family posterior-methods`.

### coef.bgmCompare() — Minor gaps
Same as above.

**Action:** Add `@examples`, `@seealso [bgmCompare()]`,
`@family posterior-methods`.

### print.bgms() / print.bgmCompare() — Minor gaps
Missing `@examples`, `@seealso`.

**Action:** Add `@examples`, `@seealso`, `@family posterior-methods`.

### summary.bgms() / summary.bgmCompare() — Minor gaps
Missing `@examples`, `@seealso`.

**Action:** Add `@examples`, `@seealso`, `@family posterior-methods`.

### Extractor functions (12 functions, 1 Rd page) — Gaps

| Issue | Fix |
|-------|-----|
| `@keywords internal` on exported functions | Remove |
| No per-function `@return` | Add return descriptions for each function |
| No `@examples` | Add examples (at least for `extract_arguments`, `extract_pairwise_interactions`, `extract_indicators`) |
| No `@seealso` | Add `@seealso [bgm()], [bgmCompare()]` |

**Action:** This is the highest-priority Tier 1 item. Split into steps:

1. Remove `@keywords internal`
2. Add `@return` to each `extract_*` function
3. Add `@examples` block
4. Add `@seealso` and `@family extractors`

---

## Tier 2: Internal R functions

### Current state
Files use consistent `# ====` / `# ----` / `# @param` pattern.
This matches the strategy. Quality varies:

| File | State | Action needed |
|------|-------|---------------|
| validate_data.R | Very good | None |
| validate_model.R | Very good | None |
| validate_sampler.R | Very good | None |
| bgm_spec.R | Very good | None |
| run_sampler.R | Good | None |
| build_output.R | Good | None |
| mcmc_summary.R | Good | None |
| compute_utils.R | Minimal | Add header blocks to undocumented helpers |
| nuts_diagnostics.R | Adequate | Add `# @param` / `Returns:` |
| simulate_predict.R | Good | Verify new GGM functions have headers |
| zzz.R | Minimal | Low priority (package hooks only) |

**Action:** Review `compute_utils.R` and `nuts_diagnostics.R` only.
Rest is compliant.

---

## Tier 3: C++ headers

### Well-documented ✅

- `omrf_model.h` — class doc, `@param` on constructor, all 21 private
  methods documented. **Gap:** ~30 public accessors/getters/setters
  lack Doxygen (see action items below).
- `bgmCompareOutput` (in `bgmCompare_output.h`) — struct with member
  shapes and purposes in a single `/** */` block. Good template for
  data-container documentation.
- `nuts.h` — `/** */` with `@param`, `@return`, `///` on struct fields
- `hmc.h` — `/** */` with `@param`, `@return`
- `leapfrog.h` — Memoizer class documented
- `metropolis.h` — `/** */` with `@param`, `@return`
- `sampler_base.h` — class doc, method docs
- `chain_runner.h` — `/** */` with `@param` on free functions
- `chain_result.h` — class doc, `@param` on reserve()
- `step_result.h` — all structs documented
- `warmup_schedule.h` — staged schedule explanation
- `edge_prior.h` — base + derived classes documented
- `common_helpers.h` — `/** */` on functions, inline on enums
- `variable_helpers.h` — block comments explaining algorithms
- `progress_manager.h` — `@brief` with feature list

### Needs Doxygen blocks

| File | What's missing | Priority |
|------|---------------|----------|
| `omrf_model.h` | ~30 public accessors lack `/** */` | High |
| `ggm_model.h` | Class doc exists but ~30 private methods undocumented, constructors lack `@param` | High |
| `cholupdate.h` | Two bare function signatures, no docs | High |
| `bgmCompare` free functions | `bgmCompare_helper.h`, `bgmCompare_logp_and_grad.h`, `bgmCompare_sampler.h` have zero Doxygen on functions | High |
| `base_model.h` | No class-level Doxygen, methods use inline `//` only | Medium |
| `sbm_edge_prior.h` | Section dividers only, no Doxygen | Medium |
| `custom_explog.h` | Bare forward declarations | Low |
| `custom_arma_explog.h` | Minimal inline comments | Low |
| `rng_utils.h` | Inline `//` only, no Doxygen | Low |
| `explog_macros.h` | Block comment exists, adequate | None |

**Action (in priority order):**

1. `ggm_model.h` — Add `/** */` blocks to all private methods.
   Use `omrf_model.h` as template.
2. `cholupdate.h` — Add algorithm description, origin note (mgcv),
   `@param` for each argument.
3. `base_model.h` — Add class-level `/** */` explaining the virtual
   interface and what subclasses must implement.
4. `sbm_edge_prior.h` — Add Doxygen blocks on class and methods.
5. `math/` files — Low priority; utility wrappers are self-documenting.

---

## Tier 4: C++ implementations

Not audited exhaustively. Spot checks show adequate inline comments
in `ggm_model.cpp` and `omrf_model.cpp`. No systemic action needed.

---

## Vignettes

| Vignette | Exists | Covers GGM | Action |
|----------|--------|------------|--------|
| intro.Rmd | ✅ | No | Update to mention `variable_type = "continuous"` |
| comparison.Rmd | ✅ | No | None (bgmCompare doesn't support GGM yet) |
| diagnostics.Rmd | ✅ | No | None (diagnostics are model-agnostic) |
| ggm.Rmd | ❌ | — | **Write new vignette** |

**Action:** Write `vignettes/ggm.Rmd` covering:
- Fitting a GGM with `bgm(x, variable_type = "continuous")`
- Interpreting precision matrix output
- Simulation with `simulate()` and `simulate_mrf()`
- Prediction with `predict()`
- Missing data imputation with `na_action = "impute"`

---

## Package metadata

### DESCRIPTION
Current `Description` field mentions only binary/ordinal variable
selection. Needs to mention continuous variables, group comparison,
simulation, and prediction.

### _pkgdown.yml
Currently has only `url:` and `template:`. Needs grouped `reference:`
sections (see strategy doc for proposed layout).

### NEWS.md
Good. No action needed except adding entries for the GGM imputation
work when this branch merges.

---

## Communication

The documentation strategy needs to reach three audiences through
three different channels:

| Channel | Audience | Purpose |
|---------|----------|---------|
| `.github/copilot-instructions.md` | AI agents (Copilot, etc.) | Compact rules auto-injected into agent prompts |
| `CONTRIBUTING.md` | Human contributors | Onboarding guide linking to strategy + styler |
| `dev/documentation/documentation_strategy.md` | Both | Full reference (already written) |

**Note:** The existing `inst/styler/bgms_style.R` handles R *code
formatting* (whitespace, `=` assignment, no space after `if`). It
cannot enforce documentation content (roxygen tags, comment structure,
Doxygen blocks). Documentation conventions must be communicated
through the files above.

---

## Execution order

Prioritised by user impact, CRAN compliance, and dependency order.
Incorporates suggestions from Review 1 and Review 2.

**All tasks land in PR #78 (`ggm_mixed`).** One commit per task,
messages reference the audit number (e.g., `docs: add @return to
extractor functions (audit #1)`).

| Priority | Task | Scope | Rationale |
|----------|------|-------|-----------|
| ~~1~~ | ~~Fix extractor functions docs~~ | ~~R/extractor_functions.R~~ | ~~Done (`3eceb74`). Removed `@keywords internal`, fixed `##'` → `#'`, added `@return`, `@examples`, `@seealso`, `@family extractors`, moved `%||%` to `zzz.R`.~~ |
| ~~2~~ | ~~Add `@examples`/`@seealso` to S3 methods~~ | ~~6 Rd pages (coef, print, summary × 2)~~ | ~~Done (`ccfceb1`). Added `@examples`, `@seealso`, `@family posterior-methods`.~~ |
| ~~3~~ | ~~Create `.github/copilot-instructions.md`~~ | ~~.github/~~ | ~~Done (`c34b674`). Includes `=` assignment, roxygen vs `#`, Doxygen rules, anti-patterns.~~ |
| ~~4~~ | ~~Create `CONTRIBUTING.md` and `dev/README.md`~~ | ~~repo root, dev/~~ | ~~Done (`fd06629`). Onboarding guide linking to strategy, `bgms_style.R`, CI, `inst/CONTRIBUTORS.md`.~~ |
| ~~5~~ | ~~Add `@family` and `@inheritParams`~~ | ~~All exported roxygen~~ | ~~Done (`8fc5fa5`). Added `@family model-fitting` to `bgm()`, `bgmCompare()`; `@family prediction` to `simulate_mrf()`, `simulate.bgms()`, `simulate.bgmCompare()`, `predict.bgms()`, `predict.bgmCompare()`. Tasks #1–#2 already added `@family extractors` and `@family posterior-methods`.~~ |
| ~~6~~ | ~~Update DESCRIPTION~~ | ~~DESCRIPTION~~ | ~~Done. Updated Title to \"Bayesian Analysis of Graphical Models\". Description now mentions continuous variables, GGM, group comparison, simulation, prediction, and missing data imputation.~~ |
| 7 | Create `.editorconfig` | repo root | CRAN-submission blocker: current text omits continuous variables, GGM, group comparison, simulation, prediction. |
| 7 | Create `.editorconfig` | repo root | 5-minute task, prevents trivial formatting diffs across editors. |
| 8 | Document `ggm_model.h` | src/models/ggm/ | High-impact C++ gap: ~30 undocumented private methods including non-trivial math. Add `@param` to constructors. Move TODO (line 33) to `dev/plans/future_tasks.md`. |
| 9 | Document `omrf_model.h` public accessors | src/models/omrf/ | ~30 public getters/setters/capability queries have no Doxygen. All private methods are documented. |
| 10 | Document `cholupdate.h` | src/models/ggm/ | Add algorithm description, origin note (mgcv), `@param` for each argument. |
| 11 | Document `bgmCompare` free functions | src/bgmCompare/ | `bgmCompare_helper.h`, `bgmCompare_logp_and_grad.h`, `bgmCompare_sampler.h` have zero Doxygen on any function. `bgmCompareOutput` struct is fine. |
| 12 | Add pkgdown reference sections | _pkgdown.yml | Now possible because `@family` tags are done (step 5). Add `articles:` section for vignettes. |
| 13 | Document `base_model.h` | src/models/ | Class-level doc explaining virtual interface and what subclasses must implement. |
| 14 | Write GGM vignette | vignettes/ggm.Rmd | Blocked on a suitable continuous dataset. Covers fitting, precision matrix output, simulation, prediction, missing data. |
| 15 | Review `compute_utils.R` / `nuts_diagnostics.R` | R/ | `nuts_diagnostics.R` has zero structured comments; add file banner and `# @param` blocks. `compute_utils.R` needs header blocks on undocumented helpers. |
| 16 | Run `bgms_style()` on assignment-inconsistent files | R/validate_data.R + others | `validate_data.R` uses `<-` throughout; should be `=` per `bgms_style.R`. |
| 17 | Document `sbm_edge_prior.h` | src/priors/ | Convert old-style `// ---` banners to Doxygen `/** */`. |
| 18 | Update intro vignette to mention GGM | vignettes/intro.Rmd | Mention `variable_type = "continuous"`. |
| 19 | Update NEWS.md for GGM imputation | NEWS.md | Add entries for GGM missing-data imputation when `ggm_mixed` merges. |
| 20 | Update README for GGM | README.md, Readme.Rmd | README only mentions binary/ordinal. Update tagline, intro paragraph, and main functions section to cover continuous variables, GGM, simulation, and prediction. Edit `Readme.Rmd` (source) and re-knit. |
