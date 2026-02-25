# Documentation Strategy & Audit Review (2026-02-25)

## 1. Completeness — Missing Audiences
- **Tests & fixtures:** Strategy never addresses documentation for the extensive `tests/testthat` suite or the `dev/fixtures` helpers, so expectations/assertions remain tribal knowledge.
- **Developer automation:** Benchmark scripts, numerical experiments, and scaffold generators under `dev/` lack guidance; new contributors cannot tell how to document their outputs or intended workflows.
- **Datasets:** Nothing states how to document the `.rda` data objects beyond their Rd files; updates to ADHD/Boredom/Wenchuan risk drifting from the style guide.
- **Error/user messages:** Rules for tone and formatting of `stop()` messages (e.g., in `R/validate_data.R`) are missing, so UX varies per author.
- **Operations & CI/CD:** Configure scripts, pkgdown deploy, and GitHub Actions are out of scope even though they need instructions for maintainers.

## 2. Tier Boundaries — Tier 1 vs Tier 2
- The Tier 1/Tier 2 split is reasonable, but core internals like `bgm_spec()`, `run_sampler()`, and `build_output()` act as supported entry points; promoting them to `#' @noRd` roxygen would make their documentation discoverable without exporting.
- Internal functions that shape user-visible results (e.g., `validate_missing_data()` and `handle_impute()` in `R/validate_data.R`) are hard to reference from tests/vignettes under the current "plain #" rule.
- Small utilities (e.g., `center_continuous_data()`) can stay Tier 2, but the blanket ban on `#' @noRd` is too rigid for complex helpers.
- Consider explicitly allowing `#' @noRd` for generics/constructors that advanced users may call directly while still marking them internal.

## 3. Practicality — Pain Points
- The mandate that every R file begins with a `# ===` banner clashes with roxygen: files like `R/extractor_functions.R` must start with `#'` or the block is ignored.
- Requiring runnable `@examples` for print/summary/coef methods forces repeated model fits during `R CMD check`, so authors will just wrap everything in `\donttest{}` against the stated goal.
- Forbidding `#' @noRd` removes help pages for complex internals; contributors must read source instead of `?function`.
- Tier 3 insists on documenting *every* method/field, which is impractical for `src/models/ggm/ggm_model.h` with dozens of short helpers; the earlier "non-obvious logic" rule is more sustainable.
- Private-field `///` comments are required, but there is no lint tooling to enforce them; reviewers must spot issues manually.

## 4. Priorities — Reordering Suggestions
- Keep extractor fixes first, but elevate communication artifacts (`.github/copilot-instructions.md`, `CONTRIBUTING.md`) so newcomers stop adding debt while Tier 1 patches land.
- Documentation for `ggm_model.h` and `cholupdate.h` should move up (before pkgdown/DESCRIPTION polish) because the new GGM feature lacks usable references.
- Pkgdown reference restructuring depends on `@family` tags and extractor fixes; schedule it after those tasks are complete.
- The GGM vignette carries user-facing impact similar to extractor docs; consider moving it immediately after Tier 1 fixes to validate workflows early.
- Low-risk internal cleanups (`compute_utils.R`, `nuts_diagnostics.R`) can remain later—they do not unblock higher-level documentation.

## 5. Communication Channels
- Besides `.github/copilot-instructions.md` and `CONTRIBUTING.md`, add a short `docs/maintenance.md` (linked from CONTRIBUTING) that summarizes the four tiers as a checklist.
- Provide tooling hooks: a pre-commit hook or CI job that runs `styler::style_pkg(style = bgms_style)` plus a roxygen/Doxygen presence check; otherwise the custom style file remains unused.
- Ship a `.lintr` or `.editorconfig` so IDEs enforce spacing rules (e.g., `if(condition)`) without manual reviewer effort.
- Add a PR template or review checklist reminding authors to update `_pkgdown.yml`, DESCRIPTION, and NEWS when they introduce exported symbols.

## 6. Additional Missing/Wrong Items
- `R/extractor_functions.R` still carries `@keywords internal` and lacks per-function `@return`, `@examples`, and `@seealso` despite being exported—contradicts Tier 1 guidance.
- Flagship roxygen blocks (e.g., `R/bgm.R`) still miss the mandated `@family model-fitting` tag.
- `src/models/ggm/ggm_model.h` only documents the class; none of the public/private helpers have `@param`/`@return`, undermining the "gold standard" narrative.
- `src/models/omrf/omrf_model.h` lacks Doxygen on several getters/setters even though the audit lists it as complete.
- `_pkgdown.yml` still lacks the referenced `reference:` sections, so pkgdown will not reflect the proposed grouping.
- The execution plan does not mention updating dataset Rd files or NEWS for the new GGM imputation support revealed by recent smoke tests.
