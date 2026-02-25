# bgms Documentation Strategy

This document defines the documentation standard for bgms—what to document,
how to write it, and at what level of detail. It applies to both human
contributors and AI agents.

## Guiding principles

1. **User docs live in roxygen; developer docs live in code comments.**
   End users read `?bgm` and vignettes. Contributors read source files.
   These are different audiences with different needs.

2. **Every exported symbol must have complete roxygen.**
   CRAN requires this. We go further: every exported function gets
   `@param` (all parameters), `@return`, `@examples`, and `@seealso`.

3. **Internal functions get structured comments, not roxygen.**
   This matches the convention in brms and rstanarm—the two most mature
   Bayesian R packages with architectures comparable to ours.

4. **C++ headers document design; implementations document mechanics.**
   Header files explain *what* a class/function does and *why*.
   Implementation files explain *how*, when the algorithm is non-obvious.

5. **Be concise.** One clear sentence beats three vague ones.
   Don't restate what the code already says.

6. **Document more rather than less.** When in doubt, write the
   comment. An over-documented codebase is easier to maintain
   than an under-documented one.

---

## Writing style and tone

bgms documentation is written by an academic research team for an
academic audience. The existing tone — visible in `bgm()` roxygen,
the intro vignette, and `validate_data.R` error messages — is the
standard. New documentation should match it.

### Rules

1. **Write in plain, direct English.** State what something does,
   then stop. Avoid filler ("it is important to note that",
   "this function serves to", "it should be noted that").
   Bad: "This function is responsible for computing the gradient
   of the log-posterior." Good: "Computes the gradient of the
   log-posterior."

2. **Use present tense, active voice, imperative mood** for
   function summaries. "Computes ...", "Validates ...",
   "Returns ...". Not "This function will compute" or "This is
   used to validate".

3. **Be precise about statistics.** Use correct terminology:
   "pseudoposterior", "precision matrix", "spike-and-slab prior",
   "Markov random field". Do not simplify these for a lay audience.

4. **Do not write AI-style prose.** Avoid:
   - Superlatives and enthusiasm ("powerful", "elegant",
     "robust solution", "seamlessly")
   - Hedging qualifiers ("arguably", "it is worth noting")
   - Transition phrases between paragraphs
     ("Furthermore", "Moreover", "Additionally")
   - Summarising what was just said or will be said
   - Bullet lists where a sentence suffices
   - Overlong parameter descriptions when a clause will do

5. **Do not write session-oriented documentation.** Avoid phrases
   like "in your R session", "you should see", "run this in your
   console". Documentation describes what a function does and
   returns, not what the reader should do interactively.

6. **Match the existing voice.** Before writing new documentation,
   read the `bgm()` roxygen block and the intro vignette. The tone
   is academic, technical, and concise — not chatty, not formal,
   not promotional. New text should be indistinguishable from the
   existing text.

7. **Error messages** follow the pattern in `validate_data.R`:
   state what went wrong, why, and what the user can do about it.
   Use complete sentences. Do not start with "Error:" — R adds
   that automatically.

8. **Code comments** are notes to the next developer, not prose.
   "Centering needed because suf_stat = X'X conflates means with
   precision" — not "We center the data here in order to ensure
   that the sufficient statistic correctly reflects the precision
   matrix structure."

---

## Tiers of documentation

bgms has four documentation tiers, each with its own rules.

### Tier 1: User-facing R (exported functions, S3 methods)

**Audience:** Package users reading `?bgm` or the pkgdown site.

**Standard:** Full roxygen2 block with the following tags (in order):

```r
#' @title   One-line imperative title (no period)
#' @description  One paragraph: what it does and when to use it.
#' @param name  Sentence fragment: what the parameter controls.
#'              Multi-line is fine for complex parameters.
#' @return  What the function returns. For complex objects, describe the
#'          structure (list elements, matrix dimensions, etc.).
#' @details  (optional) Longer explanation, math, edge cases.
#' @section Section Name:  (optional) For grouping related details.
#' @examples
#'   Minimal, runnable example. Wrap slow ones in \donttest{}.
#' @seealso  Cross-references to related functions.
#' @export
```

**Rules:**

- `@param` for *every* parameter, no exceptions. Use `@inheritParams bgm`
  for shared parameters across `bgm()`, `bgmCompare()`, and their methods.
- `@return` is mandatory. Describe the *structure*, not just the type.
  Bad: "A list." Good: "A list with elements `pairwise` (matrix),
  `thresholds` (matrix), and `indicators` (matrix)."
- `@examples` must be runnable by `R CMD check`. Use `\donttest{}` for
  anything over ~5 seconds. Prefer `\donttest{}` over `\dontrun{}`.
  Use `\dontrun{}` only when the example cannot succeed in a clean R
  session (e.g., requires credentials, network, or a pre-fitted object
  that cannot be created quickly).
- `@seealso` must link to the primary fitting function (`[bgm()]` or
  `[bgmCompare()]`) and any closely related methods.
- Use `@family` tags to group related functions:
  - `@family model-fitting` for `bgm`, `bgmCompare`
  - `@family posterior-methods` for `coef`, `print`, `summary`
  - `@family prediction` for `predict`, `simulate`, `simulate_mrf`
  - `@family extractors` for `extract_*` functions
  - `@family diagnostics` for `extract_rhat`, `extract_ess`

**Tag order:**

Preferred (not required) order: `@title`, `@description`, `@param`,
`@return`, `@details`, `@section`, `@examples`, `@seealso`, `@export`.
For very long roxygen blocks (e.g., `bgm()`) it is acceptable to place
`@seealso` earlier if it aids readability.

**S3 methods:**

- Each method gets its own Rd page (current practice).
- Use `@inheritParams bgm` for the `object` argument.
- The `@return` section should describe method-specific output, not
  repeat the generic contract.

**`@inheritParams` note:**

`@inheritParams bgm` only inherits documentation for parameters whose
names match the current function's formals. It is safe to use even when
the receiving function has fewer parameters than `bgm()`.

**Extractor functions:**

- Currently grouped on one page via `@rdname extractor_functions`.
  This is fine, but each function needs its own `@return` description.
- Remove `@keywords internal`—these are exported, user-facing functions.

### Tier 2: Internal R (non-exported helpers)

**Audience:** Developers maintaining the R layer.

**Standard:** Structured comment block (plain `#`, not `#'`).
This matches the convention in brms and rstanarm.

> **Note on `#' @noRd`:** The default for internal functions is plain
> `#` comments. However, complex internal entry points — functions with
> many parameters, non-obvious return structures, or that act as
> internal API boundaries (e.g., `bgm_spec()`, `run_sampler()`) — may
> use `#' @noRd` roxygen to gain `@inheritParams` and richer parameter
> documentation. Use `@noRd` only when the function genuinely benefits
> from roxygen tooling; do not use it on simple helpers.

```r
# ------------------------------------------------------------------------------
# function_name
# ------------------------------------------------------------------------------
#
# One-sentence purpose.
#
# @param arg1  What it is.
# @param arg2  What it is.
#
# Returns: What comes back.
# ------------------------------------------------------------------------------
function_name <- function(arg1, arg2) {
  ...
}
```

**Rules:**

- The `# @param` pseudo-tags are for developer readability—`roxygen2`
  ignores them (they use `#`, not `#'`).
- Use file-level `# ====` banners to separate sections within a file.
- Use `# ----` banners to separate individual functions.
- Explain *why* the function exists if the name doesn't make it obvious.
- Include `Returns:` when the return value isn't trivial.
- Do NOT add `#' @noRd`—we use plain `#` everywhere for internal code
  to keep a clear visual distinction from exported roxygen.
- If a function contains non-obvious logic (numerical tricks, edge
  cases, R quirks), add inline `#` comments at the relevant lines.

### Tier 2b: Test files

**Audience:** Developers writing or reviewing tests.

**Standard:** The test suite already follows a good pattern. Codify it:

- Each test file starts with a `# ====` file banner naming scope and
  test strategy (e.g., `# EXTENDS: test-methods.R`,
  `# PATTERN: fixture-based`).
- `test_that()` descriptions are imperative sentences
  (e.g., "extract pairwise interactions from a fitted model").
- Helper fixtures and setup code get `# @param` pseudo-tags.
- Integration notes explain which feature or bug a test covers.

### Error messages

Error messages produced by `stop()`, `warning()`, and `message()` should:

1. State **what** went wrong.
2. State **why** (if not obvious from the "what").
3. Suggest **what the user can do** about it.

Reference `validate_data.R` as the exemplar, e.g.:
`"You could try option na_action = 'impute'."`

### Tier 3: C++ headers (`.h` files)

**Audience:** Developers maintaining the C++ backend.

**Standard:** Doxygen-style `/** */` blocks on classes, public methods,
and free functions. Use `///` for struct fields and enum values.

```cpp
/**
 * ClassName — one-line description
 *
 * Longer explanation of the class's role in the architecture:
 * what it owns, what invariants it maintains, how it fits into
 * the BaseModel → GGMModel/OMRFModel hierarchy.
 */
class ClassName : public BaseClass {
public:
    /**
     * Brief description of this method
     *
     * @param arg1  What it is
     * @param arg2  What it is
     * @return      What comes back
     */
    ReturnType method_name(Type arg1, Type arg2);

private:
    arma::mat data_;  ///< One-line explanation of this field
};
```

**Rules:**

- **Class-level doc** is mandatory on every class. Explain the class's
  purpose and its place in the architecture.
- **Public methods** get `/** */` blocks with `@param` and `@return`.
- **Private methods** get `/** */` blocks. We prefer to document
  too much rather than too little. Even short helpers benefit from
  a one-line `/** Brief description */` block. For mathematical
  methods (e.g., `constrained_diagonal`, `compute_inv_submatrix_i`),
  also explain *which formula* they implement.
- **Private fields** get `///` one-liners when their purpose is
  non-obvious (scratch vectors, caches, flags). Self-explanatory
  fields (e.g., `size_t n_`, `size_t p_`) may omit `///` if there
  is a grouped comment or the class-level doc explains the notation.
- **Free functions** (in `algorithms/`, `math/`) get `/** */` blocks
  with `@param` and `@return`.
- **Reference examples:**
  - `bgmCompareOutput` in `bgmCompare_output.h` — well-documented
    struct: each member's shape and purpose in one `/** */` block.
  - `omrf_model.h` — all 21 private methods have `/** */` blocks.
    (Public accessors still need Doxygen — see audit.)
  - No single file is the complete gold standard yet. Use these as
    partial templates and apply the rules above consistently.
- For ported algorithms (e.g., `cholupdate.h` from mgcv), add an
  origin note: where it came from and what algorithm it implements.

### Tier 4: C++ implementations (`.cpp` files)

**Audience:** Developers debugging or modifying algorithms.

**Standard:** Inline comments only. Do not repeat the header doc.

**Rules:**

- Add `//` comments for non-obvious steps: numerical tricks, index
  arithmetic, order-of-operations dependencies.
- For long functions, use `// --- Phase 1: ... ---` section comments.
- For mathematical derivations, reference the formula by name or
  equation number rather than re-deriving it in comments.
- No Doxygen blocks—those belong in headers only.

---

## File-level organization

Every R source file should start with a file-level banner:

```r
# ==============================================================================
# Short description of what this file contains
# ==============================================================================
#
# Optional: 2-3 sentences of context. What role this file plays in the
# package architecture. What the main public/internal entry point is.
# ==============================================================================
```

**Exception:** For files whose first function is exported, the file
may start directly with the roxygen `#'` block (or with a roxygen
`NULL` documentation block like `#' @rdname ...`). The banner
convention primarily targets internal-only files and mixed files
where internal code precedes exported code.
```

Every C++ header should have the `#pragma once` directive and then
the `#include` block. The class-level Doxygen block serves as the
file-level documentation.

---

## Vignettes

**Existing:** Three vignettes (intro, comparison, diagnostics).

**Needed:** A GGM vignette showing continuous-variable workflows
(fitting, simulation, prediction, missing data imputation).

**Style:**

- Title format: "Topic with bgms" (e.g., "Gaussian Graphical Models
  with bgms")
- Structure: problem setup → data → model fitting → results →
  simulation/prediction → diagnostics
- Keep examples fast (< 30 seconds). Use pre-computed results for
  heavy MCMC if needed.
- Include references via `refs.bib`.

---

## pkgdown site

Add grouped reference sections to `_pkgdown.yml`:

```yaml
reference:
  - title: "Model Fitting"
    contents:
      - bgm
      - bgmCompare
  - title: "Model Summary"
    contents:
      - print.bgms
      - summary.bgms
      - coef.bgms
      - print.bgmCompare
      - summary.bgmCompare
      - coef.bgmCompare
  - title: "Prediction and Simulation"
    contents:
      - predict.bgms
      - predict.bgmCompare
      - simulate.bgms
      - simulate.bgmCompare
      - simulate_mrf
  - title: "Posterior Extraction"
    contents:
      - extractor_functions
  - title: "Datasets"
    contents:
      - ADHD
      - Boredom
      - Wenchuan
```

---

## DESCRIPTION

The `Description` field should reflect the full scope of the package.
Current wording mentions only variable selection and binary/ordinal
variables. Update to mention continuous variables (GGM), group
comparison, simulation, and prediction.

---

## NEWS.md

Continue the current convention:

```markdown
# bgms X.Y.Z

## New features
- Bullet per feature.

## Bug fixes
- Bullet per fix.

## Other changes
- Bullet per change.
```

Each entry should say *what changed* and *why the user cares*,
not implementation details.

---

## For AI agents

When an AI agent generates or modifies bgms code, it must follow
these rules:

1. **Exported functions:** Write full roxygen with all tags listed
   in Tier 1 above. Use `@inheritParams` when the parameter is
   already documented on `bgm()` or `bgmCompare()`.

2. **Internal functions:** Write a `# ----` header block with
   purpose, `# @param`, and `Returns:`. Do NOT use `#'` roxygen.

3. **C++ headers:** Write Doxygen `/** */` blocks on all classes
   and public methods. Include `@param` and `@return`.

4. **C++ implementations:** Add `//` inline comments only for
   non-obvious logic. Do not add Doxygen blocks.

5. **Never** add `@keywords internal` to an exported function.

6. **Never** use `\dontrun{}` when `\donttest{}` suffices.

7. When modifying an existing function's signature (adding/removing
   params), update the roxygen/comment block in the same commit.

8. When adding a new exported function, also add it to the
   appropriate `reference:` section in `_pkgdown.yml`.

9. **Use `=` for assignment** (not `<-`). The package enforces this
   via `inst/styler/bgms_style.R`. Run
   `styler::style_pkg(style = bgms_style)` before committing.

10. **No space between `if`/`for`/`while` and `(`** — also enforced
    by `bgms_style.R`.
---

## Communication channels

The documentation strategy reaches three audiences through different
channels:

| Channel | Audience | Purpose |
|---------|----------|---------|
| `.github/copilot-instructions.md` | AI agents (Copilot, etc.) | Compact rules auto-injected into agent prompts. Keep under 100 lines. |
| `CONTRIBUTING.md` | Human contributors | Onboarding: build/test, code style (link `bgms_style.R`), docs (link this file), CI. |
| This file | Both | Full reference. |

**Additional artifacts:**

- **`.editorconfig`:** Enforce indent style, trailing whitespace, final
  newline across editors. Prevents trivial formatting diffs.
- **`dev/README.md`:** One paragraph explaining that `dev/` contains
  developer-only materials excluded from the built package.
- **`inst/CONTRIBUTORS.md`:** Already exists. Reference from
  `CONTRIBUTING.md`.

**Not needed now (revisit later):**

- Pre-commit hooks for roxygen linting — the codebase is mostly
  compliant. Revisit after audit backlog clears.
- Custom CI documentation checks — `R CMD check` handles most cases.

---

## Datasets

Each `.rda` dataset in `data/` must have a corresponding `.Rd` page
(in `man/`) with `@format`, `@source`, and `@description`. When
updating a dataset, ensure the Rd page stays in sync.