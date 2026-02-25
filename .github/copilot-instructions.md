# bgms — AI Agent Instructions

Rules for AI agents (Copilot, Claude, etc.) working on this codebase.

## R code style

- Use `=` for assignment, not `<-`.
- **Exception**: inside function-call arguments where the result must
  be captured, use `<-`. With `=` R treats it as a named argument.
  Example: `expect_message(result <- foo(), "pattern")`.
- No space between `if`/`for`/`while` and `(`: write `if(`, not `if (`.
- Enforced by `inst/styler/bgms_style.R`. Run before committing:
  ```r
  source("inst/styler/bgms_style.R")
  styler::style_pkg(style = bgms_style)
  ```
  After running, check test files for `expect_*(result = ...)` and
  revert those to `result <- ...`.

## Exported R functions (Tier 1)

- Full roxygen2 block: `@title`, `@description`, `@param` (all params),
  `@return`, `@examples`, `@seealso`, `@export`.
- `@param` for every parameter, no exceptions. Use `@inheritParams bgm`
  for shared parameters.
- `@return` must describe the structure, not just the type.
- `@examples` must be runnable by `R CMD check`. Use `\donttest{}` for
  slow examples. Use `\dontrun{}` only when the example truly cannot
  run in a clean R session.
- Add `@family` tags: `model-fitting`, `posterior-methods`, `prediction`,
  `extractors`, `diagnostics`.
- Never use `@keywords internal` on exported functions.

## Internal R functions (Tier 2)

- Use plain `#` comments, not `#'` roxygen. No `@noRd` unless the
  function genuinely benefits from `@inheritParams`.
- Use `# ====` banners for file sections, `# ----` banners per function.
- Format:
  ```r
  # ------------------------------------------------------------------
  # function_name
  # ------------------------------------------------------------------
  # One-sentence purpose.
  #
  # @param arg1  What it is.
  # @param arg2  What it is.
  #
  # Returns: What comes back.
  # ------------------------------------------------------------------
  ```

## C++ headers (Tier 3)

- Doxygen `/** */` blocks on all classes, public methods, and free
  functions. Use `///` for struct fields and enum values.
- Class-level doc is mandatory: explain purpose and place in the
  `BaseModel → GGMModel/OMRFModel` hierarchy.
- Private methods get `/** */` blocks. Prefer over-documenting.
- For ported algorithms, add an origin note (source, algorithm name).
- No Doxygen in `.cpp` files — use inline `//` comments only.

## C++ implementations (Tier 4)

- Inline `//` comments for non-obvious steps only.
- Use `// --- Phase N: ... ---` section comments in long functions.
- Reference formulas by name or equation number.

## Error messages

- State what went wrong, why, and what the user can do about it.
- Use complete sentences. Do not start with "Error:" — R adds that.
- See `R/validate_data.R` as the exemplar.

## Commit conventions

- One commit per audit task. Message format:
  `docs: short description (audit #N)`.
- When modifying a function signature, update its documentation in
  the same commit.
- When adding a new exported function, add it to `_pkgdown.yml`.

## Do not

- Add `@keywords internal` to exported functions.
- Use `\dontrun{}` when `\donttest{}` suffices.
- Add Doxygen blocks in `.cpp` implementation files.
- Use `<-` for assignment in R code.
- Use AI-style prose: no superlatives, no hedging qualifiers, no
  transition phrases ("Furthermore", "Moreover", "Additionally").
- Write session-oriented documentation ("in your R session",
  "you should see", "run this in your console"). Documentation
  describes what a function does and returns, not what the reader
  should do interactively.
