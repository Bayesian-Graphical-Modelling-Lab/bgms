# Sample observations from the ordinal MRF

\`r lifecycle::badge("deprecated")\`

\`mrfSampler()\` was renamed to \[simulate_mrf()\] as of bgms 0.1.6.3 to
follow the package's naming conventions.

## Usage

``` r
mrfSampler(
  no_states,
  no_variables,
  no_categories,
  interactions,
  thresholds,
  variable_type = "ordinal",
  baseline_category,
  iter = 1000,
  seed = NULL
)
```

## Arguments

- no_states:

  The number of states of the ordinal MRF to be generated.

- no_variables:

  The number of variables in the ordinal MRF.

- no_categories:

  Either a positive integer or a vector of positive integers of length
  `no_variables`. The number of response categories on top of the base
  category: `no_categories = 1` generates binary states.

- interactions:

  A symmetric `no_variables` by `no_variables` matrix of pairwise
  interactions. Only its off-diagonal elements are used.

- thresholds:

  A `no_variables` by `max(no_categories)` matrix of category
  thresholds. The elements in row `i` indicate the thresholds of
  variable `i`. If `no_categories` is a vector, only the first
  `no_categories[i]` elements are used in row `i`. If the Blume-Capel
  model is used for the category thresholds for variable `i`, then row
  `i` requires two values (details below); the first is \\\alpha\\, the
  linear contribution of the Blume-Capel model and the second is
  \\\beta\\, the quadratic contribution.

- variable_type:

  What kind of variables are simulated? Can be a single character string
  specifying the variable type of all `p` variables at once or a vector
  of character strings of length `p` specifying the type for each
  variable separately. Currently, bgm supports “ordinal” and
  “blume-capel”. Binary variables are automatically treated as
  “ordinal”. Defaults to `variable_type = "ordinal"`.

- baseline_category:

  An integer vector of length `no_variables` specifying the
  baseline_category category that is used for the Blume-Capel model
  (details below). Can be any integer value between `0` and
  `no_categories` (or `no_categories[i]`).

- iter:

  The number of iterations used by the Gibbs sampler. The function
  provides the last state of the Gibbs sampler as output. By default set
  to `1e3`.

- seed:

  Optional integer seed for reproducibility. If `NULL`, a seed is
  generated from R's random number generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) can be used before
  calling this function).

## Value

A matrix of simulated observations (see \[simulate_mrf()\]).

## See also

\[simulate_mrf()\] for the current function.
