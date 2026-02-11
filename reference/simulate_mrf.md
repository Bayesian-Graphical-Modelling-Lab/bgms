# Simulate Observations from an Ordinal MRF

\`simulate_mrf()\` generates observations from an ordinal Markov Random
Field using Gibbs sampling with user-specified parameters.

## Usage

``` r
simulate_mrf(
  num_states,
  num_variables,
  num_categories,
  pairwise,
  main,
  variable_type = "ordinal",
  baseline_category,
  iter = 1000,
  seed = NULL
)
```

## Arguments

- num_states:

  The number of states of the ordinal MRF to be generated.

- num_variables:

  The number of variables in the ordinal MRF.

- num_categories:

  Either a positive integer or a vector of positive integers of length
  `num_variables`. The number of response categories on top of the base
  category: `num_categories = 1` generates binary states.

- pairwise:

  A symmetric `num_variables` by `num_variables` matrix of pairwise
  interactions. Only its off-diagonal elements are used.

- main:

  A `num_variables` by `max(num_categories)` matrix of category
  thresholds. The elements in row `i` indicate the thresholds of
  variable `i`. If `num_categories` is a vector, only the first
  `num_categories[i]` elements are used in row `i`. If the Blume-Capel
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

  An integer vector of length `num_variables` specifying the
  baseline_category category that is used for the Blume-Capel model
  (details below). Can be any integer value between `0` and
  `num_categories` (or `num_categories[i]`).

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

A `num_states` by `num_variables` matrix of simulated states of the
ordinal MRF.

## Details

The Gibbs sampler is initiated with random values from the response
options, after which it proceeds by simulating states for each variable
from a logistic model using the other variable states as predictor
variables.

There are two modeling options for the category thresholds. The default
option assumes that the category thresholds are free, except that the
first threshold is set to zero for identification. The user then only
needs to specify the thresholds for the remaining response categories.
This option is useful for any type of ordinal variable and gives the
user the most freedom in specifying their model.

The Blume-Capel option is specifically designed for ordinal variables
that have a special type of baseline_category category, such as the
neutral category in a Likert scale. The Blume-Capel model specifies the
following quadratic model for the threshold parameters:
\$\$\mu\_{\text{c}} = \alpha \times (\text{c} - \text{r}) + \beta \times
(\text{c} - \text{r})^2,\$\$ where \\\mu\_{\text{c}}\\ is the threshold
for category c (which now includes zero), \\\alpha\\ offers a linear
trend across categories (increasing threshold values if \\\alpha \> 0\\
and decreasing threshold values if \\\alpha \<0\\), if \\\beta \< 0\\,
it offers an increasing penalty for responding in a category further
away from the baseline_category category r, while \\\beta \> 0\\
suggests a preference for responding in the baseline_category category.

## See also

[`simulate.bgms`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/simulate.bgms.md)
for simulating from a fitted model.

## Examples

``` r
# Generate responses from a network of five binary and ordinal variables.
num_variables = 5
num_categories = sample(1:5, size = num_variables, replace = TRUE)

Pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
Pairwise[2, 1] = Pairwise[4, 1] = Pairwise[3, 2] =
  Pairwise[5, 2] = Pairwise[5, 4] = .25
Pairwise = Pairwise + t(Pairwise)
Main = matrix(0, nrow = num_variables, ncol = max(num_categories))

x = simulate_mrf(
  num_states = 1e3,
  num_variables = num_variables,
  num_categories = num_categories,
  pairwise = Pairwise,
  main = Main
)
#> Warning: The matrix ``main'' contains numeric values for variable 1 for category 
#> (categories, i.e., columns) exceding the maximum of 1. These values will 
#> be ignored.
#> Warning: The matrix ``main'' contains numeric values for variable 3 for category 
#> (categories, i.e., columns) exceding the maximum of 2. These values will 
#> be ignored.
#> Warning: The matrix ``main'' contains numeric values for variable 4 for category 
#> (categories, i.e., columns) exceding the maximum of 2. These values will 
#> be ignored.
#> Warning: The matrix ``main'' contains numeric values for variable 5 for category 
#> (categories, i.e., columns) exceding the maximum of 3. These values will 
#> be ignored.

# Generate responses from a network of 2 ordinal and 3 Blume-Capel variables.
num_variables = 5
num_categories = 4

Pairwise = matrix(0, nrow = num_variables, ncol = num_variables)
Pairwise[2, 1] = Pairwise[4, 1] = Pairwise[3, 2] =
  Pairwise[5, 2] = Pairwise[5, 4] = .25
Pairwise = Pairwise + t(Pairwise)

Main = matrix(NA, num_variables, num_categories)
Main[, 1] = -1
Main[, 2] = -1
Main[3, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
Main[5, ] = sort(-abs(rnorm(4)), decreasing = TRUE)

x = simulate_mrf(
  num_states = 1e3,
  num_variables = num_variables,
  num_categories = num_categories,
  pairwise = Pairwise,
  main = Main,
  variable_type = c("b", "b", "o", "b", "o"),
  baseline_category = 2
)
```
