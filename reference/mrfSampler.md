# Sample observations from the ordinal MRF

This function samples states from the ordinal MRF using a Gibbs sampler.
The Gibbs sampler is initiated with random values from the response
options, after which it proceeds by simulating states for each variable
from a logistic model using the other variable states as predictor
variables.

## Usage

``` r
mrfSampler(
  no_states,
  no_variables,
  no_categories,
  interactions,
  thresholds,
  variable_type = "ordinal",
  reference_category,
  iter = 1000
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
  “ordinal’’. Defaults to `variable_type = "ordinal"`.

- reference_category:

  An integer vector of length `no_variables` specifying the
  reference_category category that is used for the Blume-Capel model
  (details below). Can be any integer value between `0` and
  `no_categories` (or `no_categories[i]`).

- iter:

  The number of iterations used by the Gibbs sampler. The function
  provides the last state of the Gibbs sampler as output. By default set
  to `1e3`.

## Value

A `no_states` by `no_variables` matrix of simulated states of the
ordinal MRF.

## Details

There are two modeling options for the category thresholds. The default
option assumes that the category thresholds are free, except that the
first threshold is set to zero for identification. The user then only
needs to specify the thresholds for the remaining response categories.
This option is useful for any type of ordinal variable and gives the
user the most freedom in specifying their model.

The Blume-Capel option is specifically designed for ordinal variables
that have a special type of reference_category category, such as the
neutral category in a Likert scale. The Blume-Capel model specifies the
following quadratic model for the threshold parameters:
\$\$\mu\_{\text{c}} = \alpha \times \text{c} + \beta \times (\text{c} -
\text{r})^2,\$\$ where \\\mu\_{\text{c}}\\ is the threshold for category
c (which now includes zero), \\\alpha\\ offers a linear trend across
categories (increasing threshold values if \\\alpha \> 0\\ and
decreasing threshold values if \\\alpha \<0\\), if \\\beta \< 0\\, it
offers an increasing penalty for responding in a category further away
from the reference_category category r, while \\\beta \> 0\\ suggests a
preference for responding in the reference_category category.

## Examples

``` r
# Generate responses from a network of five binary and ordinal variables.
no_variables = 5
no_categories = sample(1:5, size = no_variables, replace = TRUE)

Interactions = matrix(0, nrow = no_variables, ncol = no_variables)
Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
  Interactions[5, 2] = Interactions[5, 4] = .25
Interactions = Interactions + t(Interactions)
Thresholds = matrix(0, nrow = no_variables, ncol = max(no_categories))

x = mrfSampler(
  no_states = 1e3,
  no_variables = no_variables,
  no_categories = no_categories,
  interactions = Interactions,
  thresholds = Thresholds
)
#> Warning: The matrix ``thresholds'' contains numeric values for variable 1 for category 
#> (categories, i.e., columns) exceding the maximum of 1. These values will 
#> be ignored.
#> Warning: The matrix ``thresholds'' contains numeric values for variable 5 for category 
#> (categories, i.e., columns) exceding the maximum of 2. These values will 
#> be ignored.

# Generate responses from a network of 2 ordinal and 3 Blume-Capel variables.
no_variables = 5
no_categories = 4

Interactions = matrix(0, nrow = no_variables, ncol = no_variables)
Interactions[2, 1] = Interactions[4, 1] = Interactions[3, 2] =
  Interactions[5, 2] = Interactions[5, 4] = .25
Interactions = Interactions + t(Interactions)

Thresholds = matrix(NA, no_variables, no_categories)
Thresholds[, 1] = -1
Thresholds[, 2] = -1
Thresholds[3, ] = sort(-abs(rnorm(4)), decreasing = TRUE)
Thresholds[5, ] = sort(-abs(rnorm(4)), decreasing = TRUE)

x = mrfSampler(
  no_states = 1e3,
  no_variables = no_variables,
  no_categories = no_categories,
  interactions = Interactions,
  thresholds = Thresholds,
  variable_type = c("b", "b", "o", "b", "o"),
  reference_category = 2
)
```
