# Sample Precision Matrices from the GGM Prior

Draws precision matrices \\K\\ from the prior of a Gaussian graphical
model using the same constrained NUTS sampler that drives
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
for continuous data. The likelihood is omitted (\\n = 0\\, \\S = 0\\),
so the chain targets the prior alone.

## Usage

``` r
sample_precision_prior(
  p,
  n_samples,
  n_warmup = 2000,
  interaction_prior = cauchy_prior(scale = 2.5),
  precision_scale_prior = gamma_prior(shape = 1, rate = 1),
  step_size = 0.1,
  max_depth = 10L,
  seed = 1L,
  verbose = TRUE,
  edge_indicators = NULL
)
```

## Arguments

- p:

  Integer. Dimension of the precision matrix (\\p \ge 2\\).

- n_samples:

  Integer. Number of post-warmup draws to keep.

- n_warmup:

  Integer. NUTS warmup iterations. Default `2000`.

- interaction_prior:

  A `bgms_parameter_prior` for the off-diagonal entries. Use
  [`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md)
  or
  [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md);
  [`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md)
  is not supported here. Default: `cauchy_prior(scale = 2.5)`.

- precision_scale_prior:

  A `bgms_scale_prior` for the diagonal entries of \\K\\. Use
  [`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md)
  or
  [`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md).
  Default: `gamma_prior(1, 1)`.

- step_size:

  Positive numeric. Initial NUTS step size used to seed dual-averaging
  adaptation. Default `0.1`.

- max_depth:

  Integer. Maximum NUTS tree depth. Default `10`.

- seed:

  Integer. RNG seed for the chain. Default `1L`.

- verbose:

  Logical. If `TRUE` (default), print a progress bar.

- edge_indicators:

  Optional integer \\p \times p\\ matrix with `1` = edge included, `0` =
  excluded. Must be symmetric with `1`s on the diagonal. Default: full
  graph (all edges included).

## Value

A list with components

- `K_offdiag`:

  Numeric matrix of size `n_samples` x `p * (p - 1) / 2` containing the
  upper-triangle off-diagonal entries of \\K\\ for each draw, in
  column-major order \\(K\_{12}, K\_{13}, K\_{23}, K\_{14}, \ldots)\\.
  Excluded edges are returned as `0`.

- `K_diag`:

  Numeric matrix of size `n_samples` x `p` containing the diagonal
  entries \\K\_{11}, \ldots, K\_{pp}\\.

- `offdiag_names`:

  Character vector of length `p * (p - 1) / 2` naming the columns of
  `K_offdiag` (e.g. `"K_1_2"`).

- `diag_names`:

  Character vector of length `p` naming the columns of `K_diag`.

- `step_size`:

  The (initial) NUTS step size used.

- `edge_indicators`:

  The integer edge-indicator matrix used (full graph if not supplied).

## Details

Off-diagonals are placed on the association scale \\K\_{yy,ij} =
-K\_{ij}/2\\ and assigned the supplied `interaction_prior`. A
`normal_prior(scale = s)` therefore constrains \\K\_{yy}\\ with standard
deviation \\s\\, equivalent to a \\\textrm{Normal}(0, 2s)\\ prior on
\\K\_{ij}\\ itself. The diagonals \\K\_{ii}\\ are drawn from the
supplied `precision_scale_prior`. When `edge_indicators` is supplied,
off-diagonals at excluded positions are constrained to zero throughout
the chain.

## See also

[`cauchy_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md),
[`normal_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md),
[`gamma_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md),
[`exponential_prior`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md),
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)

## Examples

``` r
# \donttest{
# Default Cauchy(0, 2.5) off-diagonal, Gamma(1, 1) diagonal, p = 4.
draws = sample_precision_prior(
  p = 4, n_samples = 200, n_warmup = 200,
  verbose = FALSE
)
dim(draws$K_offdiag) # 200 x 6
#> [1] 200   6
colnames(draws$K_offdiag) = draws$offdiag_names
head(draws$K_offdiag)
#>          K_1_2      K_1_3      K_1_4       K_2_3       K_2_4
#> [1,] -1.233684 -1.9146341  0.7934043  2.45108392 -1.99732141
#> [2,]  1.285579 -0.5028680  1.0132068  0.59625222 -0.17293213
#> [3,] -1.037669  0.2762245 -1.4285318 -0.64342766  0.87019167
#> [4,]  1.297038 -0.4459275  1.6381413 -0.07162156  0.65340822
#> [5,]  1.664220 -0.3653032  1.3650031  0.48766110 -0.03671729
#> [6,] -1.240221  0.3221408 -0.5824119 -0.33342054 -0.83995087
#>            K_3_4
#> [1,] -1.04512815
#> [2,] -2.55826117
#> [3,]  1.45388788
#> [4,] -0.01859332
#> [5,] -1.40563154
#> [6,] -0.79436832

# Sparser graph: drop the (1, 4) edge.
E = matrix(1L, 4, 4)
E[1, 4] = E[4, 1] = 0L
draws = sample_precision_prior(
  p = 4, n_samples = 200, n_warmup = 200,
  edge_indicators = E, verbose = FALSE
)
colnames(draws$K_offdiag) = draws$offdiag_names
all(draws$K_offdiag[, "K_1_4"] == 0) # TRUE
#> [1] TRUE
# }
```
