# Sample from the GGM (Partial-Association) Prior

Draws from the prior of a Gaussian graphical model using the same
constrained NUTS sampler that drives
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
for continuous data. The likelihood is omitted (\\n = 0\\, \\S = 0\\),
so the chain targets the prior alone.

## Usage

``` r
sample_ggm_prior(
  p,
  n_samples,
  n_warmup = 2000,
  interaction_prior = cauchy_prior(scale = 2.5),
  precision_scale_prior = gamma_prior(shape = 1, rate = 1),
  step_size = 0.1,
  max_depth = 10L,
  seed = 1L,
  verbose = TRUE,
  edge_indicators = NULL,
  delta = 0
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

  A `bgms_parameter_prior` for the partial-association off-diagonals
  \\K\_{yy,ij} = -K\_{ij}/2\\. Use
  [`cauchy_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/cauchy_prior.md)
  or
  [`normal_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/normal_prior.md);
  [`beta_prime_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/beta_prime_prior.md)
  is not supported here. Default: `cauchy_prior(scale = 2.5)` (i.e.
  \\K\_{ij}\\ has an implied \\\textrm{Cauchy}(0, 5)\\ prior).

- precision_scale_prior:

  A `bgms_scale_prior` for \\K\_{ii}/2\\. Use
  [`gamma_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/gamma_prior.md)
  or
  [`exponential_prior()`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/exponential_prior.md).
  Default: `gamma_prior(1, 1)`, which implies \\K\_{ii}/2 \sim
  \textrm{Exp}(1)\\ and therefore \\K\_{ii} \sim \textrm{Exp}(1/2)\\
  (mean \\2\\).

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

- delta:

  Non-negative numeric. Determinant-tilt exponent: multiplies the prior
  by \\\|K\|^{\delta}\\, pushing the chain away from the
  positive-definite cone boundary. `delta = 0` (default) recovers the
  untilted prior.

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

The priors are specified on the partial-association scale \\K\_{yy} =
-K/2\\: `interaction_prior` acts on \\K\_{yy,ij} = -K\_{ij}/2\\, and
`precision_scale_prior` acts on \\-K\_{yy,ii} = K\_{ii}/2\\. The same
convention is used by
[`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
and by the continuous block of the mixed-MRF model, so a prior argument
passed here means the same distribution it would mean there. Output
samples are reported as entries of \\K\\; convert with \\K\_{yy} =
-K/2\\ if you want them on the partial-association scale.

When `edge_indicators` is supplied, off-diagonals at excluded positions
are constrained to zero throughout the chain.

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
draws = sample_ggm_prior(
  p = 4, n_samples = 200, n_warmup = 200,
  verbose = FALSE
)
dim(draws$K_offdiag) # 200 x 6
#> [1] 200   6
colnames(draws$K_offdiag) = draws$offdiag_names
head(draws$K_offdiag)
#>           K_1_2        K_1_3     K_1_4       K_2_3      K_2_4
#> [1,]  5.0610611  0.008940521  3.321618 -6.62467649  0.4571407
#> [2,]  1.7664714  0.690624986  4.522583 -0.03245810 -0.8471553
#> [3,] -1.7589055 -0.171594899 -2.997448 -0.20447511  3.6198555
#> [4,]  0.2839390  0.892527750  5.298511 -0.03497855 -0.8921625
#> [5,]  0.9719396  0.116299417 -1.233758  1.49677586  1.6575632
#> [6,] -1.9725988 -0.042434279  2.600275 -0.77422461 -1.8287448
#>            K_3_4
#> [1,] -2.86041681
#> [2,] -2.96334173
#> [3,]  1.86349007
#> [4,] -0.01269069
#> [5,] -1.84625649
#> [6,]  0.67962532

# Sparser graph: drop the (1, 4) edge.
E = matrix(1L, 4, 4)
E[1, 4] = E[4, 1] = 0L
draws = sample_ggm_prior(
  p = 4, n_samples = 200, n_warmup = 200,
  edge_indicators = E, verbose = FALSE
)
colnames(draws$K_offdiag) = draws$offdiag_names
all(draws$K_offdiag[, "K_1_4"] == 0) # TRUE
#> [1] TRUE
# }
```
