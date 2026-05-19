# Sample from the GGM (Partial-Association) Prior

Draws from the prior of a Gaussian graphical model. The likelihood is
omitted (\\n = 0\\, \\S = 0\\), so the chain targets the prior alone.
Two specifications are supported via the `spec` argument:

- `"conditional"` (default): fix a graph \\\Gamma\\ and sample \\K \mid
  \Gamma\\ via the same constrained NUTS sampler that drives
  [`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  for continuous data. The chain targets \\p(K \mid \Gamma) \propto
  \mathrm{slab}(K) \cdot \mathrm{diag}(K) \cdot \|K\|^{\delta} \cdot
  \mathbf{1}\\K \in \mathcal{M}^{+}(\Gamma)\\ / Z(\Gamma)\\.

- `"joint"`: sample \\(K, \Gamma)\\ jointly from the un-normalised joint
  prior \\p(K, \Gamma) \propto \mathrm{slab}(K) \cdot \mathrm{diag}(K)
  \cdot \|K\|^{\delta} \cdot \mathbf{1}\\K \in \mathcal{M}^{+}(\Gamma)\\
  \cdot \pi(\Gamma)\\. Uses the adaptive-Metropolis MH chain from
  [`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)
  with edge selection on and the likelihood off, so the marginal on
  \\\Gamma\\ is \\\pi(\Gamma) \cdot Z(\Gamma)\\ (joint specification,
  not hierarchical). Useful for simulation-based calibration of
  [`bgm`](https://bayesian-graphical-modelling-lab.github.io/bgms/reference/bgm.md)'s
  default sampler.

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
  delta = NULL,
  spec = c("conditional", "joint"),
  edge_inclusion_prob = 0.5
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
  adaptation. Default `0.1`. Used only for `spec = "conditional"` (NUTS
  path); ignored for the `"joint"` MH path.

- max_depth:

  Integer. Maximum NUTS tree depth. Default `10`. Used only for
  `spec = "conditional"`.

- seed:

  Integer. RNG seed for the chain. Default `1L`.

- verbose:

  Logical. If `TRUE` (default), print a progress bar.

- edge_indicators:

  Optional integer \\p \times p\\ matrix with `1` = edge included, `0` =
  excluded. Must be symmetric with `1`s on the diagonal. Default: full
  graph (all edges included). Used only for `spec = "conditional"` (the
  chain samples \\K \mid \Gamma\\); ignored for `spec = "joint"`.

- delta:

  Non-negative numeric, or `NULL` for the dimension- adaptive default.
  Determinant-tilt exponent: multiplies the prior by \\\|K\|^{\delta}\\,
  softly repelling the chain from the positive-definite cone boundary.
  `delta = NULL` (default) auto-resolves to \\0.5 \log(p)\\, the simple
  form of the dimension-adaptive rule \\\delta(p) = c \log p\\ with \\c
  \in (0.3, 0.6)\\ discussed in the companion paper on
  determinant-tilted spike-and-slab priors (Marsman et al., in
  preparation). Pass `delta = 0` for the untilted prior (the
  companion-paper baseline) or a non-negative numeric to override.

- spec:

  One of `"conditional"` (default, sample \\K \mid \Gamma\\ at fixed
  \\\Gamma\\) or `"joint"` (sample \\(K, \Gamma)\\ jointly from the
  un-normalised joint prior).

- edge_inclusion_prob:

  Probability in \\(0, 1)\\ for the Bernoulli edge prior used when
  `spec = "joint"`. Default `0.5`. Ignored when `spec = "conditional"`.

## Value

A list with components

- `K_offdiag`:

  Numeric matrix of size `n_samples` x `p * (p - 1) / 2` containing the
  upper-triangle off-diagonal entries of \\K\\ for each draw, in
  column-major order \\(K\_{12}, K\_{13}, K\_{23}, K\_{14}, \ldots)\\.
  Under `spec = "conditional"`, excluded edges are returned as `0`;
  under `spec = "joint"`, off-diagonals at excluded edges are sampled at
  `0` per the inclusion indicator.

- `K_diag`:

  Numeric matrix of size `n_samples` x `p` containing the diagonal
  entries \\K\_{11}, \ldots, K\_{pp}\\.

- `offdiag_names`:

  Character vector of length `p * (p - 1) / 2` naming the columns of
  `K_offdiag` (e.g. `"K_1_2"`).

- `diag_names`:

  Character vector of length `p` naming the columns of `K_diag`.

- `edge_indicators`:

  Under `spec = "conditional"`, the `p x p` integer matrix of fixed
  inclusion indicators used (full graph if not supplied). Under
  `spec = "joint"`, an `n_samples x p(p-1)/2` integer matrix of sampled
  \\\Gamma\_{ij}\\ indicators (column order matches `K_offdiag`).

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

When `spec = "conditional"` and `edge_indicators` is supplied,
off-diagonals at excluded positions are constrained to zero throughout
the chain. `edge_indicators` is ignored when `spec = "joint"` (the chain
samples \\\Gamma\\).

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
#>            K_1_2      K_1_3      K_1_4       K_2_3      K_2_4
#> [1,] -0.99665479 -0.5067048  1.4327354  1.27206122  0.2890139
#> [2,]  1.75872847  0.5042364  0.6405899 -1.74735880  1.3637273
#> [3,] -1.85006737 -0.9258400 -0.1878773  2.41282600 -0.7835808
#> [4,] -0.06287992 -4.2151040 -1.3857240 -2.46759734  7.0181652
#> [5,]  0.70793130  4.8484312  0.8789915  1.89853902 -2.7138898
#> [6,] -2.37009472 -4.2148468  0.8678973 -0.02827887  1.8515615
#>            K_3_4
#> [1,]  0.25977772
#> [2,] -0.08441584
#> [3,] -1.37181313
#> [4,] -3.31800185
#> [5,] -0.97592686
#> [6,] -2.79041300

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
