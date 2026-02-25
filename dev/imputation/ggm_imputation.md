# Plan: Missing Data Imputation for the GGM

**Date:** 2026-02-25  
**Context:** Phase 3 Part G of the GGM cleanup (PR #78, branch `ggm_mixed`)  
**Prerequisite:** GGM centering fix (commit `564097e`) — data is centered in R before reaching C++

---

## 1. Mathematical Background

### 1.1 The GGM Likelihood

The GGM models $n$ observations $x_1, \ldots, x_n \in \mathbb{R}^p$ as i.i.d. draws from $\mathcal{N}(0, \Omega^{-1})$, where $\Omega \succ 0$ is the precision matrix. (The data is centered in R via `center_continuous_data()` before reaching C++, so the zero-mean assumption holds.)

The log-likelihood is:

$$\ell(\Omega) = \frac{n}{2} \log|\Omega| - \frac{1}{2} \text{tr}(\Omega \cdot S)$$

where $S = X'X$ is the sufficient statistic. The current C++ code stores only $S$ and discards the raw data matrix $X$.

### 1.2 Conditional Distribution for Imputation

For a missing entry $x_{iv}$ with observed values $x_{i,-v}$, the full conditional under the GGM is:

$$x_{iv} \mid x_{i,-v}, \Omega \sim \mathcal{N}\!\left(-\frac{1}{\omega_{vv}} \sum_{k \neq v} \omega_{vk} x_{ik}, \;\; \frac{1}{\omega_{vv}}\right)$$

This is the standard conditional Gaussian formula. It is already implemented in `compute_conditional_ggm()` in `src/mrf_prediction.cpp` for the prediction path.

### 1.3 Sufficient Statistic Update

When $x_{iv}$ is changed from $x_\text{old}$ to $x_\text{new}$, the sufficient statistic $S = X'X$ changes as follows. Let $\delta = x_\text{new} - x_\text{old}$. Then:

$$S'_{av} = S_{av} + \delta \cdot x_{ia} \quad \text{for all } a$$

by symmetry, $S'_{va} = S'_{av}$, and the diagonal element:

$$S'_{vv} = S_{vv} + \delta \cdot (x_{i,v,\text{new}} + x_{i,v,\text{old}}) = S_{vv} + \delta \cdot (2x_\text{old} + \delta)$$

In code, the update is a rank-1 operation on column/row $v$ of $S$. **Important:** `observations_(i, v)` must still hold $x_\text{old}$ during the loop so the diagonal gets exactly $2\delta x_\text{old}$; the $\delta^2$ correction is applied afterwards.

```
for q in 0..p:
    suf_stat_(v, q) += delta * observations_(i, q)
    suf_stat_(q, v) += delta * observations_(i, q)
// At q == v, the loop contributed 2 * delta * x_old to (v,v).
// We still need the delta^2 term:
suf_stat_(v, v) += delta * delta
// Now update the observation:
observations_(i, v) = x_new
```

Or equivalently, update `observations_` first and recompute column $v$:
```
observations_(i, v) = x_new
suf_stat_.col(v) = observations_.t() * observations_.col(v)
suf_stat_.row(v) = suf_stat_.col(v).t()
```

The first approach is $O(p)$ per imputed entry; the second is $O(np)$ per imputed entry. Since the number of missing entries is typically small relative to $n \cdot p$, the $O(p)$ incremental approach is preferred.

### 1.4 No Cholesky Recomputation Needed

The Cholesky factor and covariance matrix track the **precision matrix** $\Omega$, not the data. Changing $S$ does not invalidate these matrices. The next Metropolis-Hastings steps will automatically use the updated $S$ in their `log_density_impl_*()` evaluations. This is a key advantage: imputation has $O(p \cdot M)$ cost where $M$ is the number of missing entries, with no $O(p^2)$ or $O(p^3)$ overhead.

### 1.5 Centering Interaction

The data is centered in R (by `center_continuous_data()`) before being passed to C++. The column means used for centering are computed from the initially-imputed data (NAs filled with random draws from the observed column values). Two design considerations:

**Option A — Center once, impute on centered data:** The column means are computed on the initially-imputed data and never updated. C++ imputation works on centered data. The imputed values are centered values (mean ≈ 0), which is correct for the zero-mean GGM. Predictions and simulations use the stored `column_means` to shift back to the original scale.

**Option B — Re-center each iteration:** Recompute column means after each imputation round and re-center. This would mean `suf_stat_` always reflects the true centered cross-product. However, this adds $O(np)$ cost per iteration and complicates the pipeline (the mean becomes a moving target).

**Recommendation: Option A.** The initial centering is an approximation, but:
- For small amounts of missing data (the typical case), the initial fill is close to the true conditional mean.
- The GGM posterior on $\Omega$ does not depend on the mean (the likelihood uses $S = X'X$ with zero-mean $X$). As long as the centered imputed values have approximately zero column means, the inference is correct.
- The Gibbs conditional $\mathcal{N}(-\omega_{vv}^{-1}\omega_{v,-v}'x_{-v}, \omega_{vv}^{-1})$ has $E[\tilde{x}_{iv}] = 0$ under the zero-mean model, so column means of imputed entries do not systematically drift — they follow a zero-mean random walk with step size $O(1/n)$.
- The OMRF has the same pattern: initial random imputation, then Gibbs updates within MCMC.
- The column means used for `predict()` and `simulate()` are approximate regardless.

**Note on high missingness (>20%):** When a large fraction of a column is missing, the initial fill may bias $\hat{\mu}_v$, producing an $O(n(m - m_0)^2)$ perturbation in $S$. This is non-negligible when $\sigma_v$ is small. The plan handles this via a user-facing warning (see §3.5) rather than re-centering. Document that `column_means` is computed from the initially-filled data.

---

## 2. Design Decisions

### 2.1 C++ vs R-side Imputation

**Decision: C++ side.** 

Reasons:
- The OMRF does C++ imputation — consistency.
- The conditional Gaussian formula requires the current $\Omega$ at each iteration — this is only available in C++.
- R-side imputation would require round-tripping data and parameters between R and C++ each iteration — impractical within the current MCMC loop architecture.

### 2.2 Storing Raw Observations in GGMModel

The current `GGMModel` stores only `suf_stat_ = X'X` and discards the raw observation matrix. For imputation, we need the raw data to:
1. Compute conditional means: $\mu_{iv} = -\omega_{vv}^{-1} \sum_{k \neq v} \omega_{vk} x_{ik}$
2. Update $S$ incrementally after changing an observation.

**Decision: Add `arma::mat observations_` member variable to `GGMModel`.**

This is only populated when imputation is active. The memory cost is $n \times p \times 8$ bytes per chain — acceptable for the typical use case (e.g., $n=500, p=20$ → 80KB per chain).

### 2.3 Constructor Changes

Currently `GGMModel` has two constructors:
1. From raw observations: computes `suf_stat_ = X'X`, discards observations
2. From pre-computed sufficient statistics: no observations available

For imputation support:
- Constructor 1 needs to **conditionally store** the raw observations when imputation is active.
- Constructor 2 cannot support imputation (no raw data). This is fine — this constructor is used by `mrf_simulation.cpp` for simulation, which has no missing data.

**Decision: Add an `na_impute` flag to the from-observations constructor.** When true, store the observations matrix. When false, keep current behavior (discard observations, save memory).

Alternative considered: always store observations. Rejected because it wastes memory for the common case (no missing data).

### 2.4 Factory Function Changes

`createGGMModelFromR()` currently dispatches between the two constructors based on the presence of `"X"` vs `"n"/"suf_stat"` keys. For imputation, the factory will need to know whether to pass the `na_impute` flag. This can be inferred from the presence of `"X"` in `inputFromR` combined with the `na_impute` config flag passed from R.

**Decision: Pass `na_impute` to `createGGMModelFromR()` and through to the constructor.** The `"suf_stat"` path remains unchanged (no imputation support).

### 2.5 Imputation Within the MCMC Loop

The chain runner (`chain_runner.cpp`) already calls `model.impute_missing()` when `config.na_impute && model.has_missing_data()`. The GGM just needs to override `has_missing_data()` and `impute_missing()` — no changes to the MCMC loop.

Imputation order: **before** parameter updates, **after** `prepare_iteration()`. This is the same as the OMRF and is correct — we want the sufficient statistics to reflect the current imputed data when the MH steps evaluate the likelihood.

---

## 3. Implementation Plan

### Step G.1 — C++ `GGMModel` Changes

**File:** `src/models/ggm/ggm_model.h`

#### G.1a — New Member Variables

Add to the private section:

```cpp
// Missing data
bool has_missing_ = false;
arma::imat missing_index_;     // M × 2 matrix of 0-based (row, col) indices
arma::mat observations_;       // n × p matrix, only populated when has_missing_
```

#### G.1b — Constructor Modification

Modify the from-observations constructor to accept `na_impute`:

```cpp
GGMModel(
    const arma::mat& observations,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection = true,
    const double pairwise_scale = 2.5,
    const bool na_impute = false          // ← new parameter
) : n_(observations.n_rows),
    p_(observations.n_cols),
    dim_((p_ * (p_ + 1)) / 2),
    suf_stat_(observations.t() * observations),
    ...
    observations_(na_impute ? observations : arma::mat())   // ← conditional storage
{}
```

The from-suf_stat constructor and copy constructor need corresponding updates:
- From-suf_stat: no change needed (no observations to store)
- Copy constructor: copy `has_missing_`, `missing_index_`, `observations_`

#### G.1c — `set_missing_data()` Method

Add a public method (matching the OMRF pattern), with a guard ensuring `observations_` was retained:

```cpp
void set_missing_data(const arma::imat& missing_index) {
    BGMS_ASSERT(observations_.n_elem > 0,
        "set_missing_data() called but observations_ is empty. "
        "The model must be constructed with na_impute=true to retain observations.");
    missing_index_ = missing_index;
    has_missing_ = (missing_index.n_rows > 0 && missing_index.n_cols == 2);
}
```

#### G.1d — Override `has_missing_data()`

```cpp
bool has_missing_data() const override { return has_missing_; }
```

#### G.1e — `update_suf_stat_for_imputation()` Helper

Encapsulate the incremental `suf_stat_` update into a private helper. This localises the critical invariant that `observations_(person, variable)` must still hold $x_\text{old}$ during the update loop, preventing future refactors from accidentally breaking the ordering.

Add to `ggm_model.h` (private):

```cpp
void update_suf_stat_for_imputation(int variable, int person, double delta);
```

Add to `ggm_model.cpp`:

```cpp
void GGMModel::update_suf_stat_for_imputation(int variable, int person, double delta) {
    // INVARIANT: observations_(person, variable) must still hold x_old when
    // this function is called. The loop adds 2 * delta * x_old to the (v,v)
    // entry; the delta^2 correction completes the diagonal update.
    for (size_t q = 0; q < p_; q++) {
        suf_stat_(variable, q) += delta * observations_(person, q);
        suf_stat_(q, variable) += delta * observations_(person, q);
    }
    suf_stat_(variable, variable) += delta * delta;
}
```

#### G.1f — `impute_missing()` Implementation

Add to `ggm_model.cpp`:

```cpp
void GGMModel::impute_missing() {
    if (!has_missing_) return;
    
    const int num_missings = missing_index_.n_rows;
    
    for (int miss = 0; miss < num_missings; miss++) {
        const int person = missing_index_(miss, 0);
        const int variable = missing_index_(miss, 1);
        
        // Compute conditional mean: mu = -sum_{k != v} omega_{vk} * x_{ik} / omega_{vv}
        double conditional_mean = 0.0;
        for (size_t k = 0; k < p_; k++) {
            if (k != static_cast<size_t>(variable)) {
                conditional_mean += precision_matrix_(variable, k) * observations_(person, k);
            }
        }
        conditional_mean = -conditional_mean / precision_matrix_(variable, variable);
        
        // Conditional variance: 1 / omega_{vv}
        double conditional_sd = std::sqrt(1.0 / precision_matrix_(variable, variable));
        
        // Sample new value
        double x_new = rnorm(rng_, conditional_mean, conditional_sd);
        double x_old = observations_(person, variable);
        double delta = x_new - x_old;
        
        // Incrementally update suf_stat_ (observations_ still holds x_old)
        update_suf_stat_for_imputation(variable, person, delta);
        
        // Now update the observation
        observations_(person, variable) = x_new;
    }
    
    // Full recompute at end of sweep to eliminate floating-point drift
    // (matches OMRF pattern; cost is O(np^2), negligible for typical sizes)
    suf_stat_ = observations_.t() * observations_;
}
```

**Correctness check for the `suf_stat_` update:**

$S_{vq} = \sum_i x_{iv} x_{iq}$. After changing $x_{kv}$ from $x_\text{old}$ to $x_\text{new}$:

$S'_{vq} = S_{vq} + (x_\text{new} - x_\text{old}) \cdot x_{kq}$ for $q \neq v$

$S'_{vv} = S_{vv} + x_\text{new}^2 - x_\text{old}^2 = S_{vv} + \delta(\delta + 2x_\text{old})$

In the loop above, `suf_stat_(variable, variable)` gets two additions of `delta * observations_(person, variable)` (once from q=variable, once from the symmetric line). But `observations_(person, variable)` is still $x_\text{old}$ at this point, so the two additions contribute $2 \cdot \delta \cdot x_\text{old}$. Then the correction line adds $\delta^2$. Total change to $(v,v)$: $2\delta x_\text{old} + \delta^2 = \delta(\delta + 2x_\text{old}) = x_\text{new}^2 - x_\text{old}^2$. ✓

**Note:** The loop updates `suf_stat_(q, variable)` and `suf_stat_(variable, q)` separately (two writes per $q$). Since `suf_stat_` is symmetric, this is correct and avoids branching. An alternative is to update only the upper triangle and mirror, but the current approach is simpler and the $O(p)$ cost is negligible compared to the $O(p^4)$ MH step.

---

### Step G.2 — C++ Interface Changes

**File:** `src/sample_ggm.cpp`

#### G.2a — Add Parameters to `sample_ggm()`

Add `na_impute` and `missing_index_nullable` **at the end** of the parameter list (after all existing parameters) for clean convention:

```cpp
// [[Rcpp::export]]
Rcpp::List sample_ggm(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const int no_iter,
    const int no_warmup,
    const int no_chains,
    const bool edge_selection,
    const int seed,
    const int no_threads,
    const int progress_type,
    const std::string& edge_prior = "Bernoulli",
    const double beta_bernoulli_alpha = 1.0,
    ...
    const bool na_impute = false,                                          // ← new (at end)
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_nullable = R_NilValue   // ← new (at end)
)
```

#### G.2b — Factory Call and Missing Data Setup

```cpp
// Create model from R input (pass na_impute for conditional observation storage)
GGMModel model = createGGMModelFromR(
    inputFromR, prior_inclusion_prob, initial_edge_indicators, 
    edge_selection, 2.5, na_impute);                           // ← pass na_impute

// Set up missing data imputation (same pattern as OMRF)
if (na_impute && missing_index_nullable.isNotNull()) {
    arma::imat missing_index = Rcpp::as<arma::imat>(
        Rcpp::IntegerMatrix(missing_index_nullable.get()));
    model.set_missing_data(missing_index);
}
```

#### G.2c — SamplerConfig

```cpp
config.na_impute = na_impute;
```

#### G.2d — Factory Function Update

**File:** `src/models/ggm/ggm_model.cpp` (or `.h`)

Update `createGGMModelFromR()` to pass `na_impute` to the from-observations constructor:

```cpp
GGMModel createGGMModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection,
    const double pairwise_scale,
    const bool na_impute = false              // ← new parameter
) {
    if (inputFromR.containsElementNamed("X")) {
        arma::mat X = Rcpp::as<arma::mat>(inputFromR["X"]);
        return GGMModel(X, inclusion_probability, initial_edge_indicators, 
                        edge_selection, pairwise_scale, na_impute);  // ← pass through
    } else if (...) {
        // suf_stat path — no change, no imputation possible
    }
}
```

---

### Step G.3 — R-side Plumbing

#### G.3a — Remove the Guard and Add Validation in `validate_missing_data()`

**File:** `R/validate_data.R`

Remove the `stop()` call at lines 97–101 and add two new guards:

```r
# REMOVE:
if (is_continuous && na_action == "impute") {
    stop("Imputation is not yet supported for the Gaussian model. ",
         "Use na_action = 'listwise'.")
}

# ADD (in handle_impute or validate_missing_data):

# Guard: entire-column-missing
for (v in seq_len(ncol(x))) {
    if (all(is.na(x[, v]))) {
        stop("Variable '", colnames(x)[v], "' has no observed values. ",
             "Remove it before fitting.")
    }
}

# Warning: high missingness
pct_missing <- sum(is.na(x)) / length(x)
if (pct_missing > 0.20) {
    warning(round(pct_missing * 100), "% of entries are missing. ",
            "Imputation assumes MCAR/MAR and quality degrades with high missingness.")
}
```

#### G.3b — Pass Missing Data to C++

**File:** `R/run_sampler.R`

In `run_sampler_ggm()`, add the missing data arguments:

```r
run_sampler_ggm <- function(spec) {
  d <- spec$data
  p <- spec$prior
  s <- spec$sampler
  m <- spec$missing                           # ← new: extract missing sub-list
  
  ...
  
  out_raw <- sample_ggm(
    inputFromR               = list(X = d$x),
    ...
    edge_prior               = p$edge_prior,
    na_impute                = m$na_impute,     # ← new
    missing_index_nullable   = m$missing_index, # ← new
    beta_bernoulli_alpha     = p$beta_bernoulli_alpha,
    ...
  )
  
  out_raw
}
```

#### G.3c — Regenerate `RcppExports`

After modifying `sample_ggm.cpp`, run `Rcpp::compileAttributes()` to regenerate `R/RcppExports.R` and `src/RcppExports.cpp`.

---

### Step G.4 — Tests

**File:** `tests/testthat/test-bgm.R`

#### G.4a — Basic Imputation Smoke Test

```r
test_that("GGM with na_action='impute' runs without error", {
  set.seed(1)
  n <- 100; p <- 5
  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("V", 1:p)
  
  # Introduce ~5% missing data
  missing_mask <- sample(length(x), size = round(0.05 * length(x)))
  x[missing_mask] <- NA
  
  fit <- bgm(x, iter = 200, na_action = "impute")
  
  expect_s3_class(fit, "bgms")
  expect_true(fit$arguments$na_impute)
  expect_equal(nrow(fit$posterior_samples$pairwise), 200)
})
```

#### G.4b — Posterior Unbiasedness Test

```r
test_that("GGM imputation preserves posterior accuracy", {
  skip_on_cran()
  set.seed(42)
  
  p <- 6; n <- 300
  # Known sparse precision matrix
  Omega_true <- diag(p)
  Omega_true[1,2] <- Omega_true[2,1] <- 0.4
  Omega_true[3,4] <- Omega_true[4,3] <- -0.3
  
  Sigma <- solve(Omega_true)
  x_full <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(x_full) <- paste0("V", 1:p)
  
  # Fit on complete data
  fit_full <- bgm(x_full, iter = 2000, edge_selection = FALSE)
  
  # Introduce 5% MCAR missing data
  x_miss <- x_full
  miss_idx <- sample(length(x_miss), size = round(0.05 * length(x_miss)))
  x_miss[miss_idx] <- NA
  
  fit_miss <- bgm(x_miss, iter = 2000, edge_selection = FALSE, na_action = "impute")
  
  # Posterior means should be correlated > 0.85
  cor_pairwise <- cor(
    as.numeric(fit_full$posterior_mean_pairwise),
    as.numeric(fit_miss$posterior_mean_pairwise)
  )
  expect_gt(cor_pairwise, 0.85)
})
```

#### G.4c — Listwise vs Impute Comparison

```r
test_that("GGM imputation gives comparable results to listwise", {
  skip_on_cran()
  set.seed(123)
  
  p <- 5; n <- 200
  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("V", 1:p)
  
  # Introduce 3% missing
  miss_idx <- sample(length(x), size = round(0.03 * length(x)))
  x[miss_idx] <- NA
  
  fit_listwise <- bgm(x, iter = 500, na_action = "listwise", edge_selection = FALSE)
  fit_impute   <- bgm(x, iter = 500, na_action = "impute",   edge_selection = FALSE)
  
  # Both should produce valid output
  expect_s3_class(fit_listwise, "bgms")
  expect_s3_class(fit_impute, "bgms")
  
  # Posterior means should be in the same ballpark
  cor_val <- cor(
    as.numeric(fit_listwise$posterior_mean_pairwise),
    as.numeric(fit_impute$posterior_mean_pairwise)
  )
  expect_gt(cor_val, 0.80)
})
```

#### G.4d — Edge Cases

```r
test_that("GGM imputation handles edge cases", {
  set.seed(1)
  n <- 50; p <- 4
  x <- matrix(rnorm(n * p), n, p)
  colnames(x) <- paste0("V", 1:p)
  
  # Single missing value
  x[1, 1] <- NA
  fit <- bgm(x, iter = 100, na_action = "impute")
  expect_s3_class(fit, "bgms")
  
  # Multiple missing in same row
  x[2, 1:3] <- NA
  fit2 <- bgm(x, iter = 100, na_action = "impute")
  expect_s3_class(fit2, "bgms")
})
```

#### G.4e — Column Means Stored Correctly

```r
test_that("GGM imputation stores column_means correctly", {
  set.seed(1)
  n <- 80; p <- 4
  x <- matrix(rnorm(n * p, mean = 5), n, p)  # non-zero means
  colnames(x) <- paste0("V", 1:p)
  x[1:3, 1] <- NA
  
  fit <- bgm(x, iter = 100, na_action = "impute")
  
  # column_means should be numeric of length p
  expect_true(is.numeric(fit$arguments$column_means))
  expect_length(fit$arguments$column_means, p)
  # means should be approximately 5
  expect_true(all(abs(fit$arguments$column_means - 5) < 2))
})
```

#### G.4f — Posterior Samples Remain Finite (Indirect `suf_stat_` Consistency Check)

```r
test_that("GGM impute: posterior samples remain finite and bounded", {
  set.seed(1)
  x <- matrix(rnorm(500), 100, 5)
  colnames(x) <- paste0("V", 1:5)
  x[sample(500, 25)] <- NA
  fit <- bgm(x, iter = 500, na_action = "impute", edge_selection = FALSE)
  samples <- fit$posterior_samples$pairwise
  expect_true(all(is.finite(samples)))
  expect_true(all(abs(samples) < 100))  # no blow-ups
})
```

#### G.4g — No-NA Data with `na_action = "impute"`

```r
test_that("GGM with na_action='impute' but no NAs works transparently", {
  x <- matrix(rnorm(100), 20, 5)
  colnames(x) <- paste0("V", 1:5)
  fit <- bgm(x, iter = 100, na_action = "impute")
  expect_s3_class(fit, "bgms")
})
```

#### G.4h — Predict and Simulate After Imputed Fit

```r
test_that("predict and simulate work after GGM imputation fit", {
  set.seed(1)
  x <- matrix(rnorm(200, mean = 3), 40, 5)
  colnames(x) <- paste0("V", 1:5)
  x[1:2, 1] <- NA
  fit <- bgm(x, iter = 100, na_action = "impute")
  
  pred <- predict(fit, newdata = x[3:5, , drop = FALSE])
  expect_true(is.list(pred))
  
  sim <- simulate(fit, nsim = 10)
  expect_true(is.matrix(sim) || is.data.frame(sim))
})
```

#### G.4i — Multi-Chain with Imputation

```r
test_that("GGM imputation works with multiple chains", {
  set.seed(1)
  x <- matrix(rnorm(200), 40, 5)
  colnames(x) <- paste0("V", 1:5)
  x[1:3, 1] <- NA
  fit <- bgm(x, iter = 100, chains = 2, na_action = "impute")
  expect_s3_class(fit, "bgms")
})
```

#### G.4j — Failure Modes

```r
test_that("GGM impute: entire-column-missing gives clear error", {
  x <- matrix(rnorm(100), 20, 5)
  colnames(x) <- paste0("V", 1:5)
  x[, 3] <- NA
  expect_error(bgm(x, iter = 100, na_action = "impute"), "no observed values")
})
```

---

## 4. Files Changed Summary

| File | Change Type | Description |
|------|------------|-------------|
| `src/models/ggm/ggm_model.h` | Modify | Add `observations_`, `has_missing_`, `missing_index_` members; add `na_impute` constructor param; override `has_missing_data()`, add `set_missing_data()` declaration, `update_suf_stat_for_imputation()` (private), `impute_missing()` declaration; update copy constructor |
| `src/models/ggm/ggm_model.cpp` | Modify | Add `update_suf_stat_for_imputation()` and `impute_missing()` implementations; update `createGGMModelFromR()` signature |
| `src/sample_ggm.cpp` | Modify | Add `na_impute`, `missing_index_nullable` parameters; add missing data setup; set `config.na_impute` |
| `R/validate_data.R` | Modify | Remove the `is_continuous && na_action == "impute"` guard |
| `R/run_sampler.R` | Modify | Pass `m$na_impute` and `m$missing_index` to `sample_ggm()` |
| `R/RcppExports.R` | Auto-regen | Via `Rcpp::compileAttributes()` |
| `src/RcppExports.cpp` | Auto-regen | Via `Rcpp::compileAttributes()` |
| `tests/testthat/test-bgm.R` | Modify | Add imputation tests (G.4a–G.4j) |
| `man/bgm.Rd` (or roxygen in `R/bgm.R`) | Modify | Document center-once imputation approach in `na_action` parameter docs |

**Files NOT changed:**
- `R/bgm_spec.R` — already assembles the `missing` sub-list correctly for GGM
- `R/build_output.R` — shared builder, no imputation-specific output needed
- `R/simulate_predict.R` — `column_means` handling works as-is
- `src/mcmc/execution/chain_runner.cpp` — already calls `impute_missing()` when configured
- `src/mcmc/execution/sampler_config.h` — already has `na_impute` field

---

## 5. Execution Order

```
G.1  C++ GGMModel changes (header + implementation)
 ↓
G.2  C++ sample_ggm.cpp interface + factory function
 ↓
G.3a Remove R guard in validate_data.R + add guards/warning
G.3b Pass missing data args in run_sampler.R
G.3c Regenerate RcppExports
G.3d Document centering approach in bgm() help page
 ↓
     Compile and verify: devtools::load_all() works
 ↓
G.4  Tests (add to test-bgm.R, run devtools::test())
```

Steps G.3a, G.3b, G.3c can be done in parallel. Everything else is sequential.

---

## 6. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Sufficient statistic update has an off-by-one or double-counting bug | Medium | High | The math is verified algebraically in §1.3. Add a unit test that manually computes $X'X$ after imputation and compares against the incremental update. |
| Initial centering on imputed data introduces bias | Low | Low | For small missingness the initial fill is close to truth. The posterior is invariant to the centering mean — it only affects prediction/simulation, where the approximation is acceptable. |
| Memory overhead from storing `observations_` | Low | Low | Only stored when `na_impute = true`. For typical sizes ($n=500, p=20$), it's 80KB per chain. |
| Multiple missing values in the same row interact during imputation | Low | Low | Each missing value is imputed from its full conditional given **current** values (including already-imputed values in the same row from this MCMC step). This is standard "sequential Gibbs" and is valid. The OMRF does the same. |
| Large fraction of missing data degrades posterior quality | Medium | Medium | This is inherent to any imputation scheme. Document that imputation assumes MCAR/MAR and works best with <20% missingness. Out of scope for this implementation. |
| Numerical issues: imputed values far from zero cause large `suf_stat_` entries | Low | Medium | The conditional Gaussian samples are well-behaved (bounded variance $1/\omega_{vv}$). The centering ensures values are around zero. |

---

## 7. Design Alternatives Considered

### 7.1 Full Recomputation vs Incremental Update of `suf_stat_`

**Full recomputation:** After all missing values are imputed, recompute $S = X'X$ from scratch. Cost: $O(np^2)$.

**Incremental:** Update column $v$ of $S$ after each imputation. Cost: $O(p)$ per imputed entry, total $O(pM)$.

**Decision: Hybrid (incremental + end-of-sweep recompute).** Use incremental $O(p)$ updates during the loop (required for correct conditional means of later entries in the same sweep), then add a full `suf_stat_ = observations_.t() * observations_` recompute at the end of `impute_missing()`. This matches the OMRF pattern exactly (which does incremental residual updates per entry, then full recomputation of `pairwise_stats_` at the end). The full recompute is one $O(np^2)$ matrix multiply per iteration — negligible for typical sizes ($n=500, p=20$: ~0.2M flops) — and eliminates any theoretical concern about accumulated floating-point drift.

### 7.2 Block vs Sequential Imputation

**Block:** Sample all missing values simultaneously from the joint conditional. For the GGM, this requires computing the conditional covariance of all missing entries given all observed entries — a potentially large matrix factorization.

**Sequential (chosen):** Sample each missing entry from its univariate full conditional. Standard Gibbs strategy. Simpler, $O(p)$ per entry, and the OMRF does the same.

### 7.3 Storing Imputed Values in Output

The OMRF does not store the final imputed dataset in the output object. We follow the same pattern. If users want the imputed values, that would be a separate feature (e.g., storing imputed values at each MCMC iteration for multiple imputation). This is out of scope.

---

## 8. Validation Strategy

After implementation, verify correctness with:

1. **Manual `suf_stat_` check:** Create data with known missing values, impute manually, compute $X'X$, compare against the incremental update.

2. **Posterior recovery:** Fit on complete data, introduce MCAR missingness, fit with imputation, compare posterior means (expect correlation > 0.85 for 5% missing, > 0.70 for 15% missing).

3. **Conditional distribution check:** With a known $\Omega$ and fixed data, verify that the `impute_missing()` samples from the correct $\mathcal{N}(-\omega_{vv}^{-1}\omega_{v,-v}'x_{-v}, \omega_{vv}^{-1})$ by running many imputation steps and checking the empirical distribution.

4. **Regression against listwise:** For small missingness, imputation and listwise should give similar posteriors.

---

## 9. Review Log

### Changes Applied from Reviews

The following items were agreed upon by both reviewers and have been incorporated:

1. **§1.3 pseudocode fixed** — The diagonal correction now uses `+= delta * delta` (matching G.1e) instead of the incorrect subtraction. The dependency on `observations_` still holding $x_\text{old}$ during the loop is explicitly documented.

2. **End-of-sweep full recompute added** to `impute_missing()` (§G.1e, §7.1) — `suf_stat_ = observations_.t() * observations_` after the loop. Matches the OMRF pattern, eliminates floating-point drift concerns, negligible cost.

3. **`set_missing_data()` guard** (§G.1c) — Asserts that `observations_` was retained before enabling imputation. Prevents silent corruption if called after the sufficient-statistics constructor.

4. **Entire-column-missing guard** (§G.3a) — Explicit error before `sample()` crashes on an empty vector.

5. **High-missingness warning** (§G.3a) — Warns when >20% of entries are missing.

6. **Additional tests** (§G.4f–G.4j) added: finite/bounded posterior check, no-NA impute path, predict/simulate after imputation, multi-chain, failure modes.

7. **Option A (center once)** confirmed as the correct choice by both reviewers. The centering error is $O(1/\sqrt{n})$, does not drift systematically, and re-centering would add cost and conceptual problems.

### Decisions on Open Discussion Points

1. **D1 — Centering drift monitoring:** No runtime monitoring. The full end-of-sweep `suf_stat_` recompute and the >20% missingness warning are sufficient. The center-once approach should be documented in the `bgm()` help page under `na_action = "impute"` (see §G.3d below).

2. **D2 — Encapsulate `suf_stat_` update:** Yes. A private `update_suf_stat_for_imputation()` helper method has been added (§G.1e) to localise the ordering invariant.

3. **D3 — Parameter ordering:** `na_impute` and `missing_index_nullable` moved to the end of the `sample_ggm()` signature (§G.2a).

4. **D4 — Stochastic test stability:** G.4b iterations increased from 1000 to 2000 for more stable correlation estimates.

5. **D5 — Centering drift regression test:** Skipped. The existing posterior accuracy tests (G.4b) and the high-missingness warning cover this indirectly.

### Step G.3d — Documentation

**File:** `man/bgm.Rd` (or roxygen in `R/bgm.R`)

Add to the `na_action` parameter documentation:

> When `na_action = "impute"`, missing values are imputed within the MCMC loop using Gibbs sampling from the conditional Gaussian distribution. The data is centered once using the column means of the initially-filled data (NAs replaced by random draws from observed values in the same column). These column means are stored and used by `predict()` and `simulate()` to map back to the original scale. The centering is not updated during MCMC; this is valid because the Gibbs conditionals preserve the zero-mean property in expectation and the centering error is $O(1/\sqrt{n})$. Imputation assumes MCAR or MAR missingness and works best with less than 20% missing entries.
