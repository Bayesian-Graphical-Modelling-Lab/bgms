# Review: Missing Data Imputation for the GGM

**Reviewer:** GitHub Copilot  
**Date:** 2026-02-25  
**Document reviewed:** `dev/plans/ggm_imputation.md`  
**Branch:** `ggm_mixed` (PR #78)

---

## 1. Mathematical Correctness

### 1.1 The `suf_stat_` Incremental Update

**§G.1e is correct. §1.3 pseudocode has a bug.**

The target update for the diagonal entry is:

$$S'_{vv} = S_{vv} + x_\text{new}^2 - x_\text{old}^2 = S_{vv} + \delta(2x_\text{old} + \delta)$$

The G.1e code does:

1. Loop over all $q$: add `delta * x_old_q` to both `suf_stat_(v, q)` and `suf_stat_(q, v)`.
2. When `q == variable`, this contributes $2 \delta x_\text{old}$ to `suf_stat_(v, v)` (two writes, `observations_` still holds $x_\text{old}$).
3. The correction line `suf_stat_(v, v) += delta * delta` adds $\delta^2$.
4. Total change to $(v,v)$: $2\delta x_\text{old} + \delta^2 = x_\text{new}^2 - x_\text{old}^2$. ✓

The §1.3 pseudocode, by contrast, writes:

```
suf_stat_(v, v) -= delta * observations_(i, v)   // undo one of the two additions
```

This would leave the diagonal change at $2\delta x_\text{old} - \delta x_\text{old} = \delta x_\text{old}$, which is **wrong** (missing the $\delta^2$ term). Either fix §1.3 to match G.1e (replace the subtraction with `+= delta * delta`), or remove the §1.3 pseudocode entirely and point to G.1e as canonical.

### 1.2 Conditional Distribution

The formula in §1.2 is the standard conditional Gaussian:

$$x_{iv} \mid x_{i,-v}, \Omega \sim \mathcal{N}\!\left(-\frac{\sum_{k \neq v} \omega_{vk} x_{ik}}{\omega_{vv}}, \;\; \frac{1}{\omega_{vv}}\right)$$

This matches `compute_conditional_ggm()` in `src/mrf_prediction.cpp`, providing independent validation. Correct.

### 1.3 Sequential Gibbs Ordering

Imputing each missing entry from its univariate full conditional, given current values (including earlier imputed entries in the same sweep), is valid single-site Gibbs sampling. The precision matrix $\Omega$ does not change during the imputation sweep — only `observations_` and `suf_stat_` are updated. The next Metropolis–Hastings step then correctly uses the updated `suf_stat_` for its likelihood evaluations. This matches the OMRF pattern and is correct.

### 1.4 Cholesky Invariance

Correct: the Cholesky factor and covariance matrix track $\Omega$, not the data. Changing $S$ does not invalidate them. No recomputation needed.

---

## 2. Centering Interaction — Deep Analysis

This is the most subtle aspect of the plan. The question is whether centering once (Option A) introduces bias that accumulates during MCMC.

### 2.1 The Pipeline

The R-side pipeline is:

1. Fill NAs with random draws from observed column values (`handle_impute()`)
2. Compute column means from filled data: $\hat{\mu}_v = \bar{x}_v^{(\text{filled})}$
3. Subtract: $\tilde{x}_{iv} = x_{iv} - \hat{\mu}_v$
4. Pass $\tilde{X}$ and `column_means` $= \hat{\mu}$ to C++

C++ then operates entirely on $\tilde{X}$, with the GGM likelihood assuming $\tilde{X} \sim \mathcal{N}(0, \Omega^{-1})$.

### 2.2 Does It Matter? — No, for Two Distinct Reasons

**Reason 1: The centering error is small and static.**

The centering mean $\hat{\mu}_v$ differs from the true population mean $\mu_v$ by:

$$\hat{\mu}_v - \mu_v = \underbrace{\frac{1}{n}\sum_{i \in \text{obs}} (x_{iv} - \mu_v)}_{O(1/\sqrt{n})} + \underbrace{\frac{1}{n}\sum_{i \in \text{mis}} (x_{iv}^{\text{fill}} - \mu_v)}_{O(\sqrt{M}/n)}$$

For $n = 200$ and 5% missingness ($M = 10p$), this is $O(0.07)$ — negligible compared to the data's standard deviation ($\sigma_v \sim 1$). The effect on $S$ is:

$$S^{\text{biased}}_{vv} = S^{\text{true}}_{vv} + n(\hat{\mu}_v - \mu_v)^2 \approx S^{\text{true}}_{vv} + O(1)$$

Since $S^{\text{true}}_{vv} \sim n\sigma_v^2 \sim n$, this is an $O(1/n)$ relative perturbation — comparable to the familiar $n$ vs $n-1$ correction in sample covariance. It does not compound over iterations because $\hat{\mu}$ is fixed once.

**Reason 2: The Gibbs imputation preserves zero mean in expectation.**

The conditional draw $\tilde{x}_{kv}^{\text{new}} \sim \mathcal{N}(-\omega_{vv}^{-1}\sum_{j\neq v}\omega_{vj}\tilde{x}_{kj}, \; \omega_{vv}^{-1})$ is a linear function of centered data. Under the stationary distribution of the joint chain $(\tilde{X}_\text{mis}, \Omega)$, the marginal expectation $E[\tilde{x}_{kv}]$ is zero (by the symmetry of the zero-mean Gaussian model). So the column means of imputed entries fluctuate around zero — there is **no systematic drift**.

More precisely: at each iteration, changing one imputed value shifts the column mean by $\delta/n$. Since $E[\delta] = 0$ (the conditional distribution is centered at the correct conditional mean under the zero-mean model), the column mean follows a zero-mean random walk with step size $O(1/n)$. Over $T$ iterations, the deviation is $O(\sqrt{T}/n)$. For $T = 1000, n = 200$, this is $\sim 0.16$ — still small relative to $\sigma_v$, and it doesn't accumulate (it's a random walk, not a drift).

### 2.3 Would Re-centering Help?

No — re-centering would actually be worse:

1. **Cost:** $O(np)$ for recomputing column means + $O(np^2)$ for recomputing $S$ from scratch, every iteration.

2. **Moving-target problem:** The stored `column_means` used by `predict()` and `simulate()` to map back to the original scale would depend on the MCMC iteration. Which one do you store? The last? The average? This complicates the interface with no real benefit.

3. **Conceptual muddle:** Re-centering each iteration implicitly treats the mean $\mu$ as a parameter that is re-estimated each step, but without a proper prior or update rule. This is neither full Bayesian (no prior on $\mu$) nor consistent with the current framework (which assumes $\mu$ is known). It's a hybrid that's harder to justify theoretically than "center once."

4. **The OMRF analogy:** The OMRF fills missing values initially, then imputes within MCMC, without re-adjusting any global statistics (other than incremental `residual_matrix_` and `pairwise_stats_` updates). The "center once" approach matches this philosophy.

### 2.4 The One Edge Case Where Centering Matters

If missingness is very high (say >40% in some column), the initial random fill could substantially bias $\hat{\mu}_v$. Example: true $\mu_v = 5$, but 50% of values are missing and filled from the observed half, which has a within-sample mean of, say, 4.5. Then $\hat{\mu}_v \approx 4.75$, and all centered values for column $v$ are shifted by $0.25$. This is still absorbed into the sufficient statistic as an $O(n \cdot 0.25^2) \approx 0.06n$ perturbation — small relative to $S_{vv} \sim n\sigma_v^2$, but non-negligible if $\sigma_v$ is small.

**Recommendation:** This doesn't warrant re-centering, but it's another reason to issue a warning for high missingness (see §4.2 below).

### 2.5 Summary

**Option A (center once) is sound.** The centering error is $O(1/\sqrt{n})$, does not drift during MCMC, and its effect on the precision matrix posterior is negligible. Re-centering would add complexity, cost, and conceptual problems. The only documentation action is noting that `column_means` is computed from the initially-filled data.

---

## 3. Architectural Fit

### 3.1 What Works Well

- **Conditional `observations_` storage** is the right design. Memory is only allocated when imputation is active (the common case has no missing data).
- **`na_impute` flag threading** (constructor → factory → interface) mirrors the OMRF exactly. The chain runner already calls `impute_missing()` when `config.na_impute && model.has_missing_data()` — no loop changes needed.
- **`clone()` / copy constructor** — critical for multi-chain. Armadillo matrices deep-copy automatically, so as long as `observations_`, `missing_index_`, and `has_missing_` are in the copy constructor initializer list (plan says they will be), each chain gets independent state. ✓
- **Virtual dispatch via `BaseModel`** keeps the MCMC loop model-agnostic. No `if (model_type == "ggm")` branches anywhere.

### 3.2 Minor Suggestions

**Parameter ordering in `sample_ggm()`.** §G.2a inserts `na_impute` and `missing_index_nullable` in the middle of the parameter list. Since R calls use named arguments, this works, but adding new parameters at the **end** of the signature is cleaner convention. Consider moving them after the existing edge-prior parameters.

**Two-phase initialization.** The plan uses a two-step setup: construct model, then call `set_missing_data()`. An alternative is passing `missing_index` to the constructor when `na_impute = true`, making the object fully initialized at construction. However, the OMRF uses the same two-step pattern, so consistency argues for keeping it. Not a correctness issue.

**`createGGMModelFromR()` signature.** Adding `na_impute` as a parameter with default `false` is clean. The `suf_stat` path (used by simulation) correctly cannot support imputation since it has no raw data — good separation.

---

## 4. Missing Design Considerations

### 4.1 Entire-column-missing Guard

If all values in a column are `NA`, `handle_impute()` calls `sample(x[-mis, node], size = 1)` where `x[-mis, node]` is empty. This will error with a confusing R message. Add an explicit check:

```r
if (length(mis) == nrow(x)) {
  stop("Variable '", colnames(x)[node], "' has no observed values. ",
       "Remove it before fitting.")
}
```

### 4.2 High-missingness Warning

The plan acknowledges this risk (§6 risk table) but proposes no user-facing action. Add a warning in `handle_impute()`:

```r
pct_missing <- num_missings / length(x)
if (pct_missing > 0.20) {
  warning(round(pct_missing * 100), "% of entries are missing. ",
          "Imputation assumes MCAR/MAR and quality degrades with high missingness.")
}
```

### 4.3 Floating-point Drift in Incremental Updates

The plan discusses this in §7.1 and correctly notes the GGM has lower drift risk than the OMRF (well-scaled floating-point arithmetic, not integers). However, the OMRF includes a **full recomputation** at the end of each `impute_missing()` sweep:

```cpp
// Recompute pairwise sufficient statistics
arma::mat ps = observations_double_.t() * observations_double_;
pairwise_stats_ = arma::conv_to<arma::imat>::from(ps);
```

Adding the analogous safeguard for the GGM is cheap:

```cpp
// At end of impute_missing():
suf_stat_ = observations_.t() * observations_;
```

This is one $O(np^2)$ matrix multiply per iteration — negligible for typical sizes ($n = 500, p = 20$: $\sim 0.2$M flops). It eliminates any theoretical concern about accumulated rounding errors and makes the implementation more robust. Recommended.

### 4.4 MNAR / Non-ignorable Missing Data

Correctly out of scope. The documentation for `na_action = "impute"` should state that the procedure assumes MCAR or MAR missingness. This is a documentation-only action; no code change needed.

### 4.5 Numerical Issues with Extreme Conditional Distributions

If $\omega_{vv}$ is very small (near the positive-definiteness boundary), the conditional variance $1/\omega_{vv}$ is very large, and imputed values could be extreme. This would inflate $S_{vv}$ and potentially destabilize subsequent Metropolis–Hastings steps. In practice, the prior (Cauchy slab + positive-definiteness constraint) prevents $\omega_{vv}$ from collapsing to near-zero, so this is very unlikely. No code change needed, but worth noting as a theoretical edge case.

---

## 5. Test Coverage

### 5.1 Proposed Tests (G.4a–G.4e) — Assessment

| Test | Coverage | Verdict |
|------|----------|---------|
| G.4a — Smoke test | Runs without error, correct class and dimensions | ✓ Good |
| G.4b — Posterior unbiasedness | Compares complete vs imputed posterior means | ✓ Good, core correctness check |
| G.4c — Listwise vs impute | Checks the two approaches give similar results | ✓ Good for low missingness |
| G.4d — Edge cases | Single missing value, multiple in same row | ✓ Good |
| G.4e — Column means stored | Verifies `column_means` correct after imputation | ✓ Good for predict/simulate path |

### 5.2 Missing Tests — Recommended Additions

**Test: Deterministic `suf_stat_` verification** (HIGH PRIORITY)

This is described in §8 (validation strategy, item 1) but not included in §G.4. It's the single most valuable test for catching incremental update bugs. The idea: after fitting with imputation, verify that `suf_stat_` equals the cross-product of the current `observations_` matrix. This requires exposing internals via a C++ test helper. Alternatively, check indirectly by verifying that log-likelihood evaluations produce finite values and that the Metropolis acceptance rate is reasonable.

```r
test_that("GGM impute: suf_stat remains consistent with observations", {
  # After N iterations of imputation, if suf_stat != X'X,
  # the MH acceptance rate degrades and posteriors are wrong.
  # Indirect check: posterior samples should be finite and 
  # acceptance diagnostics should not flag issues.
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

**Test: No-NA path with `na_action = "impute"`**

When `na_action = "impute"` but the data has no NAs, `handle_impute()` returns `na_impute = FALSE` and the imputation machinery is never activated. This should work transparently:

```r
test_that("GGM with na_action='impute' but no NAs works", {
  x <- matrix(rnorm(100), 20, 5)
  colnames(x) <- paste0("V", 1:5)
  fit <- bgm(x, iter = 100, na_action = "impute")
  expect_s3_class(fit, "bgms")
})
```

**Test: predict/simulate after imputed fit**

Ensures `column_means` are correctly stored and used for the inverse mapping:

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

**Test: Multi-chain with imputation**

Validates that `clone()` correctly deep-copies missing data state:

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

### 5.3 Stochastic Test Stability

Tests G.4b and G.4c use correlation thresholds (0.85 and 0.80) to compare posterior means. Even with `skip_on_cran()` and fixed seeds, these can be fragile across platforms/compilers due to floating-point nondeterminism. Consider:

- Using a slightly lower threshold (0.75) for robustness, or
- Running more iterations (2000+) and reducing the missingness fraction to 3%, or
- Testing the sign pattern of the posterior mean matrix (which edges are positive/negative) instead of correlations — this is a weaker but more stable check.

---

## 6. Summary of Recommendations

| Priority | Item | Action |
|----------|------|--------|
| **High** | §1.3 pseudocode bug | Fix the diagonal correction or remove the sketch, point to G.1e |
| **High** | Add `suf_stat_` consistency test | Either expose internals or add indirect finite-value checks |
| **Medium** | Full `suf_stat_` recompute at sweep end | Add `suf_stat_ = observations_.t() * observations_` at end of `impute_missing()` (matches OMRF pattern, cheap insurance) |
| **Medium** | Entire-column-missing guard | Error before `sample()` crashes on empty vector |
| **Medium** | High-missingness warning | Warn user when >20% entries are missing |
| **Low** | Parameter ordering | Move `na_impute`/`missing_index` to end of `sample_ggm()` signature |
| **Low** | Additional tests | No-NA impute path, predict/simulate after imputation, multi-chain |
| **Low** | MCAR/MAR documentation | Note assumption in `bgm()` help page |

### Centering Verdict

**Option A (center once, never re-center) is the correct choice.** The centering error is $O(1/\sqrt{n})$, does not drift during MCMC (the Gibbs conditional preserves zero mean in expectation), and its impact on the posterior is negligible. Re-centering would add $O(np^2)$ cost per iteration and create conceptual problems (moving `column_means` target). No action needed beyond documenting that `column_means` is computed from the initially-filled data.

### Overall Assessment

The plan is well-designed, mathematically sound (modulo the §1.3 pseudocode inconsistency), and architecturally clean. The main risk is the incremental `suf_stat_` update, which is correctly implemented in G.1e but should be safeguarded with a full recomputation at the end of each sweep. The centering interaction is a non-issue. The test suite is reasonable but would benefit from the additions listed above.
