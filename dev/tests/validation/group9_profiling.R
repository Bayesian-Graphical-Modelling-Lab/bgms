# ===========================================================================
# GROUP 9: Performance profiling for the mixed MRF sampler
# ===========================================================================
# Benchmarks bgm() for mixed data across network sizes and sampler types.
#
# Measures:
#   1. Wall-clock time at several (p, q) grid points
#   2. Per-iteration cost (total / n_iter)
#   3. MH vs NUTS comparison
#   4. Marginal vs conditional PL comparison
#   5. Comparison with pure ordinal (OMRF) and pure continuous (GGM) models
#   6. Scaling curves: per-iteration cost vs p (fixing q) and vs q (fixing p)
#
# Usage:
#   Rscript dev/tests/validation/group9_profiling.R
# ===========================================================================

library(bgms)

source("dev/tests/validation/helpers.R")

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("===================================================================\n")
cat("  GROUP 9: Performance Profiling\n")
cat(sprintf("  Started: %s\n", Sys.time()))
cat("===================================================================\n\n")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
n_obs = 300          # observations per dataset
n_iter = 500         # post-warmup iterations
n_warmup = 200       # warmup iterations
n_chains = 1         # single chain for timing consistency
seed = 2026

# Grid of (p, q) sizes for the mixed MRF
mixed_grid = data.frame(
  p = c(3, 4, 6, 8, 10, 12),
  q = c(2, 3, 4, 6, 8,  10)
)

# ---------------------------------------------------------------------------
# Helper: time a single bgm() call
# ---------------------------------------------------------------------------
time_bgm = function(data, variable_type, pairwise_scale = 2.5,
                    iter = n_iter, warmup = n_warmup, ...) {
  gc(verbose = FALSE)
  t0 = proc.time()
  fit = bgm(
    x = data,
    variable_type = variable_type,
    pairwise_scale = pairwise_scale,
    iter = iter,
    warmup = warmup,
    chains = n_chains,
    cores = 1,
    display_progress = FALSE,
    ...
  )
  elapsed = (proc.time() - t0)["elapsed"]
  list(elapsed = unname(elapsed), fit = fit)
}

# ---------------------------------------------------------------------------
# 1. Mixed MRF: wall-clock across (p, q) grid — MH, conditional PL
# ---------------------------------------------------------------------------
cat("--- 1. Mixed MRF wall-clock (MH, conditional PL) ---\n")

mixed_results = data.frame(
  p = integer(), q = integer(), elapsed = numeric(),
  per_iter_ms = numeric(), stringsAsFactors = FALSE
)

for(i in seq_len(nrow(mixed_grid))) {
  p = mixed_grid$p[i]
  q = mixed_grid$q[i]
  n_cat = rep(2L, p)   # binary ordinal for simplicity

  net = make_network(p = p, q = q, n_cat = n_cat, density = 0.4, seed = seed)
  dat = generate_data(net, n = n_obs, seed = seed)

  vt = c(rep("ordinal", p), rep("continuous", q))

  res = time_bgm(dat, variable_type = vt,
                 iter = n_iter, warmup = n_warmup)

  per_iter_ms = res$elapsed / n_iter * 1000

  mixed_results = rbind(mixed_results, data.frame(
    p = p, q = q, elapsed = round(res$elapsed, 2),
    per_iter_ms = round(per_iter_ms, 2)
  ))

  cat(sprintf("  p=%2d q=%2d : %.2fs total, %.2f ms/iter\n",
              p, q, res$elapsed, per_iter_ms))
}

cat("\nMixed MRF timing summary (MH, conditional PL):\n")
print(mixed_results, row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# 2. Conditional vs Marginal PL comparison
# ---------------------------------------------------------------------------
cat("--- 2. Conditional vs Marginal PL comparison ---\n")

pl_grid = data.frame(p = c(4, 8, 12), q = c(3, 6, 10))
pl_results = data.frame(
  p = integer(), q = integer(),
  cond_ms = numeric(), marg_ms = numeric(),
  ratio = numeric(), stringsAsFactors = FALSE
)

for(i in seq_len(nrow(pl_grid))) {
  p = pl_grid$p[i]
  q = pl_grid$q[i]
  n_cat = rep(2L, p)

  net = make_network(p = p, q = q, n_cat = n_cat, density = 0.4, seed = seed)
  dat = generate_data(net, n = n_obs, seed = seed)
  vt = c(rep("ordinal", p), rep("continuous", q))

  res_cond = time_bgm(dat, variable_type = vt,
                      iter = n_iter, warmup = n_warmup,
                      pseudolikelihood = "conditional")
  res_marg = time_bgm(dat, variable_type = vt,
                      iter = n_iter, warmup = n_warmup,
                      pseudolikelihood = "marginal")

  cond_ms = res_cond$elapsed / n_iter * 1000
  marg_ms = res_marg$elapsed / n_iter * 1000

  pl_results = rbind(pl_results, data.frame(
    p = p, q = q,
    cond_ms = round(cond_ms, 2),
    marg_ms = round(marg_ms, 2),
    ratio = round(marg_ms / cond_ms, 2)
  ))

  cat(sprintf("  p=%2d q=%2d : cond=%.2f ms/iter, marg=%.2f ms/iter (ratio=%.2fx)\n",
              p, q, cond_ms, marg_ms, marg_ms / cond_ms))
}

cat("\nPL comparison:\n")
print(pl_results, row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# 3. MH vs NUTS (hybrid-nuts) comparison
# ---------------------------------------------------------------------------
cat("--- 3. MH vs hybrid-NUTS comparison ---\n")

nuts_grid = data.frame(p = c(4, 8), q = c(3, 6))
nuts_results = data.frame(
  p = integer(), q = integer(),
  mh_ms = numeric(), nuts_ms = numeric(),
  ratio = numeric(), stringsAsFactors = FALSE
)

for(i in seq_len(nrow(nuts_grid))) {
  p = nuts_grid$p[i]
  q = nuts_grid$q[i]
  n_cat = rep(2L, p)

  net = make_network(p = p, q = q, n_cat = n_cat, density = 0.4, seed = seed)
  dat = generate_data(net, n = n_obs, seed = seed)
  vt = c(rep("ordinal", p), rep("continuous", q))

  res_mh = time_bgm(dat, variable_type = vt,
                    iter = n_iter, warmup = n_warmup,
                    update_method = "adaptive-metropolis")
  res_nuts = time_bgm(dat, variable_type = vt,
                     iter = n_iter, warmup = n_warmup,
                     update_method = "nuts")

  mh_ms = res_mh$elapsed / n_iter * 1000
  nuts_ms = res_nuts$elapsed / n_iter * 1000

  nuts_results = rbind(nuts_results, data.frame(
    p = p, q = q,
    mh_ms = round(mh_ms, 2),
    nuts_ms = round(nuts_ms, 2),
    ratio = round(nuts_ms / mh_ms, 2)
  ))

  cat(sprintf("  p=%2d q=%2d : MH=%.2f ms/iter, NUTS=%.2f ms/iter (ratio=%.2fx)\n",
              p, q, mh_ms, nuts_ms, nuts_ms / mh_ms))
}

cat("\nMH vs NUTS:\n")
print(nuts_results, row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# 4. Cross-model comparison: Mixed vs pure OMRF vs pure GGM
# ---------------------------------------------------------------------------
cat("--- 4. Cross-model comparison (mixed vs OMRF vs GGM) ---\n")

cross_sizes = c(4, 8, 12)
cross_results = data.frame(
  size = integer(),
  mixed_ms = numeric(), omrf_ms = numeric(), ggm_ms = numeric(),
  stringsAsFactors = FALSE
)

for(sz in cross_sizes) {
  # Mixed: p = sz discrete, q = sz continuous
  net = make_network(p = sz, q = sz, n_cat = rep(2L, sz),
                     density = 0.4, seed = seed)
  dat_mixed = generate_data(net, n = n_obs, seed = seed)
  vt_mixed = c(rep("ordinal", sz), rep("continuous", sz))

  res_mixed = time_bgm(dat_mixed, variable_type = vt_mixed,
                       iter = n_iter, warmup = n_warmup)

  # Pure OMRF: 2*sz ordinal variables
  dat_omrf = dat_mixed[, 1:sz, drop = FALSE]
  vt_omrf = rep("ordinal", sz)
  res_omrf = time_bgm(dat_omrf, variable_type = vt_omrf,
                      iter = n_iter, warmup = n_warmup)

  # Pure GGM: 2*sz continuous variables
  dat_ggm = dat_mixed[, (sz + 1):(2 * sz), drop = FALSE]
  res_ggm = time_bgm(dat_ggm, variable_type = rep("continuous", sz),
                     iter = n_iter, warmup = n_warmup)

  mixed_ms = res_mixed$elapsed / n_iter * 1000
  omrf_ms = res_omrf$elapsed / n_iter * 1000
  ggm_ms = res_ggm$elapsed / n_iter * 1000

  cross_results = rbind(cross_results, data.frame(
    size = sz,
    mixed_ms = round(mixed_ms, 2),
    omrf_ms = round(omrf_ms, 2),
    ggm_ms = round(ggm_ms, 2)
  ))

  cat(sprintf("  size=%2d : mixed=%.2f, omrf=%.2f, ggm=%.2f ms/iter\n",
              sz, mixed_ms, omrf_ms, ggm_ms))
}

cat("\nCross-model comparison:\n")
print(cross_results, row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# 5. Scaling plots
# ---------------------------------------------------------------------------
cat("--- 5. Generating scaling plots ---\n")

pdf_path = file.path(out_dir, "group9_profiling.pdf")
pdf(pdf_path, width = 10, height = 8)

# Plot 1: Per-iteration cost vs total variables (p + q)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

plot(mixed_results$p + mixed_results$q, mixed_results$per_iter_ms,
     type = "b", pch = 19, col = "steelblue",
     xlab = "Total variables (p + q)", ylab = "ms / iteration",
     main = "Mixed MRF: per-iteration cost vs network size",
     log = "y")
grid()

# Plot 2: Conditional vs Marginal PL
barplot_mat = t(as.matrix(pl_results[, c("cond_ms", "marg_ms")]))
colnames(barplot_mat) = sprintf("p=%d,q=%d", pl_results$p, pl_results$q)
barplot(barplot_mat, beside = TRUE,
        col = c("steelblue", "coral"),
        ylab = "ms / iteration",
        main = "Conditional vs Marginal PL")
legend("topleft", legend = c("Conditional", "Marginal"),
       fill = c("steelblue", "coral"), bty = "n")

# Plot 3: MH vs NUTS
barplot_mat2 = t(as.matrix(nuts_results[, c("mh_ms", "nuts_ms")]))
colnames(barplot_mat2) = sprintf("p=%d,q=%d", nuts_results$p, nuts_results$q)
barplot(barplot_mat2, beside = TRUE,
        col = c("steelblue", "seagreen"),
        ylab = "ms / iteration",
        main = "MH vs hybrid-NUTS")
legend("topleft", legend = c("MH", "hybrid-NUTS"),
       fill = c("steelblue", "seagreen"), bty = "n")

# Plot 4: Cross-model comparison
barplot_mat3 = t(as.matrix(cross_results[, c("mixed_ms", "omrf_ms", "ggm_ms")]))
colnames(barplot_mat3) = sprintf("size=%d", cross_results$size)
barplot(barplot_mat3, beside = TRUE,
        col = c("steelblue", "coral", "seagreen"),
        ylab = "ms / iteration",
        main = "Mixed vs OMRF vs GGM")
legend("topleft", legend = c("Mixed", "OMRF", "GGM"),
       fill = c("steelblue", "coral", "seagreen"), bty = "n")

# Plot 5: Scaling curves — p vs cost (q fixed) and q vs cost (p fixed)
# Extract from mixed_results: fix q at smallest/middle value
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

plot(mixed_results$p, mixed_results$per_iter_ms,
     type = "b", pch = 19, col = "steelblue",
     xlab = "p (discrete)", ylab = "ms / iteration",
     main = "Per-iteration cost vs p")
grid()

plot(mixed_results$q, mixed_results$per_iter_ms,
     type = "b", pch = 19, col = "coral",
     xlab = "q (continuous)", ylab = "ms / iteration",
     main = "Per-iteration cost vs q")
grid()

dev.off()

cat(sprintf("  Plots saved to: %s\n", pdf_path))

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n===================================================================\n")
cat("  GROUP 9: Performance Profiling — COMPLETE\n")
cat(sprintf("  Finished: %s\n", Sys.time()))
cat("===================================================================\n")
