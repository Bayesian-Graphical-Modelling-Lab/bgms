# ===========================================================================
# Group 6: MCMC diagnostics
# ===========================================================================
# Evaluate chain convergence and mixing quality for the mixed MRF sampler.
# Reports R-hat, ESS, and effective samples per second. Visualises trace
# plots, density overlays per chain, and autocorrelation.
#
# Special focus on Kyy parameters (the MH component in hybrid NUTS) to
# verify that MH acceptance is reasonable and chains mix.
#
# Output: numerical summary + PDF with diagnostic panels.
# ===========================================================================

devtools::load_all(quiet = TRUE)
source(file.path("dev", "tests", "validation", "helpers.R"))

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("  Group 6: MCMC diagnostics\n")
cat("=======================================================================\n\n")

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

# Split-Rhat (Vehtari et al. 2021 simplified version)
split_rhat = function(chains) {
  # chains: list of numeric vectors (one per chain)
  n = min(sapply(chains, length))
  # Split each chain in half
  split_chains = unlist(lapply(chains, function(ch) {
    half = floor(n / 2)
    list(ch[1:half], ch[(half + 1):(2 * half)])
  }), recursive = FALSE)
  m = length(split_chains)
  n_half = length(split_chains[[1]])
  chain_means = sapply(split_chains, mean)
  chain_vars = sapply(split_chains, var)
  W = mean(chain_vars)
  B = var(chain_means) * n_half
  var_hat = (1 - 1 / n_half) * W + B / n_half
  sqrt(var_hat / W)
}

# ESS (simple formula based on autocorrelation)
simple_ess = function(x) {
  n = length(x)
  if(n < 10) return(NA)
  acf_vals = acf(x, lag.max = min(n - 1, 500), plot = FALSE)$acf[, 1, 1]
  # Geyer's initial positive sequence
  T_max = length(acf_vals)
  rho_sum = 0
  for(k in seq(2, T_max, by = 2)) {
    pair_sum = acf_vals[k] + if(k + 1 <= T_max) acf_vals[k + 1] else 0
    if(pair_sum < 0) break
    rho_sum = rho_sum + pair_sum
  }
  n / (1 + 2 * rho_sum)
}

# ------------------------------------------------------------------
# 6a. Fit with multiple chains
# ------------------------------------------------------------------
cat("--- 6a: Fitting 4-chain models -----------------------------------\n")

net = make_network(p = 4, q = 3, n_cat = c(1L, 2L, 3L, 1L), density = 0.4, seed = 102)
dat = generate_data(net, n = 3000, source = "bgms", seed = 202)
vtype = c(rep("ordinal", 4), rep("continuous", 3))

t_mh = system.time({
  fit_mh = bgm(dat, variable_type = vtype,
               pseudolikelihood = "marginal",
               update_method = "adaptive-metropolis",
               edge_selection = FALSE,
               iter = 10000, warmup = 5000, chains = 4,
               seed = 801)
})

t_nuts = system.time({
  fit_nuts = bgm(dat, variable_type = vtype,
                 pseudolikelihood = "marginal",
                 update_method = "nuts",
                 edge_selection = FALSE,
                 iter = 5000, warmup = 3000, chains = 4,
                 seed = 802)
})

cat(sprintf("  MH wall time:   %.1f sec\n", t_mh["elapsed"]))
cat(sprintf("  NUTS wall time: %.1f sec\n", t_nuts["elapsed"]))

# ------------------------------------------------------------------
# 6b. Extract per-chain diagnostics
# ------------------------------------------------------------------
cat("\n--- 6b: Per-parameter diagnostics ---------------------------------\n")

compute_diagnostics = function(fit, label, wall_time) {
  main_chains = fit$raw_samples$main
  pair_chains = fit$raw_samples$pairwise
  nchains = fit$raw_samples$nchains
  main_names = fit$raw_samples$parameter_names$main
  pair_names = fit$raw_samples$parameter_names$pairwise

  all_names = c(main_names, pair_names)
  n_main = length(main_names)
  n_pair = length(pair_names)

  rhat_vals = ess_vals = numeric(length(all_names))
  names(rhat_vals) = names(ess_vals) = all_names

  for(j in seq_len(n_main)) {
    chains_j = lapply(main_chains, function(m) m[, j])
    rhat_vals[j] = split_rhat(chains_j)
    ess_vals[j] = simple_ess(do.call(c, chains_j))
  }
  for(j in seq_len(n_pair)) {
    chains_j = lapply(pair_chains, function(m) m[, j])
    rhat_vals[n_main + j] = split_rhat(chains_j)
    ess_vals[n_main + j] = simple_ess(do.call(c, chains_j))
  }

  # Classify parameters
  block = ifelse(seq_along(all_names) <= n_main, "main", "pairwise")

  data.frame(
    parameter = all_names,
    block = block,
    rhat = round(rhat_vals, 4),
    ess = round(ess_vals, 0),
    ess_per_sec = round(ess_vals / wall_time, 1),
    method = label,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

diag_mh = compute_diagnostics(fit_mh, "MH", t_mh["elapsed"])
diag_nuts = compute_diagnostics(fit_nuts, "NUTS", t_nuts["elapsed"])
diag_all = rbind(diag_mh, diag_nuts)

# Summary statistics
cat("\n  MH diagnostics:\n")
cat(sprintf("    R-hat: median = %.3f, max = %.3f, n > 1.01 = %d\n",
            median(diag_mh$rhat), max(diag_mh$rhat), sum(diag_mh$rhat > 1.01)))
cat(sprintf("    ESS:   median = %.0f, min = %.0f, ESS/s median = %.1f\n",
            median(diag_mh$ess), min(diag_mh$ess), median(diag_mh$ess_per_sec)))

cat("\n  NUTS diagnostics:\n")
cat(sprintf("    R-hat: median = %.3f, max = %.3f, n > 1.01 = %d\n",
            median(diag_nuts$rhat), max(diag_nuts$rhat), sum(diag_nuts$rhat > 1.01)))
cat(sprintf("    ESS:   median = %.0f, min = %.0f, ESS/s median = %.1f\n",
            median(diag_nuts$ess), min(diag_nuts$ess), median(diag_nuts$ess_per_sec)))

# ------------------------------------------------------------------
# 6c. NUTS-specific diagnostics
# ------------------------------------------------------------------
cat("\n--- 6c: NUTS-specific diagnostics ---------------------------------\n")
if(!is.null(fit_nuts$nuts_diag)) {
  nuts_d = fit_nuts$nuts_diag
  if(is.list(nuts_d)) {
    # Tree depth distribution
    if(!is.null(nuts_d$tree_depth) || !is.null(nuts_d[[1]]$tree_depth)) {
      cat("  NUTS diagnostics available\n")
      # Try to extract tree depths from first chain
      td = if(is.list(nuts_d[[1]])) nuts_d[[1]]$tree_depth else nuts_d$tree_depth
      if(!is.null(td)) {
        cat("  Tree depth distribution:\n")
        print(table(td))
        cat(sprintf("  Mean tree depth: %.2f\n", mean(td)))
      }
    }
  }
} else {
  cat("  No NUTS diagnostics available in fit object.\n")
}

# ------------------------------------------------------------------
# 6d. Identify problematic parameters
# ------------------------------------------------------------------
cat("\n--- 6d: Flagged parameters ----------------------------------------\n")
flagged = diag_all[diag_all$rhat > 1.05 | diag_all$ess < 100, ]
if(nrow(flagged) > 0) {
  cat(sprintf("  %d parameters flagged (R-hat > 1.05 or ESS < 100):\n", nrow(flagged)))
  print(flagged, row.names = FALSE)
} else {
  cat("  No parameters flagged. All R-hat <= 1.05 and ESS >= 100.\n")
}

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
pdf(file.path(out_dir, "group6_diagnostics.pdf"), width = 14, height = 16)

layout(matrix(1:8, nrow = 4, byrow = TRUE))
par(mar = c(4.5, 4.5, 2.5, 1))

# Panel 1: R-hat comparison
plot(diag_mh$rhat, diag_nuts$rhat, pch = 19,
     col = adjustcolor(ifelse(diag_mh$block == "main", "steelblue", "firebrick"), 0.6),
     xlab = "R-hat (MH)", ylab = "R-hat (NUTS)",
     main = "R-hat: MH vs NUTS")
abline(h = 1.01, v = 1.01, lty = 2, col = "grey40")
abline(0, 1, lty = 3, col = "grey60")
legend("topleft", legend = c("main", "pairwise"),
       col = c("steelblue", "firebrick"), pch = 19, bty = "n")

# Panel 2: ESS comparison
plot(diag_mh$ess, diag_nuts$ess, pch = 19,
     col = adjustcolor(ifelse(diag_mh$block == "main", "steelblue", "firebrick"), 0.6),
     xlab = "ESS (MH)", ylab = "ESS (NUTS)",
     main = "ESS: MH vs NUTS", log = "xy")
abline(0, 1, lty = 2, col = "grey40")

# Panel 3: ESS/sec comparison
plot(diag_mh$ess_per_sec, diag_nuts$ess_per_sec, pch = 19,
     col = adjustcolor(ifelse(diag_mh$block == "main", "steelblue", "firebrick"), 0.6),
     xlab = "ESS/sec (MH)", ylab = "ESS/sec (NUTS)",
     main = "Efficiency: MH vs NUTS", log = "xy")
abline(0, 1, lty = 2, col = "grey40")

# Panel 4: R-hat distributions
boxplot(rhat ~ method, data = diag_all, col = c("steelblue", "firebrick"),
        main = "R-hat distribution by method",
        ylab = "R-hat")
abline(h = 1.01, lty = 2, col = "grey40")

# Panels 5-6: Trace plots for selected parameters (MH)
pair_names_mh = fit_mh$raw_samples$parameter_names$pairwise
sel = pair_names_mh[seq_len(min(2, length(pair_names_mh)))]
for(nm_idx in seq_along(sel)) {
  nm = sel[nm_idx]
  j = match(nm, pair_names_mh)
  cols_chain = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  plot(0, 0, type = "n",
       xlim = c(1, fit_mh$raw_samples$niter),
       ylim = range(sapply(fit_mh$raw_samples$pairwise, function(m) range(m[, j]))),
       xlab = "Iteration", ylab = nm,
       main = sprintf("MH trace: %s", nm))
  for(ch in seq_len(fit_mh$raw_samples$nchains)) {
    lines(fit_mh$raw_samples$pairwise[[ch]][, j],
          col = adjustcolor(cols_chain[ch], 0.4))
  }
}

# Panels 7-8: Trace plots for selected parameters (NUTS)
pair_names_nuts = fit_nuts$raw_samples$parameter_names$pairwise
sel_n = pair_names_nuts[seq_len(min(2, length(pair_names_nuts)))]
for(nm_idx in seq_along(sel_n)) {
  nm = sel_n[nm_idx]
  j = match(nm, pair_names_nuts)
  plot(0, 0, type = "n",
       xlim = c(1, fit_nuts$raw_samples$niter),
       ylim = range(sapply(fit_nuts$raw_samples$pairwise, function(m) range(m[, j]))),
       xlab = "Iteration", ylab = nm,
       main = sprintf("NUTS trace: %s", nm))
  for(ch in seq_len(fit_nuts$raw_samples$nchains)) {
    lines(fit_nuts$raw_samples$pairwise[[ch]][, j],
          col = adjustcolor(cols_chain[ch], 0.4))
  }
}

dev.off()

# --- Autocorrelation comparison ---
pdf(file.path(out_dir, "group6_autocorrelation.pdf"), width = 12, height = 10)

# Pick 4 pairwise parameters and compare autocorrelation between MH and NUTS
n_show = min(4, length(pair_names_mh))
par(mfrow = c(n_show, 2), mar = c(3, 4, 2.5, 1), mgp = c(2, 0.6, 0))

for(j in seq_len(n_show)) {
  nm = pair_names_mh[j]
  # MH autocorrelation (chain 1)
  mh_samp = fit_mh$raw_samples$pairwise[[1]][, j]
  acf(mh_samp, lag.max = 100, main = sprintf("MH: %s", nm), cex.main = 0.9)

  nm_n = pair_names_nuts[j]
  nuts_samp = fit_nuts$raw_samples$pairwise[[1]][, j]
  acf(nuts_samp, lag.max = 100, main = sprintf("NUTS: %s", nm_n), cex.main = 0.9)
}

dev.off()

cat(sprintf("\nPlots saved to %s/group6_*.pdf\n", out_dir))
cat("=== Group 6 complete =============================================\n\n")
