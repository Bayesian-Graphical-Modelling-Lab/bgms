# ===========================================================================
# Group 2: Metropolis vs NUTS posterior agreement
# ===========================================================================
# Fit the same data with adaptive-metropolis and hybrid-NUTS. Compare
# posterior means; if both samplers target the same pseudo-posterior,
# estimates should agree closely.
#
# Output: numerical summary + PDF with agreement scatter plots.
# ===========================================================================

devtools::load_all(quiet = TRUE)
source(file.path("dev", "tests", "validation", "helpers.R"))

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("  Group 2: Metropolis vs NUTS posterior agreement\n")
cat("=======================================================================\n\n")

# Use medium network from group 1 for consistency
net = make_network(p = 4, q = 3, n_cat = c(1L, 2L, 3L, 1L), density = 0.4, seed = 102)
dat = generate_data(net, n = 3000, source = "bgms", seed = 202)
vtype = c(rep("ordinal", 4), rep("continuous", 3))

true_blocks = list(
  mux = net$mux, muy = net$muy,
  Kxx = net$Kxx, Kxy = net$Kxy, Kyy = net$Kyy
)

# ------------------------------------------------------------------
# 2a. Conditional PL: Metropolis vs NUTS
# ------------------------------------------------------------------
cat("--- 2a: Conditional PL — Metropolis vs NUTS ----------------------\n")

fit_mh_cond = bgm(
  dat, variable_type = vtype,
  pseudolikelihood = "conditional",
  update_method = "adaptive-metropolis",
  edge_selection = FALSE,
  iter = 15000, warmup = 10000, chains = 2,
  seed = 401
)

fit_nuts_cond = bgm(
  dat, variable_type = vtype,
  pseudolikelihood = "conditional",
  update_method = "nuts",
  edge_selection = FALSE,
  iter = 5000, warmup = 3000, chains = 2,
  seed = 402
)

est_mh_cond = extract_bgms_blocks(fit_mh_cond, net)
est_nuts_cond = extract_bgms_blocks(fit_nuts_cond, net)

tab_mh_cond = recovery_table(true_blocks, est_mh_cond, "MH-conditional")
tab_nuts_cond = recovery_table(true_blocks, est_nuts_cond, "NUTS-conditional")

cat("  Metropolis:\n")
summarise_recovery(tab_mh_cond, "MH cond")
cat("  NUTS:\n")
summarise_recovery(tab_nuts_cond, "NUTS cond")

# Direct agreement
r_cond = cor(tab_mh_cond$estimate, tab_nuts_cond$estimate)
rmse_cond = sqrt(mean((tab_mh_cond$estimate - tab_nuts_cond$estimate)^2))
cat(sprintf("  MH vs NUTS agreement: r = %.4f, RMSE = %.4f\n", r_cond, rmse_cond))

# ------------------------------------------------------------------
# 2b. Marginal PL: Metropolis vs NUTS
# ------------------------------------------------------------------
cat("\n--- 2b: Marginal PL — Metropolis vs NUTS -------------------------\n")

fit_mh_marg = bgm(
  dat, variable_type = vtype,
  pseudolikelihood = "marginal",
  update_method = "adaptive-metropolis",
  edge_selection = FALSE,
  iter = 15000, warmup = 10000, chains = 2,
  seed = 403
)

fit_nuts_marg = bgm(
  dat, variable_type = vtype,
  pseudolikelihood = "marginal",
  update_method = "nuts",
  edge_selection = FALSE,
  iter = 5000, warmup = 3000, chains = 2,
  seed = 404
)

est_mh_marg = extract_bgms_blocks(fit_mh_marg, net)
est_nuts_marg = extract_bgms_blocks(fit_nuts_marg, net)

tab_mh_marg = recovery_table(true_blocks, est_mh_marg, "MH-marginal")
tab_nuts_marg = recovery_table(true_blocks, est_nuts_marg, "NUTS-marginal")

cat("  Metropolis:\n")
summarise_recovery(tab_mh_marg, "MH marg")
cat("  NUTS:\n")
summarise_recovery(tab_nuts_marg, "NUTS marg")

r_marg = cor(tab_mh_marg$estimate, tab_nuts_marg$estimate)
rmse_marg = sqrt(mean((tab_mh_marg$estimate - tab_nuts_marg$estimate)^2))
cat(sprintf("  MH vs NUTS agreement: r = %.4f, RMSE = %.4f\n", r_marg, rmse_marg))

# ------------------------------------------------------------------
# 2c. Mixed ordinal + Blume-Capel: Metropolis vs NUTS (marginal PL)
# ------------------------------------------------------------------
cat("\n--- 2c: BC mixed — Metropolis vs NUTS (marginal PL) --------------\n")

net_bc = make_network(
  p = 4, q = 2,
  n_cat = c(1L, 2L, 2L, 3L),
  variable_type = c("ordinal", "ordinal", "blume-capel", "blume-capel"),
  baseline_category = c(0L, 0L, 1L, 1L),
  density = 0.5, seed = 104
)
dat_bc = generate_data(net_bc, n = 3000, source = "bgms", seed = 204)

vtype_bc = c("ordinal", "ordinal", "blume-capel", "blume-capel",
             "continuous", "continuous")
bl_bc = c(0L, 0L, 1L, 1L)

true_bc = list(mux = net_bc$mux, muy = net_bc$muy,
               Kxx = net_bc$Kxx, Kxy = net_bc$Kxy, Kyy = net_bc$Kyy)

fit_mh_bc = bgm(
  dat_bc, variable_type = vtype_bc, baseline_category = bl_bc,
  pseudolikelihood = "marginal",
  update_method = "adaptive-metropolis",
  edge_selection = FALSE,
  iter = 15000, warmup = 10000, chains = 2,
  seed = 405
)

fit_nuts_bc = bgm(
  dat_bc, variable_type = vtype_bc, baseline_category = bl_bc,
  pseudolikelihood = "marginal",
  update_method = "nuts",
  edge_selection = FALSE,
  iter = 5000, warmup = 3000, chains = 2,
  seed = 406
)

est_mh_bc = extract_bgms_blocks(fit_mh_bc, net_bc)
est_nuts_bc = extract_bgms_blocks(fit_nuts_bc, net_bc)

tab_mh_bc = recovery_table(true_bc, est_mh_bc, "MH-bc")
tab_nuts_bc = recovery_table(true_bc, est_nuts_bc, "NUTS-bc")

cat("  Metropolis:\n")
summarise_recovery(tab_mh_bc, "MH BC")
cat("  NUTS:\n")
summarise_recovery(tab_nuts_bc, "NUTS BC")

r_bc = cor(tab_mh_bc$estimate, tab_nuts_bc$estimate)
rmse_bc = sqrt(mean((tab_mh_bc$estimate - tab_nuts_bc$estimate)^2))
cat(sprintf("  MH vs NUTS agreement: r = %.4f, RMSE = %.4f\n", r_bc, rmse_bc))

# ------------------------------------------------------------------
# 2d. Numerical summary table
# ------------------------------------------------------------------
cat("\n--- 2d: Summary table ---------------------------------------------\n")
summary_df = data.frame(
  comparison = c("conditional", "marginal", "BC-marginal"),
  correlation = round(c(r_cond, r_marg, r_bc), 4),
  rmse = round(c(rmse_cond, rmse_marg, rmse_bc), 4),
  max_diff = round(c(max(abs(tab_mh_cond$estimate - tab_nuts_cond$estimate)),
                      max(abs(tab_mh_marg$estimate - tab_nuts_marg$estimate)),
                      max(abs(tab_mh_bc$estimate - tab_nuts_bc$estimate))), 4)
)
print(summary_df, row.names = FALSE)

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
pdf(file.path(out_dir, "group2_mh_vs_nuts.pdf"), width = 14, height = 14)

layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, byrow = TRUE))
par(mar = c(4.5, 4.5, 2.5, 1))

# Row 1: conditional PL
recovery_scatter(tab_mh_cond, "MH conditional: recovery")
recovery_scatter(tab_nuts_cond, "NUTS conditional: recovery")
agreement_scatter(tab_mh_cond, tab_nuts_cond,
                  "MH cond", "NUTS cond",
                  sprintf("MH vs NUTS conditional (r=%.3f)", r_cond))

# Row 2: marginal PL
recovery_scatter(tab_mh_marg, "MH marginal: recovery")
recovery_scatter(tab_nuts_marg, "NUTS marginal: recovery")
agreement_scatter(tab_mh_marg, tab_nuts_marg,
                  "MH marg", "NUTS marg",
                  sprintf("MH vs NUTS marginal (r=%.3f)", r_marg))

# Row 3: BC marginal PL
recovery_scatter(tab_mh_bc, "MH BC marginal: recovery")
recovery_scatter(tab_nuts_bc, "NUTS BC marginal: recovery")
agreement_scatter(tab_mh_bc, tab_nuts_bc,
                  "MH BC", "NUTS BC",
                  sprintf("MH vs NUTS BC marginal (r=%.3f)", r_bc))

dev.off()

# --- Trace plot comparison for selected parameters ---
pdf(file.path(out_dir, "group2_traces.pdf"), width = 10, height = 14)

# Extract raw samples (combine chains)
mh_cond_pw = do.call(rbind, fit_mh_cond$raw_samples$pairwise)
nuts_cond_pw = do.call(rbind, fit_nuts_cond$raw_samples$pairwise)
mh_marg_pw = do.call(rbind, fit_mh_marg$raw_samples$pairwise)
nuts_marg_pw = do.call(rbind, fit_nuts_marg$raw_samples$pairwise)

# Assign column names
colnames(mh_cond_pw) = fit_mh_cond$raw_samples$parameter_names$pairwise
colnames(nuts_cond_pw) = fit_nuts_cond$raw_samples$parameter_names$pairwise
colnames(mh_marg_pw) = fit_mh_marg$raw_samples$parameter_names$pairwise
colnames(nuts_marg_pw) = fit_nuts_marg$raw_samples$parameter_names$pairwise

# Pick first 4 pairwise params for trace comparison
sel_params = colnames(mh_cond_pw)[seq_len(min(4, ncol(mh_cond_pw)))]

par(mfrow = c(length(sel_params), 2), mar = c(3, 3, 2, 1), mgp = c(2, 0.6, 0))
for(nm in sel_params) {
  # MH trace
  plot(mh_cond_pw[, nm], type = "l", col = adjustcolor("steelblue", 0.4),
       main = paste(nm, "— MH cond"), xlab = "Iteration", ylab = nm, cex.main = 0.9)
  # NUTS trace
  plot(nuts_cond_pw[, nm], type = "l", col = adjustcolor("firebrick", 0.4),
       main = paste(nm, "— NUTS cond"), xlab = "Iteration", ylab = nm, cex.main = 0.9)
}

dev.off()

cat(sprintf("\nPlots saved to %s/group2_mh_vs_nuts.pdf\n", out_dir))
cat(sprintf("Traces saved to %s/group2_traces.pdf\n", out_dir))
cat("=== Group 2 complete =============================================\n\n")
