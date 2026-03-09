# ===========================================================================
# Group 1: Parameter recovery on known mixed-MRF networks
# ===========================================================================
# Simulate data from networks with known parameters, fit with bgms, and
# verify that posterior means recover the true values. Tests three network
# sizes (small, medium, large) and both PL methods.
#
# Output: numerical summary + PDF with recovery scatter plots.
# ===========================================================================

devtools::load_all(quiet = TRUE)
source(file.path("dev", "tests", "validation", "helpers.R"))

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("  Group 1: Parameter recovery on known networks\n")
cat("=======================================================================\n\n")

# ------------------------------------------------------------------
# 1a. Small network: 2 ordinal + 2 continuous, dense
# ------------------------------------------------------------------
cat("--- 1a: Small network (p=2, q=2) --------------------------------\n")
net_small = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
dat_small = generate_data(net_small, n = 2000, source = "bgms", seed = 201)

true_small = list(
  mux = net_small$mux, muy = net_small$muy,
  Kxx = net_small$Kxx, Kxy = net_small$Kxy, Kyy = net_small$Kyy
)

fit_small_cond = bgm(
  dat_small,
  variable_type = c(rep("ordinal", 2), rep("continuous", 2)),
  pseudolikelihood = "conditional",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 301
)

fit_small_marg = bgm(
  dat_small,
  variable_type = c(rep("ordinal", 2), rep("continuous", 2)),
  pseudolikelihood = "marginal",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 302
)

est_small_cond = extract_bgms_blocks(fit_small_cond, net_small)
est_small_marg = extract_bgms_blocks(fit_small_marg, net_small)

tab_sc = recovery_table(true_small, est_small_cond, "conditional")
tab_sm = recovery_table(true_small, est_small_marg, "marginal")

summarise_recovery(tab_sc, "Small conditional")
summarise_recovery(tab_sm, "Small marginal")

# ------------------------------------------------------------------
# 1b. Medium network: 4 ordinal + 3 continuous, sparse
# ------------------------------------------------------------------
cat("\n--- 1b: Medium network (p=4, q=3) -------------------------------\n")
net_med = make_network(p = 4, q = 3, n_cat = c(1L, 2L, 3L, 1L), density = 0.4, seed = 102)
dat_med = generate_data(net_med, n = 3000, source = "bgms", seed = 202)

true_med = list(
  mux = net_med$mux, muy = net_med$muy,
  Kxx = net_med$Kxx, Kxy = net_med$Kxy, Kyy = net_med$Kyy
)

fit_med_cond = bgm(
  dat_med,
  variable_type = c(rep("ordinal", 4), rep("continuous", 3)),
  pseudolikelihood = "conditional",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 303
)

fit_med_marg = bgm(
  dat_med,
  variable_type = c(rep("ordinal", 4), rep("continuous", 3)),
  pseudolikelihood = "marginal",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 304
)

est_med_cond = extract_bgms_blocks(fit_med_cond, net_med)
est_med_marg = extract_bgms_blocks(fit_med_marg, net_med)

tab_mc = recovery_table(true_med, est_med_cond, "conditional")
tab_mm = recovery_table(true_med, est_med_marg, "marginal")

summarise_recovery(tab_mc, "Medium conditional")
summarise_recovery(tab_mm, "Medium marginal")

# ------------------------------------------------------------------
# 1c. Larger network: 6 ordinal + 4 continuous, moderate density
# ------------------------------------------------------------------
cat("\n--- 1c: Larger network (p=6, q=4) -------------------------------\n")
net_large = make_network(p = 6, q = 4, n_cat = c(1L, 2L, 3L, 1L, 2L, 1L),
                         density = 0.3, seed = 103)
dat_large = generate_data(net_large, n = 5000, source = "bgms", seed = 203)

true_large = list(
  mux = net_large$mux, muy = net_large$muy,
  Kxx = net_large$Kxx, Kxy = net_large$Kxy, Kyy = net_large$Kyy
)

fit_large_marg = bgm(
  dat_large,
  variable_type = c(rep("ordinal", 6), rep("continuous", 4)),
  pseudolikelihood = "marginal",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 305
)

est_large_marg = extract_bgms_blocks(fit_large_marg, net_large)
tab_lm = recovery_table(true_large, est_large_marg, "marginal")
summarise_recovery(tab_lm, "Large marginal")

# ------------------------------------------------------------------
# 1d. Mixed ordinal + Blume-Capel network: 4 discrete + 2 continuous
# ------------------------------------------------------------------
cat("\n--- 1d: Mixed ordinal + BC (p=4, q=2) ---------------------------\n")
net_bc = make_network(
  p = 4, q = 2,
  n_cat = c(1L, 2L, 2L, 3L),
  variable_type = c("ordinal", "ordinal", "blume-capel", "blume-capel"),
  baseline_category = c(0L, 0L, 1L, 1L),
  density = 0.5, seed = 104
)
dat_bc = generate_data(net_bc, n = 3000, source = "bgms", seed = 204)

true_bc = list(
  mux = net_bc$mux, muy = net_bc$muy,
  Kxx = net_bc$Kxx, Kxy = net_bc$Kxy, Kyy = net_bc$Kyy
)

fit_bc_marg = bgm(
  dat_bc,
  variable_type = c("ordinal", "ordinal", "blume-capel", "blume-capel",
                     "continuous", "continuous"),
  baseline_category = c(0L, 0L, 1L, 1L),
  pseudolikelihood = "marginal",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 306
)

est_bc_marg = extract_bgms_blocks(fit_bc_marg, net_bc)
tab_bcm = recovery_table(true_bc, est_bc_marg, "marginal")
summarise_recovery(tab_bcm, "BC marginal")

fit_bc_cond = bgm(
  dat_bc,
  variable_type = c("ordinal", "ordinal", "blume-capel", "blume-capel",
                     "continuous", "continuous"),
  baseline_category = c(0L, 0L, 1L, 1L),
  pseudolikelihood = "conditional",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 307
)

est_bc_cond = extract_bgms_blocks(fit_bc_cond, net_bc)
tab_bcc = recovery_table(true_bc, est_bc_cond, "conditional")
summarise_recovery(tab_bcc, "BC conditional")

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
pdf(file.path(out_dir, "group1_parameter_recovery.pdf"), width = 12, height = 10)
par(mfrow = c(2, 4))
recovery_scatter(tab_sc, "Small conditional PL")
recovery_scatter(tab_sm, "Small marginal PL")
recovery_scatter(tab_mc, "Medium conditional PL")
recovery_scatter(tab_mm, "Medium marginal PL")
recovery_scatter(tab_lm, "Large marginal PL")
recovery_scatter(tab_bcm, "BC marginal PL")
recovery_scatter(tab_bcc, "BC conditional PL")
plot.new()  # blank panel
dev.off()

cat(sprintf("\nPlots saved to %s/group1_parameter_recovery.pdf\n", out_dir))
cat("=== Group 1 complete =============================================\n\n")
