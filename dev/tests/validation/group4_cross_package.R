# ===========================================================================
# Group 4: Cross-package validation (bgms vs mixedGM)
# ===========================================================================
# Both packages implement the same mixed MRF pseudo-likelihood sampler.
# We fit the same data with both and compare posterior means.
# Also runs a simulate-refit cycle: fit -> simulate -> refit -> compare.
#
# Note: mixedGM uses total-categories convention (binary = 2), while
# bgms uses max-index convention (binary = 1).
# The threshold/mean prior in mixedGM must be set to log_prior_logistic
# to match bgms's default Beta(1,1)-on-logistic prior.
#
# Output: numerical summary + PDF with cross-package comparisons.
# ===========================================================================

devtools::load_all(quiet = TRUE)
library(mixedGM)
source(file.path("dev", "tests", "validation", "helpers.R"))

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("  Group 4: Cross-package validation (bgms vs mixedGM)\n")
cat("=======================================================================\n\n")

# ------------------------------------------------------------------
# 4a. Setup: shared true network + data
# ------------------------------------------------------------------
cat("--- 4a: Generating shared data -----------------------------------\n")

net = make_network(p = 3, q = 2, n_cat = c(1L, 2L, 1L), density = 0.6, seed = 601)

# Generate data via mixedGM (both packages should accept this)
set.seed(611)
sim = mixedGM::mixed_gibbs_generate(
  n = 2000,
  Kxx = net$Kxx, Kxy = net$Kxy, Kyy = net$Kyy,
  mux = net$mux, muy = net$muy,
  num_categories = net$n_cat + 1L,
  n_burnin = 1000
)

x_disc = sim$x   # n x p integer matrix
y_cont = sim$y   # n x q numeric matrix

true_blocks = list(mux = net$mux, muy = net$muy,
                   Kxx = net$Kxx, Kxy = net$Kxy, Kyy = net$Kyy)

cat(sprintf("  p = %d ordinal, q = %d continuous, n = %d\n",
            net$p, net$q, nrow(x_disc)))
cat("  n_cat (bgms convention):", net$n_cat, "\n")
cat("  n_cat (mixedGM convention):", net$n_cat + 1L, "\n\n")

# ------------------------------------------------------------------
# 4b. Fit with mixedGM (marginal PL)
# ------------------------------------------------------------------
cat("--- 4b: mixedGM marginal PL fit ----------------------------------\n")

set.seed(621)
fit_mgm_marg = mixedGM::mixed_sampler(
  x = x_disc, y = y_cont,
  num_categories = net$n_cat + 1L,      # mixedGM: total categories
  pseudolikelihood = "marginal",
  log_prior_mean = mixedGM:::log_prior_logistic,
  edge_selection = FALSE,
  n_warmup = 5000, n_samples = 10000,
  verbose = FALSE
)

est_mgm_marg = extract_mgm_blocks(fit_mgm_marg, n_cat = net$n_cat)
tab_mgm_marg = recovery_table(true_blocks, est_mgm_marg, "mixedGM-marginal")
summarise_recovery(tab_mgm_marg, "mixedGM marginal")

# ------------------------------------------------------------------
# 4c. Fit with bgms (marginal PL)
# ------------------------------------------------------------------
cat("\n--- 4c: bgms marginal PL fit -------------------------------------\n")

# Combine x and y into single data frame (bgms convention)
bgms_dat = as.data.frame(cbind(x_disc, y_cont))
names(bgms_dat) = c(paste0("X", seq_len(net$p)), paste0("Y", seq_len(net$q)))
vtype = c(rep("ordinal", net$p), rep("continuous", net$q))

fit_bgms_marg = bgm(
  bgms_dat, variable_type = vtype,
  pseudolikelihood = "marginal",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 622
)

est_bgms_marg = extract_bgms_blocks(fit_bgms_marg, net)
tab_bgms_marg = recovery_table(true_blocks, est_bgms_marg, "bgms-marginal")
summarise_recovery(tab_bgms_marg, "bgms marginal")

# ------------------------------------------------------------------
# 4d. Head-to-head comparison
# ------------------------------------------------------------------
cat("\n--- 4d: bgms vs mixedGM comparison --------------------------------\n")

r_cross = cor(tab_bgms_marg$estimate, tab_mgm_marg$estimate)
rmse_cross = sqrt(mean((tab_bgms_marg$estimate - tab_mgm_marg$estimate)^2))
max_d_cross = max(abs(tab_bgms_marg$estimate - tab_mgm_marg$estimate))
cat(sprintf("  bgms vs mixedGM: r = %.4f, RMSE = %.4f, max|diff| = %.4f\n",
            r_cross, rmse_cross, max_d_cross))

# Per-block comparison
blocks = unique(tab_bgms_marg$block)
cross_stats = do.call(rbind, lapply(blocks, function(b) {
  idx = tab_bgms_marg$block == b
  d = tab_bgms_marg$estimate[idx] - tab_mgm_marg$estimate[idx]
  data.frame(
    block = b, n = sum(idx),
    mean_diff = round(mean(d), 4),
    max_abs_diff = round(max(abs(d)), 4),
    cor = round(cor(tab_bgms_marg$estimate[idx], tab_mgm_marg$estimate[idx]), 4),
    stringsAsFactors = FALSE
  )
}))
cat("\n  Per-block bgms vs mixedGM:\n")
print(cross_stats, row.names = FALSE)

# ------------------------------------------------------------------
# 4e. Conditional PL comparison
# ------------------------------------------------------------------
cat("\n--- 4e: Conditional PL comparison ---------------------------------\n")

set.seed(631)
fit_mgm_cond = mixedGM::mixed_sampler(
  x = x_disc, y = y_cont,
  num_categories = net$n_cat + 1L,
  pseudolikelihood = "conditional",
  log_prior_mean = mixedGM:::log_prior_logistic,
  edge_selection = FALSE,
  n_warmup = 5000, n_samples = 10000,
  verbose = FALSE
)
est_mgm_cond = extract_mgm_blocks(fit_mgm_cond, n_cat = net$n_cat)
tab_mgm_cond = recovery_table(true_blocks, est_mgm_cond, "mixedGM-conditional")

fit_bgms_cond = bgm(
  bgms_dat, variable_type = vtype,
  pseudolikelihood = "conditional",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 632
)
est_bgms_cond = extract_bgms_blocks(fit_bgms_cond, net)
tab_bgms_cond = recovery_table(true_blocks, est_bgms_cond, "bgms-conditional")

r_cond = cor(tab_bgms_cond$estimate, tab_mgm_cond$estimate)
rmse_cond = sqrt(mean((tab_bgms_cond$estimate - tab_mgm_cond$estimate)^2))
cat(sprintf("  bgms vs mixedGM (conditional): r = %.4f, RMSE = %.4f\n",
            r_cond, rmse_cond))

# ------------------------------------------------------------------
# 4f. Simulate-refit cycle (bgms)
# ------------------------------------------------------------------
cat("\n--- 4f: bgms simulate-refit cycle ---------------------------------\n")

sim_bgms = simulate(fit_bgms_marg, nsim = 2000, method = "posterior-mean", seed = 641)

refit_bgms = bgm(
  sim_bgms, variable_type = vtype,
  pseudolikelihood = "marginal",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 642
)

est_orig = extract_bgms_blocks(fit_bgms_marg, net)
est_refit = extract_bgms_blocks(refit_bgms, net)

tab_orig = recovery_table(true_blocks, est_orig, "bgms-original")
tab_refit = recovery_table(true_blocks, est_refit, "bgms-refit")

r_cycle = cor(tab_orig$estimate, tab_refit$estimate)
rmse_cycle = sqrt(mean((tab_orig$estimate - tab_refit$estimate)^2))
cat(sprintf("  bgms cycle: original vs refit r = %.4f, RMSE = %.4f\n",
            r_cycle, rmse_cycle))

# ------------------------------------------------------------------
# Master summary
# ------------------------------------------------------------------
cat("\n--- Summary of all cross-package comparisons ----------------------\n")
master_summary = data.frame(
  comparison = c("bgms vs mixedGM (marginal)", "bgms vs mixedGM (conditional)",
                 "bgms cycle (marginal)"),
  correlation = round(c(r_cross, r_cond, r_cycle), 4),
  rmse = round(c(rmse_cross, rmse_cond, rmse_cycle), 4),
  max_diff = round(c(max_d_cross,
                      max(abs(tab_bgms_cond$estimate - tab_mgm_cond$estimate)),
                      max(abs(tab_orig$estimate - tab_refit$estimate))), 4)
)
print(master_summary, row.names = FALSE)

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
pdf(file.path(out_dir, "group4_cross_package.pdf"), width = 15, height = 12)
layout(matrix(1:9, nrow = 3, byrow = TRUE))
par(mar = c(4.5, 4.5, 2.5, 1))

# Row 1: Marginal PL — recovery from truth
recovery_scatter(tab_mgm_marg, "mixedGM marginal: recovery")
recovery_scatter(tab_bgms_marg, "bgms marginal: recovery")
agreement_scatter(tab_bgms_marg, tab_mgm_marg,
                  "bgms", "mixedGM",
                  sprintf("bgms vs mixedGM marginal (r=%.3f)", r_cross))

# Row 2: Conditional PL — recovery from truth
recovery_scatter(tab_mgm_cond, "mixedGM conditional: recovery")
recovery_scatter(tab_bgms_cond, "bgms conditional: recovery")
agreement_scatter(tab_bgms_cond, tab_mgm_cond,
                  "bgms", "mixedGM",
                  sprintf("bgms vs mixedGM conditional (r=%.3f)", r_cond))

# Row 3: bgms cycle
recovery_scatter(tab_orig, "bgms original fit")
recovery_scatter(tab_refit, "bgms refit (simulated data)")
agreement_scatter(tab_orig, tab_refit,
                  "Original", "Refit",
                  sprintf("bgms cycle (r=%.3f)", r_cycle))

dev.off()

cat(sprintf("\nPlots saved to %s/group4_cross_package.pdf\n", out_dir))
cat("=== Group 4 complete =============================================\n\n")
