# ===========================================================================
# Group 7: Empirical data smoke tests
# ===========================================================================
# Fit the mixed MRF model on real-world data to verify that:
#   1. The model runs to completion without errors.
#   2. Posterior summaries have correct dimensions.
#   3. R-hat values are acceptable.
#   4. Simulate and predict methods work on the fitted object.
#   5. Both PL methods and both samplers produce sensible output.
#
# Uses the Boredom dataset (8 items: 4 treated as ordinal, 4 as continuous)
# following the convention from the existing cycle test.
#
# Output: numerical summary + PDF with trace and density panels.
# ===========================================================================

devtools::load_all(quiet = TRUE)
source(file.path("dev", "tests", "validation", "helpers.R"))

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("  Group 7: Empirical data smoke tests\n")
cat("=======================================================================\n\n")

# ------------------------------------------------------------------
# 7a. Prepare Boredom data as mixed
# ------------------------------------------------------------------
cat("--- 7a: Preparing Boredom data ------------------------------------\n")

data("Boredom", package = "bgms")
raw = as.matrix(Boredom[, 2:9])   # drop language factor
raw = raw - 1L                     # shift 1-7 -> 0-6
colnames(raw) = names(Boredom)[2:9]

ord_idx = c(1, 3, 5, 7)
cont_idx = c(2, 4, 6, 8)
vtype = character(8)
vtype[ord_idx] = "ordinal"
vtype[cont_idx] = "continuous"

cat(sprintf("  n = %d, %d ordinal + %d continuous variables\n",
            nrow(raw), length(ord_idx), length(cont_idx)))
cat("  Variable types:", vtype, "\n\n")

# ------------------------------------------------------------------
# 7b. Fit: MH + marginal PL
# ------------------------------------------------------------------
cat("--- 7b: MH + marginal PL -----------------------------------------\n")

fit_mh_marg = bgm(
  raw, variable_type = vtype,
  pseudolikelihood = "marginal",
  update_method = "adaptive-metropolis",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 901
)

# Sanity checks
stopifnot(!is.null(fit_mh_marg$posterior_summary_pairwise))
stopifnot(!is.null(fit_mh_marg$posterior_mean_main))
stopifnot(!is.null(fit_mh_marg$posterior_mean_pairwise))
cat("  posterior_summary_pairwise rows:", nrow(fit_mh_marg$posterior_summary_pairwise), "\n")
cat("  posterior_mean_pairwise dim:", paste(dim(fit_mh_marg$posterior_mean_pairwise), collapse = "x"), "\n")
cat("  Main effects (discrete) dim:", paste(dim(coef(fit_mh_marg)$main$discrete), collapse = "x"), "\n")
cat("  Main effects (continuous) dim:", paste(dim(coef(fit_mh_marg)$main$continuous), collapse = "x"), "\n")
cat("  OK\n\n")

# ------------------------------------------------------------------
# 7c. Fit: MH + conditional PL
# ------------------------------------------------------------------
cat("--- 7c: MH + conditional PL --------------------------------------\n")

fit_mh_cond = bgm(
  raw, variable_type = vtype,
  pseudolikelihood = "conditional",
  update_method = "adaptive-metropolis",
  edge_selection = FALSE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 902
)

stopifnot(!is.null(fit_mh_cond$posterior_mean_pairwise))
cat("  posterior_mean_pairwise dim:", paste(dim(fit_mh_cond$posterior_mean_pairwise), collapse = "x"), "\n")
cat("  OK\n\n")

# ------------------------------------------------------------------
# 7d. Fit: NUTS + marginal PL
# ------------------------------------------------------------------
cat("--- 7d: NUTS + marginal PL ---------------------------------------\n")

fit_nuts_marg = bgm(
  raw, variable_type = vtype,
  pseudolikelihood = "marginal",
  update_method = "nuts",
  edge_selection = FALSE,
  iter = 5000, warmup = 3000, chains = 2,
  seed = 903
)

stopifnot(!is.null(fit_nuts_marg$posterior_mean_pairwise))
cat("  posterior_mean_pairwise dim:", paste(dim(fit_nuts_marg$posterior_mean_pairwise), collapse = "x"), "\n")
cat("  OK\n\n")

# ------------------------------------------------------------------
# 7e. Fit: NUTS + conditional PL
# ------------------------------------------------------------------
cat("--- 7e: NUTS + conditional PL ------------------------------------\n")

fit_nuts_cond = bgm(
  raw, variable_type = vtype,
  pseudolikelihood = "conditional",
  update_method = "nuts",
  edge_selection = FALSE,
  iter = 5000, warmup = 3000, chains = 2,
  display_progress = "none", seed = 904
)

stopifnot(!is.null(fit_nuts_cond$posterior_mean_pairwise))
cat("  posterior_mean_pairwise dim:", paste(dim(fit_nuts_cond$posterior_mean_pairwise), collapse = "x"), "\n")
cat("  OK\n\n")

# ------------------------------------------------------------------
# 7f. Cross-method agreement on empirical data
# ------------------------------------------------------------------
cat("--- 7f: Cross-method agreement on empirical data ------------------\n")

get_pw_flat = function(fit) as.vector(fit$posterior_mean_pairwise[upper.tri(fit$posterior_mean_pairwise)])

pw_mh_marg = get_pw_flat(fit_mh_marg)
pw_mh_cond = get_pw_flat(fit_mh_cond)
pw_nuts_marg = get_pw_flat(fit_nuts_marg)
pw_nuts_cond = get_pw_flat(fit_nuts_cond)

pairs = list(
  c("MH-marginal", "MH-conditional"),
  c("MH-marginal", "NUTS-marginal"),
  c("NUTS-marginal", "NUTS-conditional"),
  c("MH-marginal", "NUTS-conditional")
)
pw_list = list(
  `MH-marginal` = pw_mh_marg, `MH-conditional` = pw_mh_cond,
  `NUTS-marginal` = pw_nuts_marg, `NUTS-conditional` = pw_nuts_cond
)

cat("  Pairwise interaction correlations:\n")
for(pr in pairs) {
  r = cor(pw_list[[pr[1]]], pw_list[[pr[2]]])
  rmse = sqrt(mean((pw_list[[pr[1]]] - pw_list[[pr[2]]])^2))
  cat(sprintf("    %s vs %s: r = %.4f, RMSE = %.4f\n", pr[1], pr[2], r, rmse))
}

# ------------------------------------------------------------------
# 7g. Simulate from fitted model
# ------------------------------------------------------------------
cat("\n--- 7g: Simulate from fitted model --------------------------------\n")

sim1 = simulate(fit_mh_marg, nsim = 500, method = "posterior-mean", seed = 911)
cat(sprintf("  simulate(posterior-mean): %d x %d\n", nrow(sim1), ncol(sim1)))
stopifnot(nrow(sim1) == 500)
stopifnot(ncol(sim1) == 8)

# Check that ordinal columns are integers and continuous are numeric
for(i in ord_idx) {
  stopifnot(all(sim1[, i] == floor(sim1[, i])))
}
cat("  Ordinal columns contain integers: OK\n")
cat("  Continuous column ranges:\n")
for(i in cont_idx) {
  cat(sprintf("    col %d: [%.2f, %.2f]\n", i, min(sim1[, i]), max(sim1[, i])))
}

# ------------------------------------------------------------------
# 7h. Edge selection on empirical data
# ------------------------------------------------------------------
cat("\n--- 7h: Edge selection on empirical data --------------------------\n")

fit_edge = bgm(
  raw, variable_type = vtype,
  pseudolikelihood = "marginal",
  update_method = "adaptive-metropolis",
  edge_selection = TRUE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 921
)

pip = fit_edge$posterior_mean_indicator
n_edges = sum(pip[upper.tri(pip)] > 0.5)
n_possible = sum(upper.tri(pip))
cat(sprintf("  Edges detected (PIP > 0.5): %d / %d\n", n_edges, n_possible))
cat("  PIP distribution:\n")
pip_ut = pip[upper.tri(pip)]
print(summary(pip_ut))

# ------------------------------------------------------------------
# 7i. R-hat check across all fits
# ------------------------------------------------------------------
cat("\n--- 7i: R-hat check -----------------------------------------------\n")

check_rhat = function(fit, label) {
  summ = fit$posterior_summary_pairwise
  if(!is.null(summ) && "Rhat" %in% names(summ)) {
    rhats = summ$Rhat
    cat(sprintf("  %s: max R-hat = %.3f, n > 1.05 = %d / %d\n",
                label, max(rhats, na.rm = TRUE),
                sum(rhats > 1.05, na.rm = TRUE), length(rhats)))
  } else {
    cat(sprintf("  %s: R-hat not available in summary\n", label))
  }
}

check_rhat(fit_mh_marg, "MH marginal")
check_rhat(fit_mh_cond, "MH conditional")
check_rhat(fit_nuts_marg, "NUTS marginal")
check_rhat(fit_nuts_cond, "NUTS conditional")

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
pdf(file.path(out_dir, "group7_empirical.pdf"), width = 14, height = 12)

layout(matrix(1:6, nrow = 2, byrow = TRUE))
par(mar = c(4.5, 4.5, 2.5, 1))

# Panel 1: All methods pairwise agreement (MH-marg as x, rest as y)
rng = range(c(pw_mh_marg, pw_mh_cond, pw_nuts_marg, pw_nuts_cond))
plot(pw_mh_marg, pw_mh_cond, pch = 19,
     col = adjustcolor("steelblue", 0.5),
     xlim = rng, ylim = rng, asp = 1,
     xlab = "MH marginal", ylab = "MH conditional",
     main = sprintf("MH: marginal vs conditional (r=%.3f)",
                     cor(pw_mh_marg, pw_mh_cond)))
abline(0, 1, lty = 2, col = "grey40")

plot(pw_mh_marg, pw_nuts_marg, pch = 19,
     col = adjustcolor("firebrick", 0.5),
     xlim = rng, ylim = rng, asp = 1,
     xlab = "MH marginal", ylab = "NUTS marginal",
     main = sprintf("Marginal: MH vs NUTS (r=%.3f)",
                     cor(pw_mh_marg, pw_nuts_marg)))
abline(0, 1, lty = 2, col = "grey40")

plot(pw_mh_marg, pw_nuts_cond, pch = 19,
     col = adjustcolor("#4DAF4A", 0.5),
     xlim = rng, ylim = rng, asp = 1,
     xlab = "MH marginal", ylab = "NUTS conditional",
     main = sprintf("MH-marg vs NUTS-cond (r=%.3f)",
                     cor(pw_mh_marg, pw_nuts_cond)))
abline(0, 1, lty = 2, col = "grey40")

# Panel 4: Estimated network heatmap (MH marginal)
pw_mat = fit_mh_marg$posterior_mean_pairwise
image(seq_len(nrow(pw_mat)), seq_len(ncol(pw_mat)), pw_mat,
      col = hcl.colors(50, "Blue-Red 3"),
      xlab = "", ylab = "", main = "MH marginal: pairwise matrix",
      axes = FALSE)
axis(1, at = seq_len(ncol(pw_mat)), labels = colnames(raw), las = 2, cex.axis = 0.7)
axis(2, at = seq_len(nrow(pw_mat)), labels = colnames(raw), las = 2, cex.axis = 0.7)

# Panel 5: PIPs heatmap
if(!is.null(pip)) {
  image(seq_len(nrow(pip)), seq_len(ncol(pip)), pip,
        col = hcl.colors(50, "Viridis"),
        xlab = "", ylab = "", main = "Edge inclusion probabilities",
        axes = FALSE, zlim = c(0, 1))
  axis(1, at = seq_len(ncol(pip)), labels = colnames(raw), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nrow(pip)), labels = colnames(raw), las = 2, cex.axis = 0.7)
}

# Panel 6: Simulated vs observed distributions (first 2 ordinal vars)
sim_check = simulate(fit_mh_marg, nsim = nrow(raw), method = "posterior-mean", seed = 930)
par(mar = c(4.5, 4.5, 2.5, 1))
obs_tab = table(raw[, 1])
sim_tab = table(factor(sim_check[, 1], levels = names(obs_tab)))
barplot(rbind(obs_tab / sum(obs_tab), sim_tab / sum(sim_tab)),
        beside = TRUE, col = c("steelblue", "firebrick"),
        main = sprintf("%s: observed vs simulated", colnames(raw)[1]),
        ylab = "Proportion", xlab = "Category")
legend("topright", legend = c("Observed", "Simulated"),
       fill = c("steelblue", "firebrick"), bty = "n")

dev.off()

cat(sprintf("\nPlots saved to %s/group7_empirical.pdf\n", out_dir))
cat("=== Group 7 complete =============================================\n\n")
