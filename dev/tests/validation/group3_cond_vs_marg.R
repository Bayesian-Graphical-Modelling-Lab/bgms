# ===========================================================================
# Group 3: Conditional vs Marginal pseudo-likelihood agreement
# ===========================================================================
# Both PL methods target the same model; posterior means should converge
# to the same values with enough data. This test compares them on the
# same dataset and checks where they diverge, broken down by block.
#
# Output: numerical summary + PDF with agreement scatter and difference
#         distributions.
# ===========================================================================

devtools::load_all(quiet = TRUE)
source(file.path("dev", "tests", "validation", "helpers.R"))

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("  Group 3: Conditional vs Marginal PL agreement\n")
cat("=======================================================================\n\n")

# ------------------------------------------------------------------
# 3a. Small network (p=2, q=2)
# ------------------------------------------------------------------
cat("--- 3a: Small network (p=2, q=2) ---------------------------------\n")
net_s = make_network(p = 2, q = 2, n_cat = c(1L, 2L), density = 1.0, seed = 101)
dat_s = generate_data(net_s, n = 2000, source = "bgms", seed = 201)
vtype_s = c(rep("ordinal", 2), rep("continuous", 2))

true_s = list(mux = net_s$mux, muy = net_s$muy,
              Kxx = net_s$Kxx, Kxy = net_s$Kxy, Kyy = net_s$Kyy)

fit_cond_s = bgm(dat_s, variable_type = vtype_s,
                 pseudolikelihood = "conditional", edge_selection = FALSE,
                 iter = 10000, warmup = 5000, chains = 2,
                 seed = 501)

fit_marg_s = bgm(dat_s, variable_type = vtype_s,
                 pseudolikelihood = "marginal", edge_selection = FALSE,
                 iter = 10000, warmup = 5000, chains = 2,
                 seed = 502)

est_cond_s = extract_bgms_blocks(fit_cond_s, net_s)
est_marg_s = extract_bgms_blocks(fit_marg_s, net_s)

tab_cond_s = recovery_table(true_s, est_cond_s, "conditional")
tab_marg_s = recovery_table(true_s, est_marg_s, "marginal")

summarise_recovery(tab_cond_s, "Small conditional")
summarise_recovery(tab_marg_s, "Small marginal")

r_s = cor(tab_cond_s$estimate, tab_marg_s$estimate)
max_d_s = max(abs(tab_cond_s$estimate - tab_marg_s$estimate))
cat(sprintf("  Cond vs Marg: r = %.4f, max|diff| = %.4f\n\n", r_s, max_d_s))

# ------------------------------------------------------------------
# 3b. Medium network (p=4, q=3)
# ------------------------------------------------------------------
cat("--- 3b: Medium network (p=4, q=3) --------------------------------\n")
net_m = make_network(p = 4, q = 3, n_cat = c(1L, 2L, 3L, 1L), density = 0.4, seed = 102)
dat_m = generate_data(net_m, n = 3000, source = "bgms", seed = 202)
vtype_m = c(rep("ordinal", 4), rep("continuous", 3))

true_m = list(mux = net_m$mux, muy = net_m$muy,
              Kxx = net_m$Kxx, Kxy = net_m$Kxy, Kyy = net_m$Kyy)

fit_cond_m = bgm(dat_m, variable_type = vtype_m,
                 pseudolikelihood = "conditional", edge_selection = FALSE,
                 iter = 10000, warmup = 5000, chains = 2,
                 seed = 503)

fit_marg_m = bgm(dat_m, variable_type = vtype_m,
                 pseudolikelihood = "marginal", edge_selection = FALSE,
                 iter = 10000, warmup = 5000, chains = 2,
                 seed = 504)

est_cond_m = extract_bgms_blocks(fit_cond_m, net_m)
est_marg_m = extract_bgms_blocks(fit_marg_m, net_m)

tab_cond_m = recovery_table(true_m, est_cond_m, "conditional")
tab_marg_m = recovery_table(true_m, est_marg_m, "marginal")

summarise_recovery(tab_cond_m, "Medium conditional")
summarise_recovery(tab_marg_m, "Medium marginal")

r_m = cor(tab_cond_m$estimate, tab_marg_m$estimate)
max_d_m = max(abs(tab_cond_m$estimate - tab_marg_m$estimate))
cat(sprintf("  Cond vs Marg: r = %.4f, max|diff| = %.4f\n\n", r_m, max_d_m))

# ------------------------------------------------------------------
# 3c. Mixed ordinal + Blume-Capel (p=4, q=2)
# ------------------------------------------------------------------
cat("--- 3c: Mixed ordinal + BC (p=4, q=2) ----------------------------\n")
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

fit_cond_bc = bgm(dat_bc, variable_type = vtype_bc, baseline_category = bl_bc,
                  pseudolikelihood = "conditional", edge_selection = FALSE,
                  iter = 10000, warmup = 5000, chains = 2,
                  seed = 505)

fit_marg_bc = bgm(dat_bc, variable_type = vtype_bc, baseline_category = bl_bc,
                  pseudolikelihood = "marginal", edge_selection = FALSE,
                  iter = 10000, warmup = 5000, chains = 2,
                  seed = 506)

est_cond_bc = extract_bgms_blocks(fit_cond_bc, net_bc)
est_marg_bc = extract_bgms_blocks(fit_marg_bc, net_bc)

tab_cond_bc = recovery_table(true_bc, est_cond_bc, "conditional")
tab_marg_bc = recovery_table(true_bc, est_marg_bc, "marginal")

summarise_recovery(tab_cond_bc, "BC conditional")
summarise_recovery(tab_marg_bc, "BC marginal")

r_bc = cor(tab_cond_bc$estimate, tab_marg_bc$estimate)
max_d_bc = max(abs(tab_cond_bc$estimate - tab_marg_bc$estimate))
cat(sprintf("  Cond vs Marg: r = %.4f, max|diff| = %.4f\n\n", r_bc, max_d_bc))

# ------------------------------------------------------------------
# 3d. Per-block divergence analysis
# ------------------------------------------------------------------
cat("--- 3d: Per-block divergence (medium network) --------------------\n")
blocks = unique(tab_cond_m$block)
block_stats = do.call(rbind, lapply(blocks, function(b) {
  idx = tab_cond_m$block == b
  d = tab_cond_m$estimate[idx] - tab_marg_m$estimate[idx]
  data.frame(
    block = b,
    n = sum(idx),
    mean_diff = mean(d),
    sd_diff = sd(d),
    max_abs_diff = max(abs(d)),
    cor = cor(tab_cond_m$estimate[idx], tab_marg_m$estimate[idx]),
    stringsAsFactors = FALSE
  )
}))
cat("Per-block conditional vs marginal divergence:\n")
print(block_stats, row.names = FALSE)

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
pdf(file.path(out_dir, "group3_cond_vs_marg.pdf"), width = 14, height = 14)

layout(matrix(1:9, nrow = 3, byrow = TRUE))
par(mar = c(4.5, 4.5, 2.5, 1))

# Row 1: Small network
recovery_scatter(tab_cond_s, "Small: conditional recovery")
recovery_scatter(tab_marg_s, "Small: marginal recovery")
agreement_scatter(tab_cond_s, tab_marg_s,
                  "Conditional", "Marginal",
                  sprintf("Small: cond vs marg (r=%.3f)", r_s))

# Row 2: Medium network
recovery_scatter(tab_cond_m, "Medium: conditional recovery")
recovery_scatter(tab_marg_m, "Medium: marginal recovery")
agreement_scatter(tab_cond_m, tab_marg_m,
                  "Conditional", "Marginal",
                  sprintf("Medium: cond vs marg (r=%.3f)", r_m))

# Row 3: BC network
recovery_scatter(tab_cond_bc, "BC: conditional recovery")
recovery_scatter(tab_marg_bc, "BC: marginal recovery")
agreement_scatter(tab_cond_bc, tab_marg_bc,
                  "Conditional", "Marginal",
                  sprintf("BC: cond vs marg (r=%.3f)", r_bc))
dev.off()

# Difference distribution plots
pdf(file.path(out_dir, "group3_diff_distributions.pdf"), width = 10, height = 8)
par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))

block_cols = c(mux = "#E41A1C", muy = "#377EB8", Kxx = "#4DAF4A",
               Kxy = "#984EA3", Kyy = "#FF7F00")

for(b in blocks) {
  idx = tab_cond_m$block == b
  d = tab_cond_m$estimate[idx] - tab_marg_m$estimate[idx]
  if(length(d) > 2) {
    hist(d, breaks = 15, col = adjustcolor(block_cols[b], 0.5),
         border = block_cols[b],
         main = sprintf("%s diffs (n=%d)", b, sum(idx)),
         xlab = "Conditional - Marginal")
    abline(v = 0, lty = 2, col = "grey30")
  } else {
    barplot(d, names.arg = tab_cond_m$parameter[idx],
            col = adjustcolor(block_cols[b], 0.5),
            main = sprintf("%s diffs (n=%d)", b, sum(idx)),
            ylab = "Conditional - Marginal")
    abline(h = 0, lty = 2, col = "grey30")
  }
}

# Overall
all_d = tab_cond_m$estimate - tab_marg_m$estimate
hist(all_d, breaks = 20, col = adjustcolor("grey60", 0.5),
     main = sprintf("All diffs (n=%d)", length(all_d)),
     xlab = "Conditional - Marginal")
abline(v = 0, lty = 2, col = "red")

dev.off()

cat(sprintf("\nPlots saved to %s/group3_*.pdf\n", out_dir))
cat("=== Group 3 complete =============================================\n\n")
