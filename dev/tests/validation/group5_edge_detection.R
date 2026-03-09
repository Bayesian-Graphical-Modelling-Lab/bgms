# ===========================================================================
# Group 5: Edge detection accuracy
# ===========================================================================
# Tests Bayesian edge selection on known sparse networks. With edge_selection
# TRUE, the sampler reports posterior inclusion probabilities (PIPs). We
# evaluate sensitivity, specificity, precision, F1, and MCC at a 0.5 PIP
# threshold and also via ROC-like analysis.
#
# Conditions: weak edges (small J), medium edges, strong edges,
#             both PL methods, both samplers.
#
# Output: numerical summary + PDF with ROC-like curves and confusion matrices.
# ===========================================================================

devtools::load_all(quiet = TRUE)
source(file.path("dev", "tests", "validation", "helpers.R"))

out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")
cat("  Group 5: Edge detection accuracy\n")
cat("=======================================================================\n\n")

# ------------------------------------------------------------------
# Helper: compute edge detection metrics from PIPs and truth
# ------------------------------------------------------------------
# @param pip       Matrix of posterior inclusion probabilities (p+q x p+q).
# @param true_adj  Binary adjacency matrix (p+q x p+q).
# @param threshold PIP threshold for declaring an edge present.
#
# Returns: named list with TP, FP, FN, TN, sensitivity, specificity,
#          precision, F1, MCC.
# ------------------------------------------------------------------
edge_metrics = function(pip, true_adj, threshold = 0.5) {
  # Upper triangle only (avoid double counting)
  ut = upper.tri(pip)
  pred = ifelse(pip[ut] >= threshold, 1, 0)
  truth = ifelse(true_adj[ut] != 0, 1, 0)

  tp = sum(pred == 1 & truth == 1)
  fp = sum(pred == 1 & truth == 0)
  fn = sum(pred == 0 & truth == 1)
  tn = sum(pred == 0 & truth == 0)

  sens = if(tp + fn > 0) tp / (tp + fn) else NA
  spec = if(tn + fp > 0) tn / (tn + fp) else NA
  prec = if(tp + fp > 0) tp / (tp + fp) else NA
  f1 = if(!is.na(prec) && !is.na(sens) && (prec + sens) > 0)
    2 * prec * sens / (prec + sens) else NA
  denom = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc = if(denom > 0) (tp * tn - fp * fn) / denom else NA

  list(TP = tp, FP = fp, FN = fn, TN = tn,
       sensitivity = sens, specificity = spec,
       precision = prec, F1 = f1, MCC = mcc)
}

# ------------------------------------------------------------------
# Helper: ROC-like curve from PIPs
# ------------------------------------------------------------------
roc_curve = function(pip, true_adj) {
  ut = upper.tri(pip)
  scores = pip[ut]
  truth = ifelse(true_adj[ut] != 0, 1, 0)
  thresholds = sort(unique(c(0, scores, 1)))

  tpr = fpr = numeric(length(thresholds))
  for(i in seq_along(thresholds)) {
    pred = ifelse(scores >= thresholds[i], 1, 0)
    tp = sum(pred == 1 & truth == 1)
    fp = sum(pred == 1 & truth == 0)
    fn = sum(pred == 0 & truth == 1)
    tn = sum(pred == 0 & truth == 0)
    tpr[i] = if(tp + fn > 0) tp / (tp + fn) else 0
    fpr[i] = if(fp + tn > 0) fp / (fp + tn) else 0
  }

  # AUC by trapezoidal rule
  ord = order(fpr, tpr)
  fpr_s = fpr[ord]; tpr_s = tpr[ord]
  auc = sum(diff(fpr_s) * (tpr_s[-1] + tpr_s[-length(tpr_s)]) / 2)

  list(fpr = fpr_s, tpr = tpr_s, auc = auc)
}

# ------------------------------------------------------------------
# Helper: create true adjacency matrix from network
# ------------------------------------------------------------------
true_adjacency = function(net) {
  p = net$p; q = net$q; tot = p + q
  adj = matrix(0, tot, tot)
  adj[seq_len(p), seq_len(p)] = (net$Kxx != 0) * 1
  adj[p + seq_len(q), p + seq_len(q)] = (net$Kyy != 0) * 1
  diag(adj[p + seq_len(q), p + seq_len(q)]) = 0  # diag doesn't count
  adj[seq_len(p), p + seq_len(q)] = (net$Kxy != 0) * 1
  adj[p + seq_len(q), seq_len(p)] = t(adj[seq_len(p), p + seq_len(q)])
  adj
}

# ------------------------------------------------------------------
# 5a. Define edge-strength scenarios
# ------------------------------------------------------------------
scenarios = list(
  weak = list(density = 0.4, scale = 0.3, seed = 701, label = "weak (J~0.3)"),
  medium = list(density = 0.4, scale = 0.6, seed = 702, label = "medium (J~0.6)"),
  strong = list(density = 0.4, scale = 1.0, seed = 703, label = "strong (J~1.0)")
)

# Modified make_network with controllable edge strength
make_network_scaled = function(p, q, n_cat, density, scale, seed,
                               variable_type = rep("ordinal", p),
                               baseline_category = rep(0L, p)) {
  net = make_network(p, q, n_cat, variable_type = variable_type,
                     baseline_category = baseline_category, density = density,
                     seed = seed)
  # Rescale non-zero edges
  net$Kxx = net$Kxx * scale / max(abs(net$Kxx[net$Kxx != 0]), 1)
  net$Kxy = net$Kxy * scale / max(abs(net$Kxy[net$Kxy != 0]), 1)
  off_diag_yy = net$Kyy; diag(off_diag_yy) = 0
  if(any(off_diag_yy != 0)) {
    net$Kyy[upper.tri(net$Kyy)] = net$Kyy[upper.tri(net$Kyy)] * scale / max(abs(off_diag_yy[off_diag_yy != 0]), 1)
    net$Kyy[lower.tri(net$Kyy)] = t(net$Kyy)[lower.tri(net$Kyy)]
  }
  # Re-enforce diagonal dominance for Kyy
  diag(net$Kyy) = abs(rowSums(net$Kyy - diag(diag(net$Kyy)))) + runif(q, 1.0, 2.0)
  net
}

results = list()
roc_data = list()

for(sc_name in names(scenarios)) {
  sc = scenarios[[sc_name]]
  cat(sprintf("--- 5_%s: %s ---\n", sc_name, sc$label))

  net = make_network_scaled(p = 4, q = 3, n_cat = c(1L, 2L, 3L, 1L),
                            density = sc$density, scale = sc$scale, seed = sc$seed)
  dat = generate_data(net, n = 3000, source = "bgms", seed = sc$seed + 100)
  vtype = c(rep("ordinal", 4), rep("continuous", 3))
  adj = true_adjacency(net)

  n_true_edges = sum(adj[upper.tri(adj)])
  n_possible = sum(upper.tri(adj))
  cat(sprintf("  True edges: %d / %d possible\n", n_true_edges, n_possible))

  # Fit with edge selection
  fit = bgm(
    dat, variable_type = vtype,
    pseudolikelihood = "marginal",
    edge_selection = TRUE,
    iter = 10000, warmup = 5000, chains = 2,
    seed = sc$seed + 200
  )

  pip = fit$posterior_mean_indicator
  met = edge_metrics(pip, adj, threshold = 0.5)
  roc = roc_curve(pip, adj)

  results[[sc_name]] = data.frame(
    scenario = sc$label,
    n_true_edges = n_true_edges,
    TP = met$TP, FP = met$FP, FN = met$FN, TN = met$TN,
    sensitivity = round(met$sensitivity, 3),
    specificity = round(met$specificity, 3),
    precision = round(met$precision, 3),
    F1 = round(met$F1, 3),
    MCC = round(met$MCC, 3),
    AUC = round(roc$auc, 3),
    stringsAsFactors = FALSE
  )
  roc_data[[sc_name]] = roc

  cat(sprintf("  Sens = %.3f | Spec = %.3f | F1 = %.3f | MCC = %.3f | AUC = %.3f\n",
              met$sensitivity, met$specificity, met$F1, met$MCC, roc$auc))
}

# ------------------------------------------------------------------
# 5b. Blume-Capel edge detection (medium signal strength)
# ------------------------------------------------------------------
cat("\n--- 5b: BC mixed — edge detection (medium signal) ----------------\n")

net_bc = make_network_scaled(
  p = 4, q = 2,
  n_cat = c(1L, 2L, 2L, 3L),
  density = 0.4, scale = 0.6, seed = 710,
  variable_type = c("ordinal", "ordinal", "blume-capel", "blume-capel"),
  baseline_category = c(0L, 0L, 1L, 1L)
)
dat_bc = generate_data(net_bc, n = 3000, source = "bgms", seed = 810)

vtype_bc = c("ordinal", "ordinal", "blume-capel", "blume-capel",
             "continuous", "continuous")
adj_bc = true_adjacency(net_bc)

n_true_bc = sum(adj_bc[upper.tri(adj_bc)])
n_poss_bc = sum(upper.tri(adj_bc))
cat(sprintf("  True edges: %d / %d possible\n", n_true_bc, n_poss_bc))

fit_bc = bgm(
  dat_bc, variable_type = vtype_bc,
  baseline_category = c(0L, 0L, 1L, 1L),
  pseudolikelihood = "marginal",
  edge_selection = TRUE,
  iter = 10000, warmup = 5000, chains = 2,
  seed = 910
)

pip_bc = fit_bc$posterior_mean_indicator
met_bc = edge_metrics(pip_bc, adj_bc, threshold = 0.5)
roc_bc = roc_curve(pip_bc, adj_bc)

results[["bc_medium"]] = data.frame(
  scenario = "BC medium (J~0.6)",
  n_true_edges = n_true_bc,
  TP = met_bc$TP, FP = met_bc$FP, FN = met_bc$FN, TN = met_bc$TN,
  sensitivity = round(met_bc$sensitivity, 3),
  specificity = round(met_bc$specificity, 3),
  precision = round(met_bc$precision, 3),
  F1 = round(met_bc$F1, 3),
  MCC = round(met_bc$MCC, 3),
  AUC = round(roc_bc$auc, 3),
  stringsAsFactors = FALSE
)
roc_data[["bc_medium"]] = roc_bc

cat(sprintf("  Sens = %.3f | Spec = %.3f | F1 = %.3f | MCC = %.3f | AUC = %.3f\n",
            met_bc$sensitivity, met_bc$specificity, met_bc$F1, met_bc$MCC, roc_bc$auc))

# ------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------
cat("\n--- Edge detection summary ----------------------------------------\n")
summary_df = do.call(rbind, results)
print(summary_df, row.names = FALSE)

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
pdf(file.path(out_dir, "group5_edge_detection.pdf"), width = 14, height = 14)

layout(matrix(1:8, nrow = 2, byrow = TRUE))
par(mar = c(4.5, 4.5, 2.5, 1))

# Row 1: ROC curves + PIP panels for ordinal scenarios
cols = c(weak = "#E41A1C", medium = "#377EB8", strong = "#4DAF4A",
         bc_medium = "#984EA3")

# Combined ROC (all scenarios including BC)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = "ROC curves by scenario", asp = 1)
abline(0, 1, lty = 2, col = "grey60")
for(sc_name in names(roc_data)) {
  roc = roc_data[[sc_name]]
  lines(roc$fpr, roc$tpr, col = cols[sc_name], lwd = 2)
}
legend("bottomright",
       legend = paste0(names(roc_data), sprintf(" (AUC=%.2f)", sapply(roc_data, `[[`, "auc"))),
       col = cols[names(roc_data)], lwd = 2, bty = "n")

# PIP panels for ordinal scenarios
for(sc_name in names(scenarios)) {
  sc = scenarios[[sc_name]]
  net = make_network_scaled(p = 4, q = 3, n_cat = c(1L, 2L, 3L, 1L),
                            density = sc$density, scale = sc$scale, seed = sc$seed)
  dat = generate_data(net, n = 3000, source = "bgms", seed = sc$seed + 100)
  vtype = c(rep("ordinal", 4), rep("continuous", 3))

  fit = bgm(dat, variable_type = vtype,
            pseudolikelihood = "marginal", edge_selection = TRUE,
            iter = 10000, warmup = 5000, chains = 2,
            display_progress = "none", seed = sc$seed + 200)

  pip = fit$posterior_mean_indicator
  adj = true_adjacency(net)

  ut_pip = pip[upper.tri(pip)]
  ut_adj = adj[upper.tri(adj)]
  plot(jitter(ut_adj, 0.1), ut_pip, pch = 19,
       col = adjustcolor(ifelse(ut_adj == 1, "firebrick", "steelblue"), 0.6),
       xlab = "True edge (0/1)", ylab = "Posterior inclusion prob",
       main = sprintf("%s: PIPs", sc$label),
       xlim = c(-0.3, 1.3), ylim = c(0, 1))
  abline(h = 0.5, lty = 2, col = "grey40")
  legend("center",
         legend = c("True edge", "True null"),
         col = c("firebrick", "steelblue"), pch = 19, bty = "n")
}

# Row 2: BC PIP panel + barplot
# BC PIP panel
ut_pip_bc = pip_bc[upper.tri(pip_bc)]
ut_adj_bc = adj_bc[upper.tri(adj_bc)]
plot(jitter(ut_adj_bc, 0.1), ut_pip_bc, pch = 19,
     col = adjustcolor(ifelse(ut_adj_bc == 1, "firebrick", "steelblue"), 0.6),
     xlab = "True edge (0/1)", ylab = "Posterior inclusion prob",
     main = "BC medium: PIPs",
     xlim = c(-0.3, 1.3), ylim = c(0, 1))
abline(h = 0.5, lty = 2, col = "grey40")
legend("center",
       legend = c("True edge", "True null"),
       col = c("firebrick", "steelblue"), pch = 19, bty = "n")

# Barplot of metrics (all scenarios)
all_cols = cols[names(results)]
metrics_mat = t(sapply(results, function(r)
  c(Sensitivity = r$sensitivity, Specificity = r$specificity,
    F1 = r$F1, MCC = r$MCC)))
barplot(metrics_mat, beside = TRUE, col = all_cols,
        main = "Edge detection metrics by scenario",
        ylab = "Score", ylim = c(0, 1.1))
legend("topright", legend = names(results), fill = all_cols, bty = "n")
abline(h = 1, lty = 2, col = "grey60")

# Blank panels
plot.new()
plot.new()

dev.off()

cat(sprintf("\nPlots saved to %s/group5_edge_detection.pdf\n", out_dir))
cat("=== Group 5 complete =============================================\n\n")
