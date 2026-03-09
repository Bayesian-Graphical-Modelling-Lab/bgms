# ==============================================================================
# Compare bgms vs mixedGM marginal PL cycle test
# ==============================================================================
# Uses Boredom data (8 items: 4 ordinal + 4 continuous).
# Fits marginal PL with both packages, simulates, refits, compares main effects.
# ==============================================================================

devtools::load_all(quiet = TRUE)            # bgms
library(mixedGM)

data("Boredom", package = "bgms")
raw = as.matrix(Boredom[, 2:9])
raw = raw - 1L                              # shift 1-7 to 0-6
colnames(raw) = names(Boredom)[2:9]
n_obs = nrow(raw)

# Assign variable types: alternating ordinal / continuous
ord_idx = c(1, 3, 5, 7)
cont_idx = c(2, 4, 6, 8)

x_disc = raw[, ord_idx, drop = FALSE]       # n x 4 integer ordinal (0-6)
y_cont = raw[, cont_idx, drop = FALSE]      # n x 4 integer-valued, treated as continuous
storage.mode(y_cont) = "double"

num_cats_disc = apply(x_disc, 2, max) + 1L  # total categories (mixedGM convention)
p = ncol(x_disc)
q = ncol(y_cont)
seed = 44321

cat(sprintf("p = %d ordinal, q = %d continuous, n = %d\n", p, q, n_obs))
cat("num_categories (max index):", num_cats_disc, "\n\n")

# ==============================================================================
# mixedGM marginal PL cycle
# ==============================================================================
cat("=== mixedGM: marginal PL ===\n")
set.seed(seed)
fit_mg = mixed_sampler(x = x_disc, y = y_cont,
                       num_categories = num_cats_disc,
                       pseudolikelihood = "marginal",
                       log_prior_mean = mixedGM:::log_prior_logistic,
                       edge_selection = FALSE,
                       n_warmup = 5000, n_samples = 5000,
                       verbose = FALSE)

mux_orig_mg = apply(fit_mg$samples$mux, c(2, 3), mean)
muy_orig_mg = colMeans(fit_mg$samples$muy)
Kxx_orig_mg = apply(fit_mg$samples$Kxx, c(2, 3), mean)
Kyy_orig_mg = apply(fit_mg$samples$Kyy, c(2, 3), mean)
Kxy_orig_mg = apply(fit_mg$samples$Kxy, c(2, 3), mean)

cat("Original mux:\n"); print(round(mux_orig_mg, 3))
cat("Original muy:", round(muy_orig_mg, 3), "\n\n")

# Simulate from posterior means
sim_mg = mixed_gibbs_generate(
  n = n_obs,
  Kxx = Kxx_orig_mg, Kxy = Kxy_orig_mg, Kyy = Kyy_orig_mg,
  mux = mux_orig_mg, muy = muy_orig_mg,
  num_categories = num_cats_disc,
  n_burnin = 1000
)

# Refit
set.seed(seed)
refit_mg = mixed_sampler(x = sim_mg$x, y = sim_mg$y,
                         num_categories = num_cats_disc,
                         pseudolikelihood = "marginal",
                         log_prior_mean = mixedGM:::log_prior_logistic,
                         edge_selection = FALSE,
                         n_warmup = 5000, n_samples = 5000,
                         verbose = FALSE)

mux_refit_mg = apply(refit_mg$samples$mux, c(2, 3), mean)
muy_refit_mg = colMeans(refit_mg$samples$muy)

cat("Refit mux:\n"); print(round(mux_refit_mg, 3))
cat("Refit muy:", round(muy_refit_mg, 3), "\n\n")

# Flat comparison
orig_flat_mg = c(as.vector(t(mux_orig_mg)), muy_orig_mg)
refit_flat_mg = c(as.vector(t(mux_refit_mg)), muy_refit_mg)
keep_mg = !is.na(orig_flat_mg) & !is.na(refit_flat_mg)

cat("=== mixedGM main-effect comparison ===\n")
df_mg = data.frame(original = round(orig_flat_mg[keep_mg], 3),
                   refit = round(refit_flat_mg[keep_mg], 3),
                   diff = round(refit_flat_mg[keep_mg] - orig_flat_mg[keep_mg], 3))
print(df_mg)
cat("Correlation:", round(cor(orig_flat_mg[keep_mg], refit_flat_mg[keep_mg]), 3), "\n")
cat("Max |diff|:", round(max(abs(df_mg$diff)), 3), "\n\n")

# ==============================================================================
# bgms marginal PL cycle (for comparison)
# ==============================================================================
cat("=== bgms: marginal PL ===\n")
vtype = c("ordinal", "continuous", "ordinal", "continuous",
          "ordinal", "continuous", "ordinal", "continuous")

fit_bg = bgm(raw, variable_type = vtype, pseudolikelihood = "marginal",
             edge_selection = FALSE, seed = seed,
             iter = 5000, warmup = 5000, chains = 1,
             display_progress = "none")

set.seed(seed)
sim_bg = simulate(fit_bg, nsim = n_obs, method = "posterior-mean", seed = seed)

refit_bg = bgm(sim_bg, variable_type = vtype, pseudolikelihood = "marginal",
               edge_selection = FALSE, seed = seed,
               iter = 5000, warmup = 5000, chains = 1,
               display_progress = "none")

extract_main_flat = function(fit) {
  pm = fit$posterior_mean_main
  c(as.vector(t(pm$discrete)), as.vector(t(pm$continuous)))
}
orig_flat_bg = extract_main_flat(fit_bg)
refit_flat_bg = extract_main_flat(refit_bg)
keep_bg = !is.na(orig_flat_bg) & !is.na(refit_flat_bg)

cat("=== bgms main-effect comparison ===\n")
df_bg = data.frame(original = round(orig_flat_bg[keep_bg], 3),
                   refit = round(refit_flat_bg[keep_bg], 3),
                   diff = round(refit_flat_bg[keep_bg] - orig_flat_bg[keep_bg], 3))
print(df_bg)
cat("Correlation:", round(cor(orig_flat_bg[keep_bg], refit_flat_bg[keep_bg]), 3), "\n")
cat("Max |diff|:", round(max(abs(df_bg$diff)), 3), "\n\n")

# ==============================================================================
# Also compare original estimates between bgms and mixedGM
# ==============================================================================
cat("=== Original estimates: bgms vs mixedGM ===\n")
cat("mux comparison (bgms discrete thresholds vs mixedGM mux):\n")
bgms_mux = fit_bg$posterior_mean_main$discrete
cat("bgms:\n"); print(round(bgms_mux, 3))
cat("mixedGM:\n"); print(round(mux_orig_mg, 3))
cat("Max |diff|:", round(max(abs(bgms_mux[!is.na(bgms_mux)] - mux_orig_mg[!is.na(mux_orig_mg)])), 3), "\n\n")

cat("muy comparison:\n")
bgms_muy = fit_bg$posterior_mean_main$continuous[, "mean"]
cat("bgms:", round(bgms_muy, 3), "\n")
cat("mixedGM:", round(muy_orig_mg, 3), "\n")
cat("Max |diff|:", round(max(abs(bgms_muy - muy_orig_mg)), 3), "\n")

# ==============================================================================
# Scatterplots: 2x2 grid (mixedGM / bgms) x (pairwise / main)
# ==============================================================================
# Also extract pairwise from mixedGM
Kxx_refit_mg = apply(refit_mg$samples$Kxx, c(2, 3), mean)
Kyy_refit_mg = apply(refit_mg$samples$Kyy, c(2, 3), mean)
Kxy_refit_mg = apply(refit_mg$samples$Kxy, c(2, 3), mean)

# Flatten upper-tri pairwise for mixedGM
flatten_pw_mg = function(Kxx, Kyy, Kxy) {
  c(Kxx[upper.tri(Kxx)], Kyy[upper.tri(Kyy)], as.vector(Kxy))
}
orig_pw_mg = flatten_pw_mg(Kxx_orig_mg, Kyy_orig_mg, Kxy_orig_mg)
refit_pw_mg = flatten_pw_mg(Kxx_refit_mg, Kyy_refit_mg, Kxy_refit_mg)

# bgms pairwise
orig_pw_bg = colMeans(extract_pairwise_interactions(fit_bg))
refit_pw_bg = colMeans(extract_pairwise_interactions(refit_bg))

cor_pw_mg = cor(orig_pw_mg, refit_pw_mg)
cor_pw_bg = cor(orig_pw_bg, refit_pw_bg)
cor_main_mg = cor(orig_flat_mg[keep_mg], refit_flat_mg[keep_mg])
cor_main_bg = cor(orig_flat_bg[keep_bg], refit_flat_bg[keep_bg])

out_file = file.path("dev", "tests", "cycle_scatterplots.pdf")
pdf(out_file, width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))

# Row 1: pairwise
rng_pw = range(c(orig_pw_mg, refit_pw_mg, orig_pw_bg, refit_pw_bg))

plot(orig_pw_mg, refit_pw_mg,
     xlab = "Original pairwise estimates",
     ylab = "Refit pairwise estimates",
     main = sprintf("mixedGM marginal PL - pairwise (r = %.3f)", cor_pw_mg),
     pch = 19, col = adjustcolor("firebrick", 0.7),
     xlim = rng_pw, ylim = rng_pw, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

plot(orig_pw_bg, refit_pw_bg,
     xlab = "Original pairwise estimates",
     ylab = "Refit pairwise estimates",
     main = sprintf("bgms marginal PL - pairwise (r = %.3f)", cor_pw_bg),
     pch = 19, col = adjustcolor("steelblue", 0.7),
     xlim = rng_pw, ylim = rng_pw, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

# Row 2: main effects
rng_main = range(c(orig_flat_mg[keep_mg], refit_flat_mg[keep_mg],
                    orig_flat_bg[keep_bg], refit_flat_bg[keep_bg]))

plot(orig_flat_mg[keep_mg], refit_flat_mg[keep_mg],
     xlab = "Original main-effect estimates",
     ylab = "Refit main-effect estimates",
     main = sprintf("mixedGM marginal PL - main (r = %.3f)", cor_main_mg),
     pch = 19, col = adjustcolor("firebrick", 0.7),
     xlim = rng_main, ylim = rng_main, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

plot(orig_flat_bg[keep_bg], refit_flat_bg[keep_bg],
     xlab = "Original main-effect estimates",
     ylab = "Refit main-effect estimates",
     main = sprintf("bgms marginal PL - main (r = %.3f)", cor_main_bg),
     pch = 19, col = adjustcolor("steelblue", 0.7),
     xlim = rng_main, ylim = rng_main, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

dev.off()
cat(sprintf("\nSaved to %s\n", out_file))
