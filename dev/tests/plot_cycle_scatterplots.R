# ==============================================================================
# Scatterplots for mixed MRF simulation-recovery cycle tests
# ==============================================================================
# Fits mixed MRF with conditional and marginal PL on the full Boredom data
# (8 items, 4 ordinal + 4 continuous), simulates from each fit, refits on
# simulated data, and plots original vs refit estimates for both pairwise
# interactions and main effects.
# ==============================================================================

devtools::load_all(quiet = TRUE)

data("Boredom", package = "bgms")
x = as.matrix(Boredom[, 2:9])
x = x - 1L
colnames(x) = names(Boredom)[2:9]
n_obs = nrow(x)
vtype = c("ordinal", "continuous", "ordinal", "continuous",
          "ordinal", "continuous", "ordinal", "continuous")
seed = 44321

# Helper: extract a flat main-effects vector from a mixed fit
extract_main_flat = function(fit) {
  pm = fit$posterior_mean_main
  c(as.vector(t(pm$discrete)),  # discrete thresholds (row-major)
    as.vector(t(pm$continuous))) # continuous mean + precision (row-major)
}

# ==============================================================================
# Conditional PL
# ==============================================================================
cat("Fitting conditional PL...\n")
fit_cond = bgm(x, variable_type = vtype, pseudolikelihood = "conditional",
               edge_selection = FALSE, seed = seed,
               iter = 1000, warmup = 1000, chains = 1,
               display_progress = "none")

orig_pw_cond = colMeans(extract_pairwise_interactions(fit_cond))
orig_main_cond = extract_main_flat(fit_cond)

cat("Simulating from conditional PL fit...\n")
set.seed(seed)
sim_cond = simulate(fit_cond, nsim = n_obs, method = "posterior-mean",
                    seed = seed)

cat("Refitting on conditional PL simulated data...\n")
refit_cond = bgm(sim_cond, variable_type = vtype,
                 pseudolikelihood = "conditional",
                 edge_selection = FALSE, seed = seed,
                 iter = 1000, warmup = 1000, chains = 1,
                 display_progress = "none")

refit_pw_cond = colMeans(extract_pairwise_interactions(refit_cond))
refit_main_cond = extract_main_flat(refit_cond)
cor_pw_cond = cor(orig_pw_cond, refit_pw_cond)
cor_main_cond = cor(orig_main_cond, refit_main_cond, use = "complete.obs")
cat(sprintf("Conditional PL — pairwise r = %.3f, main r = %.3f\n",
            cor_pw_cond, cor_main_cond))

# ==============================================================================
# Marginal PL
# ==============================================================================
cat("Fitting marginal PL...\n")
fit_marg = bgm(x, variable_type = vtype, pseudolikelihood = "marginal",
               edge_selection = FALSE, seed = seed,
               iter = 1000, warmup = 1000, chains = 1,
               display_progress = "none")

orig_pw_marg = colMeans(extract_pairwise_interactions(fit_marg))
orig_main_marg = extract_main_flat(fit_marg)

cat("Simulating from marginal PL fit...\n")
set.seed(seed)
sim_marg = simulate(fit_marg, nsim = n_obs, method = "posterior-mean",
                    seed = seed)

cat("Refitting on marginal PL simulated data...\n")
refit_marg = bgm(sim_marg, variable_type = vtype,
                 pseudolikelihood = "marginal",
                 edge_selection = FALSE, seed = seed,
                 iter = 1000, warmup = 1000, chains = 1,
                 display_progress = "none")

refit_pw_marg = colMeans(extract_pairwise_interactions(refit_marg))
refit_main_marg = extract_main_flat(refit_marg)
cor_pw_marg = cor(orig_pw_marg, refit_pw_marg)
cor_main_marg = cor(orig_main_marg, refit_main_marg, use = "complete.obs")
cat(sprintf("Marginal PL — pairwise r = %.3f, main r = %.3f\n",
            cor_pw_marg, cor_main_marg))

# ==============================================================================
# Plot: 2 x 2 grid (pairwise / main) x (conditional / marginal)
# ==============================================================================
out_file = file.path("dev", "tests", "cycle_scatterplots.pdf")
pdf(out_file, width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))

# Row 1: pairwise
rng_pw = range(c(orig_pw_cond, refit_pw_cond, orig_pw_marg, refit_pw_marg))

plot(orig_pw_cond, refit_pw_cond,
     xlab = "Original pairwise estimates",
     ylab = "Refit pairwise estimates",
     main = sprintf("Conditional PL — pairwise  (r = %.3f)", cor_pw_cond),
     pch = 19, col = adjustcolor("steelblue", 0.7),
     xlim = rng_pw, ylim = rng_pw, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

plot(orig_pw_marg, refit_pw_marg,
     xlab = "Original pairwise estimates",
     ylab = "Refit pairwise estimates",
     main = sprintf("Marginal PL — pairwise  (r = %.3f)", cor_pw_marg),
     pch = 19, col = adjustcolor("firebrick", 0.7),
     xlim = rng_pw, ylim = rng_pw, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

# Row 2: main effects (drop NAs from padding)
keep_cond = !is.na(orig_main_cond) & !is.na(refit_main_cond)
keep_marg = !is.na(orig_main_marg) & !is.na(refit_main_marg)
rng_main = range(c(orig_main_cond[keep_cond], refit_main_cond[keep_cond],
                    orig_main_marg[keep_marg], refit_main_marg[keep_marg]))

plot(orig_main_cond[keep_cond], refit_main_cond[keep_cond],
     xlab = "Original main-effect estimates",
     ylab = "Refit main-effect estimates",
     main = sprintf("Conditional PL — main  (r = %.3f)", cor_main_cond),
     pch = 19, col = adjustcolor("steelblue", 0.7),
     xlim = rng_main, ylim = rng_main, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

plot(orig_main_marg[keep_marg], refit_main_marg[keep_marg],
     xlab = "Original main-effect estimates",
     ylab = "Refit main-effect estimates",
     main = sprintf("Marginal PL — main  (r = %.3f)", cor_main_marg),
     pch = 19, col = adjustcolor("firebrick", 0.7),
     xlim = rng_main, ylim = rng_main, asp = 1)
abline(0, 1, lty = 2, col = "grey40")

dev.off()
cat(sprintf("Saved to %s\n", out_file))
