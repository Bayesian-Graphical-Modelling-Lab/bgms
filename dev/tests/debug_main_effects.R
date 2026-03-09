devtools::load_all(quiet = TRUE)

data("Boredom", package = "bgms")
x = as.matrix(Boredom[, 2:9])
x = x - 1L
colnames(x) = names(Boredom)[2:9]
n_obs = nrow(x)
vtype = c("ordinal", "continuous", "ordinal", "continuous",
          "ordinal", "continuous", "ordinal", "continuous")
seed = 44321

extract_main_flat = function(fit) {
  pm = fit$posterior_mean_main
  c(as.vector(t(pm$discrete)), as.vector(t(pm$continuous)))
}

# --- Marginal PL ---
cat("Fitting marginal PL...\n")
fit_marg = bgm(x, variable_type = vtype, pseudolikelihood = "marginal",
               edge_selection = FALSE, seed = seed,
               iter = 1000, warmup = 1000, chains = 1,
               display_progress = "none")

cat("=== Original marginal PL main effects ===\n")
cat("Discrete:\n"); print(round(fit_marg$posterior_mean_main$discrete, 3))
cat("Continuous:\n"); print(round(fit_marg$posterior_mean_main$continuous, 3))

set.seed(seed)
sim_marg = simulate(fit_marg, nsim = n_obs, method = "posterior-mean", seed = seed)

cat("\nRefitting on simulated data...\n")
refit_marg = bgm(sim_marg, variable_type = vtype, pseudolikelihood = "marginal",
                 edge_selection = FALSE, seed = seed,
                 iter = 1000, warmup = 1000, chains = 1,
                 display_progress = "none")

cat("=== Refit marginal PL main effects ===\n")
cat("Discrete:\n"); print(round(refit_marg$posterior_mean_main$discrete, 3))
cat("Continuous:\n"); print(round(refit_marg$posterior_mean_main$continuous, 3))

o = extract_main_flat(fit_marg)
r = extract_main_flat(refit_marg)
keep = !is.na(o) & !is.na(r)
cat("\n=== Flat comparison (non-NA) ===\n")
df = data.frame(original = round(o[keep], 3), refit = round(r[keep], 3),
                diff = round(r[keep] - o[keep], 3))
print(df)
cat("\nCorrelation:", round(cor(o[keep], r[keep]), 3), "\n")
cat("Original range:", round(range(o[keep]), 3), "\n")
cat("Refit range:   ", round(range(r[keep]), 3), "\n")

# --- Also compare with conditional PL for reference ---
cat("\n\n--- Conditional PL for comparison ---\n")
fit_cond = bgm(x, variable_type = vtype, pseudolikelihood = "conditional",
               edge_selection = FALSE, seed = seed,
               iter = 1000, warmup = 1000, chains = 1,
               display_progress = "none")

set.seed(seed)
sim_cond = simulate(fit_cond, nsim = n_obs, method = "posterior-mean", seed = seed)

refit_cond = bgm(sim_cond, variable_type = vtype, pseudolikelihood = "conditional",
                 edge_selection = FALSE, seed = seed,
                 iter = 1000, warmup = 1000, chains = 1,
                 display_progress = "none")

oc = extract_main_flat(fit_cond)
rc = extract_main_flat(refit_cond)
keepc = !is.na(oc) & !is.na(rc)
cat("\n=== Conditional PL flat comparison (non-NA) ===\n")
dfc = data.frame(original = round(oc[keepc], 3), refit = round(rc[keepc], 3),
                 diff = round(rc[keepc] - oc[keepc], 3))
print(dfc)
cat("\nCorrelation:", round(cor(oc[keepc], rc[keepc]), 3), "\n")
