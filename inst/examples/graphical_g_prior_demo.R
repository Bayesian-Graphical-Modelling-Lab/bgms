# ============================================================================
# Demo: Graphical G-prior + inclusion Bayes factors
# ============================================================================
#
# A walk-through of the new GG-prior machinery on the feat/graphical-g-prior
# branch. Sections are independent; you can step through with Cmd-Enter /
# Ctrl-Enter or source the whole file.
#
# Sections:
#   0. Setup
#   1. Fit with the default g-hyperprior (conjugate Gamma)
#   2. Inspect the g-trace
#   3. The inclusion_bayes_factor() extractor
#   4. Compare the six g-hyperprior variants on the same data
#   5. Compare GG-prior vs. plain Cauchy(0, 2.5) on equal footing
#   6. Customised tCCH (the umbrella family)
# ============================================================================


# ---- 0. Setup ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(bgms)        # load the package (or devtools::load_all(".")
  library(MASS)        # mvrnorm for data generation
})

set.seed(42)
p <- 5L
n <- 250L

# Truth: one strong edge between V1 and V2; another weaker edge V3-V4.
K_true        <- diag(p)
K_true[1, 2]  <- K_true[2, 1] <- -0.6   # partial association ~ 0.6
K_true[3, 4]  <- K_true[4, 3] <- -0.3   # partial association ~ 0.3
Sigma         <- solve(K_true)
X             <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
colnames(X)   <- paste0("V", seq_len(p))


# ---- 1. Fit with the default g-hyperprior -----------------------------------

# graphical_g_prior() defaults to the conjugate-Gamma hyperprior on
# t = sqrt(g) with Gamma(a0 = 1, b0 = 1). The slab variance is
# V_ij = n / (4 * S_ii_bar * S_jj_bar) — data-adaptive per edge.

fit_cg <- bgm(
  x = X,
  variable_type = "continuous",
  interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                        a0 = 1, b0 = 1, g_init = 1),
  edge_selection = TRUE,
  iter = 2000L, warmup = 1000L, chains = 2L,
  display_progress = "none"
)

# Standard outputs: posterior PIPs, posterior partial-association means.
round(fit_cg$posterior_mean_indicator, 3)
round(fit_cg$posterior_mean_pairwise,  3)


# ---- 2. Inspect the g-trace -------------------------------------------------

# fit$gg_diagnostics$g[[c]] is the per-iteration trajectory of
# g = t^2 for chain c. Use it to diagnose mixing.
str(fit_cg$gg_diagnostics, max.level = 2)

# Quick eyeball:
op <- par(mfrow = c(1, 2))
plot(fit_cg$gg_diagnostics$g[[1]], type = "l",
     xlab = "iter", ylab = "g", main = "Chain 1: g trace")
plot(fit_cg$gg_diagnostics$g[[2]], type = "l",
     xlab = "iter", ylab = "g", main = "Chain 2: g trace")
par(op)


# ---- 3. Prior and posterior inclusion probability extractors ---------------

# Posterior PIPs: p x p symmetric matrix (already in the fit, but the
# extractor gives you a consistent named matrix).
post_pip <- extract_posterior_inclusion_probabilities(fit_cg)
round(post_pip, 3)

# Prior PIPs: same shape. For graphical_g_prior the matrix is estimated
# by a one-time prior-only chain on the same V_ij(Y) — the result is
# cached on the fit. For other priors (cauchy, normal, ...) the matrix
# is filled analytically from the edge prior (here Bernoulli(0.5)).
prior_pip <- extract_prior_inclusion_probabilities(fit_cg)
round(prior_pip, 3)

# Inclusion BF is a one-liner from the two:
post_odds  <- post_pip  / (1 - post_pip)
prior_odds <- prior_pip / (1 - prior_pip)
bf_cg <- post_odds / prior_odds
round(bf_cg, 3)

# NOTE: under heavy-tailed g-hyperpriors (Zellner-Siow, hyper-g, tCCH,
# ...) the simulated prior PIPs can saturate near 1 because the chain
# explores large g where the joint prior favours all edges — the
# Lindley-Bartlett paradox. Expect a Z(γ)-based analytic replacement
# once the mNLO method lands.


# ---- 4. Compare the six g-hyperprior variants -------------------------------

# All five live g-hyperprior settings (plus fixed) on the same data.
# Compare g-trajectories and BFs side by side.

hyperpriors <- list(
  fixed_g1        = graphical_g_prior(g_hyperprior = "fixed", g_fixed = 1),
  conjugate_gamma = graphical_g_prior(g_hyperprior = "conjugate_gamma"),
  zellner_siow    = graphical_g_prior(g_hyperprior = "zellner_siow", b = 1),
  hyper_g         = graphical_g_prior(g_hyperprior = "hyper_g",      a = 3),
  hyper_g_over_n  = graphical_g_prior(g_hyperprior = "hyper_g_over_n", a = 3),
  tCCH_default    = graphical_g_prior(g_hyperprior = "tCCH")
)

run_one <- function(prior_obj, seed = 1L) {
  fit <- bgm(
    x = X, variable_type = "continuous",
    interaction_prior = prior_obj,
    edge_selection = TRUE,
    iter = 1500L, warmup = 750L, chains = 1L,
    seed = seed, display_progress = "none"
  )
  list(
    g_trace = fit$gg_diagnostics$g[[1]],
    post    = extract_posterior_inclusion_probabilities(fit),
    prior   = extract_prior_inclusion_probabilities(fit)
  )
}

results <- lapply(hyperpriors, run_one)

# Summary table: posterior PIP on V1-V2 (strong) and V3-V4 (weak)
# across hyperpriors.
pip_table <- data.frame(
  hyperprior  = names(results),
  mean_g      = vapply(results, function(r) mean(r$g_trace), numeric(1)),
  post_V1_V2  = vapply(results, function(r) r$post[1, 2],    numeric(1)),
  post_V3_V4  = vapply(results, function(r) r$post[3, 4],    numeric(1)),
  prior_V1_V2 = vapply(results, function(r) r$prior[1, 2],   numeric(1))
)
print(pip_table, row.names = FALSE)


# ---- 5. GG-prior vs. plain Cauchy(0, 2.5) ----------------------------------

# Both extractors work for any continuous prior, so you can put
# data-adaptive (GG-prior) and data-independent (Cauchy) priors side
# by side. For the Cauchy fit the prior PIP matrix is constant at the
# edge-prior marginal (no simulation needed).

fit_cauchy <- bgm(
  x = X, variable_type = "continuous",
  interaction_prior = cauchy_prior(scale = 2.5),
  edge_selection = TRUE,
  iter = 1500L, warmup = 750L, chains = 1L,
  seed = 1L, display_progress = "none"
)
post_cauchy  <- extract_posterior_inclusion_probabilities(fit_cauchy)
prior_cauchy <- extract_prior_inclusion_probabilities(fit_cauchy)

cat("\nCauchy(0, 2.5):\n")
cat(sprintf("  posterior PIP[V1-V2] = %.3f   prior PIP[V1-V2] = %.3f\n",
            post_cauchy[1, 2], prior_cauchy[1, 2]))
cat("\nConjugateGamma GG-prior:\n")
cat(sprintf("  posterior PIP[V1-V2] = %.3f   prior PIP[V1-V2] = %.3f\n",
            post_pip[1, 2], prior_pip[1, 2]))


# ---- 6. Customised tCCH (the umbrella family) -------------------------------

# tCCH is parameterised as
#   log pi(g) = (b-1) log g - r log(1 + g/u) - (a+b) log(1+g) - s g
# with a, b > 0 and r, s >= 0, u > 0. Defaults reduce to (1+g)^(-2).
# Special cases (BAS conventions):
#   hyper_g(alpha):     a = alpha - 2, b = 1, r = 0, s = 0, u = 1
#   beta_prime(alpha,beta): a = beta, b = alpha, r = 0, s = 0, u = 1
#   "robust prior":     b = 1, r = 1.5, s = 0, u = n/(p+1)

fit_tcch <- bgm(
  x = X, variable_type = "continuous",
  interaction_prior = graphical_g_prior(
    g_hyperprior = "tCCH",
    tcch = list(a = 1, b = 1, r = 1.5, s = 0, u = n / (p + 1))   # robust-style
  ),
  edge_selection = TRUE,
  iter = 1500L, warmup = 750L, chains = 1L,
  seed = 1L, display_progress = "none"
)
extract_posterior_inclusion_probabilities(fit_tcch)
extract_prior_inclusion_probabilities(fit_tcch)


# ============================================================================
# That's the tour. The main hooks for the simulation study are
#   - graphical_g_prior(g_hyperprior = ..., a = ..., b = ..., tcch = ...)
#   - fit$gg_diagnostics$g  for shrinkage-scale diagnostics
#   - extract_posterior_inclusion_probabilities(fit)  posterior PIP (p x p)
#   - extract_prior_inclusion_probabilities(fit)      prior PIP (p x p)
# BF is just (post_odds / prior_odds) elementwise on those matrices.
# ============================================================================
