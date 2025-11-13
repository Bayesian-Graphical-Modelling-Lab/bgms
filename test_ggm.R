library(bgms)

# Dimension and true precision
p <- 10

adj <- matrix(0, nrow = p, ncol = p)
adj[lower.tri(adj)] <- rbinom(p * (p - 1) / 2, size = 1, prob = 0.3)
adj <- adj + t(adj)
# qgraph::qgraph(adj)
Omega <- BDgraph::rgwish(1, adj = adj, b = p + sample(0:p, 1), D = diag(p))
Sigma <- solve(Omega)
zapsmall(Omega)

# Data
n <- 1e3
x <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sigma)


# ---- Run MCMC with warmup and sampling ------------------------------------

# debugonce(mbgms:::bgm_gaussian)
sampling_results <- bgms:::sample_ggm(
  X = x,
  prior_inclusion_prob = matrix(.5, p, p),
  initial_edge_indicators = adj,
  no_iter = 500,
  no_warmup = 500,
  no_chains = 3,
  edge_selection = FALSE,
  no_threads = 1,
  seed = 123,
  progress_type = 1
)

true_values     <- zapsmall(Omega[upper.tri(Omega, TRUE)])
posterior_means <- rowMeans(sampling_results[[2]]$samples)
cbind(true_values, posterior_means)

plot(true_values, posterior_means)
abline(0, 1)

sampling_results2 <- bgms:::sample_ggm(
  X = x,
  prior_inclusion_prob = matrix(.5, p, p),
  initial_edge_indicators = adj,
  no_iter = 500,
  no_warmup = 500,
  no_chains = 3,
  edge_selection = TRUE,
  no_threads = 1,
  seed = 123,
  progress_type = 1
)

true_values     <- zapsmall(Omega[upper.tri(Omega, TRUE)])
posterior_means <- rowMeans(sampling_results2[[2]]$samples)

plot(true_values, posterior_means)
abline(0, 1)

plot(posterior_means, rowMeans(sampling_results2[[2]]$samples != 0))


mmm <- matrix(c(
  1.6735,   0,        0,        0,        0,
  0,   1.0000,        0,        0,  -3.4524,
  0,        0,   1.0000,        0,        0,
  0,        0,        0,   1.0000,        0,
  0,  -3.4524,        0,        0,   9.6674
), p, p)
mmm
chol(mmm)
base::isSymmetric(mmm)
eigen(mmm)

profvis::profvis({
  sampling_results <- bgm_gaussian(
    x = x,
    n = n,
    n_iter = 400,
    n_warmup = 400,
    n_phases = 10
  )
})

# Extract results
aveOmega <- sampling_results$aveOmega
aveGamma <- sampling_results$aveGamma
aOmega <- sampling_results$aOmega
aGamma <- sampling_results$aGamma
prob <- sampling_results$prob
proposal_sd <- sampling_results$proposal_sd

library(patchwork)
library(ggplot2)
df <- data.frame(
  true = aveOmega[lower.tri(aveOmega)],
  Omega[lower.tri(Omega)],
  estimated = aveOmega[lower.tri(aveOmega)],
  p_inclusion = aveGamma[lower.tri(aveGamma)]
)
p1 <- ggplot(df, aes(x = true, y = estimated)) +
  geom_point(size = 5, alpha = 0.8, shape = 21, fill = "grey") +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  labs(x = "True Values Omega", y = "Estimated Values Omega (Posterior Mean)")
p2 <- ggplot(df, aes(x = estimated, y = p_inclusion)) +
  geom_point(size = 5, alpha = 0.8, shape = 21, fill = "grey") +
  labs(
    x = "Estimated Values Omega (Posterior Mean)",
    y = "Estimated Inclusion Probabilities"
  )
(p1 + p2) + plot_layout(ncol = 1) & theme_bw(base_size = 20)

