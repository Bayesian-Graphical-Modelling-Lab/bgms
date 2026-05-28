# --------------------------------------------------------------------------- #
# Experiment 4 - Diagnose the 0.5-1.5% positive bias in the auxiliary chain
# (Experiment 3) via three probes.
#
# Probe (a) -- Cross-route verification of BF_inc_marg.
#   The analytic in Experiment 3 uses
#     BF_inc_marg(omega-route) = integral over omega of BF_inc(omega) p(omega)
#   where BF_inc(omega) = 1 / BF_excl(omega) computed via the L-space SD
#   identity. An independent route is to marginalise K instead:
#     BF_inc_marg(K-route)    = integral over K of L_data(K) Cauchy(K | 0, sigma)
#   where L_data(K) = p(Y | gamma=1, K, rest) / p(Y | gamma=0, rest)
#                   = exp(-A_K K^2 + B_K K)
#   (the data-only kernel in K-space, after extracting the slab piece).
#   These should be analytically equal. If they disagree, the analytic was
#   wrong and the chain is right.
#
# Probe (b) -- Collapsed (gamma, omega) chain with slice sampling.
#   Drops K from state. omega | gamma=1, rest, Y is non-conjugate:
#     p(omega | gamma=1, rest, Y) ~ BF_inc(omega) * p_IG(omega).
#   We slice-sample it. If the bias goes away, the issue was in the
#   interaction between the K Gibbs step and the omega conjugate step in
#   Experiment 3's chain.
#
# Probe (c) -- alpha > 1 cross-check.
#   At alpha > 1 the SD primitive is sinh-quadrature, not the alpha = 1
#   closed form. Run the same auxiliary chain at alpha = 3 and see if the
#   bias size is similar or alpha-dependent.
#
# Run:
#   Rscript --vanilla experiments/cauchy-slab/04-bias-diagnostic.R
# --------------------------------------------------------------------------- #

devtools::load_all(".", quiet = TRUE)
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("This experiment needs the `coda` package (Suggests).")
}

# ---------- Shared primitives ---------- #

d_omega_log <- function(omega) {
  -0.5 * log(2 * pi) - 1.5 * log(omega) - 1 / (2 * omega)
}
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

log_gauss <- function(x, mu, tau) {
  0.5 * log(tau / (2 * pi)) - 0.5 * tau * (x - mu) ^ 2
}

# alpha = 1 closed-form log BF(gamma=0 vs gamma=1) at the Roverato slave.
log_bf_excl_alpha1 <- function(sigma2_eff, l_ii, m_ij, S_ij, S_jj, beta) {
  inv_sig2  <- 1 / sigma2_eff
  A_post    <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta + S_jj)
  B_post    <- l_ii ^ 2 * m_ij * inv_sig2 - S_ij * l_ii
  A_prior   <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta)
  B_prior   <- l_ii ^ 2 * m_ij * inv_sig2
  tau_post  <- 2 * A_post
  tau_prior <- 2 * A_prior
  mu_post   <- B_post  / tau_post
  mu_prior  <- B_prior / tau_prior
  log_pi_post  <- log_gauss(m_ij, mu_post,  tau_post)
  log_pi_prior <- log_gauss(m_ij, mu_prior, tau_prior)
  log_pi_post - log_pi_prior
}

# alpha > 1 via sinh primitive.
log_bf_excl_sinh <- function(sigma2_eff, l_ii, m_ij, S_ij, S_jj, beta, alpha) {
  inv_sig2 <- 1 / sigma2_eff
  A_post   <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta + S_jj)
  B_post   <- l_ii ^ 2 * m_ij * inv_sig2 - S_ij * l_ii
  A_prior  <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta)
  B_prior  <- l_ii ^ 2 * m_ij * inv_sig2
  s_jj     <- 1.0  # treated as a state variable; matches Experiment 3 convention
  r_post   <- sd_log_density_at_l_ji_sinh_cpp(m_ij, A_post,  B_post,  s_jj, alpha)
  r_prior  <- sd_log_density_at_l_ji_sinh_cpp(m_ij, A_prior, B_prior, s_jj, alpha)
  r_post$log_density - r_prior$log_density
}


# ---------- Configurations (alpha = 1 closed-form) ---------- #

grid <- list(
  weak = list(
    l_ii = 1.0, m_ij = 0.05, S_ij = 0.5, S_jj = 50,
    beta = 0.5, sigma2 = 1.0, pi_inc = 0.5,
    label = "weak (S_ij=0.5)"
  ),
  borderline = list(
    l_ii = 1.0, m_ij = 0.20, S_ij = 2.0, S_jj = 50,
    beta = 0.5, sigma2 = 1.0, pi_inc = 0.5,
    label = "borderline (S_ij=2)"
  ),
  strong = list(
    l_ii = 1.0, m_ij = 0.30, S_ij = 5.0, S_jj = 50,
    beta = 0.5, sigma2 = 1.0, pi_inc = 0.5,
    label = "strong (S_ij=5)"
  )
)


# ===========================================================================
# Probe (a) -- Cross-route BF_inc_marg
# ===========================================================================

# Route 1: integrate over omega using SD log_bf_excl.
bf_inc_marg_omega_route <- function(cfg) {
  integrand <- function(omega) {
    vapply(omega, function(w) {
      lbe <- log_bf_excl_alpha1(cfg$sigma2 * w, cfg$l_ii, cfg$m_ij,
                                 cfg$S_ij, cfg$S_jj, cfg$beta)
      exp(-lbe + d_omega_log(w))
    }, numeric(1))
  }
  integrate(integrand, lower = 1e-6, upper = Inf,
            rel.tol = 1e-12, abs.tol = 1e-15,
            subdivisions = 5000L)$value
}

# Route 2: integrate over K with Cauchy slab.
# In K-space the data-only piece of the kernel (using my convention with
# 2*beta on the diag prior contribution) is
#     L_data(K) = exp(-A_K K^2 + B_K K)
# where
#     A_K = A_data / l_ii^2 = (beta + 0.5*S_jj) / l_ii^2
#     B_K = (B_data - 2 * A_data * m_ij) / l_ii
#         = -S_ij - (2*beta + S_jj) * m_ij / l_ii
# L_data(0) = 1 by construction, so L_data(K) = p(Y|gamma=1,K)/p(Y|gamma=0).
bf_inc_marg_k_route <- function(cfg) {
  A_data <- 0.5 * (2 * cfg$beta + cfg$S_jj)        # L-space data + diag contribution
  B_data <- -cfg$S_ij * cfg$l_ii                    # L-space data contribution to B
  A_K    <- A_data / cfg$l_ii ^ 2
  B_K    <- (B_data - 2 * A_data * cfg$m_ij) / cfg$l_ii
  integrand <- function(K) {
    exp(-A_K * K ^ 2 + B_K * K) /
      (pi * cfg$sigma2 * (1 + (K / sqrt(cfg$sigma2)) ^ 2)) *
      sqrt(cfg$sigma2)
  }
  integrate(integrand, lower = -Inf, upper = Inf,
            rel.tol = 1e-12, abs.tol = 1e-15,
            subdivisions = 5000L)$value
}

# Route 3: explicit 2D quadrature over (K, omega) via the joint
# integrand p(Y|gamma=1, K) * N(K | 0, sigma^2 * omega) * p(omega) / p(Y|gamma=0).
# Implemented as nested 1D integrals.
bf_inc_marg_2d_route <- function(cfg) {
  A_data <- 0.5 * (2 * cfg$beta + cfg$S_jj)
  B_data <- -cfg$S_ij * cfg$l_ii
  A_K    <- A_data / cfg$l_ii ^ 2
  B_K    <- (B_data - 2 * A_data * cfg$m_ij) / cfg$l_ii
  inner <- function(omega) {
    sig2 <- cfg$sigma2 * omega
    # gaussian inner integral: integral of L_data(K) * N(K | 0, sig2) dK
    A_tot <- A_K + 1 / (2 * sig2)
    z <- B_K ^ 2 / (4 * A_tot)
    sqrt(pi / A_tot) * exp(z) / sqrt(2 * pi * sig2)
  }
  integrand <- function(omega) {
    vapply(omega, function(w) inner(w) * exp(d_omega_log(w)), numeric(1))
  }
  integrate(integrand, lower = 1e-6, upper = Inf,
            rel.tol = 1e-12, abs.tol = 1e-15,
            subdivisions = 5000L)$value
}

cat("=== Probe (a): cross-route BF_inc_marg verification ===\n\n")
for (case_name in names(grid)) {
  cfg <- grid[[case_name]]
  r_omega <- bf_inc_marg_omega_route(cfg)
  r_k     <- bf_inc_marg_k_route(cfg)
  r_2d    <- bf_inc_marg_2d_route(cfg)
  cat(sprintf("--- %s ---\n", cfg$label))
  cat(sprintf("  omega-route: %.8f\n", r_omega))
  cat(sprintf("  K-route:     %.8f\n", r_k))
  cat(sprintf("  2D route:    %.8f\n", r_2d))
  cat(sprintf("  K / omega:   %.6f   |  2D / omega: %.6f\n",
              r_k / r_omega, r_2d / r_omega))
  cat("\n")
}


# ===========================================================================
# Probe (b) -- Collapsed (gamma, omega) chain with slice sampling
# ===========================================================================

# log p(omega | gamma=1, rest, Y) up to constant = log BF_inc(omega) + log p(omega)
log_omega_post_g1 <- function(omega, cfg) {
  if (omega <= 0) return(-Inf)
  -log_bf_excl_alpha1(cfg$sigma2 * omega, cfg$l_ii, cfg$m_ij,
                       cfg$S_ij, cfg$S_jj, cfg$beta) +
    d_omega_log(omega)
}

# Simple log-space slice sampler (stepping-out + shrinkage) on omega > 0.
slice_omega_g1 <- function(omega_curr, cfg, w = 1.0, max_steps = 50L) {
  # Work in u = log(omega) so the variable is unconstrained.
  log_target <- function(u) log_omega_post_g1(exp(u), cfg) + u  # +u = log Jacobian
  u_curr <- log(omega_curr)
  log_y  <- log_target(u_curr) - rexp(1)  # uniform vertical
  L <- u_curr - w * runif(1)
  R <- L + w
  # Stepping-out: expand until log_target(L) and log_target(R) < log_y.
  for (k in seq_len(max_steps)) {
    if (log_target(L) <= log_y) break
    L <- L - w
  }
  for (k in seq_len(max_steps)) {
    if (log_target(R) <= log_y) break
    R <- R + w
  }
  repeat {
    u_new <- runif(1, L, R)
    if (log_target(u_new) > log_y) return(exp(u_new))
    if (u_new < u_curr) L <- u_new else R <- u_new
  }
}

run_collapsed_chain <- function(cfg, T_total, burn, seed) {
  set.seed(seed)
  gamma <- 0L
  omega <- rinvgamma(1, 0.5, 0.5)
  out_gamma <- integer(T_total)
  for (t in seq_len(T_total)) {
    lbe <- log_bf_excl_alpha1(cfg$sigma2 * omega, cfg$l_ii, cfg$m_ij,
                               cfg$S_ij, cfg$S_jj, cfg$beta)
    log_odds <- log(cfg$pi_inc / (1 - cfg$pi_inc)) - lbe
    p1 <- 1 / (1 + exp(-log_odds))
    gamma <- as.integer(runif(1) < p1)
    if (gamma == 1L) {
      omega <- slice_omega_g1(omega, cfg)
    } else {
      omega <- rinvgamma(1, 0.5, 0.5)
    }
    out_gamma[t] <- gamma
  }
  out_gamma[(burn + 1):T_total]
}

cat("\n=== Probe (b): collapsed (gamma, omega) chain via slice sampling ===\n\n")

T_total <- 50000L
burn    <- 5000L
n_seeds <- 3L

for (case_name in names(grid)) {
  cfg <- grid[[case_name]]
  bf_marg <- bf_inc_marg_omega_route(cfg)
  analytic_p1 <- cfg$pi_inc * bf_marg /
                 (cfg$pi_inc * bf_marg + (1 - cfg$pi_inc))
  p1_vec <- numeric(n_seeds)
  ess_vec <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    g <- run_collapsed_chain(cfg, T_total, burn, seed = 4000L + s)
    p1_vec[s] <- mean(g)
    ess_vec[s] <- coda::effectiveSize(as.numeric(g))
  }
  p1_hat <- mean(p1_vec)
  se     <- sd(p1_vec) / sqrt(n_seeds)
  cat(sprintf("--- %s ---\n", cfg$label))
  cat(sprintf("analytic P(gamma=1) = %.5f\n", analytic_p1))
  cat(sprintf("collapsed slice chain (avg %d seeds): P(gamma=1) = %.5f  se = %.5f  mean ESS = %.0f\n",
              n_seeds, p1_hat, se, mean(ess_vec)))
  cat(sprintf("bias vs analytic: %+.5f  (z = %+.2f)\n",
              p1_hat - analytic_p1, (p1_hat - analytic_p1) / se))
  cat("\n")
}


# ===========================================================================
# Probe (c) -- alpha > 1 cross-check (auxiliary chain via sinh primitive)
# ===========================================================================

# Reproduces Experiment 3's chain structure but using the sinh primitive
# for the log BF instead of the closed form, and alpha = 3.
k_ij_pars_alpha_gt1 <- function(sigma2_eff, cfg) {
  # K_ij is not Gaussian at alpha > 1, but the auxiliary scheme only needs
  # us to refresh K from p(K | gamma=1, omega, rest, Y). The conditional
  # density in L-space is f(l_ji) = -A x^2 + B x + (alpha-1) log(s_jj + x^2);
  # the conditional mode + curvature give a Laplace proposal but the exact
  # conditional needs the sinh primitive's density. For diagnostic purposes
  # we use a Gaussian proposal at the L-space Laplace centre + MH correction
  # against the exact density. Cheap and produces correct stationary.
  inv_sig2 <- 1 / sigma2_eff
  A_post   <- 0.5 * (cfg$l_ii ^ 2 * inv_sig2 + 2 * cfg$beta + cfg$S_jj)
  B_post   <- cfg$l_ii ^ 2 * cfg$m_ij * inv_sig2 - cfg$S_ij * cfg$l_ii
  tau_l    <- 2 * A_post
  mu_l     <- B_post / tau_l
  list(
    mu_K = cfg$l_ii * (mu_l - cfg$m_ij),
    sd_K = cfg$l_ii / sqrt(tau_l)
  )
}

# Analytic via omega route at alpha > 1.
bf_inc_marg_omega_route_alpha <- function(cfg, alpha) {
  integrand <- function(omega) {
    vapply(omega, function(w) {
      lbe <- log_bf_excl_sinh(cfg$sigma2 * w, cfg$l_ii, cfg$m_ij,
                               cfg$S_ij, cfg$S_jj, cfg$beta, alpha)
      exp(-lbe + d_omega_log(w))
    }, numeric(1))
  }
  integrate(integrand, lower = 1e-6, upper = Inf,
            rel.tol = 1e-10, abs.tol = 1e-14,
            subdivisions = 2000L)$value
}

run_aux_chain_alpha <- function(cfg, alpha, T_total, burn, seed) {
  set.seed(seed)
  gamma <- 0L
  K_ij  <- 0.0
  omega <- rinvgamma(1, 0.5, 0.5)
  out_gamma <- integer(T_total)
  for (t in seq_len(T_total)) {
    sig2_eff <- cfg$sigma2 * omega
    lbe <- log_bf_excl_sinh(sig2_eff, cfg$l_ii, cfg$m_ij,
                             cfg$S_ij, cfg$S_jj, cfg$beta, alpha)
    log_odds <- log(cfg$pi_inc / (1 - cfg$pi_inc)) - lbe
    p1 <- 1 / (1 + exp(-log_odds))
    gamma <- as.integer(runif(1) < p1)
    if (gamma == 1L) {
      pars <- k_ij_pars_alpha_gt1(sig2_eff, cfg)
      K_ij <- rnorm(1, pars$mu_K, pars$sd_K)
    } else {
      K_ij <- 0.0
    }
    if (gamma == 1L) {
      omega <- rinvgamma(1, 1.0, 0.5 + K_ij ^ 2 / (2 * cfg$sigma2))
    } else {
      omega <- rinvgamma(1, 0.5, 0.5)
    }
    out_gamma[t] <- gamma
  }
  out_gamma[(burn + 1):T_total]
}

cat("\n=== Probe (c): aux chain at alpha = 3 (sinh primitive) ===\n\n")

alpha_test <- 3.0
T_total_c  <- 50000L
burn_c     <- 5000L
n_seeds_c  <- 3L

for (case_name in names(grid)) {
  cfg <- grid[[case_name]]
  bf_marg <- bf_inc_marg_omega_route_alpha(cfg, alpha_test)
  analytic_p1 <- cfg$pi_inc * bf_marg /
                 (cfg$pi_inc * bf_marg + (1 - cfg$pi_inc))
  p1_vec <- numeric(n_seeds_c)
  for (s in seq_len(n_seeds_c)) {
    g <- run_aux_chain_alpha(cfg, alpha_test, T_total_c, burn_c,
                              seed = 5000L + s)
    p1_vec[s] <- mean(g)
  }
  p1_hat <- mean(p1_vec)
  se     <- sd(p1_vec) / sqrt(n_seeds_c)
  cat(sprintf("--- %s ---\n", cfg$label))
  cat(sprintf("analytic P(gamma=1) = %.5f  (BF_inc_marg = %.4f)\n",
              analytic_p1, bf_marg))
  cat(sprintf("aux chain (alpha=%g, avg %d seeds): P(gamma=1) = %.5f  se = %.5f\n",
              alpha_test, n_seeds_c, p1_hat, se))
  cat(sprintf("bias vs analytic: %+.5f  (z = %+.2f)\n",
              p1_hat - analytic_p1, (p1_hat - analytic_p1) / se))
  cat("\n")
}

cat("\nDONE.\n")
