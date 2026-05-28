# --------------------------------------------------------------------------- #
# Experiment 6 - Localise the bias to the diag-prior / omega-conjugate
# interaction.
#
# Hypothesis (from Experiments 3 + 5):
#   The aux chain's omega-conjugate update uses pure-slab conjugacy:
#     omega | K, gamma=1 ~ InvGamma(1, 1/2 + K^2 / (2*sigma^2))
#   This is correct when the prior on K_ij | gamma=1, omega is just the
#   slab N(0, sigma^2 omega). But under the L-space framework that bgms's
#   SD formula targets, the effective prior on K_ij also includes the
#   diag-prior contribution on K_jj (via Roverato slaving), introducing
#   an omega-dependent normalisation Z(omega) that my conjugate update
#   silently ignores.
#
# Test: re-run Experiment 3's aux chain with beta = 0 (no diag-prior
# contribution to A). If the 1% bias vanishes, the diag-prior / omega
# interaction is the source. The fix in production would be slice / rejection
# sampling for the full omega conditional, replacing the pure-slab conjugate.
#
# Run:
#   Rscript --vanilla experiments/cauchy-slab/06-bias-localization-probe.R
# --------------------------------------------------------------------------- #

devtools::load_all(".", quiet = TRUE)
if (!requireNamespace("coda", quietly = TRUE)) stop("needs `coda`.")

d_omega_log <- function(omega) {
  -0.5 * log(2 * pi) - 1.5 * log(omega) - 1 / (2 * omega)
}
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

log_gauss <- function(x, mu, tau) {
  0.5 * log(tau / (2 * pi)) - 0.5 * tau * (x - mu) ^ 2
}

# log BF(gamma=0 vs gamma=1) at the Roverato slave m_ij, alpha = 1.
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

# K_ij conditional Gaussian (in K-space) for gamma = 1.
k_ij_pars_K <- function(sigma2_eff, l_ii, m_ij, S_ij, S_jj, beta) {
  inv_sig2  <- 1 / sigma2_eff
  A_post    <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta + S_jj)
  B_post    <- l_ii ^ 2 * m_ij * inv_sig2 - S_ij * l_ii
  tau_l     <- 2 * A_post
  mu_l      <- B_post / tau_l
  list(mu_K = l_ii * (mu_l - m_ij), sd_K = l_ii / sqrt(tau_l))
}

# BF_inc_marg via 1D quadrature over omega prior.
bf_inc_marg <- function(cfg) {
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

# Aux chain with pure-slab omega conjugate.
run_aux <- function(cfg, T_total, burn, seed) {
  set.seed(seed)
  gamma <- 0L
  K_ij  <- 0.0
  omega <- rinvgamma(1, 0.5, 0.5)
  out <- integer(T_total)
  for (t in seq_len(T_total)) {
    sig2_eff <- cfg$sigma2 * omega
    lbe <- log_bf_excl_alpha1(sig2_eff, cfg$l_ii, cfg$m_ij,
                               cfg$S_ij, cfg$S_jj, cfg$beta)
    log_odds <- log(cfg$pi_inc / (1 - cfg$pi_inc)) - lbe
    p1 <- 1 / (1 + exp(-log_odds))
    gamma <- as.integer(runif(1) < p1)
    if (gamma == 1L) {
      pars <- k_ij_pars_K(sig2_eff, cfg$l_ii, cfg$m_ij, cfg$S_ij,
                          cfg$S_jj, cfg$beta)
      K_ij <- rnorm(1, pars$mu_K, pars$sd_K)
    } else {
      K_ij <- 0.0
    }
    if (gamma == 1L) {
      omega <- rinvgamma(1, 1.0, 0.5 + K_ij ^ 2 / (2 * cfg$sigma2))
    } else {
      omega <- rinvgamma(1, 0.5, 0.5)
    }
    out[t] <- gamma
  }
  out[(burn + 1):T_total]
}


# ---------- Grid: same 3 configs as Experiment 3, plus beta=0 variants ---- #

base_configs <- list(
  weak       = list(l_ii=1.0, m_ij=0.05, S_ij=0.5, S_jj=50, sigma2=1.0,
                    pi_inc=0.5, label="weak (S_ij=0.5)"),
  borderline = list(l_ii=1.0, m_ij=0.20, S_ij=2.0, S_jj=50, sigma2=1.0,
                    pi_inc=0.5, label="borderline (S_ij=2)"),
  strong     = list(l_ii=1.0, m_ij=0.30, S_ij=5.0, S_jj=50, sigma2=1.0,
                    pi_inc=0.5, label="strong (S_ij=5)")
)

# Two variants per config: beta = 0.5 (Exp 3 default) and beta = 0.
betas <- c(beta_05 = 0.5, beta_00 = 0.0)

T_total <- 100000L
burn    <- 10000L
n_seeds <- 5L

cat("=== Experiment 6: bias localisation via beta sweep ===\n")
cat(sprintf("T = %d post-burn x %d seeds\n\n", T_total - burn, n_seeds))

for (case in names(base_configs)) {
  cfg_base <- base_configs[[case]]
  cat(sprintf("--- %s ---\n", cfg_base$label))
  for (bname in names(betas)) {
    cfg <- c(cfg_base, list(beta = betas[[bname]]))
    bfm <- bf_inc_marg(cfg)
    analytic_p1 <- cfg$pi_inc * bfm /
                   (cfg$pi_inc * bfm + (1 - cfg$pi_inc))
    p1_vec <- numeric(n_seeds)
    for (s in seq_len(n_seeds)) {
      p1_vec[s] <- mean(run_aux(cfg, T_total, burn, seed = 7000L + s))
    }
    p1_hat <- mean(p1_vec)
    se     <- sd(p1_vec) / sqrt(n_seeds)
    cat(sprintf("  beta = %.2f:  analytic = %.5f,  aux = %.5f  se = %.5f  bias = %+.5f  z = %+.2f\n",
                cfg$beta, analytic_p1, p1_hat, se,
                p1_hat - analytic_p1,
                (p1_hat - analytic_p1) / se))
  }
  cat("\n")
}

cat("DONE.\n")
