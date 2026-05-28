# --------------------------------------------------------------------------- #
# Experiment 3 - Mixing of the (gamma, K_ij, omega) auxiliary-variable chain
# at fixed rest-conditional state.
#
# Setup. The full GGM SD between-step targets p(gamma, K | Y, rest). For a
# single edge with rest fixed, the joint posterior of the auxiliary triple
# (gamma_ij, K_ij, omega_ij) is
#
#   p(gamma=1, K_ij, omega | rest, Y) prop  p(Y | gamma=1, K_ij, rest)
#                                         * N(K_ij | 0, sigma^2 omega)
#                                         * InvGamma(omega | 1/2, 1/2)
#                                         * pi_inc
#   p(gamma=0, K_ij=0, omega | rest, Y) prop p(Y | gamma=0, rest)
#                                          * InvGamma(omega | 1/2, 1/2)
#                                          * (1 - pi_inc)
#
# When K_ij is integrated out, the marginal posterior on gamma is
#   P(gamma=1) = pi_inc * BF_Cauchy / (pi_inc * BF_Cauchy + (1 - pi_inc))
# where BF_Cauchy = E_omega[BF_aux(sigma^2 omega)] (= the Experiment 2
# marginal quadrature result).
#
# This experiment runs the 3-step Gibbs scan:
#   (1) gamma | omega, rest, Y -- exact Gibbs using BF_aux at current omega
#   (2) K_ij  | gamma, omega, rest, Y -- Gaussian conditional when gamma=1
#   (3) omega | gamma, K_ij    -- InvGamma conjugate
# and compares its stationary frequency P(gamma=1) against the analytic
# marginal, also reporting effective sample size for the gamma chain.
#
# Reference chain: marginalised Gibbs at the same fixed rest, sampling
# gamma directly from the marginal P(gamma=1). This is the best possible
# mixing for a Gibbs on gamma alone -- the upper bound the auxiliary chain
# is competing against.
#
# CURRENT FINDING (open thread):
#   The auxiliary (gamma, K, omega) chain shows a small but reproducible
#   positive bias in P(gamma=1) of order 0.5%-1.5% across all three
#   configurations (z = 10 to 25 across 5 seeds * 200k iterations). The
#   marginal Gibbs chain matches the analytic to MC error, confirming the
#   analytic value is correct. The bias is largest in the borderline regime
#   where gamma flips most often, suggesting a subtle issue at the MoMS
#   boundary. Each conditional update has been re-derived (collapsed Gibbs
#   for gamma, conjugate IG for omega, K-space Gaussian for K_ij with the
#   l_ii * (mu_l - m_ij) transform) and looks correct. Worth deeper
#   investigation before committing to the C++ implementation; the
#   algorithm appears "essentially correct" but isn't bit-exact against
#   the analytic.
#
# Run:
#   Rscript --vanilla experiments/cauchy-slab/03-aux-chain-mixing.R
# --------------------------------------------------------------------------- #

devtools::load_all(".", quiet = TRUE)
if (!requireNamespace("coda", quietly = TRUE)) {
  stop("This experiment needs the `coda` package (Suggests).")
}

set.seed(2026L)

d_omega_log <- function(omega) {
  -0.5 * log(2 * pi) - 1.5 * log(omega) - 1 / (2 * omega)
}

# Closed-form alpha = 1 log-BF(gamma=0 vs gamma=1) -- the Savage-Dickey
# *exclusion* BF.  By SD, BF(M_0 vs M_1) = post density at K=0 / prior
# density at K=0; this is evaluated in L-space at the Roverato slave
# point m_ij (post density / prior density at m_ij). High BF here means
# the spike (gamma=0) is preferred.
log_bf_excl_alpha1 <- function(sigma2_eff, l_ii, m_ij, S_ij, S_jj, beta) {
  inv_sig2 <- 1 / sigma2_eff
  A_post   <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta + S_jj)
  B_post   <- l_ii ^ 2 * m_ij * inv_sig2 - S_ij * l_ii
  A_prior  <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta)
  B_prior  <- l_ii ^ 2 * m_ij * inv_sig2
  tau_post  <- 2 * A_post
  tau_prior <- 2 * A_prior
  mu_post   <- B_post  / tau_post
  mu_prior  <- B_prior / tau_prior
  log_pi_post <- 0.5 * log(tau_post  / (2 * pi)) -
                 0.5 * tau_post  * (m_ij - mu_post)  ^ 2
  log_pi_prior <- 0.5 * log(tau_prior / (2 * pi)) -
                  0.5 * tau_prior * (m_ij - mu_prior) ^ 2
  log_pi_post - log_pi_prior
}

# Conditional Gaussian for K_ij | gamma = 1, omega, rest, Y, returned in
# K-space.  The kernel f(x) = -A x^2 + B x lives in L-space (x = l_ji),
# so the closed-form conditional gives mu_l, tau_l in L-space.  K-space
# follows via K_ij = l_ii * (l_ji - m_ij) (the Roverato relation).
k_ij_conditional_params_K <- function(sigma2_eff, l_ii, m_ij, S_ij, S_jj,
                                       beta) {
  inv_sig2  <- 1 / sigma2_eff
  A_post    <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta + S_jj)
  B_post    <- l_ii ^ 2 * m_ij * inv_sig2 - S_ij * l_ii
  tau_l     <- 2 * A_post
  mu_l      <- B_post / tau_l
  list(
    mu_K = l_ii * (mu_l - m_ij),
    sd_K = l_ii / sqrt(tau_l)
  )
}

# Inverse-Gamma sampling: omega = 1 / Gamma(shape, rate).
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

# BF_Cauchy marginal (inclusion direction) via 1D quadrature.
#
# Marginal likelihood ratio under the Cauchy slab:
#   BF_inc_Cauchy = p(Y | gamma=1, rest) / p(Y | gamma=0, rest)
#                 = integral_omega [p(Y | gamma=1, omega) / p(Y | gamma=0)] p(omega) d omega
#                 = integral_omega BF_inc(omega) p(omega) d omega
#                 = E_omega[1 / BF_excl(omega)]
#
# Note: this is NOT 1/E_omega[BF_excl(omega)] in general (Jensen). The
# "average BF_excl" is the wrong quantity for the marginal slab.
bf_cauchy_inc_marginal <- function(sigma2_nominal, l_ii, m_ij, S_ij, S_jj,
                                    beta, log_bf_excl_fn) {
  integrand <- function(omega) {
    vapply(omega,
           function(w) exp(-log_bf_excl_fn(sigma2_nominal * w, l_ii, m_ij,
                                            S_ij, S_jj, beta) +
                            d_omega_log(w)),
           numeric(1))
  }
  integrate(integrand, lower = 1e-6, upper = Inf,
            rel.tol = 1e-10, abs.tol = 1e-14,
            subdivisions = 2000L)$value
}


# ---------- Configurations ---------- #
# alpha = 1 only here -- closed-form log_BF lets us isolate the auxiliary
# mixing question from any quadrature numerics in the sinh primitive.

grid <- list(
  weak = list(
    l_ii = 1.0, m_ij = 0.05, S_ij = 0.5, S_jj = 50,
    beta = 0.5, sigma2 = 1.0, pi_inc = 0.5,
    label = "weak data (S_ij=0.5, pi_inc=0.5)"
  ),
  strong = list(
    l_ii = 1.0, m_ij = 0.30, S_ij = 5.0, S_jj = 50,
    beta = 0.5, sigma2 = 1.0, pi_inc = 0.5,
    label = "strong data (S_ij=5,   pi_inc=0.5)"
  ),
  borderline = list(
    l_ii = 1.0, m_ij = 0.20, S_ij = 2.0, S_jj = 50,
    beta = 0.5, sigma2 = 1.0, pi_inc = 0.5,
    label = "borderline (S_ij=2,    pi_inc=0.5)"
  )
)


# ---------- Auxiliary chain ---------- #

run_aux_chain <- function(cfg, T_total, burn) {
  # State.
  gamma <- 0L
  K_ij  <- 0.0
  omega <- rinvgamma(1, shape = 0.5, rate = 0.5)

  gamma_chain <- integer(T_total)
  K_chain     <- numeric(T_total)
  omega_chain <- numeric(T_total)
  log_bf_chain <- numeric(T_total)

  for (t in seq_len(T_total)) {
    sigma2_eff <- cfg$sigma2 * omega
    log_bf_excl <- log_bf_excl_alpha1(sigma2_eff, cfg$l_ii, cfg$m_ij,
                                       cfg$S_ij, cfg$S_jj, cfg$beta)

    # (1) gamma Gibbs.  log_odds(gamma=1) = log(pi/(1-pi)) - log_bf_excl,
    # since BF_excl > 1 favours the spike.
    log_odds <- log(cfg$pi_inc / (1 - cfg$pi_inc)) - log_bf_excl
    p1 <- 1 / (1 + exp(-log_odds))
    gamma <- as.integer(runif(1) < p1)

    # (2) K_ij Gibbs in K-space.
    if (gamma == 1L) {
      pars <- k_ij_conditional_params_K(sigma2_eff, cfg$l_ii, cfg$m_ij,
                                         cfg$S_ij, cfg$S_jj, cfg$beta)
      K_ij <- rnorm(1, pars$mu_K, pars$sd_K)
    } else {
      K_ij <- 0.0
    }

    # (3) omega conjugate update.
    # Prior:           InvGamma(1/2, 1/2)
    # Likelihood on K: gamma==1 -> N(0, sigma^2 omega); gamma==0 -> none
    # Posterior:       InvGamma(1/2 + gamma/2,  1/2 + gamma*K_ij^2/(2 sigma^2))
    if (gamma == 1L) {
      omega <- rinvgamma(1, shape = 1.0,
                         rate = 0.5 + K_ij ^ 2 / (2 * cfg$sigma2))
    } else {
      omega <- rinvgamma(1, shape = 0.5, rate = 0.5)
    }

    gamma_chain[t]    <- gamma
    K_chain[t]        <- K_ij
    omega_chain[t]    <- omega
    log_bf_chain[t]   <- log_bf_excl
  }

  keep <- (burn + 1):T_total
  list(
    gamma = gamma_chain[keep],
    K     = K_chain[keep],
    omega = omega_chain[keep],
    log_bf = log_bf_chain[keep]
  )
}


# ---------- Marginalised reference chain ---------- #
# Sample gamma directly from the marginal P(gamma=1). This is the upper
# bound for gamma mixing -- iid samples, no autocorrelation.

run_marginal_chain <- function(cfg, T_total, burn, bf_inc_marg) {
  log_odds <- log(cfg$pi_inc / (1 - cfg$pi_inc)) + log(bf_inc_marg)
  p1 <- 1 / (1 + exp(-log_odds))
  gamma <- as.integer(runif(T_total) < p1)
  list(gamma = gamma[(burn + 1):T_total], p1 = p1)
}


# ---------- Run + report ---------- #

T_total <- 200000L
burn    <- 20000L
n_seeds <- 5L

cat("=== Experiment 3: auxiliary (gamma, K_ij, omega) chain at fixed rest ===\n")
cat(sprintf("T = %d (post-burn %d) across %d seeds\n\n",
            T_total, T_total - burn, n_seeds))

for (case in names(grid)) {
  cfg <- grid[[case]]

  bf_inc_marg <- bf_cauchy_inc_marginal(cfg$sigma2, cfg$l_ii, cfg$m_ij,
                                         cfg$S_ij, cfg$S_jj, cfg$beta,
                                         log_bf_excl_alpha1)
  analytic_p1 <- cfg$pi_inc * bf_inc_marg /
                 (cfg$pi_inc * bf_inc_marg + (1 - cfg$pi_inc))

  aux_p1_vec <- numeric(n_seeds)
  mar_p1_vec <- numeric(n_seeds)
  ess_aux_vec <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    set.seed(2026L + s)
    aux <- run_aux_chain(cfg, T_total, burn)
    mar <- run_marginal_chain(cfg, T_total, burn, bf_inc_marg)
    aux_p1_vec[s]  <- mean(aux$gamma)
    mar_p1_vec[s]  <- mean(mar$gamma)
    ess_aux_vec[s] <- coda::effectiveSize(as.numeric(aux$gamma))
  }
  aux_p1   <- mean(aux_p1_vec)
  aux_se   <- sd(aux_p1_vec) / sqrt(n_seeds)
  mar_p1   <- mean(mar_p1_vec)
  mar_se   <- sd(mar_p1_vec) / sqrt(n_seeds)
  ess_mean <- mean(ess_aux_vec)

  cat("--- ", cfg$label, " ---\n", sep = "")
  cat(sprintf("BF_Cauchy_inc(marg) = %.6f,  analytic P(gamma=1) = %.6f\n",
              bf_inc_marg, analytic_p1))
  cat(sprintf("Auxiliary chain (avg over %d seeds):  P(gamma=1) = %.5f (se across seeds = %.5f),  mean ESS = %.0f\n",
              n_seeds, aux_p1, aux_se, ess_mean))
  cat(sprintf("Marginal  chain (avg over %d seeds):  P(gamma=1) = %.5f (se across seeds = %.5f)\n",
              n_seeds, mar_p1, mar_se))
  cat(sprintf("Aux bias vs analytic: %+.5f  (z = %+.2f)\n",
              aux_p1 - analytic_p1, (aux_p1 - analytic_p1) / aux_se))
  cat(sprintf("Mar bias vs analytic: %+.5f  (z = %+.2f)\n",
              mar_p1 - analytic_p1, (mar_p1 - analytic_p1) / mar_se))
  cat(sprintf("Per-seed aux estimates: %s\n",
              paste(sprintf("%.5f", aux_p1_vec), collapse = ", ")))
  cat("\n")
}

cat("\nDONE.\n")
