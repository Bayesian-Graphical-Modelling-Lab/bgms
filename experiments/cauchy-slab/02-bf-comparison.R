# --------------------------------------------------------------------------- #
# Experiment 2 — Compare BF_Normal, BF_Cauchy_aux(omega), and BF_Cauchy_marginal
# on a small grid of (l_ii, m_ij, S_ij_data, S_jj_data) configurations.
#
# Goals:
#  (a) Verify E_{omega ~ IG(1/2, 1/2)} [BF_aux(sigma^2 * omega)] = BF_marginal
#      (both numerically via integrate() and by Monte Carlo over omega).
#  (b) Characterise the integrand BF_aux(sigma^2 * omega) * p_IG(omega) as a
#      function of omega -- where does the mass live? Tails OK? Integrable?
#  (c) Decide between auxiliary-variable (A) and marginalised (B) SD designs.
#      The integrand's shape and the quadrature cost are the deciding factors.
#
# Notation matches src/models/ggm/ggm_model.cpp::update_edge_indicator_*_pair:
#   A_post  = 0.5 * (l_ii^2 / sigma^2_eff + 2*beta + S_jj_data)
#   B_post  = l_ii^2 * m_ij / sigma^2_eff  -  S_ij_data * l_ii
#   A_prior = 0.5 * (l_ii^2 / sigma^2_eff + 2*beta)
#   B_prior = l_ii^2 * m_ij / sigma^2_eff
# sigma^2_eff = sigma^2 (Normal slab) or sigma^2 * omega (Cauchy via scale mix).
#
# At alpha = 1 the kernel f(x) = -A x^2 + B x is exact Gaussian: tau = 2A,
# mu = B/(2A), log_density(m_ij) computable in closed form. We use this
# closed form as the alpha = 1 cross-check on the C++ sinh primitive.
#
# Run:
#   Rscript --vanilla experiments/cauchy-slab/02-bf-comparison.R
# --------------------------------------------------------------------------- #

devtools::load_all(".", quiet = TRUE)

# Inverse-Gamma(1/2, 1/2) log-density. p(omega) = (1/sqrt(2 pi)) * omega^(-3/2) * exp(-1/(2 omega)).
d_omega_log <- function(omega) {
  -0.5 * log(2 * pi) - 1.5 * log(omega) - 1 / (2 * omega)
}
d_omega <- function(omega) exp(d_omega_log(omega))

# Closed-form alpha = 1 Gaussian log density at x_eval.
log_gauss <- function(x, mu, tau) {
  0.5 * log(tau / (2 * pi)) - 0.5 * tau * (x - mu) ^ 2
}

# log BF (state 1 over state 0) at the Roverato slave m_ij, evaluated at a
# given sigma^2_eff (the per-call effective slab variance).
log_bf_at_sigma2 <- function(sigma2_eff,
                             l_ii, m_ij, S_ij, S_jj, alpha, beta, s_jj_state,
                             route = c("alpha1_closed", "sinh")) {
  route <- match.arg(route)
  inv_sig2 <- 1 / sigma2_eff
  A_post   <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta + S_jj)
  B_post   <- l_ii ^ 2 * m_ij * inv_sig2 - S_ij * l_ii
  A_prior  <- 0.5 * (l_ii ^ 2 * inv_sig2 + 2 * beta)
  B_prior  <- l_ii ^ 2 * m_ij * inv_sig2

  if (route == "alpha1_closed") {
    stopifnot(alpha == 1)
    tau_post  <- 2 * A_post
    tau_prior <- 2 * A_prior
    mu_post   <- B_post  / tau_post
    mu_prior  <- B_prior / tau_prior
    log_pi_post  <- log_gauss(m_ij, mu_post,  tau_post)
    log_pi_prior <- log_gauss(m_ij, mu_prior, tau_prior)
    return(log_pi_post - log_pi_prior)
  }

  # sinh-quadrature primitive (the production path).
  r_post  <- sd_log_density_at_l_ji_sinh_cpp(
    m_ij, A_post,  B_post,  s_jj_state, alpha
  )
  r_prior <- sd_log_density_at_l_ji_sinh_cpp(
    m_ij, A_prior, B_prior, s_jj_state, alpha
  )
  if (r_post$status != 0 || r_prior$status != 0) {
    return(NA_real_)
  }
  r_post$log_density - r_prior$log_density
}


# ---------- Grid of SD-step states to evaluate ---------- #

grid <- list(
  null_alpha1 = list(
    l_ii = 1.0, m_ij = 0.05, S_ij = 0.5,  S_jj = 50, alpha = 1.0,
    beta = 0.5, s_jj_state = 1.0,
    label = "alpha=1, weak data (S_ij=0.5)"
  ),
  signal_alpha1 = list(
    l_ii = 1.0, m_ij = 0.30, S_ij = 5.0,  S_jj = 50, alpha = 1.0,
    beta = 0.5, s_jj_state = 1.0,
    label = "alpha=1, strong data (S_ij=5)"
  ),
  null_alpha3 = list(
    l_ii = 1.0, m_ij = 0.05, S_ij = 0.5,  S_jj = 50, alpha = 3.0,
    beta = 0.5, s_jj_state = 1.0,
    label = "alpha=3, weak data (S_ij=0.5)"
  ),
  signal_alpha3 = list(
    l_ii = 1.0, m_ij = 0.30, S_ij = 5.0,  S_jj = 50, alpha = 3.0,
    beta = 0.5, s_jj_state = 1.0,
    label = "alpha=3, strong data (S_ij=5)"
  )
)

sigma2_nominal <- 1.0  # nominal Cauchy scale^2 (sigma = 1.0)


# ---------- Per-grid analysis ---------- #

cat("=== Experiment 2: BF_Normal vs BF_Cauchy_aux vs BF_Cauchy_marginal ===\n")
cat(sprintf("Cauchy slab nominal scale sigma = %.3f\n\n", sqrt(sigma2_nominal)))

for (case_name in names(grid)) {
  cfg <- grid[[case_name]]
  route <- if (cfg$alpha == 1) "alpha1_closed" else "sinh"

  log_bf_fn <- function(sigma2_eff) {
    log_bf_at_sigma2(
      sigma2_eff,
      l_ii = cfg$l_ii, m_ij = cfg$m_ij, S_ij = cfg$S_ij, S_jj = cfg$S_jj,
      alpha = cfg$alpha, beta = cfg$beta, s_jj_state = cfg$s_jj_state,
      route = route
    )
  }

  bf_normal <- exp(log_bf_fn(sigma2_nominal))

  # Sweep BF_aux(sigma^2 * omega) across omega.
  omega_grid <- exp(seq(log(1e-4), log(1e+4), length.out = 41))
  log_bf_aux <- vapply(omega_grid, function(w) log_bf_fn(sigma2_nominal * w),
                       numeric(1))

  # Quadrature: BF_marginal = integral of BF_aux(omega) * p_IG(omega) over omega.
  bf_marginal_integrand <- function(omega) {
    vapply(omega, function(w) exp(log_bf_fn(sigma2_nominal * w) +
                                  d_omega_log(w)),
           numeric(1))
  }
  q_res <- tryCatch(
    integrate(
      bf_marginal_integrand,
      lower = 1e-6, upper = Inf,
      rel.tol = 1e-8, abs.tol = 1e-12,
      subdivisions = 1000L
    ),
    error = function(e) list(value = NA_real_, abs.error = NA_real_,
                             message = conditionMessage(e))
  )
  bf_marginal_quad <- q_res$value

  # Monte Carlo: BF_marginal ~= mean(BF_aux(sigma^2 * omega_i)),
  #              omega_i ~ IG(1/2, 1/2).
  N_mc <- 2e4
  omega_mc <- 1 / rgamma(N_mc, shape = 0.5, rate = 0.5)
  log_bf_mc <- vapply(omega_mc, function(w) log_bf_fn(sigma2_nominal * w),
                      numeric(1))
  bf_marginal_mc <- mean(exp(log_bf_mc), na.rm = TRUE)
  bf_marginal_mc_se <- sd(exp(log_bf_mc), na.rm = TRUE) / sqrt(N_mc)

  # Where does the integrand mass live? Find the median (in omega) of
  # the integrand mass.
  integrand_vals <- bf_marginal_integrand(omega_grid)
  integrand_mass <- integrand_vals * c(diff(omega_grid),
                                        diff(omega_grid)[length(omega_grid) - 1])
  cdf <- cumsum(integrand_mass) / sum(integrand_mass)
  med_omega <- omega_grid[which.min(abs(cdf - 0.5))]
  q05_omega <- omega_grid[which.min(abs(cdf - 0.05))]
  q95_omega <- omega_grid[which.min(abs(cdf - 0.95))]

  cat("--- ", cfg$label, " ---\n", sep = "")
  cat(sprintf("BF_Normal(sigma^2 = %.2f):      %12.6f\n",
              sigma2_nominal, bf_normal))
  cat(sprintf("BF_Cauchy_marginal (integrate): %12.6f  (abs.err = %.2e)\n",
              bf_marginal_quad, q_res$abs.error))
  cat(sprintf("BF_Cauchy_marginal (MC N=%d):  %12.6f  (se = %.2e)\n",
              N_mc, bf_marginal_mc, bf_marginal_mc_se))
  cat(sprintf("  -> agreement quad/MC: ratio = %.4f\n",
              bf_marginal_quad / bf_marginal_mc))
  cat(sprintf("Integrand mass (omega): q05 = %.4f, median = %.4f, q95 = %.4f\n",
              q05_omega, med_omega, q95_omega))
  cat(sprintf("BF_aux range over omega in [1e-4, 1e4]: log-BF = [%.3f, %.3f]\n",
              min(log_bf_aux, na.rm = TRUE), max(log_bf_aux, na.rm = TRUE)))
  cat("\n")
}


# ---------- Diagnostic for the alpha=3 sinh primitive: status codes? ---------- #

cat("--- alpha=3 sinh primitive status sweep ---\n")
cfg <- grid$signal_alpha3
omega_test <- c(1e-4, 1e-2, 0.1, 1, 10, 100, 1e4)
for (w in omega_test) {
  sigma2_eff <- sigma2_nominal * w
  inv_sig2 <- 1 / sigma2_eff
  A_post  <- 0.5 * (cfg$l_ii ^ 2 * inv_sig2 + 2 * cfg$beta + cfg$S_jj)
  B_post  <- cfg$l_ii ^ 2 * cfg$m_ij * inv_sig2 - cfg$S_ij * cfg$l_ii
  r <- sd_log_density_at_l_ji_sinh_cpp(
    cfg$m_ij, A_post, B_post, cfg$s_jj_state, cfg$alpha
  )
  cat(sprintf("  omega = %.0e:  A_post = %.3e, B_post = %.3e, status = %d, log_Z = %.3e\n",
              w, A_post, B_post, r$status, r$log_Z))
}

cat("\nDONE.\n")
