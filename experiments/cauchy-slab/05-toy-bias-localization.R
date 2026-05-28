# --------------------------------------------------------------------------- #
# Experiment 5 - Localise the 1% aux-chain bias on a minimal toy.
#
# Strip everything specific to bgms (L-space SD, Roverato slaving, sinh
# quadrature, GGM Cholesky structure) and run the same (gamma, K, omega)
# auxiliary chain on the smallest possible spike-and-Cauchy-slab problem:
#
#     y       ~ N(K, 1)              (single observation)
#     gamma   ~ Bern(pi_inc)
#     K | gamma=0 = 0                 (spike)
#     K | gamma=1, omega ~ N(0, omega) (slab; sigma^2 = 1)
#     omega    ~ InvGamma(1/2, 1/2)
#
# Marginally K | gamma=1 ~ Cauchy(0, 1), so this is exactly the same
# scale-mixture-of-normals construction as the GGM extension would use.
#
# Closed forms:
#   posterior K | gamma=1, omega, y ~ N( y*omega/(omega+1), omega/(omega+1) )
#   marginal likelihood: p(y | gamma=1, omega) = N(y | 0, omega+1)
#   BF_inc(omega) = (omega+1)^(-1/2) * exp(y^2 * omega / (2*(omega+1)))
#   BF_inc_marg   = integral BF_inc(omega) * IG(omega | 1/2, 1/2) d omega
#   analytic P(gamma=1) = pi_inc * BF_inc_marg / (pi_inc BF_inc_marg + (1-pi_inc))
#
# If the aux chain (conjugate omega update from K) shows the same ~1% bias
# here as in Experiment 3, the bias is intrinsic to the (gamma, K, omega)
# Gibbs-scan structure -- bgms-irrelevant, and the C++ design has to use
# slice/rejection sampling for omega instead. If the aux chain is bias-free
# in this toy, the bias in Experiment 3 lives in the GGM-specific code
# (L-space densities, Roverato transform), not in the auxiliary scheme.
#
# Run:
#   Rscript --vanilla experiments/cauchy-slab/05-toy-bias-localization.R
# --------------------------------------------------------------------------- #

if (!requireNamespace("coda", quietly = TRUE)) {
  stop("This experiment needs the `coda` package (Suggests).")
}

d_omega_log <- function(omega) {
  -0.5 * log(2 * pi) - 1.5 * log(omega) - 1 / (2 * omega)
}
rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

# BF_inc(omega) = p(y | gamma=1, omega) / p(y | gamma=0)
#               = N(y | 0, omega+1) / N(y | 0, 1)
log_bf_inc <- function(omega, y) {
  -0.5 * log(omega + 1) + y ^ 2 * omega / (2 * (omega + 1))
}

# Marginal BF integrated against the IG(1/2, 1/2) prior on omega.
bf_inc_marg <- function(y) {
  integrand <- function(w) {
    vapply(w, function(wi) exp(log_bf_inc(wi, y) + d_omega_log(wi)),
           numeric(1))
  }
  integrate(integrand, lower = 1e-6, upper = Inf,
            rel.tol = 1e-12, abs.tol = 1e-15,
            subdivisions = 5000L)$value
}

# ---------------- aux chain (suspected biased) ---------------- #
# Gibbs scan:
#   gamma | omega : collapsed Gibbs using BF_inc(omega)
#   K     | gamma, omega, y : Gaussian conjugate (or 0 if gamma=0)
#   omega | gamma, K : InvGamma conjugate (or prior if gamma=0)
run_aux <- function(y, pi_inc, T_total, burn, seed) {
  set.seed(seed)
  gamma <- 0L
  K     <- 0
  omega <- rinvgamma(1, 0.5, 0.5)
  out   <- integer(T_total)
  for (t in seq_len(T_total)) {
    lbf <- log_bf_inc(omega, y)
    log_odds <- log(pi_inc / (1 - pi_inc)) + lbf
    p1 <- 1 / (1 + exp(-log_odds))
    gamma <- as.integer(runif(1) < p1)
    if (gamma == 1L) {
      s2 <- omega / (omega + 1)
      mu <- y * s2
      K  <- rnorm(1, mu, sqrt(s2))
    } else {
      K <- 0
    }
    if (gamma == 1L) {
      omega <- rinvgamma(1, 1.0, 0.5 + K ^ 2 / 2)
    } else {
      omega <- rinvgamma(1, 0.5, 0.5)
    }
    out[t] <- gamma
  }
  out[(burn + 1):T_total]
}

# ---------------- slice chain (collapsed, bias-free reference) -------- #
# Gibbs scan:
#   gamma | omega : same as aux
#   omega | gamma, y : non-conjugate. For gamma=1 use log-space slice
#                       sampling on p(omega | gamma=1, y) prop BF_inc(omega) p(omega).
slice_omega_g1 <- function(omega_curr, y, w_step = 1.0, max_steps = 50L) {
  log_target <- function(u) {
    w <- exp(u)
    log_bf_inc(w, y) + d_omega_log(w) + u   # +u = log Jacobian for u = log omega
  }
  u_curr <- log(omega_curr)
  log_y  <- log_target(u_curr) - rexp(1)
  L <- u_curr - w_step * runif(1)
  R <- L + w_step
  for (k in seq_len(max_steps)) {
    if (log_target(L) <= log_y) break
    L <- L - w_step
  }
  for (k in seq_len(max_steps)) {
    if (log_target(R) <= log_y) break
    R <- R + w_step
  }
  repeat {
    u_new <- runif(1, L, R)
    if (log_target(u_new) > log_y) return(exp(u_new))
    if (u_new < u_curr) L <- u_new else R <- u_new
  }
}

run_slice <- function(y, pi_inc, T_total, burn, seed) {
  set.seed(seed)
  gamma <- 0L
  omega <- rinvgamma(1, 0.5, 0.5)
  out <- integer(T_total)
  for (t in seq_len(T_total)) {
    lbf <- log_bf_inc(omega, y)
    log_odds <- log(pi_inc / (1 - pi_inc)) + lbf
    p1 <- 1 / (1 + exp(-log_odds))
    gamma <- as.integer(runif(1) < p1)
    if (gamma == 1L) {
      omega <- slice_omega_g1(omega, y)
    } else {
      omega <- rinvgamma(1, 0.5, 0.5)
    }
    out[t] <- gamma
  }
  out[(burn + 1):T_total]
}

# ---------------- Independent MC verification of the analytic ---------- #
mc_p1 <- function(y, pi_inc, N) {
  omega_draws <- rinvgamma(N, 0.5, 0.5)
  # P(gamma=1) = pi_inc * E_omega[BF_inc(omega)] / (pi_inc E_omega[...] + (1-pi_inc))
  bf_inc_draws <- exp(log_bf_inc(omega_draws, y))
  bf_mc <- mean(bf_inc_draws)
  bf_mc_se <- sd(bf_inc_draws) / sqrt(N)
  list(
    bf_marg = bf_mc,
    bf_marg_se = bf_mc_se,
    p1 = pi_inc * bf_mc / (pi_inc * bf_mc + (1 - pi_inc))
  )
}

# ---------------- Run + report ---------------- #
y_values <- c(0.5, 1.5, 3.0)   # weak / borderline / strong evidence
pi_inc   <- 0.5
T_total  <- 100000L
burn     <- 10000L
n_seeds  <- 5L

cat("=== Experiment 5: toy spike-and-Cauchy-slab bias localization ===\n")
cat(sprintf("y ~ N(K, 1);  pi_inc = %.2f;  T = %d post-burn x %d seeds\n\n",
            pi_inc, T_total - burn, n_seeds))

for (y in y_values) {
  bf_marg_quad <- bf_inc_marg(y)
  mc <- mc_p1(y, pi_inc, N = 2e6)
  analytic_p1  <- pi_inc * bf_marg_quad /
                  (pi_inc * bf_marg_quad + (1 - pi_inc))

  aux_p1   <- numeric(n_seeds)
  slice_p1 <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    aux_p1[s]   <- mean(run_aux  (y, pi_inc, T_total, burn, seed = 5000L + s))
    slice_p1[s] <- mean(run_slice(y, pi_inc, T_total, burn, seed = 6000L + s))
  }

  aux_mean   <- mean(aux_p1);   aux_se   <- sd(aux_p1)   / sqrt(n_seeds)
  slice_mean <- mean(slice_p1); slice_se <- sd(slice_p1) / sqrt(n_seeds)

  cat(sprintf("--- y = %.2f ---\n", y))
  cat(sprintf("BF_inc_marg (quad): %.6f\n", bf_marg_quad))
  cat(sprintf("BF_inc_marg (MC N=2e6): %.6f +- %.6f  (ratio quad/MC = %.6f)\n",
              mc$bf_marg, mc$bf_marg_se, bf_marg_quad / mc$bf_marg))
  cat(sprintf("analytic P(gamma=1):     %.5f\n", analytic_p1))
  cat(sprintf("aux chain  (5 seeds):    %.5f  se = %.5f   bias = %+.5f  z = %+.2f\n",
              aux_mean, aux_se, aux_mean - analytic_p1,
              (aux_mean - analytic_p1) / aux_se))
  cat(sprintf("slice chain (5 seeds):    %.5f  se = %.5f   bias = %+.5f  z = %+.2f\n",
              slice_mean, slice_se, slice_mean - analytic_p1,
              (slice_mean - analytic_p1) / slice_se))
  cat("\n")
}

cat("DONE.\n")
