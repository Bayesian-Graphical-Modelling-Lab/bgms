# --------------------------------------------------------------------------- #
# Experiment 1 — Scale-mixture-of-normals identity for the Cauchy slab.
#
# Goal: numerically verify that
#     omega ~ InvGamma(1/2, 1/2)
#     X | omega ~ N(0, sigma^2 * omega)
# induces  X ~ Cauchy(0, sigma).
#
# If this fails, the entire auxiliary-variable design (Option A) is wrong
# regardless of how we wire it into the SD between-step. This is the
# foundational sanity check.
#
# Run:
#   Rscript --vanilla experiments/cauchy-slab/01-identity-check.R
# --------------------------------------------------------------------------- #

set.seed(2026L)

sigma <- 1.7   # arbitrary non-unit scale to flush sigma-handling bugs
N     <- 5e5

# omega ~ InvGamma(1/2, 1/2) via 1 / Gamma(1/2, 1/2)
omega <- 1 / rgamma(N, shape = 0.5, rate = 0.5)

# X | omega ~ N(0, sigma^2 * omega)
x_mixture <- rnorm(N, mean = 0, sd = sigma * sqrt(omega))

# Reference: direct Cauchy samples at the same scale.
x_direct <- rcauchy(N, location = 0, scale = sigma)

# Truncated KS to avoid Cauchy-tail dominance of one or two outliers; the
# practical agreement we care about is in the bulk where the SD primitives
# actually evaluate the kernel.
trim <- 200 * sigma
ks <- suppressWarnings(ks.test(
  x_mixture[abs(x_mixture) < trim],
  x_direct  [abs(x_direct)  < trim]
))

# Quantile agreement at coverage-relevant probabilities.
probs <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
q_mix <- quantile(x_mixture, probs)
q_ref <- qcauchy(probs, location = 0, scale = sigma)
q_dir <- quantile(x_direct,  probs)

cat("=== Experiment 1: scale-mixture identity ===\n")
cat(sprintf("sigma = %.3f,  N = %d\n", sigma, N))
cat(sprintf("trim threshold for KS = %.1f * sigma = %.1f\n", trim / sigma, trim))
cat("\nKS (mixture vs. direct Cauchy, both trimmed):\n")
cat(sprintf("  D = %.4f,  p = %.4f  (trimmed N_mix=%d, N_dir=%d)\n",
            ks$statistic, ks$p.value,
            sum(abs(x_mixture) < trim), sum(abs(x_direct) < trim)))

cat("\nQuantiles (mixture vs. qcauchy vs. direct sample):\n")
tab <- data.frame(
  prob      = probs,
  mixture   = q_mix,
  qcauchy   = q_ref,
  direct    = q_dir,
  err_vs_q  = q_mix - q_ref
)
print(tab, row.names = FALSE, digits = 4)

# Tail probability check at a moderate point. For X ~ Cauchy(0, sigma),
# P(|X| > t) = 1 - 2/pi * atan(t/sigma) approximately = 2/pi * sigma / t
# for t >> sigma.
for (t in c(3, 10, 30) * sigma) {
  p_mix    <- mean(abs(x_mixture) > t)
  p_direct <- mean(abs(x_direct)  > t)
  p_exact  <- 2 * pcauchy(-t, scale = sigma)
  cat(sprintf("P(|X| > %.1f):  mixture = %.5f,  direct = %.5f,  exact = %.5f\n",
              t, p_mix, p_direct, p_exact))
}

# Diagnostic: omega distribution — verify it really is InvGamma(1/2, 1/2).
cat("\nomega ~ InvGamma(1/2, 1/2) diagnostics (median, IQR, P(omega > 10)):\n")
cat(sprintf("  median(omega) = %.4f  (exact = %.4f)\n",
            median(omega), 1 / qgamma(0.5, 0.5, 0.5)))
cat(sprintf("  IQR(omega)    = %.4f\n", IQR(omega)))
cat(sprintf("  P(omega > 10) = %.4f  (exact = %.4f)\n",
            mean(omega > 10),
            pgamma(1/10, 0.5, 0.5)))  # P(1/omega < 1/10) = P(v < 1/10)

cat("\nDONE.\n")
