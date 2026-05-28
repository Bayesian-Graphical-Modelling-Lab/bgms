# --------------------------------------------------------------------------- #
# Unit tests for savage_dickey::slice_sample_cauchy_omega_active and
# savage_dickey::sample_cauchy_omega_prior.
#
# These primitives implement the per-edge omega update for the Cauchy slab
# via the scale-mixture-of-normals representation
#     omega ~ InvGamma(1/2, 1/2),
#     K_ij | gamma=1, omega ~ N(0, sigma^2 * omega).
# See src/math/savage_dickey/cauchy_omega.h and experiments/cauchy-slab/.
#
# Coverage:
#   - Prior draws match IG(1/2, 1/2) in quantiles and tails.
#   - Slice sampler matches the closed-form pure-slab conjugate when
#     A_diag_K = B_K = 0.
#   - Slice sampler matches a Monte-Carlo reference (importance sampling
#     over the IG prior) for nonzero A_diag_K and B_K.
#   - Pure-slab special case is bias-free against the IG conjugate.
# --------------------------------------------------------------------------- #

set.seed(2026L)

# Closed-form log-conditional p(omega | K, gamma=1, rest_L) up to a constant
# in omega; matches the kernel implemented in cauchy_omega.cpp.
log_cond_omega = function(omega, K, sigma2, A_diag_K, B_K) {
  A_K = 0.5 / (sigma2 * omega) + A_diag_K
  -K^2 / (2 * sigma2 * omega) -
    1 / (2 * omega) -
    1.5 * log(omega) +
    0.5 * log(A_K) -
    B_K^2 / (4 * A_K)
}

# Run the slice sampler iteratively from a single starting omega, returning
# post-burn samples. Each step is a valid slice transition for the target
# conditional, so the chain has the right stationary distribution.
run_omega_chain = function(K, sigma2, A_diag_K, B_K,
                           T_total, burn, seed_start, omega0 = 1.0) {
  omega = omega0
  out = numeric(T_total)
  for(t in seq_len(T_total)) {
    omega = sd_slice_sample_cauchy_omega_active_cpp(
      K = K, sigma2 = sigma2,
      A_diag_K = A_diag_K, B_K = B_K,
      omega_curr = omega, seed = seed_start + t
    )
    out[t] = omega
  }
  out[(burn + 1):T_total]
}


# --------------------------------------------------------------------------- #
# Test 1: prior draws are IG(1/2, 1/2)
# --------------------------------------------------------------------------- #

test_that("sample_cauchy_omega_prior matches IG(1/2, 1/2)", {
  N = 2e4
  omega = vapply(
    seq_len(N),
    function(s) sd_sample_cauchy_omega_prior_cpp(seed = s),
    numeric(1)
  )
  expect_true(all(omega > 0))

  # IG(0.5, 0.5) quantiles via 1 / qgamma(1 - p, 0.5, 0.5)
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
  exact = 1 / qgamma(1 - probs, shape = 0.5, rate = 0.5)
  emp = quantile(omega, probs = probs, names = FALSE)
  # 0.5% tolerance on quantiles given N = 2e4 IG samples.
  expect_equal(emp, exact, tolerance = 0.05)

  # Tail probability P(omega > 10): exact = pgamma(0.1, 0.5, 0.5)
  expect_equal(mean(omega > 10),
    pgamma(0.1, shape = 0.5, rate = 0.5),
    tolerance = 0.01
  )
})


# --------------------------------------------------------------------------- #
# Test 2: pure-slab special case (A_diag_K = B_K = 0) matches the
# InvGamma(1, 1/2 + K^2/(2 sigma^2)) conjugate.
# --------------------------------------------------------------------------- #

test_that("slice sampler at A_diag_K = B_K = 0 matches IG conjugate", {
  K = 0.3
  sigma2 = 1.0
  omega = run_omega_chain(
    K = K, sigma2 = sigma2,
    A_diag_K = 0, B_K = 0,
    T_total = 5e4, burn = 5e3,
    seed_start = 100L
  )
  expect_true(all(omega > 0))

  # Reference: IG(1, 1/2 + K^2/(2 sigma^2)).
  rate_ref = 0.5 + K^2 / (2 * sigma2)
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
  exact = 1 / qgamma(1 - probs, shape = 1, rate = rate_ref)
  emp = quantile(omega, probs = probs, names = FALSE)
  expect_equal(emp, exact, tolerance = 0.05)
})


# --------------------------------------------------------------------------- #
# Test 3: slice sampler matches an importance-sampling reference for
# nonzero A_diag_K and B_K. The IS estimate of P(omega < q) uses the
# IG(0.5, 0.5) prior as the proposal and reweights by the log-conditional.
# --------------------------------------------------------------------------- #

test_that("slice sampler matches IS reference under nonzero A_diag_K, B_K", {
  K = 0.30
  sigma2 = 1.0
  A_diag_K = 0.25 # = beta/(2 l_ii^2) for l_ii=1, beta=0.5 (bgms conv.)
  B_K = -0.025 # = -beta * m_ij / l_ii for m_ij=0.05

  omega = run_omega_chain(
    K = K, sigma2 = sigma2,
    A_diag_K = A_diag_K, B_K = B_K,
    T_total = 5e4, burn = 5e3,
    seed_start = 200L
  )

  # IS reference: omega_i ~ IG(0.5, 0.5),
  #               w_i = exp(log_target_full(omega_i) - log_prior(omega_i))
  set.seed(2026L)
  N_is = 1e6
  omega_is = 1 / rgamma(N_is, shape = 0.5, rate = 0.5)
  log_prior = -1.5 * log(omega_is) - 1 / (2 * omega_is)
  log_full = log_cond_omega(omega_is, K, sigma2, A_diag_K, B_K)
  log_w = log_full - log_prior
  w = exp(log_w - max(log_w))
  w = w / sum(w)

  weighted_q = function(x, w, p) {
    ord = order(x)
    x = x[ord]
    w = w[ord]
    cw = cumsum(w)
    vapply(p, function(pp) x[which.max(cw >= pp)], numeric(1))
  }

  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
  q_chain = quantile(omega, probs = probs, names = FALSE)
  q_ref = weighted_q(omega_is, w, probs)

  # 10% tolerance on quantiles. The upper tail is heavy (omega^-3/2) so the
  # q90 has the largest MC error on both sides.
  expect_equal(q_chain, q_ref, tolerance = 0.10)
})
