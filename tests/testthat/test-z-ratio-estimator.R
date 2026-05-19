# --------------------------------------------------------------------------- #
# Phase 3: Russian-Roulette V(Γ, U) estimator of 1/Z(Γ).
#
# Exercises:
#   degord_V_at_Gamma_pi_cpp   - signed V at fixed (G_pi, c, rho)
#   degord_draw_U_rr_cpp       - fresh (K_depth, pools_t) draw
#
# Acceptance per the Phase 3 plan section:
#   - Mean V across many independent (K_depth, pools) draws sits within MC
#     CI of 1/Z(Γ).
#   - K_depth has the Geom(1 - rho) distribution.
#   - Outputs are deterministic in the SafeRNG seed.
# --------------------------------------------------------------------------- #


# ---- Helpers -----------------------------------------------------------------

draw_random_graph <- function(q, seed, p_edge = 0.5) {
  set.seed(seed)
  G <- matrix(0L, q, q)
  if (q < 2) return(G)
  for (i in 1:(q - 1)) for (j in (i + 1):q)
    if (runif(1) < p_edge) { G[i, j] <- 1L; G[j, i] <- 1L }
  G
}


# ---- V is unbiased for 1/Z within MC noise ---------------------------------

test_that("mean V across independent pools converges to 1/Z within MC noise", {
  # Test V's pointwise unbiasedness: E[V(Γ, U)] = 1/Z(Γ).
  # "Truth" is built from a high-precision log_Zhat at very large M_inner,
  # averaged over many independent pools to dampen MC noise to << test
  # tolerance.
  q <- 5L
  G_pi <- draw_random_graph(q, seed = 1L)
  alpha <- 1.0; beta <- 1.0; sigma <- 1.0; delta <- 0.5

  # Truth proxy: M = 5000, n_truth = 20 → ~22000 effective samples; the SE on
  # log_Zhat is dominated by Z_truth's MC noise which is ~ sd/sqrt(n_truth).
  M_truth <- 5000L
  n_truth <- 20L
  log_Zhat_truth <- vapply(seq_len(n_truth), function(k) {
    pool_t <- degord_draw_bartlett_pool_cpp(q, M_truth, seed = 9000L + k)
    degord_log_Zhat_pi_from_pool_cpp(
      pool_t, G_pi, alpha, beta, sigma, delta, 0L
    )
  }, double(1))
  Z_truth    <- exp(mean(log_Zhat_truth))
  invZ_truth <- 1 / Z_truth

  # V estimator
  kappa   <- 1.5
  rho     <- 0.5
  M_inner <- 100L
  c_val   <- kappa * Z_truth
  n_outer <- 1000L
  V_samples <- vapply(seq_len(n_outer), function(m) {
    U <- degord_draw_U_rr_cpp(M_inner, q, rho, seed = 100000L + m)
    degord_V_at_Gamma_pi_cpp(
      U$K_depth, U$pools_t, G_pi,
      alpha, beta, sigma, delta, c_val, rho, 0L
    )
  }, double(1))

  mean_V  <- mean(V_samples)
  se_mean <- sd(V_samples) / sqrt(n_outer)
  z_score <- (mean_V - invZ_truth) / se_mean

  # Allow up to 4 SD - tail probability < 1e-4, generous enough to absorb
  # Z_truth's MC residual without flaking on CI.
  expect_lt(abs(z_score), 4.0,
            label = sprintf("|z|=%.2f (mean_V=%.4f vs 1/Z=%.4f, SE=%.4g)",
                            abs(z_score), mean_V, invZ_truth, se_mean))
})


# ---- K_depth follows Geom(1 - rho) -----------------------------------------

test_that("K_depth ~ Geom(1 - rho) under SafeRNG", {
  rho <- 0.5
  n_draws <- 10000L
  K_samples <- vapply(seq_len(n_draws), function(m) {
    U <- degord_draw_U_rr_cpp(M_inner = 1L, q = 3L, rho = rho,
                              seed = 50000L + m)
    U$K_depth
  }, integer(1))
  # Geom(1 - rho) has mean rho / (1 - rho) and P(K=0) = (1 - rho).
  # Tolerances are absolute (testthat's `tolerance` is relative for numeric);
  # we compute the absolute gap directly to avoid relative-tolerance traps on
  # the small tail probability.
  expect_lt(abs(mean(K_samples) - rho / (1 - rho)), 0.05)
  expect_lt(abs(mean(K_samples == 0L) - (1 - rho)), 0.02)
  # Tail: P(K >= 5) = rho^5 = 0.03125. SE at n=10000 ~ 0.0017, allow 3*SE.
  expect_lt(abs(mean(K_samples >= 5L) - rho^5), 0.005)
})


# ---- K_depth = 0 short-circuits to V = 1/c ---------------------------------

test_that("V at K_depth = 0 equals 1/c exactly", {
  q <- 5L
  G_pi <- draw_random_graph(q, seed = 17L)
  empty_pools <- list()
  for (c_val in c(0.5, 1.0, 3.7, 12.0)) {
    v <- degord_V_at_Gamma_pi_cpp(
      0L, empty_pools, G_pi,
      1.0, 1.0, 1.0, 0.5, c_val, 0.5, 0L
    )
    expect_equal(v, 1 / c_val, tolerance = 1e-12,
                 info = sprintf("c_val=%g", c_val))
  }
})


# ---- draw_U_degord_rr is deterministic in seed -----------------------------

test_that("draw_U_degord_rr produces identical output under identical seed", {
  U_a <- degord_draw_U_rr_cpp(M_inner = 30L, q = 5L, rho = 0.5, seed = 42L)
  U_b <- degord_draw_U_rr_cpp(M_inner = 30L, q = 5L, rho = 0.5, seed = 42L)
  expect_equal(U_a$K_depth, U_b$K_depth)
  for (n in seq_along(U_a$pools_t))
    expect_identical(U_a$pools_t[[n]], U_b$pools_t[[n]])

  # Different seed → different output (with high prob).
  U_c <- degord_draw_U_rr_cpp(M_inner = 30L, q = 5L, rho = 0.5, seed = 43L)
  has_pool <- length(U_c$pools_t) > 0L && length(U_a$pools_t) > 0L &&
              U_a$K_depth == U_c$K_depth
  if (has_pool) {
    diffs <- abs(U_a$pools_t[[1]] - U_c$pools_t[[1]])
    expect_gt(max(diffs), 0)
  }
})


# ---- Signed V can flip sign for K_depth >= 1 -------------------------------

test_that("V tracks the alternating-series sign", {
  # When Zhat_n - c is consistently small relative to c, V stays positive;
  # when Zhat_n - c can be larger than c (or negative), V can change sign.
  # Force an extreme case: choose c far away from 1/Z so the (Zhat - c)/c
  # products grow and the alternating series can produce negative V.
  q <- 5L
  G_pi <- draw_random_graph(q, seed = 23L)
  alpha <- 1.0; beta <- 1.0; sigma <- 1.0; delta <- 0.5
  rho   <- 0.6
  M_inner <- 50L
  # Tiny c (far below 1/Z) inflates each (Zhat - c)/c factor; K_depth >= 1
  # samples will fluctuate widely in sign.
  c_val <- 1e-4
  n_outer <- 200L
  V_samples <- vapply(seq_len(n_outer), function(m) {
    U <- degord_draw_U_rr_cpp(M_inner, q, rho, seed = 200000L + m)
    if (U$K_depth == 0L) return(NA_real_)
    degord_V_at_Gamma_pi_cpp(
      U$K_depth, U$pools_t, G_pi,
      alpha, beta, sigma, delta, c_val, rho, 0L
    )
  }, double(1))
  V_samples <- V_samples[!is.na(V_samples) & is.finite(V_samples)]
  # We should see SOME negative values when c is small (alternating series).
  expect_true(any(V_samples < 0),
              label = "no negative V samples under small c_val (alternating sign expected)")
})
