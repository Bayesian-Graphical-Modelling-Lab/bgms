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
    if (runif(1) < p_edge) {
      G[i, j] <- 1L
      G[j, i] <- 1L
    }
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
  alpha <- 1.0
  beta  <- 1.0
  sigma <- 1.0
  delta <- 0.5

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


# ---- Log-space V matches linear V at small p, survives underflow at large p --

test_that("V_log_at_Gamma_pi matches V_at_Gamma_pi in the safe regime", {
  # F5 bit-equality smoke. log-space and linear forms operate on the same
  # pools_t / G_pi / chain_aux and should produce identical |V| (modulo
  # FP reordering in the signed log-sum-exp).
  alpha <- 1.0
  beta  <- 1.0
  sigma <- 1.0
  delta <- 0.5
  kappa <- 1.0
  rho   <- 0.5
  for (q in c(5L, 10L, 20L)) {
    G_pi <- draw_random_graph(q, seed = 31L + q)
    # Centre c near 1/Z by approximating log_Z with a single Bartlett pool
    # at large M. Tolerance is set on differences, so the centring need not
    # be perfect; we just need the linear form to stay finite.
    pool_truth <- degord_draw_bartlett_pool_cpp(q, 2000L, seed = 41L + q)
    log_Z <- degord_log_Zhat_pi_from_pool_cpp(
      pool_truth, G_pi, alpha, beta, sigma, delta, 0L
    )
    c_val <- kappa * exp(log_Z)
    log_c <- log(kappa) + log_Z
    M_inner <- 100L
    n_samples <- 50L
    diffs <- numeric(n_samples)
    n_finite <- 0L
    for (m in seq_len(n_samples)) {
      U <- degord_draw_U_rr_cpp(M_inner, q, rho, seed = 300000L + 1000L * q + m)
      v_lin <- degord_V_at_Gamma_pi_cpp(
        U$K_depth, U$pools_t, G_pi,
        alpha, beta, sigma, delta, c_val, rho, 0L
      )
      v_log <- degord_V_log_at_Gamma_pi_cpp(
        U$K_depth, U$pools_t, G_pi,
        alpha, beta, sigma, delta, log_c, rho, 0L
      )
      if (!is.finite(v_lin) || v_lin == 0 ||
          !is.finite(v_log$log_abs) || v_log$sign == 0) next
      lhs <- log(abs(v_lin))
      rhs <- v_log$log_abs
      n_finite <- n_finite + 1L
      diffs[n_finite] <- abs(lhs - rhs)
      # Sign agreement is exact (no FP reordering can flip sign).
      expect_equal(sign(v_lin), v_log$sign,
                   info = sprintf("q=%d, m=%d, K=%d", q, m, U$K_depth))
    }
    diffs <- diffs[seq_len(n_finite)]
    # At q = 5 the series barely truncates; tolerance ~ 1e-12. At q = 20 the
    # alternating series can accumulate more cancellation; loosen to 1e-9.
    tol <- if (q <= 5L) 1e-12 else if (q <= 10L) 1e-10 else 1e-9
    expect_true(max(diffs) < tol,
                info = sprintf("q=%d, max|log_abs gap|=%.3g (tol=%.3g, n_finite=%d)",
                               q, max(diffs), tol, n_finite))
  }
})


test_that("V_log_at_Gamma_pi stays finite where V_at_Gamma_pi underflows", {
  # F5 motivating regime: when log_c is far below 0 in magnitude (c = exp(log_c)
  # flushes to 0 in double precision), the linear form returns Inf / NaN. The
  # log-space form must remain finite.
  q <- 5L
  G_pi <- draw_random_graph(q, seed = 7L)
  alpha <- 1.0
  beta  <- 1.0
  sigma <- 1.0
  delta <- 0.5
  rho   <- 0.5
  # log_c = -3500 mirrors what log_Z_NLO + log_kappa can reach at p = 100;
  # exp(-3500) is 0 in double precision.
  log_c <- -3500
  c_val <- exp(log_c)
  expect_equal(c_val, 0)
  M_inner <- 50L
  n_samples <- 20L
  n_log_finite <- 0L
  n_lin_finite <- 0L
  for (m in seq_len(n_samples)) {
    U <- degord_draw_U_rr_cpp(M_inner, q, rho, seed = 400000L + m)
    v_lin <- degord_V_at_Gamma_pi_cpp(
      U$K_depth, U$pools_t, G_pi,
      alpha, beta, sigma, delta, c_val, rho, 0L
    )
    v_log <- degord_V_log_at_Gamma_pi_cpp(
      U$K_depth, U$pools_t, G_pi,
      alpha, beta, sigma, delta, log_c, rho, 0L
    )
    if (is.finite(v_lin) && v_lin != 0) n_lin_finite <- n_lin_finite + 1L
    if (is.finite(v_log$log_abs) && v_log$sign != 0) n_log_finite <- n_log_finite + 1L
  }
  # The linear form should mostly explode here.
  expect_lt(n_lin_finite, n_samples / 2L)
  # The log-space form should mostly stay finite.
  expect_gt(n_log_finite, n_samples / 2L)
})


# ---- F6: paired V_log with within-pool cache reuse ------------------------

test_that("paired V_log matches two independent V_log calls bit-equal", {
  # F6 bit-equality. The paired call shares the inner Phi-build across
  # G_pi_curr / G_pi_star by caching (rw_head, S_trail) under a_curr.
  # The single-graph call rebuilds Phi from scratch for each graph. Both
  # paths should produce identical (log|V|, sign(V)) for both curr and
  # star to FP-reordering tolerance.
  alpha <- 1.0
  beta  <- 1.0
  sigma <- 1.0
  delta <- 0.5
  rho   <- 0.5
  cases <- list(
    list(q = 5L,  K = 0L, seed = 11L),
    list(q = 5L,  K = 2L, seed = 12L),
    list(q = 5L,  K = 5L, seed = 13L),
    list(q = 10L, K = 0L, seed = 21L),
    list(q = 10L, K = 2L, seed = 22L),
    list(q = 10L, K = 5L, seed = 23L),
    list(q = 20L, K = 0L, seed = 31L),
    list(q = 20L, K = 2L, seed = 32L),
    list(q = 20L, K = 5L, seed = 33L)
  )
  for (case in cases) {
    q <- case$q
    K <- case$K
    seed <- case$seed
    # G_pi_curr is a random symmetric 0/1 matrix; G_pi_star toggles only
    # the trailing slot (q-2, q-1). This mirrors what the chain hot path
    # passes after degord_permutation places the toggled edge at (q-2, q-1).
    G_pi_curr <- draw_random_graph(q, seed = seed, p_edge = 0.5)
    G_pi_star <- G_pi_curr
    G_pi_star[q - 1L, q] <- 1L - G_pi_curr[q - 1L, q]
    G_pi_star[q, q - 1L] <- G_pi_star[q - 1L, q]
    # Build a K-deep pool list deterministically.
    M_inner <- 50L
    pools_t <- vector("list", K)
    for (n in seq_len(K)) {
      pools_t[[n]] <- degord_draw_bartlett_pool_cpp(
        q, M_inner, seed = 700000L + 1000L * seed + n)
    }
    # log_c values: choose so c ~ exp(log_Z) is near 1/Z (any reasonable
    # value lets V be finite).
    if (K > 0L) {
      log_Z_proxy <- degord_log_Zhat_pi_from_pool_cpp(
        pools_t[[1L]], G_pi_curr, alpha, beta, sigma, delta, 0L
      )
    } else {
      log_Z_proxy <- 0
    }
    log_c_curr <- log_Z_proxy
    log_c_star <- log_Z_proxy + 0.1   # small offset; both must stay finite
    # Independent path.
    v_curr_ind <- degord_V_log_at_Gamma_pi_cpp(
      K, pools_t, G_pi_curr,
      alpha, beta, sigma, delta, log_c_curr, rho, 0L
    )
    v_star_ind <- degord_V_log_at_Gamma_pi_cpp(
      K, pools_t, G_pi_star,
      alpha, beta, sigma, delta, log_c_star, rho, 0L
    )
    # Paired path.
    v_pair <- degord_V_log_pair_at_Gamma_curr_star_cpp(
      K, pools_t, G_pi_curr, G_pi_star,
      alpha, beta, sigma, delta, log_c_curr, log_c_star, rho, 0L
    )
    # Sign agreement is exact.
    expect_equal(v_pair$curr$sign, v_curr_ind$sign,
                 info = sprintf("q=%d, K=%d (curr)", q, K))
    expect_equal(v_pair$star$sign, v_star_ind$sign,
                 info = sprintf("q=%d, K=%d (star)", q, K))
    # log_abs agreement modulo FP reordering inside the inner kernel.
    # Both paths execute the same per-sample row_logw arithmetic; only the
    # row q-2 evaluation differs (cache reuses S_trail). Tolerance must
    # absorb the order-of-additions difference in the log-sum-exp wrap-up
    # at large q.
    tol <- if (q <= 10L) 1e-12 else 1e-9
    if (is.finite(v_pair$curr$log_abs) && is.finite(v_curr_ind$log_abs)) {
      expect_lt(abs(v_pair$curr$log_abs - v_curr_ind$log_abs), tol,
                label = sprintf("q=%d, K=%d, curr log_abs gap=%.3g",
                                q, K,
                                abs(v_pair$curr$log_abs - v_curr_ind$log_abs)))
    }
    if (is.finite(v_pair$star$log_abs) && is.finite(v_star_ind$log_abs)) {
      expect_lt(abs(v_pair$star$log_abs - v_star_ind$log_abs), tol,
                label = sprintf("q=%d, K=%d, star log_abs gap=%.3g",
                                q, K,
                                abs(v_pair$star$log_abs - v_star_ind$log_abs)))
    }
  }
})


test_that("log_Zhat_star_from_cache matches fresh log_Zhat at G_pi_star", {
  # Direct check on the cache adapter. log_Zhat_star_from_cache must
  # produce the same value as a fresh log_Zhat_pi_from_pool on the star
  # graph, to FP-reordering tolerance (the cache path reuses row_logw[r]
  # for r != q-2 instead of recomputing).
  alpha <- 1.0
  beta  <- 1.0
  sigma <- 1.0
  delta <- 0.5
  for (q in c(5L, 10L, 20L)) {
    G_pi_curr <- draw_random_graph(q, seed = 41L + q, p_edge = 0.5)
    G_pi_star <- G_pi_curr
    G_pi_star[q - 1L, q] <- 1L - G_pi_curr[q - 1L, q]
    G_pi_star[q, q - 1L] <- G_pi_star[q - 1L, q]
    pool <- degord_draw_bartlett_pool_cpp(
      q, M_inner = 200L, seed = 800000L + q)
    log_Zhat_star_cache <- degord_log_Zhat_star_from_cache_cpp(
      pool, G_pi_curr, G_pi_star,
      alpha, beta, sigma, delta, 0L
    )
    log_Zhat_star_direct <- degord_log_Zhat_pi_from_pool_cpp(
      pool, G_pi_star, alpha, beta, sigma, delta, 0L
    )
    tol <- if (q <= 10L) 1e-12 else 1e-9
    expect_lt(abs(log_Zhat_star_cache - log_Zhat_star_direct), tol,
              label = sprintf("q=%d, cache=%.6f, direct=%.6f, gap=%.3g",
                              q, log_Zhat_star_cache, log_Zhat_star_direct,
                              abs(log_Zhat_star_cache - log_Zhat_star_direct)))
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
  alpha <- 1.0
  beta  <- 1.0
  sigma <- 1.0
  delta <- 0.5
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
