# --------------------------------------------------------------------------- #
# DEGORD sampler (Phase 2): port of the Bartlett-Cholesky importance
# sampler for log Zhat(G) from ~/SV/Z/R/src/degord_sampler.h (v4).
#
# Exercises:
#   degord_chain_aux_cpp                     - ChainAux constants
#   degord_pi_aux_cpp                        - PiAux: nu_pi, E_count, log_C0
#   degord_permute_graph_cpp                 - DEGORD permutation -> (q-2, q-1)
#   degord_log_Zhat_pi_from_pool_cpp         - main Zhat-from-pool entry
#   degord_delta_log_Zhat_pi_toggle_cpp      - cache-based delta under trailing toggle
#   degord_draw_bartlett_pool_cpp            - SafeRNG-based standard-normal pool
#
# Ground truth is a fixture
# (tests/testthat/fixtures/degord_sampler_reference.rds) pre-generated from
# the z reference (see dev/numerical_analyses/generate_degord_fixture.R).
# --------------------------------------------------------------------------- #


# ---- log_Zhat_pi_from_pool: bit-parity against z on shared noise pool -------

test_that("log_Zhat_pi_from_pool matches the z reference bit-exact", {
  fixture_path <- testthat::test_path("fixtures", "degord_sampler_reference.rds")
  cases <- readRDS(fixture_path)
  for (case in cases) {
    for (alpha in c(1.0, 2.0)) {
      for (delta in c(0.0, 0.5, 1.0)) {
        for (tilt_mode in c(0L, 1L)) {
          key <- sprintf("a%g_d%g_t%d", alpha, delta, tilt_mode)
          ours <- degord_log_Zhat_pi_from_pool_cpp(
            case$pool_t, case$G_pi, alpha, 1.0, 1.0, delta, tilt_mode
          )
          ref <- case$values[[key]]
          expect_equal(
            ours, ref,
            tolerance = 1e-12,
            info = sprintf("q=%d rep=%d %s", case$q, case$rep, key)
          )
        }
      }
    }
  }
})


# ---- DEGORD permutation correctness ----------------------------------------

test_that("DEGORD permutation sends (i, j) to (q-2, q-1)", {
  set.seed(31)
  for (q in c(3L, 5L, 7L)) {
    G <- matrix(0L, q, q)
    for (i in 1:(q - 1)) for (j in (i + 1):q)
      if (runif(1) < 0.5) {
        G[i, j] <- 1L
        G[j, i] <- 1L
      }
    for (i0 in 0:(q - 2)) {
      for (j0 in (i0 + 1):(q - 1)) {
        G_pi <- degord_permute_graph_cpp(G, i0, j0)
        # The DEGORD permutation sends (i, j) -> (q-2, q-1). So
        # G_pi[q-1, q] (1-based) must equal G[i+1, j+1] (1-based).
        expect_equal(
          G_pi[q - 1, q], G[i0 + 1, j0 + 1],
          info = sprintf("q=%d (i,j)=(%d,%d)", q, i0, j0)
        )
        # G_pi must be symmetric, integer-valued.
        expect_true(isSymmetric(G_pi))
        expect_true(all(G_pi %in% c(0L, 1L)))
      }
    }
  }
})


# ---- ChainAux constants: bgms vs z (interactive, requires z wrapper) -------

# This is more thorough but requires the z source to be available, so it is
# gated behind the wrapper. The fixture above is sufficient for CI parity;
# this block runs locally when the developer has the z project on disk.

test_that("ChainAux nu-tables match the z reference (local-only)", {
  wrapper <- testthat::test_path("..", "..", "dev", "numerical_analyses",
                                  "z_degord_wrapper.cpp")
  skip_if_not(file.exists(wrapper),
              "z DEGORD wrapper not present (regenerate fixture instead).")
  z_src <- "~/Dropbox/Projecten/sv/z/R/src/degord_sampler.h"
  skip_if_not(file.exists(z_src), "z reference header not present.")
  suppressMessages(Rcpp::sourceCpp(wrapper))

  fields <- c("sigma_diag", "nu_chi_df", "nu_mu_l", "nu_H_e_saddle",
              "nu_lgamma_half_k", "nu_diag_const", "nu_slab_const_saddle",
              "nu_per_vertex")
  for (q in c(3L, 5L, 10L)) {
    for (alpha in c(1.0, 2.0)) {
      for (delta in c(0.0, 0.5, 1.0)) {
        z <- z_degord_chain_aux(q, alpha, 1.0, 1.0, delta)
        b <- degord_chain_aux_cpp(q, alpha, 1.0, 1.0, delta)
        for (f in fields) {
          expect_equal(
            b[[f]], z[[f]],
            tolerance = 1e-12,
            info = sprintf("q=%d alpha=%g delta=%g %s", q, alpha, delta, f)
          )
        }
      }
    }
  }
})


# ---- delta_log_Zhat_pi_toggle: cache-delta vs direct full-recompute ---------

# This used to be a known-discrepancy note (~0.17 nat gap between cached
# delta and the direct log_Zhat(star) - log_Zhat(curr)). Two bugs in the
# v4 cache trick caused it:
#
#   1. rw_head spans rows 0..q-3 but row q-1's diagonal log-weight is
#      invariant across (curr, star) AND sample-dependent; omitting it
#      from the star aggregation left a sample-shifted bias.
#      Fix: rw_head extended to include row q-1 (Option 2; companion
#      Z-side proposal 2026-05-19, applied at degord_sampler.cpp:302).
#
#   2. delta_log_Zhat_pi_toggle passed z_trail = 0.0 to row_qm2_logw_from_S
#      under slab_tilt_mode == 1, dropping the saddle-shifted slab
#      innovation noise[q + edge_offset(q-2, q-1)] = noise[q + (q-2)(q+1)/2].
#      Fix: pass the actual slab slot via noise_pool.colptr(slab_idx)
#      (applied at degord_sampler.cpp:398-404).
#
# With both fixes the cache-delta matches the direct full-recompute
# difference at machine precision under both slab_tilt_modes. This is
# the SBC-relevant invariant — without it the chain's delta_log_Zhat
# acceptance contribution silently drifts from log Zhat(star) /
# log Zhat(curr) by a tilt-amplified bias.

test_that("delta_log_Zhat_pi_toggle equals direct full-recompute at machine precision", {
  # CI-portable regression net for the cache-fix + z_trail-fix pair: the
  # cache-trick delta MUST equal the direct log_Zhat(star) - log_Zhat(curr)
  # at machine precision under BOTH slab_tilt_modes. Without this assertion
  # silent regressions in the cache aggregation (rw_head, S_trail) or in
  # the z_trail slot index would re-introduce a tilt-amplified bias that
  # only shows up under SBC stress at high delta.
  draw_G <- function(q, seed) {
    set.seed(seed)
    G <- matrix(0L, q, q)
    for (i in 1:(q - 1)) for (j in (i + 1):q)
      if (runif(1) < 0.5) {
        G[i, j] <- 1L
        G[j, i] <- 1L
      }
    G
  }

  # Explicit regression row: q=10, delta=0, slab_tilt_mode=1, toggle (3, 9).
  # This is the cell that surfaced the z_trail bug during the bgms port.
  {
    set.seed(2026)
    q <- 10L
    G <- draw_G(q, q + 13L)
    dim_pool <- q + q * (q - 1) / 2
    M <- 100L
    pool   <- matrix(rnorm(M * dim_pool), M, dim_pool)
    pool_t <- t(pool)
    d_cache <- degord_delta_log_Zhat_pi_toggle_cpp(
      pool, pool_t, G, 3L, 9L, 1.0, 1.0, 1.0, 0.0, 1L
    )
    G_pi_curr <- degord_permute_graph_cpp(G, 3L, 9L)
    G_pi_star <- G_pi_curr
    G_pi_star[q - 1, q] <- 1L - G_pi_curr[q - 1, q]
    G_pi_star[q, q - 1] <- G_pi_star[q - 1, q]
    lZ_curr <- degord_log_Zhat_pi_from_pool_cpp(
      pool_t, G_pi_curr, 1.0, 1.0, 1.0, 0.0, 1L
    )
    lZ_star <- degord_log_Zhat_pi_from_pool_cpp(
      pool_t, G_pi_star, 1.0, 1.0, 1.0, 0.0, 1L
    )
    expect_equal(
      d_cache, lZ_star - lZ_curr, tolerance = 1e-10,
      info = "z_trail-fix regression row: q=10, delta=0, tilt=1, toggle (3, 9)"
    )
  }

  # Full sweep across (q, alpha, delta, tilt_mode, toggles).
  for (q in c(5L, 10L)) {
    set.seed(2026 + q)
    G <- draw_G(q, q + 13L)
    dim_pool <- q + q * (q - 1) / 2
    M <- 50L
    pool   <- matrix(rnorm(M * dim_pool), M, dim_pool)
    pool_t <- t(pool)
    for (alpha in c(1.0, 2.0)) {
      for (delta in c(0.0, 0.5, 1.0, 2.0)) {
        for (tilt_mode in c(0L, 1L)) {
          for (i0 in 0:(q - 2)) {
            for (j0 in (i0 + 1):(q - 1)) {
              d_cache <- degord_delta_log_Zhat_pi_toggle_cpp(
                pool, pool_t, G, i0, j0, alpha, 1.0, 1.0, delta, tilt_mode
              )
              G_pi_curr <- degord_permute_graph_cpp(G, i0, j0)
              G_pi_star <- G_pi_curr
              G_pi_star[q - 1, q] <- 1L - G_pi_curr[q - 1, q]
              G_pi_star[q, q - 1] <- G_pi_star[q - 1, q]
              lZ_curr <- degord_log_Zhat_pi_from_pool_cpp(
                pool_t, G_pi_curr, alpha, 1.0, 1.0, delta, tilt_mode
              )
              lZ_star <- degord_log_Zhat_pi_from_pool_cpp(
                pool_t, G_pi_star, alpha, 1.0, 1.0, delta, tilt_mode
              )
              expect_equal(
                d_cache, lZ_star - lZ_curr, tolerance = 1e-10,
                info = sprintf("q=%d alpha=%g delta=%g tilt=%d (%d,%d)",
                               q, alpha, delta, tilt_mode, i0, j0)
              )
            }
          }
        }
      }
    }
  }
})


test_that("delta_log_Zhat_pi_toggle matches the z reference bit-exact", {
  wrapper <- testthat::test_path("..", "..", "dev", "numerical_analyses",
                                  "z_degord_wrapper.cpp")
  skip_if_not(file.exists(wrapper),
              "z DEGORD wrapper not present (delta-toggle z parity is local-only).")
  z_src <- "~/Dropbox/Projecten/sv/z/R/src/degord_sampler.h"
  skip_if_not(file.exists(z_src), "z reference header not present.")
  suppressMessages(Rcpp::sourceCpp(wrapper))

  set.seed(99)
  for (q in c(3L, 5L, 7L)) {
    G <- matrix(0L, q, q)
    for (i in 1:(q - 1)) for (j in (i + 1):q)
      if (runif(1) < 0.5) {
        G[i, j] <- 1L
        G[j, i] <- 1L
      }
    dim_pool <- q + q * (q - 1) / 2
    M <- 50L
    pool   <- matrix(rnorm(M * dim_pool), M, dim_pool)
    pool_t <- t(pool)

    for (alpha in c(1.0, 2.0)) {
      for (delta in c(0.0, 0.5)) {
        for (tilt_mode in c(0L, 1L)) {
          for (i0 in 0:(q - 2)) {
            for (j0 in (i0 + 1):(q - 1)) {
              d_bgms <- degord_delta_log_Zhat_pi_toggle_cpp(
                pool, pool_t, G, i0, j0, alpha, 1.0, 1.0, delta, tilt_mode
              )
              d_z <- z_degord_delta_log_Zhat_pi_toggle(
                pool, pool_t, G, i0, j0, alpha, 1.0, 1.0, delta, tilt_mode
              )
              expect_equal(
                d_bgms, d_z, tolerance = 1e-12,
                info = sprintf("q=%d (i,j)=(%d,%d) alpha=%g delta=%g tilt=%d",
                               q, i0, j0, alpha, delta, tilt_mode)
              )
            }
          }
        }
      }
    }
  }
})


# ---- Variance of log Zhat scales as 1/M_inner ------------------------------

test_that("variance of log Zhat scales as 1/M_inner (Phase 2 acceptance)", {
  q <- 5L
  set.seed(42)
  G_pi <- matrix(0L, q, q)
  for (i in 1:(q - 1)) for (j in (i + 1):q)
    if (runif(1) < 0.5) {
      G_pi[i, j] <- 1L
      G_pi[j, i] <- 1L
    }
  alpha <- 2.0
  beta  <- 1.0
  sigma <- 1.0
  delta <- 0.5

  M_grid <- c(30L, 1000L)
  n_reps <- 200L
  vars <- numeric(length(M_grid))
  for (k in seq_along(M_grid)) {
    M <- M_grid[k]
    vals <- vapply(seq_len(n_reps), function(r) {
      pool_t <- degord_draw_bartlett_pool_cpp(q, M, seed = r + 1000L * k)
      degord_log_Zhat_pi_from_pool_cpp(
        pool_t, G_pi, alpha, beta, sigma, delta
      )
    }, double(1))
    vars[k] <- var(vals)
  }
  # Under 1/M scaling, vars[1] / vars[2] approx M_grid[2] / M_grid[1] = 33.3.
  # MC noise on the variance estimate at n_reps=200 is large; allow a
  # 2x window around the theoretical ratio.
  ratio <- vars[1] / vars[2]
  expected_ratio <- M_grid[2] / M_grid[1]
  expect_gt(ratio, 0.5 * expected_ratio)
  expect_lt(ratio, 2.0 * expected_ratio)
})


# ---- SafeRNG-based Bartlett pool is the right shape ------------------------

test_that("draw_bartlett_pool returns a (dim x M) standard-normal matrix", {
  for (q in c(3L, 5L, 10L)) {
    M <- 50L
    pool_t <- degord_draw_bartlett_pool_cpp(q, M, seed = 1L)
    dim_expected <- q + q * (q - 1) / 2
    expect_equal(nrow(pool_t), dim_expected, info = sprintf("q=%d", q))
    expect_equal(ncol(pool_t), M)
    # Mean ~ 0, sd ~ 1 over many samples. Use a wide window for MC noise.
    expect_lt(abs(mean(pool_t)), 0.1)
    expect_gt(sd(as.numeric(pool_t)), 0.85)
    expect_lt(sd(as.numeric(pool_t)), 1.15)
  }
})


# ---- Same SafeRNG seed produces the same pool ------------------------------

test_that("draw_bartlett_pool is deterministic in seed", {
  pool_a <- degord_draw_bartlett_pool_cpp(5L, 30L, seed = 42L)
  pool_b <- degord_draw_bartlett_pool_cpp(5L, 30L, seed = 42L)
  expect_identical(pool_a, pool_b)
  pool_c <- degord_draw_bartlett_pool_cpp(5L, 30L, seed = 43L)
  expect_false(identical(pool_a, pool_c))
})
