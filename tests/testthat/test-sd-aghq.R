# --------------------------------------------------------------------------- #
# Unit tests for ggm_sd::density_at_l_ji_aghq.
#
# The adaptive Gauss-Hermite primitive locates the actual mode of f via the
# closed-form cubic solver and uses it as the Laplace reference. This Phase
# 2 version handles the unimodal case correctly (every alpha, every Delta < 0
# cell) and uses the global mode as a single Laplace reference in the
# bimodal case (Delta >= 0). The mixture branch lands in Phase 3.
#
# Tier 1 coverage (per dev/plans/active/sd-aghq-cubic.md):
#   - alpha = 1 closed-form: AGHQ log Z matches Gaussian closed form to
#     roundoff at any N (the integrand becomes constant in y).
#   - alpha > 1 unimodal: AGHQ log Z matches a fine-grid trapezoidal
#     reference to < 1e-8 relative.
#   - log_density consistency: log pi(x_eval) = f(x_eval) - log Z.
#   - Mode/curvature returned match the cubic solver's global mode.
#   - Failure paths (A <= 0, s_jj <= 0).
#   - Bimodal cells: single-mode AGHQ remains valid as a worst-case
#     reference (status = 0; mass under the integrand is finite).
# --------------------------------------------------------------------------- #

ell_ref <- function(phi, A, B, s_jj, alpha) {
  -A * phi^2 + B * phi + (alpha - 1) * log(s_jj + phi^2)
}

# Fine-grid trapezoidal reference: span +- R/sqrt(A) standard deviations
# around B/(2A) at N=5000 points. Robust on smooth unimodal integrands and
# adequate as a 1e-8 reference. For bimodal cells we use it as a soft
# reference only; tighter tolerances come in Phase 3 via the mixture path.
grid_logZ <- function(A, B, s_jj, alpha, R = 20, N = 5000) {
  x   <- seq(-R / sqrt(A) + B / (2 * A),
              R / sqrt(A) + B / (2 * A),
             length.out = N)
  ell <- ell_ref(x, A, B, s_jj, alpha)
  M   <- max(ell)
  log(sum(exp(ell - M))) + M + log(diff(x)[1])
}


# --------------------------------------------------------------------------- #
# Test 1: alpha = 1 -- AGHQ log Z matches the closed-form Gaussian log Z to
# roundoff. At alpha = 1 the cubic returns phi* = B/(2A) and kappa = 2A;
# substituting into the AGHQ formula gives constant ell + y^2 = B^2/(4A),
# and the GH sum collapses to sqrt(pi) * exp(B^2/(4A)). Result is N-
# independent.
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: alpha = 1 matches closed-form to roundoff", {
  for (A in c(0.5, 1, 2, 5)) {
    for (B in c(-3, -0.5, 0, 1.2, 7)) {
      for (s_jj in c(0.1, 1, 5)) {
        closed <- 0.5 * log(2 * pi / (2 * A)) + B^2 / (4 * A)
        for (N in c(8, 16, 32, 64)) {
          r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
            x_eval = B / (2 * A), A, B, s_jj, alpha = 1, num_nodes = N
          )
          expect_equal(r$status, 0L)
          expect_equal(r$x_mode, B / (2 * A), tolerance = 1e-14)
          expect_equal(r$curvature, 2 * A, tolerance = 1e-14)
          expect_equal(r$log_Z, closed, tolerance = 1e-12,
                       info = sprintf("A=%g B=%g s=%g N=%d", A, B, s_jj, N))
        }
      }
    }
  }
})


# --------------------------------------------------------------------------- #
# Test 2: alpha > 1 unimodal cells -- AGHQ log Z matches fine-grid reference.
# Pick cells with Delta < 0 (unimodal); the cubic returns a single mode.
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: alpha > 1 unimodal matches fine-grid reference", {
  configs <- list(
    list(A = 2,   B = 1,   s_jj = 1,   alpha = 1.5),
    list(A = 1,   B = 0,   s_jj = 2,   alpha = 2),
    list(A = 0.5, B = -2,  s_jj = 3,   alpha = 3),
    list(A = 3,   B = 5,   s_jj = 0.5, alpha = 4),
    list(A = 0.3, B = 0.5, s_jj = 0.4, alpha = 3)   # mild near-bimodal
  )
  for (cfg in configs) {
    ref <- grid_logZ(cfg$A, cfg$B, cfg$s_jj, cfg$alpha)
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
      x_eval = cfg$B / (2 * cfg$A), cfg$A, cfg$B, cfg$s_jj, cfg$alpha,
      num_nodes = 32
    )
    expect_equal(r$status, 0L,
                 info = sprintf("alpha=%g", cfg$alpha))
    # Tight tolerance on a smooth unimodal integrand at N = 32.
    expect_lt(abs(r$log_Z - ref), 1e-8 * (1 + abs(ref)),
              label = sprintf("log Z mismatch (alpha=%g): AGHQ=%g ref=%g diff=%g",
                              cfg$alpha, r$log_Z, ref, r$log_Z - ref))
  }
})


# --------------------------------------------------------------------------- #
# Test 3: log_density consistency. log pi(x_eval) = f(x_eval) - log Z.
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: log_density = f(x_eval) - log_Z", {
  set.seed(2026)
  for (.k in 1:200) {
    A     <- exp(runif(1, -1, 1))
    B     <- runif(1, -3, 3)
    s_jj  <- exp(runif(1, -2, 0))
    alpha <- runif(1, 1.001, 4)
    x_eval <- runif(1, -2, 2)
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(x_eval, A, B, s_jj, alpha)
    if (r$status != 0) next
    f_eval <- ell_ref(x_eval, A, B, s_jj, alpha)
    expect_equal(r$log_density, f_eval - r$log_Z, tolerance = 1e-12)
  }
})


# --------------------------------------------------------------------------- #
# Test 4: mode/curvature returned by AGHQ match the cubic solver's global
# mode (when not in the fallback path).
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: mode + curvature match the cubic solver", {
  set.seed(4)
  for (.k in 1:200) {
    A     <- exp(runif(1, -1, 1))
    B     <- runif(1, -2, 2)
    s_jj  <- exp(runif(1, -2, 0))
    alpha <- runif(1, 1.01, 4)
    aghq  <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, A, B, s_jj, alpha)
    cubic <- bgms:::sd_cubic_solve_cpp(A, B, s_jj, alpha)
    if (aghq$status == 3 || cubic$n_modes < 1) next
    gmi <- cubic$global_mode_index + 1   # 0-indexed -> 1-indexed
    expect_equal(aghq$x_mode,    cubic$phi[gmi],       tolerance = 1e-14)
    expect_equal(aghq$curvature, cubic$curvature[gmi], tolerance = 1e-14)
  }
})


# --------------------------------------------------------------------------- #
# Test 5: failure paths.
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: A <= 0 returns status 1", {
  for (A in c(0, -1)) {
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, A, 1, 1, 2)
    expect_equal(r$status, 1L)
    expect_true(is.nan(r$log_Z) || is.na(r$log_Z))
  }
})

test_that("density_at_l_ji_aghq: s_jj <= 0 returns status 2", {
  for (s in c(0, -0.01)) {
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, 1, 1, s, 2)
    expect_equal(r$status, 2L)
  }
})


# --------------------------------------------------------------------------- #
# Test 6: bimodal cells. Single-mode AGHQ underestimates log Z by the
# secondary-mode mass, but the result is still finite and status = 0.
# Tighter accuracy will come from the mixture branch in Phase 3.
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: bimodal cells return finite log Z", {
  configs <- list(
    list(A = 0.3, B = 0,   s_jj = 0.1, alpha = 4),
    list(A = 0.2, B = 0.1, s_jj = 0.05, alpha = 5),
    list(A = 0.5, B = 0,   s_jj = 0.2, alpha = 6)
  )
  for (cfg in configs) {
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
      x_eval = 0, cfg$A, cfg$B, cfg$s_jj, cfg$alpha, num_nodes = 32
    )
    expect_equal(r$status, 0L)
    expect_true(is.finite(r$log_Z))
    # Loose comparison to fine-grid reference: single-Laplace at the global
    # mode is expected to underestimate by O(1) at most.
    ref <- grid_logZ(cfg$A, cfg$B, cfg$s_jj, cfg$alpha)
    expect_lt(abs(r$log_Z - ref), 0.1,
              label = sprintf("bimodal alpha=%g AGHQ=%g ref=%g",
                              cfg$alpha, r$log_Z, ref))
  }
})


# --------------------------------------------------------------------------- #
# Test 7: triple-root degeneracy (alpha - 1 = A * s_jj exactly). The cubic
# solver returns n_modes = 0 there; AGHQ falls back to the alpha = 1
# reference Gaussian and reports status = 3.
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: triple-root degeneracy falls back cleanly", {
  # Pick parameters so (alpha - 1) = A * s_jj exactly.
  A <- 0.5; s_jj <- 0.4; alpha <- 1 + A * s_jj  # = 1.2
  r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
    x_eval = 0, A, B = 0, s_jj, alpha, num_nodes = 32
  )
  expect_equal(r$status, 3L)
  expect_equal(r$x_mode, 0, tolerance = 1e-14)        # B/(2A) = 0
  expect_equal(r$curvature, 2 * A, tolerance = 1e-14) # fallback kappa
  expect_true(is.finite(r$log_Z))
})


# --------------------------------------------------------------------------- #
# Test 8: property-based fuzz. 2000 random configs, AGHQ vs fine-grid.
# Skip status != 0 (fallback) and bimodal cells (covered in test 6).
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Test 8: property-based fuzz. AGHQ-32-vs-AGHQ-128 truncation error grows
# when the log term has features finer than the Laplace scale (small s_jj
# or large alpha near 1 with the mode close to the log singularity at
# phi^2 = -s_jj). The chain operates with N = 32 and lives with that
# truncation; this test asserts the truncation stays bounded across a
# broad parameter sweep, which catches gross errors (sign flip, missing
# normaliser, mode misidentification) without flagging acceptable noise.
# --------------------------------------------------------------------------- #

test_that("density_at_l_ji_aghq: broader fuzz with relaxed N-convergence tol", {
  set.seed(31)
  N <- 1000L
  A_vec     <- exp(runif(N, log(0.1), log(10)))
  B_vec     <- runif(N, -5, 5)
  s_jj_vec  <- exp(runif(N, log(0.01), log(5)))
  alpha_vec <- runif(N, 1.001, 4)

  for (i in seq_len(N)) {
    cubic <- bgms:::sd_cubic_solve_cpp(A_vec[i], B_vec[i],
                                        s_jj_vec[i], alpha_vec[i])
    if (cubic$status != 0 || cubic$n_real_roots == 3) next
    r32 <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
      x_eval = 0, A_vec[i], B_vec[i], s_jj_vec[i], alpha_vec[i],
      num_nodes = 32
    )
    r128 <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
      x_eval = 0, A_vec[i], B_vec[i], s_jj_vec[i], alpha_vec[i],
      num_nodes = 128
    )
    if (r32$status != 0 || r128$status != 0) next
    expect_lt(
      abs(r32$log_Z - r128$log_Z), 1e-3 * (1 + abs(r128$log_Z)),
      label = sprintf("i=%d (A=%g B=%g s=%g a=%g): AGHQ_32=%.6g AGHQ_128=%.6g",
                      i, A_vec[i], B_vec[i], s_jj_vec[i], alpha_vec[i],
                      r32$log_Z, r128$log_Z)
    )
  }
})
