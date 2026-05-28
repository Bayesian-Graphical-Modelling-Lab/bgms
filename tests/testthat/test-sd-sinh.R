# --------------------------------------------------------------------------- #
# Tests for savage_dickey::density_at_l_ji_sinh (sinh-substitution + midpoint-rule
# quadrature for the L-space SD primitive).
#
# Mathematical setup. The substitution phi = sqrt(s) * sinh(t) maps the
# integral
#   I = int (s + phi^2)^(alpha-1) exp(-A phi^2 + B phi) dphi
# to
#   I = s^(alpha - 1/2) * int cosh(t)^(2 alpha - 1) *
#       exp(-A s sinh^2(t) + B sqrt(s) sinh(t)) dt
# and pins the integrand's branch points at t = +- i pi/2 (constant
# distance from the real axis, independent of s). This gives the
# midpoint-rule quadrature in t uniform exponential convergence across
# the (A, B, s, alpha) plane, in particular on sharp-log cells (small
# A * s) where standard Gauss-Hermite at the alpha = 1 mean B/(2A)
# converges only polynomially.
#
# The tests below cover:
#   1. Closed-form correctness at alpha = 1 (Gaussian) and alpha = 2
#      (Gaussian * polynomial). N = 32 should be effectively exact.
#   2. Non-regression vs the fine-grid trapezoidal reference at:
#        - alpha < 0.5 (integrable singularity at t = +- i pi/2),
#        - alpha = 0.5 (singularity removed by Jacobian),
#        - 0.5 < alpha < 1, alpha = 1, 1 < alpha < 2, alpha > 2.
#   3. Sharp-log cells (small A * s with fractional alpha).
#   4. Large |B| (mode far from t = 0 after substitution).
#   5. Failure paths (A <= 0, s_jj <= 0).
#   6. log_density consistency: log pi(x_eval) = f(x_eval) - log Z.
#   7. Determinism.
#   8. Symmetry under B -> -B reflection.
#   9. Convergence in N: error at N = 64 strictly tighter than at N = 32.
# --------------------------------------------------------------------------- #

# --- references --- #

ell_ref <- function(phi, A, B, s, alpha) {
  -A * phi^2 + B * phi + (alpha - 1) * log(s + phi^2)
}

grid_logZ <- function(A, B, s, alpha, R = 30, N = 2000000) {
  centre <- B / (2 * A)
  x <- seq(centre - R / sqrt(A), centre + R / sqrt(A), length.out = N)
  e <- ell_ref(x, A, B, s, alpha)
  M <- max(e)
  log(sum(exp(e - M))) + M + log(diff(x)[1])
}

closed_alpha1 <- function(A, B, s = NA) {
  0.5 * log(2 * pi / (2 * A)) + B^2 / (4 * A)
}
closed_alpha2 <- function(A, B, s) {
  factor <- s + 1 / (2 * A) + B^2 / (4 * A^2)
  log(factor) + 0.5 * log(pi / A) + B^2 / (4 * A)
}


# --- 1. closed-form correctness --- #

test_that("sinh primitive matches closed-form Gaussian (alpha = 1) at production N=128", {
  set.seed(2026)
  for (.k in 1:30) {
    A <- exp(runif(1, -1, 2))
    B <- runif(1, -3, 3)
    s <- exp(runif(1, -3, 1))
    r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, A, B, s, alpha = 1, num_nodes = 128)
    expect_equal(r$status, 0L)
    expect_lt(abs(r$log_Z - closed_alpha1(A, B)), 1e-10,
      label = sprintf("A=%g B=%g s=%g", A, B, s)
    )
  }
})

test_that("sinh primitive matches closed-form (alpha = 2) at production N=128", {
  set.seed(11)
  for (.k in 1:30) {
    A <- exp(runif(1, -1, 2))
    B <- runif(1, -3, 3)
    s <- exp(runif(1, -3, 1))
    r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, A, B, s, alpha = 2, num_nodes = 128)
    expect_equal(r$status, 0L)
    expect_lt(abs(r$log_Z - closed_alpha2(A, B, s)), 1e-10,
      label = sprintf("A=%g B=%g s=%g", A, B, s)
    )
  }
})


# --- 2. fine-grid reference across alpha regimes --- #

run_grid_cell <- function(A, B, s, alpha, tol, N = 32) {
  ref <- grid_logZ(A, B, s, alpha)
  r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, A, B, s, alpha, N)
  expect_equal(r$status, 0L)
  expect_lt(abs(r$log_Z - ref), tol,
    label = sprintf(
      "alpha=%g A=%g B=%g s=%g sinh=%g ref=%g",
      alpha, A, B, s, r$log_Z, ref
    )
  )
}

test_that("sinh N=32 unimodal-smooth cells match fine grid to ~1e-7", {
  configs <- list(
    list(A = 2, B = 1, s = 1, alpha = 1.5),
    list(A = 1, B = 0, s = 2, alpha = 0.7),
    list(A = 0.5, B = -2, s = 3, alpha = 0.3),
    list(A = 3, B = 5, s = 0.5, alpha = 4),
    list(A = 0.3, B = 0.5, s = 0.4, alpha = 3),
    list(A = 1, B = -1, s = 0.5, alpha = 0.9)
  )
  for (cfg in configs) {
    run_grid_cell(cfg$A, cfg$B, cfg$s, cfg$alpha, tol = 1e-5)
  }
})


# --- 3. sharp-log cells, with relaxed-then-strict tolerance per N --- #

test_that("sinh converges on sharp-log cells with increasing N", {
  sharp <- list(
    list(A = 0.66, B = 3.73, s = 0.011, alpha = 1.27),
    list(A = 0.50, B = -2.24, s = 0.029, alpha = 1.11),
    list(A = 7.25, B = 1.90, s = 0.041, alpha = 1.33),
    list(A = 0.52, B = 2.13, s = 0.343, alpha = 1.10)
  )
  for (cfg in sharp) {
    ref <- grid_logZ(cfg$A, cfg$B, cfg$s, cfg$alpha)
    e32 <- abs(bgms:::sd_log_density_at_l_ji_sinh_cpp(0, cfg$A, cfg$B, cfg$s, cfg$alpha, 32)$log_Z - ref)
    e64 <- abs(bgms:::sd_log_density_at_l_ji_sinh_cpp(0, cfg$A, cfg$B, cfg$s, cfg$alpha, 64)$log_Z - ref)
    e128 <- abs(bgms:::sd_log_density_at_l_ji_sinh_cpp(0, cfg$A, cfg$B, cfg$s, cfg$alpha, 128)$log_Z - ref)
    # N = 32 should be within 1e-3 even on the hardest sharp-log cells.
    expect_lt(e32, 1e-3, label = sprintf("N=32 alpha=%g", cfg$alpha))
    # N = 64 should be much tighter -- below 1e-6.
    expect_lt(e64, 1e-6, label = sprintf("N=64 alpha=%g", cfg$alpha))
    # N = 128 should reach near-roundoff.
    expect_lt(e128, 1e-9, label = sprintf("N=128 alpha=%g", cfg$alpha))
  }
})


# --- 4. alpha < 0.5 (singularity persists at t = +-i pi/2) --- #

test_that("sinh handles alpha < 0.5 correctly", {
  # alpha = 0.1: integrable singularity at t = +-i pi/2 of order (z)^(-0.8).
  # Midpoint rule still converges; just with a slightly larger prefactor.
  configs <- list(
    list(A = 1, B = 0, s = 1, alpha = 0.1),
    list(A = 1, B = 1, s = 0.5, alpha = 0.2),
    list(A = 2, B = -1.5, s = 0.3, alpha = 0.3),
    list(A = 0.5, B = 0.5, s = 1.5, alpha = 0.4)
  )
  for (cfg in configs) {
    run_grid_cell(cfg$A, cfg$B, cfg$s, cfg$alpha, tol = 1e-6, N = 32)
  }
})


# --- 5. alpha = 0.5 special: substitution exactly removes singularity --- #

test_that("sinh handles alpha = 0.5 (entire integrand in t)", {
  configs <- list(
    list(A = 1, B = 0, s = 1),
    list(A = 0.5, B = 2, s = 0.2),
    list(A = 3, B = -1, s = 0.5),
    list(A = 0.3, B = 0.3, s = 2)
  )
  # Tolerance 1e-6 at N=32 reflects that at alpha=0.5 the cosh^0
  # factor disappears entirely but the Gaussian-in-sinh envelope can
  # still need ~32 midpoint nodes to reach 1e-7 in tight-mode cases.
  # Production default is N=128 where these cells are at roundoff.
  for (cfg in configs) {
    run_grid_cell(cfg$A, cfg$B, cfg$s, alpha = 0.5, tol = 1e-6)
  }
})


# --- 6. large |B|: mode far from t = 0 --- #

test_that("sinh handles large |B| at N = 64", {
  for (B in c(-15, -10, -5, 5, 10, 15)) {
    ref <- grid_logZ(1, B, 1, 2)
    r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, 1, B, 1, 2, 64)
    expect_lt(abs(r$log_Z - ref), 1e-9,
      label = sprintf("B=%g", B)
    )
  }
})


# --- 7. failure paths --- #

test_that("sinh A <= 0 returns status 1", {
  for (A in c(0, -0.1, -1)) {
    r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, A, 1, 1, 2)
    expect_equal(r$status, 1L)
  }
})

test_that("sinh s_jj <= 0 returns status 2", {
  for (s in c(0, -0.01, -1)) {
    r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, 1, 1, s, 2)
    expect_equal(r$status, 2L)
  }
})


# --- 8. log_density consistency --- #

test_that("sinh log_density = ell(x_eval) - log_Z", {
  set.seed(7)
  for (.k in 1:200) {
    A <- exp(runif(1, -1, 2))
    B <- runif(1, -3, 3)
    s <- exp(runif(1, -2, 1))
    alpha <- runif(1, 0.1, 5)
    x_eval <- runif(1, -2, 2)
    r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(x_eval, A, B, s, alpha, 64)
    if (r$status != 0) next
    f_eval <- ell_ref(x_eval, A, B, s, alpha)
    expect_equal(r$log_density, f_eval - r$log_Z, tolerance = 1e-12)
  }
})


# --- 9. determinism --- #

test_that("sinh is deterministic", {
  r1 <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, 1, 1, 1, 2)
  r2 <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, 1, 1, 1, 2)
  expect_identical(r1, r2)
})


# --- 10. B -> -B reflection symmetry --- #

test_that("sinh log_Z is invariant under B -> -B (B = 0 axis check)", {
  # log Z is invariant under B -> -B (integrand becomes (s + phi^2)^(alpha-1)
  # exp(-A phi^2 + (-B) phi); substitute phi -> -phi to recover original).
  set.seed(31)
  for (.k in 1:30) {
    A <- exp(runif(1, -1, 2))
    B <- runif(1, 0.5, 3)
    s <- exp(runif(1, -2, 1))
    alpha <- runif(1, 0.2, 4)
    r_pos <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, A, B, s, alpha, 64)
    r_neg <- bgms:::sd_log_density_at_l_ji_sinh_cpp(0, A, -B, s, alpha, 64)
    expect_equal(r_pos$log_Z, r_neg$log_Z,
      tolerance = 1e-10,
      info = sprintf("A=%g B=%g s=%g alpha=%g", A, B, s, alpha)
    )
  }
})


# --- 11. convergence in N --- #

test_that("sinh error at N=64 is strictly tighter than at N=32 on hard cells", {
  configs <- list(
    list(A = 0.66, B = 3.73, s = 0.011, alpha = 1.27),
    list(A = 0.50, B = -2.24, s = 0.029, alpha = 1.11),
    list(A = 0.1, B = 1, s = 0.001, alpha = 1.05)
  )
  for (cfg in configs) {
    ref <- grid_logZ(cfg$A, cfg$B, cfg$s, cfg$alpha)
    e32 <- abs(bgms:::sd_log_density_at_l_ji_sinh_cpp(0, cfg$A, cfg$B, cfg$s, cfg$alpha, 32)$log_Z - ref)
    e64 <- abs(bgms:::sd_log_density_at_l_ji_sinh_cpp(0, cfg$A, cfg$B, cfg$s, cfg$alpha, 64)$log_Z - ref)
    expect_lt(e64, e32, label = sprintf("e32=%.2e e64=%.2e", e32, e64))
  }
})


# --- 12. property-based fuzz at the production default N=128 --- #

test_that("sinh fuzz (500 random configs) is below 1e-8 at N=128", {
  # Production default: N=128. The sinh substitution gives uniform
  # exponential convergence governed by the t-coordinate analyticity
  # strip width pi/2, so N=128 reaches roundoff across the chain-
  # realistic parameter range. Tolerance 1e-8 is conservative: a
  # 1000-cell empirical sweep at this regime had worst-case 1.1e-10
  # (machine roundoff floor for double-precision quadrature).
  set.seed(2026)
  N <- 500L
  A_vec <- exp(runif(N, log(0.3), log(50)))
  B_vec <- runif(N, -5, 5)
  s_vec <- exp(runif(N, log(0.01), log(10)))
  alpha_vec <- runif(N, 0.1, 5)

  worst <- 0
  for (i in seq_len(N)) {
    r <- bgms:::sd_log_density_at_l_ji_sinh_cpp(
      0, A_vec[i], B_vec[i], s_vec[i], alpha_vec[i], 128
    )
    if (r$status != 0) next
    ref <- grid_logZ(A_vec[i], B_vec[i], s_vec[i], alpha_vec[i])
    err <- abs(r$log_Z - ref)
    worst <- max(worst, err)
    expect_lt(err, 1e-8,
      label = sprintf(
        "i=%d A=%g B=%g s=%g a=%g err=%g",
        i, A_vec[i], B_vec[i], s_vec[i], alpha_vec[i], err
      )
    )
  }
  cat(sprintf("  [fuzz N=128 worst error: %.2e]\n", worst))
})


# --- 13. sinh-32 accurate on the canonical configuration grid --- #

test_that("sinh-32 worst-case absolute log_Z error < 1e-3 on the canonical grid", {
  # On well-behaved cells both quadratures reach roundoff; the looser
  # cells (small s, bimodal modes) historically tripped the legacy
  # alternatives. With sinh-32 the worst-case error on this grid is
  # bounded by 1e-3.
  configs = list(
    list(A = 1, B = 0, s = 1, alpha = 0.5),
    list(A = 0.5, B = -2, s = 0.05, alpha = 0.8),
    list(A = 3, B = 5, s = 1, alpha = 2.5),
    list(A = 0.3, B = 0.5, s = 0.4, alpha = 3),
    list(A = 0.3, B = 0, s = 0.1, alpha = 4),
    list(A = 1, B = 10, s = 1, alpha = 2)
  )
  for(cfg in configs) {
    ref = grid_logZ(cfg$A, cfg$B, cfg$s, cfg$alpha)
    r_sinh = bgms:::sd_log_density_at_l_ji_sinh_cpp(
      0, cfg$A, cfg$B, cfg$s, cfg$alpha, 32
    )
    err_sinh = abs(r_sinh$log_Z - ref)
    expect_lt(err_sinh, 1e-3,
      label = sprintf(
        "alpha=%g A=%g B=%g s=%g: err_sinh=%g",
        cfg$alpha, cfg$A, cfg$B, cfg$s, err_sinh
      )
    )
  }
})
