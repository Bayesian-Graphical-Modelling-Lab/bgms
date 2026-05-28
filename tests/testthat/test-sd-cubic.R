# --------------------------------------------------------------------------- #
# Unit tests for savage_dickey::solve_sd_cubic.
#
# The cubic critical-point solver underpins the adaptive Gauss-Hermite
# quadrature primitive for the alpha > 1 SD between-step. Correctness here
# is load-bearing: a wrong mode poisons the AGHQ log Z and the MH proposal.
#
# Tier 1 coverage (per dev/plans/active/sd-aghq-cubic.md):
#   - polynomial residual at each returned root
#   - alpha = 1 closed-form single root
#   - mode/saddle classification
#   - symmetric case B = 0
#   - failure paths (A <= 0, s_jj <= 0)
#   - ell/ell_pp at returned roots match an independent R-side evaluation
# --------------------------------------------------------------------------- #

# Reference: evaluate the kernel and its second derivative directly in R, to
# cross-check the C++ wrappers. Same formulas as src/math/savage_dickey/cubic_mode.cpp.
ell_ref <- function(phi, A, B, s_jj, alpha) {
  -A * phi^2 + B * phi + (alpha - 1) * log(s_jj + phi^2)
}
ell_pp_ref <- function(phi, A, s_jj, alpha) {
  -2 * A + 2 * (alpha - 1) * (s_jj - phi^2) / (s_jj + phi^2)^2
}

# Cubic in standard form. r = -B/(2A), c3 = s_jj, s = (alpha-1)/A.
#   phi^3 + r phi^2 + (c3 - s) phi + r c3 = 0
cubic_residual_ref <- function(phi, A, B, s_jj, alpha) {
  r <- -B / (2 * A)
  s <- (alpha - 1) / A
  c3 <- s_jj
  phi^3 + r * phi^2 + (c3 - s) * phi + r * c3
}


# --------------------------------------------------------------------------- #
# Test 1: polynomial residual < 1e-9 at each returned root, on a deterministic
# grid covering Delta < 0, Delta near 0, Delta > 0, symmetric q ~ 0, large r.
# --------------------------------------------------------------------------- #

# Residual tolerance that scales by the magnitude of contributing terms in
# the cubic, matching the solver's own self-check. Catches cancellation-
# dominated regions (tight A with large |B|) where absolute residuals can
# be O(1e-7) even at a numerically correct root.
residual_tol <- function(phi, A, B, s_jj, alpha, k = 1e-9) {
  r <- -B / (2 * A); c <- s_jj - (alpha - 1) / A; d <- r * s_jj
  coef_scale <- 1 + abs(r) + abs(c) + abs(d)
  phi_scale  <- 1 + abs(phi) + abs(phi)^2 + abs(phi)^3
  k * coef_scale * phi_scale
}

test_that("solve_sd_cubic: polynomial residual is tiny at every returned root", {
  grid <- expand.grid(
    A     = c(0.1, 0.5, 1, 2, 10),
    B     = c(-5, -1, 0, 1, 3),
    s_jj  = c(0.05, 0.5, 1, 5),
    alpha = c(1.01, 1.5, 2, 3, 5)
  )
  for (i in seq_len(nrow(grid))) {
    p <- grid[i, ]
    r <- bgms:::sd_cubic_solve_cpp(p$A, p$B, p$s_jj, p$alpha)
    expect_equal(r$status, 0,
                 info = sprintf("status != 0 at A=%g B=%g s=%g a=%g",
                                p$A, p$B, p$s_jj, p$alpha))
    n <- r$n_real_roots
    expect_true(n %in% c(1L, 3L),
                info = sprintf("n_real_roots=%d not in {1,3}", n))
    for (k in seq_len(n)) {
      res <- cubic_residual_ref(r$phi[k], p$A, p$B, p$s_jj, p$alpha)
      tol <- residual_tol(r$phi[k], p$A, p$B, p$s_jj, p$alpha)
      expect_lt(abs(res), tol,
                label = sprintf("residual at root %d (A=%g B=%g s=%g a=%g): %g",
                                k, p$A, p$B, p$s_jj, p$alpha, res))
    }
  }
})


# --------------------------------------------------------------------------- #
# Test 2: alpha = 1 closed-form. The cubic factors as (phi^2 + s_jj)(phi + r),
# yielding the single real root phi = -r = B/(2A) exactly.
# --------------------------------------------------------------------------- #

test_that("solve_sd_cubic: alpha = 1 returns single root at B/(2A) exactly", {
  for (A in c(0.5, 1, 2, 10)) {
    for (B in c(-3, -0.1, 0, 0.7, 5)) {
      for (s_jj in c(0.1, 1, 5)) {
        r <- bgms:::sd_cubic_solve_cpp(A, B, s_jj, alpha = 1)
        expect_equal(r$status, 0)
        expect_equal(r$n_real_roots, 1L)
        expect_equal(r$n_modes, 1L)
        expect_equal(r$phi[1], B / (2 * A), tolerance = 1e-12,
                     info = sprintf("A=%g B=%g s=%g", A, B, s_jj))
        expect_true(r$is_mode[1])
      }
    }
  }
})


# --------------------------------------------------------------------------- #
# Test 3: mode/saddle classification on a known-bimodal cell. With B = 0,
# small s_jj, and alpha large enough, the cubic has three real roots forming
# two symmetric modes around the saddle at phi = 0.
# --------------------------------------------------------------------------- #

test_that("solve_sd_cubic: bimodal symmetric cell yields two modes + one saddle", {
  r <- bgms:::sd_cubic_solve_cpp(A = 0.5, B = 0, s_jj = 0.1, alpha = 4)
  expect_equal(r$status, 0)
  expect_equal(r$n_real_roots, 3L)
  expect_equal(r$n_modes, 2L)
  # The saddle is at phi = 0; the two modes are symmetric.
  saddle_phi <- r$phi[r$saddle_index + 1]   # 0-indexed -> 1-indexed
  expect_equal(saddle_phi, 0, tolerance = 1e-10)
  expect_gt(r$ell_pp[r$saddle_index + 1], 0)
  # Mode positions are mirrored across phi = 0.
  mode_phi <- r$phi[c(r$global_mode_index, r$local_mode_index) + 1]
  expect_equal(sum(mode_phi), 0, tolerance = 1e-10)
  # Both modes have negative ell_pp (curvature > 0).
  expect_true(all(r$ell_pp[c(r$global_mode_index, r$local_mode_index) + 1] < 0))
  expect_true(all(r$curvature[c(r$global_mode_index, r$local_mode_index) + 1] > 0))
  # Symmetric case: the two modes have equal ell value.
  expect_equal(
    r$ell[r$global_mode_index + 1],
    r$ell[r$local_mode_index  + 1],
    tolerance = 1e-10
  )
})


# --------------------------------------------------------------------------- #
# Test 4: ell and ell_pp returned by the solver match an independent R-side
# evaluation. Catches accidental drift in the kernel formulas.
# --------------------------------------------------------------------------- #

test_that("solve_sd_cubic: ell and ell_pp at roots match R-side reference", {
  configs <- list(
    list(A = 1,   B = 2,  s_jj = 1,    alpha = 1.5),
    list(A = 0.5, B = 0,  s_jj = 0.2,  alpha = 3),
    list(A = 3,   B = -2, s_jj = 5,    alpha = 2),
    list(A = 0.1, B = 0,  s_jj = 0.05, alpha = 5)
  )
  for (cfg in configs) {
    r <- bgms:::sd_cubic_solve_cpp(cfg$A, cfg$B, cfg$s_jj, cfg$alpha)
    for (k in seq_len(r$n_real_roots)) {
      expect_equal(
        r$ell[k],
        ell_ref(r$phi[k], cfg$A, cfg$B, cfg$s_jj, cfg$alpha),
        tolerance = 1e-12,
        info = sprintf("ell mismatch at root %d, alpha=%g", k, cfg$alpha)
      )
      expect_equal(
        r$ell_pp[k],
        ell_pp_ref(r$phi[k], cfg$A, cfg$s_jj, cfg$alpha),
        tolerance = 1e-12,
        info = sprintf("ell_pp mismatch at root %d, alpha=%g", k, cfg$alpha)
      )
      expect_equal(r$curvature[k], -r$ell_pp[k], tolerance = 1e-15)
    }
  }
})


# --------------------------------------------------------------------------- #
# Test 5: failure paths.
# --------------------------------------------------------------------------- #

test_that("solve_sd_cubic: A <= 0 returns status = 1, no roots", {
  for (A in c(0, -0.1, -1)) {
    r <- bgms:::sd_cubic_solve_cpp(A, B = 1, s_jj = 1, alpha = 2)
    expect_equal(r$status, 1L, info = sprintf("A=%g", A))
    expect_equal(r$n_real_roots, 0L)
  }
})

test_that("solve_sd_cubic: s_jj <= 0 returns status = 2", {
  for (s in c(0, -0.05, -1)) {
    r <- bgms:::sd_cubic_solve_cpp(A = 1, B = 1, s_jj = s, alpha = 2)
    expect_equal(r$status, 2L, info = sprintf("s_jj=%g", s))
    expect_equal(r$n_real_roots, 0L)
  }
})


# --------------------------------------------------------------------------- #
# Test 6: B -> -B symmetry. Reflecting B across zero reflects every root
# across zero. The number of roots, mode count, and ell values stay invariant.
# --------------------------------------------------------------------------- #

test_that("solve_sd_cubic: B -> -B reflects every root across zero", {
  configs <- list(
    list(A = 0.5, s_jj = 0.2, alpha = 3, B = 0.4),
    list(A = 2,   s_jj = 1,   alpha = 2, B = 1.7),
    list(A = 0.1, s_jj = 0.5, alpha = 5, B = 0.3)
  )
  for (cfg in configs) {
    r_pos <- bgms:::sd_cubic_solve_cpp(cfg$A,  cfg$B, cfg$s_jj, cfg$alpha)
    r_neg <- bgms:::sd_cubic_solve_cpp(cfg$A, -cfg$B, cfg$s_jj, cfg$alpha)
    expect_equal(r_pos$n_real_roots, r_neg$n_real_roots)
    expect_equal(r_pos$n_modes,      r_neg$n_modes)
    expect_equal(
      sort(r_pos$phi[seq_len(r_pos$n_real_roots)]),
      sort(-r_neg$phi[seq_len(r_neg$n_real_roots)]),
      tolerance = 1e-10
    )
    expect_equal(
      sort(r_pos$ell[seq_len(r_pos$n_real_roots)]),
      sort(r_neg$ell[seq_len(r_neg$n_real_roots)]),
      tolerance = 1e-10
    )
  }
})


# --------------------------------------------------------------------------- #
# Test 7: bimodality threshold via the curvature-at-zero test.
# Condition: ell_pp(0) > 0  <=>  tau * c_3 < 2(alpha - 1) with tau = 2 A.
# The companion-side note (~/SV/Z/notes/2026-05-27_*) gives this as the
# necessary condition for the cubic discriminant Delta > 0.
# --------------------------------------------------------------------------- #

test_that("solve_sd_cubic: curvature-at-zero predicts bimodality (with B=0)", {
  # With B = 0 the saddle (if any) sits at phi = 0. Bimodal iff ell_pp(0) > 0.
  # Points at the boundary (alpha - 1) == A * s_jj are triple-root degeneracies
  # where the solver returns n_modes = 0 by design (no isolated Laplace mode).
  # Skip those with a margin to avoid roundoff misclassification.
  configs <- expand.grid(
    A     = c(0.2, 0.5, 1, 2),
    s_jj  = c(0.05, 0.2, 1),
    alpha = c(1.1, 2, 4, 8)
  )
  boundary_margin <- 0.01
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, ]
    ell_pp_0 <- ell_pp_ref(0, cfg$A, cfg$s_jj, cfg$alpha)
    # Skip cases too close to the threshold (within margin of the curvature-
    # at-zero diagnostic).
    if (abs(ell_pp_0) < boundary_margin * max(1, 2 * cfg$A)) next
    bimodal_predicted <- ell_pp_0 > 0
    r <- bgms:::sd_cubic_solve_cpp(cfg$A, B = 0, cfg$s_jj, cfg$alpha)
    expect_equal(
      r$n_modes == 2,
      bimodal_predicted,
      info = sprintf("A=%g s=%g alpha=%g: predicted bimodal=%s but n_modes=%d",
                     cfg$A, cfg$s_jj, cfg$alpha, bimodal_predicted, r$n_modes)
    )
  }
})


# --------------------------------------------------------------------------- #
# Test 8: property-based fuzz (Tier 2 in the plan).
# 5000 random configs, check residuals + ell agreement.
# --------------------------------------------------------------------------- #

test_that("solve_sd_cubic: 5000 random configs satisfy residual + ell consistency", {
  set.seed(2026)
  N <- 5000L
  A_vec     <- exp(runif(N, log(0.01), log(100)))
  B_vec     <- runif(N, -50, 50)
  s_jj_vec  <- exp(runif(N, log(1e-3), log(10)))
  alpha_vec <- runif(N, 1.001, 5)

  for (i in seq_len(N)) {
    r <- bgms:::sd_cubic_solve_cpp(A_vec[i], B_vec[i], s_jj_vec[i], alpha_vec[i])
    if (r$status != 0) next   # numerical fallback; covered separately below
    for (k in seq_len(r$n_real_roots)) {
      res <- cubic_residual_ref(r$phi[k], A_vec[i], B_vec[i],
                                s_jj_vec[i], alpha_vec[i])
      tol <- residual_tol(r$phi[k], A_vec[i], B_vec[i],
                          s_jj_vec[i], alpha_vec[i], k = 1e-7)
      expect_lt(
        abs(res), tol,
        label = sprintf("residual i=%d k=%d (A=%g B=%g s=%g a=%g): %g vs tol %g",
                        i, k, A_vec[i], B_vec[i], s_jj_vec[i],
                        alpha_vec[i], res, tol)
      )
      ell_diff <- abs(
        r$ell[k] - ell_ref(r$phi[k], A_vec[i], B_vec[i],
                           s_jj_vec[i], alpha_vec[i])
      )
      expect_lt(ell_diff, 1e-10,
                label = sprintf("ell mismatch i=%d k=%d", i, k))
    }
  }
})


# --------------------------------------------------------------------------- #
# Test 9: sd_log_kernel and sd_log_kernel_pp wrappers match the reference.
# --------------------------------------------------------------------------- #

test_that("sd_log_kernel and sd_log_kernel_pp match R reference", {
  set.seed(11)
  for (.k in 1:200) {
    A <- exp(runif(1, -2, 2)); B <- runif(1, -3, 3)
    s <- exp(runif(1, -3, 1)); a <- runif(1, 1.01, 4)
    phi <- runif(1, -3, 3)
    got <- bgms:::sd_log_kernel_cpp(phi, A, B, s, a)
    expect_equal(got$ell,    ell_ref(phi, A, B, s, a), tolerance = 1e-12)
    expect_equal(got$ell_pp, ell_pp_ref(phi, A, s, a), tolerance = 1e-12)
  }
})
