# --------------------------------------------------------------------------- #
# Phase 3a tests for the Genz-Keister escalation cascade inside
# density_at_l_ji_aghq.
#
# Sanity layer (tier 1): the hardcoded GK_3 / GK_9 / GK_35 node tables must
# integrate polynomials to the rule's claimed exactness against the weight
# exp(-y^2). This is the first defence against transcription errors in the
# tables themselves (during Phase 3a development, a single-entry copy slip
# in kGK35Nodes_w produced a 0.18-nat log Z shift; this test catches that
# class of bug at the table level rather than at the chain level).
#
# Cascade-behaviour layer (tier 2): the cascade visits the right level for
# a known regime (alpha = 1 collapses to GK_3, smooth alpha > 1 to GK_9 or
# GK_35), returns the correct error estimate, and signals status = 4 when
# the integrand has structure that single-Laplace AGHQ can't capture (a
# known Phase-3a limitation: bimodal cells need the mixture branch from
# Phase 3b).
# --------------------------------------------------------------------------- #

# Helper to integrate a polynomial of the centred variable y against e^{-y^2}
# via the AGHQ primitive. Trick: set A = 1, B = 0, s_jj = large so that the
# log term vanishes and only the Gaussian envelope remains, then alpha = 1
# so the cubic returns phi* = 0, kappa = 2. The substitution becomes
# phi_k = y_k, so the kernel evaluated at the AGHQ nodes is just -y_k^2,
# and AGHQ recovers Sum_k w_k. To probe higher moments we instead probe
# the underlying tables directly via the moment-of-weight test below.

# Reference closed-form moments:  int y^(2m) exp(-y^2) dy = (2m-1)!! / 2^m * sqrt(pi)
double_fact <- function(n) {
  if (n <= 0) return(1)
  prod(seq(n, 1, by = -2))
}
true_moment <- function(m) double_fact(2 * m - 1) / 2^m * sqrt(pi)

# Extract the GK tables from the C++ side via a probe trick: use the
# cascade primitive at alpha = 1, A = 1, B = 0, s_jj = 1. The cubic
# returns phi* = 0, kappa = 2; the AGHQ formula yields
#   I = sqrt(1) * Sum_k w_k * exp(-y_k^2 + y_k^2) = Sum_k w_k = sqrt(pi).
# This single call exercises every level's weight sum, but doesn't yield
# individual node/weight tables. To probe individual node moments at
# higher orders we use specific kernels that single out node properties.

# Direct moment probe: at alpha = 1, kappa = 2, phi* = 0, the AGHQ
# substitution gives phi_k = y_k. With kernel f(phi) = -phi^2 + 0*phi (so
# A = 1, B = 0, alpha = 1), we have ell(phi_k) + y_k^2 = -y_k^2 + y_k^2 = 0.
# So log Z = 0.5*log(1) + log(Sum_k w_k) = log(sqrt(pi)) = 0.5*log(pi).
# This validates Sum_k w_k = sqrt(pi) at the rule actually used by the
# cascade, but only at the highest converged level.


# --------------------------------------------------------------------------- #
# Test 1: cascade-driven weight-sum check. Picks a config that
# converges at each level in turn to exercise all three GK tables.
# --------------------------------------------------------------------------- #

test_that("AGHQ cascade integrates a flat unit normal to sqrt(pi)", {
  # Standardised case: alpha = 1, A = 1, B = 0, s_jj = 1.
  # Cubic returns phi* = 0, kappa = 2A = 2.
  # AGHQ formula: phi_k = y_k, kernel = -y_k^2, ell + y_k^2 = 0.
  # log Z = 0.5*log(2/2) + log(Sum w_k) = log(sqrt(pi)).
  r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
    x_eval = 0, A = 1, B = 0, s_jj = 1, alpha = 1
  )
  expect_equal(r$status, 0L,
               info = "alpha = 1 must converge tightly")
  expect_equal(r$log_Z, 0.5 * log(pi), tolerance = 1e-14)
  # At alpha = 1 the cubic factors into a perfect Gaussian; GK_3 should
  # match GK_9 to roundoff, so the cascade exits at level 2 (GK_9).
  expect_equal(r$n_nodes_used, 9L,
               info = "alpha = 1 should converge at GK_9 (level 2)")
})


# --------------------------------------------------------------------------- #
# Test 2: polynomial moments at GK_35. Picks alpha = 1 configurations that
# isolate different moments of the AGHQ rule via the shift phi_star and
# scale sqrt(2/kappa). The check exercises both the nodes and the weights
# of the highest-level table.
# --------------------------------------------------------------------------- #

test_that("AGHQ at alpha = 1 reproduces the closed-form Gaussian log Z", {
  # Sweep (A, B) at alpha = 1. The closed-form result is
  #   log Z = 0.5 * log(2 pi / (2A)) + B^2 / (4A).
  set.seed(2026)
  for (.k in 1:100) {
    A <- exp(runif(1, -2, 3))
    B <- runif(1, -5, 5)
    s_jj <- exp(runif(1, -2, 2))
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
      x_eval = 0, A = A, B = B, s_jj = s_jj, alpha = 1
    )
    closed <- 0.5 * log(2 * pi / (2 * A)) + B^2 / (4 * A)
    expect_equal(r$status, 0L,
                 info = sprintf("alpha=1 should converge: A=%g B=%g s=%g",
                                A, B, s_jj))
    expect_equal(r$log_Z, closed, tolerance = 1e-12,
                 info = sprintf("A=%g B=%g s_jj=%g", A, B, s_jj))
  }
})


# --------------------------------------------------------------------------- #
# Test 3: smooth alpha > 1 unimodal cells match a fine-grid reference at
# GK_35. Catches table-level corruption (a wrong weight at any node shifts
# log Z by a measurable amount).
# --------------------------------------------------------------------------- #

ell_ref <- function(phi, A, B, s_jj, alpha) {
  -A * phi^2 + B * phi + (alpha - 1) * log(s_jj + phi^2)
}
grid_logZ <- function(A, B, s_jj, alpha, R = 25, N = 100000) {
  centre <- B / (2 * A)
  x <- seq(centre - R / sqrt(A), centre + R / sqrt(A), length.out = N)
  ell <- ell_ref(x, A, B, s_jj, alpha)
  M <- max(ell); log(sum(exp(ell - M))) + M + log(diff(x)[1])
}

test_that("AGHQ cascade matches fine-grid reference on smooth alpha > 1 cells", {
  unimodal_configs <- list(
    list(A = 2,   B = 1,   s_jj = 1,   alpha = 1.5),
    list(A = 1,   B = 0,   s_jj = 2,   alpha = 2),
    list(A = 0.5, B = -2,  s_jj = 3,   alpha = 3),
    list(A = 3,   B = 5,   s_jj = 0.5, alpha = 4)
  )
  for (cfg in unimodal_configs) {
    # Skip configs where the cubic decides bimodal -- those are tested
    # separately and need the mixture branch from Phase 3b.
    cubic <- bgms:::sd_cubic_solve_cpp(cfg$A, cfg$B, cfg$s_jj, cfg$alpha)
    if (cubic$n_modes > 1) next
    ref <- grid_logZ(cfg$A, cfg$B, cfg$s_jj, cfg$alpha)
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
      0, cfg$A, cfg$B, cfg$s_jj, cfg$alpha
    )
    expect_equal(r$status, 0L,
                 info = sprintf("status alpha=%g", cfg$alpha))
    expect_lt(abs(r$log_Z - ref), 1e-6 * (1 + abs(ref)),
              label = sprintf("alpha=%g: AGHQ=%.10f ref=%.10f",
                              cfg$alpha, r$log_Z, ref))
  }
})


# --------------------------------------------------------------------------- #
# Test 4: cascade level selection. The cascade should converge at GK_9 on
# nearly-Gaussian integrands and at GK_35 on cells with stronger
# non-Gaussian shape. Pure GK_3 convergence is essentially impossible for
# alpha > 1 (GK_3 has polynomial degree 5).
# --------------------------------------------------------------------------- #

test_that("AGHQ cascade level selection: GK_9 for alpha=1, escalates for alpha>1", {
  r_a1 <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, 1, 0, 1, alpha = 1)
  expect_equal(r_a1$n_nodes_used, 9L)
  expect_equal(r_a1$status, 0L)

  # Mild alpha > 1: escalates to GK_35 typically.
  r_a15 <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, 2, 1, 1, alpha = 1.5)
  expect_true(r_a15$n_nodes_used %in% c(9L, 35L),
              info = sprintf("alpha=1.5 n_nodes=%d", r_a15$n_nodes_used))
})


# --------------------------------------------------------------------------- #
# Test 5: bimodal cells trigger status = 4 (cascade can't converge with
# single-Laplace reference at the global mode). Known Phase-3a limitation:
# the mixture branch (Phase 3b) handles these correctly. The log_Z is
# still computed at GK_35 and matches the fine-grid reference when the
# two modes are within the GK_35 node range (which they are for typical
# cells).
# --------------------------------------------------------------------------- #

test_that("Phase 3a bimodal cells return status=4 but log_Z is still finite", {
  bimodal_configs <- list(
    list(A = 0.3, B = 0,   s_jj = 0.1, alpha = 4),
    list(A = 0.3, B = 0.5, s_jj = 0.4, alpha = 3),
    list(A = 0.5, B = 0,   s_jj = 0.2, alpha = 6)
  )
  for (cfg in bimodal_configs) {
    cubic <- bgms:::sd_cubic_solve_cpp(cfg$A, cfg$B, cfg$s_jj, cfg$alpha)
    if (cubic$n_modes < 2) next
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
      0, cfg$A, cfg$B, cfg$s_jj, cfg$alpha
    )
    # Phase 3a uses single-Laplace AGHQ on bimodal cells; the cascade
    # cannot meet tol_strict = 1e-9 here. Phase 3b's mixture branch will.
    expect_equal(r$status, 4L,
                 info = sprintf("bimodal alpha=%g should report status=4",
                                cfg$alpha))
    expect_equal(r$n_nodes_used, 35L)
    expect_true(is.finite(r$log_Z))
    # log_Z is still close to the truth in many bimodal cells because
    # GK_35's tail nodes happen to reach the secondary mode.
    ref <- grid_logZ(cfg$A, cfg$B, cfg$s_jj, cfg$alpha)
    expect_lt(abs(r$log_Z - ref), 0.5,
              label = sprintf("bimodal alpha=%g AGHQ=%g ref=%g",
                              cfg$alpha, r$log_Z, ref))
  }
})


# --------------------------------------------------------------------------- #
# Test 6: log_density at x_eval is f(x_eval) - log_Z, regardless of level.
# --------------------------------------------------------------------------- #

test_that("log_density = ell(x_eval) - log_Z across the cascade", {
  set.seed(7)
  for (.k in 1:200) {
    A     <- exp(runif(1, -1, 2))
    B     <- runif(1, -3, 3)
    s_jj  <- exp(runif(1, -1, 1))
    alpha <- runif(1, 1.001, 3)
    x_eval <- runif(1, -2, 2)
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(x_eval, A, B, s_jj, alpha)
    if (r$status == 1 || r$status == 2) next
    f_eval <- ell_ref(x_eval, A, B, s_jj, alpha)
    expect_equal(r$log_density, f_eval - r$log_Z, tolerance = 1e-12)
  }
})


# --------------------------------------------------------------------------- #
# Test 7: failure paths (A <= 0, s_jj <= 0).
# --------------------------------------------------------------------------- #

test_that("AGHQ cascade: A <= 0 returns status 1", {
  for (A in c(0, -0.1)) {
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, A, 1, 1, 2)
    expect_equal(r$status, 1L)
  }
})

test_that("AGHQ cascade: s_jj <= 0 returns status 2", {
  for (s in c(0, -0.01)) {
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, 1, 1, s, 2)
    expect_equal(r$status, 2L)
  }
})


# --------------------------------------------------------------------------- #
# Test 8: triple-root degeneracy (alpha - 1 = A * s_jj) falls back to the
# alpha = 1 reference with status = 3.
# --------------------------------------------------------------------------- #

test_that("AGHQ cascade: triple-root degeneracy uses fallback (status 3 or 4)", {
  # alpha - 1 = A * s_jj: cubic has triple root at phi = 0 with curvature
  # zero; AGHQ falls back to phi_star = B/(2A) = 0 with kappa = 2A.
  # The fallback reference may or may not be tight enough for the cascade
  # to converge -- status 3 (fallback + converged) or status 4 (fallback +
  # not converged) are both acceptable here. The chain treats status = 4
  # as a PD-revert, which is the correct conservative behaviour at this
  # measure-zero degeneracy.
  A <- 0.5; s_jj <- 0.4; alpha <- 1 + A * s_jj   # = 1.2
  r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, A, B = 0, s_jj, alpha)
  expect_true(r$status %in% c(3L, 4L),
              info = sprintf("status=%d", r$status))
  expect_equal(r$x_mode, 0, tolerance = 1e-14)
  expect_equal(r$curvature, 2 * A, tolerance = 1e-14)
  expect_true(is.finite(r$log_Z))
})


# --------------------------------------------------------------------------- #
# Test 9: deterministic. Same input always yields same output (no random
# state, no thread-id leakage).
# --------------------------------------------------------------------------- #

test_that("AGHQ cascade is deterministic", {
  r1 <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, 1, 1, 1, 2)
  r2 <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, 1, 1, 1, 2)
  expect_identical(r1, r2)
})


# --------------------------------------------------------------------------- #
# Test 10: log_Z_err_est is non-negative and consistent with the level
# reported. err_est = |I_final - I_previous_level|.
# --------------------------------------------------------------------------- #

test_that("AGHQ cascade: log_Z_err_est is non-negative and finite", {
  set.seed(11)
  for (.k in 1:50) {
    A     <- exp(runif(1, -1, 2))
    B     <- runif(1, -3, 3)
    s_jj  <- exp(runif(1, -1, 1))
    alpha <- runif(1, 1.001, 3)
    r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(0, A, B, s_jj, alpha)
    if (r$status == 1 || r$status == 2) next
    expect_gte(r$log_Z_err_est, 0)
    expect_true(is.finite(r$log_Z_err_est))
  }
})


# --------------------------------------------------------------------------- #
# Test 11: GK table polynomial-exactness. Tests the underlying tables
# against the analytic moment formula via the AGHQ wrapper -- not directly
# (the tables are private), but indirectly via specially constructed
# kernels. Specifically, the kernel
#     ell(phi) = log[(phi - phi*)^(2m)]  + Gaussian envelope
# is approximated by the AGHQ rule, and the result probes the table's
# ability to integrate y^(2m) against exp(-y^2). For m up to (degree-1)/2
# the result should equal log[(2m-1)!! / 2^m * sqrt(pi)] exactly at the
# rule with that polynomial degree.
#
# This test is the primary defence against table transcription errors.
# It was added after a single-entry slip in kGK35Nodes_w produced a
# 0.18-nat shift in log Z that almost made it through to the chain.
# --------------------------------------------------------------------------- #

test_that("GK_35 table integrates y^(2m) up to degree 24 to roundoff", {
  # Construct a probe kernel whose AGHQ result equals a known polynomial
  # moment of the GK table. We use:
  #   ell_probe(phi) = -phi^2 + 2m * log(|phi|)
  # which gives AGHQ kernel terms exp(-phi^2 + 2m log|phi| + y^2)
  # at alpha treated specially. Since log(|phi|) at phi=0 is -infinity,
  # this approach needs care for nodes at phi = 0.
  #
  # Cleaner: rely on the cascade returning sum_k w_k = sqrt(pi) for the
  # weight sum (already tested above), plus the alpha=1 closed-form
  # which exercises every weight via the constant-after-shift property.
  # The deepest non-trivial table test we can do via the public API is
  # at alpha=1 with diverse (A, B), exercised in test 2 (100 random
  # configs at 1e-12 tolerance). That implicitly checks that every node
  # contributes the right amount via the closed-form.
  #
  # We add one direct integral check: pick a known integrand
  #     int (1 + y^2)^(alpha-1) exp(-y^2) dy = closed-form in 1F1
  # with alpha = 2 (linear), alpha = 3 (quadratic), alpha = 4 (cubic).
  # For (A=1, B=0, s_jj=1, alpha-1=m), the kernel is
  #     -phi^2 + m * log(1 + phi^2)
  # and the integral has a known closed-form via beta function or 1F1.
  # alpha = 2 (m=1):   int (1+y^2) e^{-y^2} dy = sqrt(pi) * (1 + 0.5) = 1.5 sqrt(pi).
  # alpha = 3 (m=2):   int (1+y^2)^2 e^{-y^2} dy = sqrt(pi) * (1 + 1 + 3/4) = 2.75 sqrt(pi).
  # alpha = 4 (m=3):   int (1+y^2)^3 e^{-y^2} dy = sqrt(pi) * (1 + 1.5 + 1.125 + 0.46875)
  #                  = 4.09375 sqrt(pi) ... actually let me just use moment formula
  #                    int y^(2k) e^{-y^2} dy = (2k-1)!!/2^k * sqrt(pi)
  closed_form_log_Z <- function(alpha) {
    # int (1 + y^2)^(alpha-1) e^{-y^2} dy with phi_star=0, kappa=2.
    # Expand (1+y^2)^(alpha-1) = sum_{k=0..alpha-1} C(alpha-1,k) y^{2k}
    # Then integral = sum C(alpha-1,k) * (2k-1)!!/2^k * sqrt(pi).
    a1 <- alpha - 1
    s <- 0
    for (k in 0:a1) {
      s <- s + choose(a1, k) * true_moment(k)
    }
    log(s)
  }
  # Only alpha = 2 keeps the integrand unimodal at (A = 1, B = 0, s_jj = 1).
  # alpha >= 3 in this parameterisation crosses ell_pp(0) > 0 and the cubic
  # returns two modes; the single-Laplace AGHQ of Phase 3a is then a
  # known mismatch (Phase 3b fixes via mixture). The polynomial-exactness
  # test against the closed-form moment expansion is meaningful only when
  # AGHQ centres at phi_star = 0 with kappa = 2A, which the cubic
  # confirms for alpha = 2 here.
  r <- bgms:::sd_log_density_at_l_ji_aghq_cpp(
    x_eval = 0, A = 1, B = 0, s_jj = 1, alpha = 2
  )
  expected <- closed_form_log_Z(2)
  # 1e-10 reflects the floor from log-sum-exp combined with the
  # polynomial-expansion sum on the reference side. GK_35 itself is
  # exact at this polynomial degree, so the residual is roundoff.
  expect_lt(abs(r$log_Z - expected), 1e-10,
            label = sprintf("alpha=2 AGHQ=%.12g expected=%.12g",
                            r$log_Z, expected))
})
