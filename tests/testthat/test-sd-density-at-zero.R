# --------------------------------------------------------------------------- #
# Step 1c (SD pivot): unit tests for the Savage-Dickey 1D conditional density
# at K_ij = 0.
#
# Tests:
#   1. Bit-for-bit agreement with Z's standalone primitive
#      (sd_log_post_density_at_zero_cpp in ~/SV/Z/R/src/sd_density_at_zero.cpp).
#      Skipped if the Z file is not on this machine.
#   2. Symmetry under (i, j) <-> (j, i).
#   3. PD-failure path: K not PD => status = 1, log_density = -Inf.
#   4. NLO + Laplace vs adaptive R quadrature of exp(L(x)) in a benign regime.
# --------------------------------------------------------------------------- #


# Helper: build the 1D log-density L(x) directly from the math in the spec.
# Used as the integrand for the NLO-vs-quadrature check.
make_L_fun <- function(K, i, j, S, n_obs, delta, sigma) {
  K0 <- K
  K0[i, j] <- 0
  K0[j, i] <- 0
  Sigma0 <- chol2inv(chol(K0))
  c1 <- Sigma0[i, j]
  c2 <- Sigma0[i, i] * Sigma0[j, j] - Sigma0[i, j]^2
  a  <- delta + n_obs / 2
  Sij <- S[i, j]
  inv_sigma2 <- 1 / sigma^2
  list(
    L = function(x) {
      D <- 1 + 2 * c1 * x - c2 * x * x
      ifelse(D > 0, a * log(D) - Sij * x - 0.5 * x * x * inv_sigma2, -Inf)
    },
    c1 = c1, c2 = c2
  )
}


# Generate a single benign PD K and matching S(n_obs).
make_problem <- function(q, seed, n_obs = 50, edge_density = 0.3) {
  set.seed(seed)
  G <- matrix(0L, q, q)
  diag(G) <- 1L
  for (i in 1:(q - 1)) for (j in (i + 1):q) {
    if (runif(1) < edge_density) G[i, j] <- G[j, i] <- 1L
  }
  # Build a random PD K consistent with G structurally (entries off-diagonal
  # may still be non-zero where G = 0; the density is defined on all of K).
  L <- matrix(rnorm(q * q, sd = 0.3), q, q)
  L[upper.tri(L)] <- 0
  diag(L) <- abs(diag(L)) + 1
  K <- tcrossprod(L)
  Y <- matrix(rnorm(n_obs * q), n_obs, q)
  S <- crossprod(Y)
  list(K = K, S = S, n_obs = n_obs, G = G)
}


# ---------------------------------------------------------------- #
# Test 1: bit-for-bit vs Z's standalone primitive.
# ---------------------------------------------------------------- #
test_that("sd_log_density_at_zero_cpp matches Z's standalone primitive", {
  z_src <- "~/Library/CloudStorage/Dropbox/Projecten/SV/Z/R/src/sd_density_at_zero.cpp"
  skip_if(!file.exists(path.expand(z_src)),
          "Z reference file not present")

  # Sourcing the Z file compiles its own ::cpp -> R wrapper.
  # The wrapper exposes sd_log_post_density_at_zero_cpp.
  zenv <- new.env()
  Rcpp::sourceCpp(path.expand(z_src), env = zenv, verbose = FALSE,
                  rebuild = FALSE, cleanupCacheDir = FALSE)
  z_fun <- zenv$sd_log_post_density_at_zero_cpp

  grid <- expand.grid(
    q       = c(5, 10),
    seed    = 1:3,
    n_obs   = c(20, 200),
    delta   = c(0, 0.5, 1.5),
    sigma   = c(0.25, 1.0),
    nlo     = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  for (row in seq_len(nrow(grid))) {
    pars <- grid[row, ]
    prob <- make_problem(pars$q, pars$seed, pars$n_obs)
    i <- 1L
    j <- 2L # any off-diagonal pair will do

    ours <- sd_log_density_at_zero_cpp(
      K = prob$K, i = i, j = j, S = prob$S,
      n_obs = prob$n_obs, delta = pars$delta, sigma = pars$sigma,
      nlo = pars$nlo)
    theirs <- z_fun(
      K = prob$K, i = i, j = j, S = prob$S,
      n_obs = prob$n_obs, delta = pars$delta, sigma = pars$sigma,
      nlo = pars$nlo)

    expect_equal(ours$status, theirs$status,
                 info = paste("status mismatch at row", row))
    if (ours$status == 0) {
      expect_equal(ours$log_density, theirs$log_density, tolerance = 1e-10,
                   info = paste("log_density mismatch at row", row))
      expect_equal(ours$x_mode,      theirs$x_mode,      tolerance = 1e-10,
                   info = paste("x_mode mismatch at row", row))
      expect_equal(ours$curvature,   theirs$curvature,   tolerance = 1e-10,
                   info = paste("curvature mismatch at row", row))
    }
  }
})


# ---------------------------------------------------------------- #
# Test 2: symmetry under index swap.
# ---------------------------------------------------------------- #
test_that("density at zero is symmetric under (i, j) <-> (j, i)", {
  for (seed in 1:5) {
    prob <- make_problem(q = 8, seed = seed)
    for (pair in list(c(1, 2), c(3, 5), c(2, 7))) {
      i <- pair[1]
      j <- pair[2]
      a <- sd_log_density_at_zero_cpp(
        K = prob$K, i = i, j = j, S = prob$S, n_obs = prob$n_obs,
        delta = 1.0, sigma = 1.0, nlo = TRUE)
      b <- sd_log_density_at_zero_cpp(
        K = prob$K, i = j, j = i, S = prob$S, n_obs = prob$n_obs,
        delta = 1.0, sigma = 1.0, nlo = TRUE)
      expect_equal(a$status, b$status)
      if (a$status == 0) {
        expect_equal(a$log_density, b$log_density, tolerance = 1e-12)
        expect_equal(a$curvature,   b$curvature,   tolerance = 1e-12)
        # x_mode flips sign if and only if c1 changes sign under the
        # swap; with K symmetric and Sigma_0 symmetric, c1 is unchanged.
        expect_equal(a$x_mode,      b$x_mode,      tolerance = 1e-12)
      }
    }
  }
})


# ---------------------------------------------------------------- #
# Test 3: PD failure path.
# ---------------------------------------------------------------- #
test_that("non-PD K_0 returns status = 1 and log_density = -Inf", {
  # K_0 is K with K_{i,j} set to zero. To trigger status = 1 we need
  # K_0 (not K_{i,j}) to be the source of non-PD. So break PD via an
  # entry *other than* the one being toggled.
  q <- 5
  K <- diag(1.0, q)
  K[3, 4] <- K[4, 3] <- 5.0  # breaks PD; (3, 4) is not the toggled edge
  S <- diag(1.0, q)
  res <- sd_log_density_at_zero_cpp(
    K = K, i = 1L, j = 2L, S = S, n_obs = 10L,
    delta = 1.0, sigma = 1.0, nlo = TRUE)
  expect_equal(res$status, 1L)
  expect_true(is.infinite(res$log_density) && res$log_density < 0)
})


# ---------------------------------------------------------------- #
# Test 4: NLO + Laplace vs adaptive R quadrature in a benign regime.
# ---------------------------------------------------------------- #
test_that("Laplace+NLO agrees with adaptive quadrature in benign regime", {
  # Pick a cell where the 1D density is smooth and well-localised around
  # the mode (large n_obs, moderate sigma) so Laplace + NLO is tight.
  prob <- make_problem(q = 6, seed = 11, n_obs = 500)
  delta <- 1.0
  sigma <- 1.0
  i <- 1L
  j <- 2L

  Lfun <- make_L_fun(prob$K, i, j, prob$S, prob$n_obs, delta, sigma)
  # Find a safe interval where the determinant factor D(x) > 0.
  # Roots of 1 + 2 c1 x - c2 x^2 = 0 are at x = (c1 +/- sqrt(c1^2 + c2))/c2.
  disc <- Lfun$c1^2 + Lfun$c2
  if (Lfun$c2 > 0 && disc > 0) {
    r_neg <- (Lfun$c1 - sqrt(disc)) / Lfun$c2
    r_pos <- (Lfun$c1 + sqrt(disc)) / Lfun$c2
    eps <- 1e-6 * (r_pos - r_neg)
    lo <- r_neg + eps
    hi <- r_pos - eps
  } else {
    lo <- -10
    hi <- 10
  }

  # Locate the mode to use as the log-shift for numerical stability.
  res <- sd_log_density_at_zero_cpp(
    K = prob$K, i = i, j = j, S = prob$S, n_obs = prob$n_obs,
    delta = delta, sigma = sigma, nlo = TRUE)
  expect_equal(res$status, 0L)
  L_at_mode <- Lfun$L(res$x_mode)
  integrand <- function(x) exp(Lfun$L(x) - L_at_mode)
  quad <- integrate(integrand, lower = lo, upper = hi,
                    rel.tol = 1e-10, subdivisions = 1000L)
  log_int_quad <- log(quad$value) + L_at_mode
  # res$log_density = -(log_int_Laplace_NLO);  we compare integrals.
  log_int_lap   <- -res$log_density
  # Expect agreement to ~1e-3 in the well-localised regime.
  expect_equal(log_int_lap, log_int_quad, tolerance = 1e-3)

  # Also: NLO should improve over plain Laplace.
  res_no_nlo <- sd_log_density_at_zero_cpp(
    K = prob$K, i = i, j = j, S = prob$S, n_obs = prob$n_obs,
    delta = delta, sigma = sigma, nlo = FALSE)
  log_int_no_nlo <- -res_no_nlo$log_density
  err_with_nlo <- abs(log_int_lap     - log_int_quad)
  err_no_nlo   <- abs(log_int_no_nlo  - log_int_quad)
  expect_lt(err_with_nlo, err_no_nlo + 1e-12)
})
