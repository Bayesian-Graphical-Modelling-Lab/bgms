# --------------------------------------------------------------------------- #
# Step (a3)-1: unit tests for cholesky_helpers::perm_to_trailing_2x2.
#
# The function brings positions (i, j) to the trailing 2×2 of the permuted
# lower-triangular Cholesky factor via Bunch-style adjacent symmetric swaps
# + Givens rotations, never going through Σ_BB = K⁻¹_BB. This avoids the
# 1/det amplification that troubled the prior Σ-corner extraction path.
#
# Tests:
#   1. Bit-for-bit (to roundoff) match against chol(P K Pᵀ) for random PD K
#      across p ∈ {3..20}, all edges (i, j), multiple seeds.
#   2. Singular K handling: function reports ok = FALSE without crashing.
# --------------------------------------------------------------------------- #

bunch_truth <- function(K, i, j) {
  p <- nrow(K)
  perm <- c(setdiff(seq_len(p), c(i, j)), i, j)
  Lp <- t(chol(K[perm, perm]))
  list(
    l_ii = Lp[p - 1, p - 1],
    l_ji = Lp[p,     p - 1],
    l_jj = Lp[p,     p]
  )
}

rand_pd <- function(p, seed) {
  set.seed(seed)
  A <- matrix(rnorm(p * p), p, p)
  crossprod(A) + diag(p) * 0.5
}

test_that("perm_to_trailing_2x2 matches chol(P K P^T) for random PD K", {
  skip_if_not_installed("bgms")
  seeds <- c(1, 7, 42, 101, 999)
  for (p in c(3, 4, 5, 7, 10, 15, 20)) {
    for (seed in seeds) {
      K <- rand_pd(p, seed)
      for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
          truth <- bunch_truth(K, i, j)
          got   <- bgms:::chol_perm_trailing_2x2_cpp(K, i, j)
          expect_true(got$ok,
                      info = sprintf("p=%d seed=%d i=%d j=%d",
                                     p, seed, i, j))
          expect_equal(got$l_ii, truth$l_ii, tolerance = 1e-10,
                       info = sprintf("l_ii mismatch p=%d seed=%d i=%d j=%d",
                                      p, seed, i, j))
          expect_equal(got$l_ji, truth$l_ji, tolerance = 1e-10,
                       info = sprintf("l_ji mismatch p=%d seed=%d i=%d j=%d",
                                      p, seed, i, j))
          expect_equal(got$l_jj, truth$l_jj, tolerance = 1e-10,
                       info = sprintf("l_jj mismatch p=%d seed=%d i=%d j=%d",
                                      p, seed, i, j))
        }
      }
    }
  }
})

test_that("perm_to_trailing_2x2 handles non-PD K via ok = FALSE", {
  # K with a zero eigenvalue: chol() inside the test interface will fail and
  # the function returns ok = FALSE without throwing.
  p <- 5
  v <- c(1, -1, 0, 0, 0)
  K <- v %o% v   # rank-1, not PD
  got <- bgms:::chol_perm_trailing_2x2_cpp(K, 1, 2)
  expect_false(got$ok)
})

test_that("perm_to_trailing_2x2 is invariant under permutation of the rest", {
  # Schur(K_perm)_trailing = K_BB − K_BA K_AA⁻¹ K_AB depends on A only through
  # the symmetric sum K_BA K_AA⁻¹ K_AB, which is invariant under permutations
  # of the "rest" set A. So the full trailing 2×2 (l_ii, l_ji, l_jj) is
  # invariant. Test by comparing Bunch against a reverse-rest ground truth.
  reverse_rest_truth <- function(K, i, j) {
    p <- nrow(K)
    perm <- c(rev(setdiff(seq_len(p), c(i, j))), i, j)
    Lp <- t(chol(K[perm, perm]))
    list(l_ii = Lp[p - 1, p - 1],
         l_ji = Lp[p,     p - 1],
         l_jj = Lp[p,     p])
  }
  K <- rand_pd(8, 2026)
  for (i in 1:7) for (j in (i + 1):8) {
    a <- bgms:::chol_perm_trailing_2x2_cpp(K, i, j)
    b <- reverse_rest_truth(K, i, j)
    expect_equal(a$l_ii, b$l_ii, tolerance = 1e-10)
    expect_equal(a$l_ji, b$l_ji, tolerance = 1e-10)
    expect_equal(a$l_jj, b$l_jj, tolerance = 1e-10)
  }
})
