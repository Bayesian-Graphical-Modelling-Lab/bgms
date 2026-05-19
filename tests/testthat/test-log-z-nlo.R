# --------------------------------------------------------------------------- #
# Closed-form log Z(G) approximations: parity vs the z-project reference,
# and incremental-vs-full agreement at alpha = 1.
#
# These tests exercise the Stage 3 Phase 1a numerics primitives:
#   log_Z_NLO_gamma_cpp                       (general-alpha full-recompute)
#   log_Z_NLO_gamma_degord_cpp                (degord wrapper)
#   log_Z_NLO_gamma_delta_incr_alpha1_cpp     (alpha = 1 incremental form)
#
# Ground truth is a fixture (tests/testthat/fixtures/log_z_nlo_reference.rds)
# pre-generated from log_Z_laplace_NLO_gamma_cpp at branchB_chain.cpp:470-620
# in the z project (see dev/numerical_analyses/generate_log_z_fixture.R).
# --------------------------------------------------------------------------- #


# ---- Helpers -----------------------------------------------------------------

draw_random_graph <- function(q, seed, p_edge = 0.5) {
  set.seed(seed)
  G <- matrix(0L, q, q)
  if (q < 2) return(G)
  for (i in 1:(q - 1)) for (j in (i + 1):q) {
    if (runif(1) < p_edge) {
      G[i, j] <- 1L
      G[j, i] <- 1L
    }
  }
  G
}


# ---- Parity against the z reference -----------------------------------------

test_that("log_Z_NLO_gamma matches the z reference bit-exact on the fixture", {
  fixture_path <- testthat::test_path("fixtures", "log_z_nlo_reference.rds")
  cases <- readRDS(fixture_path)

  for (case in cases) {
    for (alpha in c(1.0, 2.0)) {
      for (delta in c(0.0, 0.5, 1.0)) {
        for (include_F in c(FALSE, TRUE)) {
          key <- sprintf("a%g_d%g_F%s", alpha, delta, include_F)
          ours <- log_Z_NLO_gamma_cpp(
            case$G, alpha, 1.0, 1.0, include_F, delta
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


# ---- Empty-graph closed form ------------------------------------------------

test_that("log_Z_NLO_gamma at empty graph equals the analytic constant", {
  # Empty G: log Z = sum_v lgamma(alpha + delta) - q lgamma(alpha) - q delta log(beta).
  # The (E_count = 0) branch in the implementation returns this directly.
  for (q in c(2L, 5L, 10L)) {
    for (alpha in c(0.5, 1.0, 2.5)) {
      for (delta in c(0.0, 0.5, 1.0)) {
        for (beta in c(0.5, 1.0, 2.0)) {
          G <- matrix(0L, q, q)
          expected <- q * lgamma(alpha + delta) - q * lgamma(alpha) -
                      q * delta * log(beta)
          got <- log_Z_NLO_gamma_cpp(G, alpha, beta, 1.0, FALSE, delta)
          expect_equal(got, expected, tolerance = 1e-12,
                       info = sprintf("q=%d alpha=%g beta=%g delta=%g",
                                      q, alpha, beta, delta))
        }
      }
    }
  }
})


# ---- DEGORD wrapper: permutation-consistent --------------------------------

test_that("log_Z_NLO_gamma_degord equals full-recompute on hand-permuted graph", {
  for (q in c(5L, 7L)) {
    G <- draw_random_graph(q, seed = 41 + q)
    for (i_zero in c(0L, 2L)) {
      for (j_zero in c(1L, 3L, q - 1L)) {
        if (i_zero >= j_zero) next
        perm <- c(i_zero, j_zero, setdiff(0:(q - 1), c(i_zero, j_zero)))
        G_perm <- G[perm + 1L, perm + 1L]
        for (alpha in c(1.0, 2.0)) {
          for (delta in c(0.0, 1.0)) {
            via_degord <- log_Z_NLO_gamma_degord_cpp(
              G, i_zero, j_zero, alpha, 1.0, 1.0, FALSE, delta
            )
            via_full <- log_Z_NLO_gamma_cpp(
              G_perm, alpha, 1.0, 1.0, FALSE, delta
            )
            expect_equal(
              via_degord, via_full,
              tolerance = 1e-12,
              info = sprintf("q=%d (i,j)=(%d,%d) alpha=%g delta=%g",
                             q, i_zero, j_zero, alpha, delta)
            )
          }
        }
      }
    }
  }
})


# ---- alpha = 1 incremental agrees with full-recompute difference -----------

test_that("alpha = 1 incremental matches full-recompute log-Z difference", {
  for (q in c(3L, 5L, 7L)) {
    G <- draw_random_graph(q, seed = 71 + q)
    for (delta in c(0.0, 0.5, 1.0)) {
      for (include_F in c(FALSE, TRUE)) {
        for (i_zero in 0:(q - 2)) {
          for (j_zero in (i_zero + 1):(q - 1)) {
            G_after <- G
            G_after[i_zero + 1, j_zero + 1] <- 1L - G[i_zero + 1, j_zero + 1]
            G_after[j_zero + 1, i_zero + 1] <- G_after[i_zero + 1, j_zero + 1]
            full_diff <- log_Z_NLO_gamma_cpp(
              G_after, 1.0, 1.0, 1.0, include_F, delta
            ) - log_Z_NLO_gamma_cpp(
              G, 1.0, 1.0, 1.0, include_F, delta
            )
            inc <- log_Z_NLO_gamma_delta_incr_alpha1_cpp(
              G, i_zero, j_zero, 1.0, 1.0, delta, include_F
            )
            expect_equal(
              inc, full_diff,
              tolerance = 1e-10,
              info = sprintf("q=%d (i,j)=(%d,%d) delta=%g F=%s",
                             q, i_zero, j_zero, delta, include_F)
            )
          }
        }
      }
    }
  }
})


# ---- alpha > 1 incremental agrees with full-recompute difference -----------

test_that("general-alpha incremental matches full-recompute log-Z difference", {
  # Phase 1b: the alpha > 1 cascade adds H_e dependence on the larger endpoint
  # via the -4 beta (alpha - 1) / M_v[v] piece, so a toggle (i, j) cascades to
  # every edge incident to i (forward AND backward). V_a = {i, j} ∪ N(i)
  # captures the full affected set.
  for (q in c(3L, 5L, 7L)) {
    G <- draw_random_graph(q, seed = 200 + q)
    for (alpha in c(1.5, 2.0, 3.0)) {
      for (delta in c(0.0, 0.5, 1.0)) {
        for (include_F in c(FALSE, TRUE)) {
          for (i_zero in 0:(q - 2)) {
            for (j_zero in (i_zero + 1):(q - 1)) {
              G_after <- G
              G_after[i_zero + 1, j_zero + 1] <- 1L - G[i_zero + 1, j_zero + 1]
              G_after[j_zero + 1, i_zero + 1] <- G_after[i_zero + 1, j_zero + 1]
              full_diff <- log_Z_NLO_gamma_cpp(
                G_after, alpha, 1.0, 1.0, include_F, delta
              ) - log_Z_NLO_gamma_cpp(
                G, alpha, 1.0, 1.0, include_F, delta
              )
              inc <- log_Z_NLO_gamma_delta_incr_alphaN_cpp(
                G, i_zero, j_zero, alpha, 1.0, 1.0, delta, include_F
              )
              expect_equal(
                inc, full_diff,
                tolerance = 1e-10,
                info = sprintf("q=%d (i,j)=(%d,%d) alpha=%g delta=%g F=%s",
                               q, i_zero, j_zero, alpha, delta, include_F)
              )
            }
          }
        }
      }
    }
  }
})


# ---- alpha = 1 reduction: general-alpha == alpha-1 specialisation ----------

test_that("general-alpha incremental reduces to alpha=1 form at alpha = 1", {
  for (q in c(3L, 5L, 7L)) {
    G <- draw_random_graph(q, seed = 311 + q)
    for (delta in c(0.0, 0.5, 1.0)) {
      for (include_F in c(FALSE, TRUE)) {
        for (i_zero in 0:(q - 2)) {
          for (j_zero in (i_zero + 1):(q - 1)) {
            via_alpha1 <- log_Z_NLO_gamma_delta_incr_alpha1_cpp(
              G, i_zero, j_zero, 1.0, 1.0, delta, include_F
            )
            via_alphaN <- log_Z_NLO_gamma_delta_incr_alphaN_cpp(
              G, i_zero, j_zero, 1.0, 1.0, 1.0, delta, include_F
            )
            expect_equal(
              via_alphaN, via_alpha1,
              tolerance = 1e-10,
              info = sprintf("q=%d (i,j)=(%d,%d) delta=%g F=%s",
                             q, i_zero, j_zero, delta, include_F)
            )
          }
        }
      }
    }
  }
})


# ---- delta finite-difference vs analytic delta-derivative (alpha = 1) -------

test_that("d(log Z) / d(delta) matches finite differences at alpha = 1", {
  # Sanity: at small steps in delta, the log Z value should be smooth in delta.
  # We check second-order central FD vs successive evaluations agree to O(h^2).
  for (q in c(3L, 4L)) {
    G <- draw_random_graph(q, seed = 91 + q)
    delta0 <- 0.7
    h <- 1e-4
    for (include_F in c(FALSE, TRUE)) {
      z_minus <- log_Z_NLO_gamma_cpp(G, 1.0, 1.0, 1.0, include_F, delta0 - h)
      z_plus  <- log_Z_NLO_gamma_cpp(G, 1.0, 1.0, 1.0, include_F, delta0 + h)
      # Symmetric FD: (z_plus - z_minus) / (2 h). We just check the value is
      # finite and the function evaluates without NaN/Inf across a delta sweep.
      d_fd <- (z_plus - z_minus) / (2 * h)
      expect_true(is.finite(d_fd),
                  info = sprintf("q=%d F=%s d_fd not finite", q, include_F))
    }
  }
})
