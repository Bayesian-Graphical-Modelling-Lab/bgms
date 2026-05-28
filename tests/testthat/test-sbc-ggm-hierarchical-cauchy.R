# --------------------------------------------------------------------------- #
# SBC for the hierarchical-spec GGM with a Cauchy slab.
#
# Companion to test-sbc-ggm-hierarchical.R (Normal slab). The Cauchy slab
# routes through the scale-mixture-of-normals representation: per-edge
# omega_ij ~ IG(1/2, 1/2) is slice-sampled by
# savage_dickey::slice_sample_cauchy_omega_active in addition to the
# usual SD between-step. See math/savage_dickey/cauchy_omega.h.
#
# Truth: sample_ggm_prior() at fixed Gamma_true with cauchy_prior +
# gamma_prior (NUTS at fixed Gamma, independent of the hierarchical-spec
# SD chain). Inference: bgm() with the same encompassing prior and
# prior_factorization = "hierarchical".
#
# Rank uniformity tested by KS at alpha = 0.01 per K_ii per delta.
# Gamma marginal calibration as a sanity check.
#
# Gated behind BGMS_RUN_SLOW_TESTS. Smaller R (= 200) than the Normal
# test to keep total runtime under an hour; ~3-4x faster per rep than
# the Normal hierarchical SBC because the smaller delta sweep avoids
# the costly high-tilt cells.
# --------------------------------------------------------------------------- #


# ---- Skip gate -------------------------------------------------------------

skip_if_not(
  identical(tolower(Sys.getenv("BGMS_RUN_SLOW_TESTS")), "true"),
  "Set BGMS_RUN_SLOW_TESTS=true to run hierarchical-spec SBC tests."
)
skip_if_not_installed("MASS")


# ---- Test parameters -------------------------------------------------------

p              <- 5L
n_obs          <- 100L
R              <- 200L
delta_grid     <- c(0.0, 0.5, 1.0)
p_inc          <- 0.5
iter_post      <- 400L
warmup_post    <- 100L
n_warmup_truth <- 500L
slab_scale     <- 1.0


# ---- One-rep helper --------------------------------------------------------

one_rep_cauchy <- function(r, delta) {
  set.seed(40000L + r * 100L + as.integer(delta * 100))
  G_true <- matrix(0L, p, p)
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    if (runif(1) < p_inc) {
      G_true[i, j] <- 1L
      G_true[j, i] <- 1L
    }
  }
  G_full <- G_true
  diag(G_full) <- 1L

  truth <- sample_ggm_prior(
    p             = p,
    n_samples     = 1L,
    n_warmup      = n_warmup_truth,
    interaction_prior     = cauchy_prior(scale = slab_scale),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta            = delta,
    edge_indicators  = G_full,
    seed             = 50000L + r,
    verbose          = FALSE
  )
  K_diag_true    <- as.numeric(truth$K_diag[1L, ])

  K_true <- diag(K_diag_true)
  for (k in seq_along(truth$offdiag_names)) {
    ij <- strsplit(truth$offdiag_names[k], "_")[[1L]]
    i <- as.integer(ij[2L])
    j <- as.integer(ij[3L])
    K_true[i, j] <- truth$K_offdiag[1L, k]
    K_true[j, i] <- K_true[i, j]
  }

  Sigma_true <- solve(K_true)
  Y <- MASS::mvrnorm(n_obs, mu = rep(0, p), Sigma = Sigma_true)

  fit <- bgm(
    Y,
    variable_type         = "continuous",
    interaction_prior     = cauchy_prior(scale = slab_scale),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta            = delta,
    prior_factorization = "hierarchical",
    iter             = iter_post,
    warmup           = warmup_post,
    update_method    = "adaptive-metropolis",
    chains           = 1L, cores = 1L,
    seed             = 60000L + r,
    display_progress = "none",
    verbose          = FALSE
  )
  raw      <- S7::prop(fit, "raw_samples")
  main_chn <- raw$main[[1L]]
  ind_chn  <- raw$indicator[[1L]]
  K_ii_rank <- vapply(seq_len(p), function(i) {
    sum(main_chn[, i] < K_diag_true[i])
  }, integer(1L))

  gamma_true_vec  <- G_true[upper.tri(G_true)]
  gamma_post_mean <- colMeans(ind_chn)

  list(
    K_ii_rank       = K_ii_rank,
    gamma_true_vec  = gamma_true_vec,
    gamma_post_mean = gamma_post_mean,
    n_iter_post     = nrow(main_chn)
  )
}


# ---- Run + uniformity tests across the delta sweep -------------------------

for (delta in delta_grid) {
  results <- vector("list", R)
  for (r in seq_len(R)) results[[r]] <- one_rep_cauchy(r, delta)

  K_ii_ranks     <- do.call(rbind, lapply(results, `[[`, "K_ii_rank"))
  gamma_true_mat <- do.call(rbind, lapply(results, `[[`, "gamma_true_vec"))
  gamma_post_mat <- do.call(rbind, lapply(results, `[[`, "gamma_post_mean"))
  n_iter <- results[[1L]]$n_iter_post

  ks_p <- vapply(seq_len(p), function(i) {
    u <- (K_ii_ranks[, i] + 0.5) / (n_iter + 1L)
    suppressWarnings(stats::ks.test(u, "punif")$p.value)
  }, double(1L))

  for (i in seq_len(p)) {
    test_that(sprintf("SBC: K_ii[%d] ranks uniform under Cauchy + hierarchical (delta=%g)",
                      i, delta), {
      expect_gt(ks_p[i], 0.01,
                label = sprintf("ks_p=%.3g (delta=%g, K_ii[%d])",
                                ks_p[i], delta, i))
    })
  }

  test_that(sprintf("Cauchy + hierarchical: gamma marginal tracks p_inc (delta=%g)",
                    delta), {
    mean_post_gamma <- mean(gamma_post_mat)
    mean_true_gamma <- mean(gamma_true_mat)
    se_gap <- sqrt(p_inc * (1 - p_inc) / (R * ncol(gamma_post_mat)))
    expect_lt(abs(mean_post_gamma - mean_true_gamma), 4 * se_gap,
              label = sprintf("|post-mean(gamma) - true-mean(gamma)|=%.3g, SE=%.3g",
                              abs(mean_post_gamma - mean_true_gamma), se_gap))
  })
}
