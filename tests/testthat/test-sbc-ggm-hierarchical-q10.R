# --------------------------------------------------------------------------- #
# SBC for the hierarchical-spec GGM with edge selection at q = 10.
#
# Sister of test-sbc-ggm-hierarchical.R (which fixes q = 5). Truth is
# the same NUTS-at-fixed-Gamma generator from sample_ggm_prior() with a
# longer warmup; we trust it at q = 10 under the operational alpha = 1
# Gamma diagonal prior. AKM-J is NOT used here (its perfect-sampler
# cap-saturation at q = 10, delta = 2 would corrupt the truth pool).
#
# Truth generation:
#   1. Draw Gamma_true ~ Bern(p_inc) on the upper triangle.
#   2. Draw one well-mixed K_true via sample_ggm_prior() at fixed
#      Gamma_true (constrained NUTS, n_warmup_truth iterations).
#   3. Generate Y | K_true ~ N(0, K_true^{-1}), n_obs cases.
#
# Inference:
#   4. Fit bgm() with graph_prior_spec = "hierarchical" on Y.
#
# Tested cells:
#   delta ∈ {0, 1, 2}, alpha = 1, R = 500 reps each.
#
# Uniformity tested by KS at alpha = 0.01 per K_ii per delta. Edge-
# indicator marginal calibration (P(gamma_ij = 1 | Y) tracks observed
# Gamma_true frequency) is a sanity check.
#
# Gated behind BGMS_RUN_SLOW_TESTS. Expected runtime is multi-hour on
# M5 Pro; this is the validation sweep meant to be run overnight.
# --------------------------------------------------------------------------- #


# ---- Skip gate -------------------------------------------------------------

skip_if_not(
  identical(tolower(Sys.getenv("BGMS_RUN_SLOW_TESTS")), "true"),
  "Set BGMS_RUN_SLOW_TESTS=true to run hierarchical-spec SBC tests."
)


# ---- Test parameters -------------------------------------------------------

p          <- 10L
n_obs      <- 200L
R          <- 500L
delta_grid <- c(0.0, 1.0, 2.0)
p_inc      <- 0.5
# Z-matched chain config (see ~/Dropbox/Projecten/SV/Z/R/scripts/
# sbc_degord_hier_paired_sweep.R). Z's q=10 SBC passes at delta ∈ {1, 2}
# only with these heavier chains; bgms's earlier 100 + 400 was 20× too
# short. delta = 0 is OUT of the operational regime for the V/RR
# estimator and is expected to fail (kept in the sweep as documentation).
iter_post  <- 2000L
warmup_post <- 2000L
n_warmup_truth <- 2000L

# Across-rep parallelism. Each rep is single-threaded (bgm(cores = 1L)); the
# fan-out happens here via parallel::mclapply. Capped at 12 to match the user-
# specified core budget; on machines with fewer cores the cap drops to
# detectCores().
N_CORES <- min(parallel::detectCores(), 12L)


# ---- One-rep helper --------------------------------------------------------

one_rep <- function(r, delta) {
  # Truth: Gamma_true ~ Bern(p_inc) on upper triangle, K_true | Gamma_true.
  set.seed(10000L + r * 100L + as.integer(delta * 100))
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
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta            = delta,
    edge_indicators  = G_full,
    seed             = 20000L + r,
    verbose          = FALSE
  )
  K_diag_true    <- as.numeric(truth$K_diag[1L, ])
  K_offdiag_true <- as.numeric(truth$K_offdiag[1L, ])

  # Reconstruct K_true (p x p). offdiag_names is row-major; parse the
  # names rather than trust the natural ordering.
  K_true <- diag(K_diag_true)
  for (k in seq_along(truth$offdiag_names)) {
    ij <- strsplit(truth$offdiag_names[k], "_")[[1L]]
    i <- as.integer(ij[2L])
    j <- as.integer(ij[3L])
    K_true[i, j] <- truth$K_offdiag[1L, k]
    K_true[j, i] <- K_true[i, j]
  }

  # Data: Y ~ N(0, K^{-1}).
  Sigma_true <- solve(K_true)
  Y <- MASS::mvrnorm(n_obs, mu = rep(0, p), Sigma = Sigma_true)

  # Inference under the hierarchical-spec chain.
  fit <- bgm(
    Y,
    variable_type         = "continuous",
    interaction_prior     = normal_prior(scale = 1),
    precision_scale_prior = gamma_prior(shape = 1, rate = 1),
    delta            = delta,
    graph_prior_spec = "hierarchical",
    # Z-matched: M_inner = 50 (half), kappa = 2 (doubled for bigger V-series
    # convergence margin at q = 10).
    z_ratio_tuning   = list(M_inner = 50L, kappa = 2.0, rho = 0.5),
    iter             = iter_post,
    warmup           = warmup_post,
    update_method    = "adaptive-metropolis",
    chains           = 1L, cores = 1L,
    seed             = 30000L + r,
    display_progress = "none",
    verbose          = FALSE
  )
  raw      <- S7::prop(fit, "raw_samples")
  main_chn <- raw$main[[1L]]            # iter x p, bgms convention: K_yy_ii = K_ii / 2
  ind_chn  <- raw$indicator[[1L]]       # iter x p(p-1)/2, 0/1
  # Rank K_ii on the K_yy partial-association scale (both truth and bgm).
  K_ii_rank <- vapply(seq_len(p), function(i) {
    sum(main_chn[, i] < K_diag_true[i])
  }, integer(1L))

  # Gamma_true upper triangle, edge-order matched to bgm indicator slots.
  gamma_true_vec <- G_true[upper.tri(G_true)]
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
  results <- parallel::mclapply(
    seq_len(R),
    function(r) one_rep(r, delta),
    mc.cores = N_CORES,
    mc.preschedule = TRUE
  )
  # Surface any worker errors before the assertion stage so the failure
  # message references the bad rep rather than a downstream stack-rbind.
  err_idx <- which(vapply(results, inherits, logical(1L), what = "try-error"))
  if (length(err_idx) > 0L) {
    stop(sprintf(
      "SBC q=10: %d/%d reps errored at delta=%g; first failure at r=%d:\n%s",
      length(err_idx), R, delta, err_idx[1L], conditionMessage(attr(results[[err_idx[1L]]], "condition"))
    ))
  }

  K_ii_ranks <- do.call(rbind, lapply(results, `[[`, "K_ii_rank"))
  gamma_true_mat <- do.call(rbind, lapply(results, `[[`, "gamma_true_vec"))
  gamma_post_mat <- do.call(rbind, lapply(results, `[[`, "gamma_post_mean"))
  n_iter <- results[[1L]]$n_iter_post

  # Per-K_ii KS test of normalised ranks ~ Uniform(0, 1).
  ks_p <- vapply(seq_len(p), function(i) {
    u <- (K_ii_ranks[, i] + 0.5) / (n_iter + 1L)
    suppressWarnings(stats::ks.test(u, "punif")$p.value)
  }, double(1L))

  for (i in seq_len(p)) {
    test_that(sprintf("SBC q=10: K_ii[%d] ranks uniform under hierarchical (delta=%g)",
                      i, delta), {
      expect_gt(ks_p[i], 0.01,
                label = sprintf("ks_p=%.3g (delta=%g, K_ii[%d])",
                                ks_p[i], delta, i))
    })
  }

  # Gamma marginal calibration: posterior P(gamma_ij = 1 | Y) averaged across
  # reps should track the prior inclusion probability p_inc.
  test_that(sprintf("SBC q=10: hierarchical-spec gamma marginal tracks p_inc (delta=%g)",
                    delta), {
    mean_post_gamma <- mean(gamma_post_mat)
    mean_true_gamma <- mean(gamma_true_mat)
    se_gap <- sqrt(p_inc * (1 - p_inc) / (R * ncol(gamma_post_mat)))
    expect_lt(abs(mean_post_gamma - mean_true_gamma), 4 * se_gap,
              label = sprintf("|post-mean(gamma) - true-mean(gamma)|=%.3g, SE=%.3g",
                              abs(mean_post_gamma - mean_true_gamma), se_gap))
  })
}
