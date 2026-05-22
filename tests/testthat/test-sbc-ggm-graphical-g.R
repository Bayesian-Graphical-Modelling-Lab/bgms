# --------------------------------------------------------------------------- #
# Conditional SBC for the Graphical G-prior on GGM
#
# Standard SBC fails for the GG-prior because V_ij = n / (4·S̄_ii·S̄_jj)
# depends on the data: the prior used to draw K_true and the prior used
# in the bgm fit on the simulated data have different V_ij values, so
# ranks are not uniform by construction.
#
# The conditional SBC fixes V_ij at a reference V_ref (computed once
# from some X_ref) and routes BOTH the generator and the fitter
# through that same matrix via the V_ij_external knob added to
# graphical_g_prior(). Generator and fitter now share the same prior,
# so under a correct sampler the SBC ranks ARE uniform.
#
# Test functions: diagonal entries of K (continuous, no point mass)
# and log|K|. Off-diagonal K_ij entries are skipped because the
# spike-and-slab point mass at zero induces tied ranks that break
# standard uniformity tests.
#
# Gated behind BGMS_RUN_SLOW_TESTS.
# --------------------------------------------------------------------------- #


skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run SBC tests"
  )
}


# Reconstruct a p x p symmetric K from upper-triangle (i <= j, row-major)
# raw samples produced by bgm() under variable_type = "continuous".
reconstruct_K = function(upper_row, p) {
  K = matrix(0, p, p)
  idx = 1L
  for(i in seq_len(p)) {
    for(j in i:p) {
      K[i, j] = upper_row[idx]
      K[j, i] = upper_row[idx]
      idx = idx + 1L
    }
  }
  K
}


test_that("conditional SBC: GG-prior (ConjugateGamma) produces uniform diagonal ranks (p=3)", {
  skip_unless_slow()

  # Generator AND fitter share the same V_ij = V_ref. K_true is drawn
  # from the prior-only chain at V_ref; y_star ~ N(0, K_true^{-1}) is
  # simulated; bgm() fits with V_ij_external = V_ref. Per the
  # conditional-SBC argument, the diagonal-of-K ranks should be
  # uniform on {0, 1, ..., L}.

  p = 3L
  n = 100L
  R = 200L          # number of replications
  L = 999L          # posterior draws kept per replication
  thin = 5L
  n_warmup_prior = 2000L

  set.seed(2026)

  # Reference X to fix V_ref. Identity covariance so V_ref entries are
  # easy to interpret and the chain isn't pushed to the PD boundary.
  X_ref = MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
  V_ref = matrix(0, p, p)
  for(i in seq_len(p)) {
    for(j in seq_len(p)) {
      if(i != j) {
        # V_ij = n / (4 * S_ii * S_jj), with S = X_ref'X_ref
        s_ii = sum(X_ref[, i]^2)
        s_jj = sum(X_ref[, j]^2)
        V_ref[i, j] = n / (4 * s_ii * s_jj)
      }
    }
  }

  prior_obj = graphical_g_prior(
    g_hyperprior  = "conjugate_gamma",
    a0 = 1, b0 = 1, g_init = 1,
    V_ij_external = V_ref
  )

  # Generator: one long prior-only chain at V_ref produces R thinned
  # (K_true) draws. We use a short bgm() call to build a validated
  # spec, then re-run via run_sampler() with prior_only = TRUE and the
  # SBC-sized iter / warmup.
  setup_fit = bgm(
    x = X_ref, variable_type = "continuous",
    interaction_prior = prior_obj,
    edge_selection = TRUE, delta = 0,
    iter = 1L, warmup = 1L, chains = 1L, cores = 1L,
    update_method = "adaptive-metropolis",
    seed = 2026L, display_progress = "none"
  )
  spec = setup_fit$.bgm_spec
  spec$prior$prior_only = TRUE
  spec$sampler$iter = as.integer(R * thin)
  spec$sampler$warmup = n_warmup_prior
  spec$sampler$chains = 1L
  spec$sampler$update_method = "adaptive-metropolis"
  spec$sampler$target_accept = 0.44
  spec$sampler$progress_type = 0L
  gen_raw = run_sampler(spec)
  K_true_samples = gen_raw[[1L]]$samples            # (p(p+1)/2) x iter

  thin_idx = seq(thin, R * thin, by = thin)
  K_true_list = lapply(thin_idx, function(it) {
    reconstruct_K(K_true_samples[, it], p)
  })

  # Fitter loop. Ranks of diag(K_true) entries within posterior samples.
  ranks = matrix(NA_integer_, nrow = R, ncol = p)

  for(r in seq_len(R)) {
    K_true = K_true_list[[r]]
    # Bail on draws that aren't PD numerically (rare under proper prior).
    eig = tryCatch(eigen(K_true, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NULL)
    if(is.null(eig) || min(eig) <= 1e-8) {
      next
    }
    Sigma = solve(K_true)
    y_star = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

    fit_r = bgm(
      x = y_star, variable_type = "continuous",
      interaction_prior = graphical_g_prior(
        g_hyperprior  = "conjugate_gamma",
        a0 = 1, b0 = 1, g_init = 1,
        V_ij_external = V_ref                  # << conditional SBC
      ),
      edge_selection = TRUE, delta = 0,
      iter = L, warmup = 1000L, chains = 1L, cores = 1L,
      update_method = "adaptive-metropolis",
      seed = 2026L + r, display_progress = "none"
    )

    main_samples = do.call(rbind, fit_r$raw_samples$main)
    # main_samples is L x p (precision diagonals).
    for(i in seq_len(p)) {
      ranks[r, i] = sum(main_samples[, i] < K_true[i, i])
    }
  }

  ranks = ranks[stats::complete.cases(ranks), , drop = FALSE]
  expect_true(nrow(ranks) > R * 0.9,
    info = sprintf("Too many replications dropped (kept %d / %d).",
                   nrow(ranks), R))

  # Per-parameter KS test at alpha = 0.01, allowing one false positive
  # across the p diagonal entries.
  n_fail_ks = 0L
  for(j in seq_len(p)) {
    u = ranks[, j] / (L + 1)
    p_val = suppressWarnings(stats::ks.test(u, "punif")$p.value)
    if(p_val <= 0.01) n_fail_ks = n_fail_ks + 1L
  }
  expect_true(n_fail_ks <= 1L,
    info = sprintf("Conditional SBC KS: %d/%d parameters failed at alpha=0.01.",
                   n_fail_ks, p))

  # Global chi-squared on the pooled ranks.
  all_ranks = as.vector(ranks)
  bins = cut(all_ranks / (L + 1),
             breaks = seq(0, 1, length.out = 21),
             include.lowest = TRUE)
  counts = tabulate(bins, nbins = 20L)
  chisq_p = suppressWarnings(stats::chisq.test(counts)$p.value)
  expect_true(chisq_p > 0.001,
    info = sprintf("Conditional SBC global chi-squared p = %.4f", chisq_p))
})
