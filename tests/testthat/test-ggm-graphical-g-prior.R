# --------------------------------------------------------------------------- #
# Operational validation for the Graphical G-prior (GGM, joint spec)
#
# These tests exercise the end-to-end GG-prior path through bgm() and the
# inclusion-Bayes-factor surface in summary(fit)$indicator. They are
# operational ("does the BF land where it should under known scenarios")
# rather than calibration ("are the posterior ranks uniform"): full
# joint-spec SBC for the GG-prior is impeded by the data-dependence of
# the per-edge slab variance V_ij = n / (4 * S̄_ii * S̄_jj), which means
# the prior used to draw K_true and the prior used in the bgm fit on the
# simulated data don't share the same V_ij. The PNAS-style smoke tests
# here therefore focus on whether BFs discriminate signal from null
# correctly, which is what matters for the simulation study downstream.
#
# Fast tests (always run):
#   - constructor + bgm() integration smoke
#   - structural properties of summary()$indicator under graphical_g_prior
#   - BF on a strong-signal edge is large; BF on null edges is moderate
#
# Slow tests (BGMS_RUN_SLOW_TESTS=true):
#   - BF calibration under a null DGP (median BF ~ 1, not systematically > 1)
#   - cross-prior consistency: Fixed vs ConjugateGamma agree qualitatively
# --------------------------------------------------------------------------- #


skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run slow GG-prior validation"
  )
}


# ---- Fast tests --------------------------------------------------------------

test_that("graphical_g_prior() returns a proper bgms_interaction_prior", {
  gp = graphical_g_prior()
  expect_s3_class(gp, "bgms_graphical_g_prior")
  expect_s3_class(gp, "bgms_interaction_prior")
  expect_equal(gp$hyper.parameters$g_hyperprior, "conjugate_gamma")
  expect_equal(gp$hyper.parameters$a0, 1)
  expect_equal(gp$hyper.parameters$b0, 1)
})


test_that("bgm() runs end-to-end with interaction_prior = graphical_g_prior()", {
  set.seed(1)
  p = 4L; n = 80L
  X = matrix(stats::rnorm(n * p), n, p)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                          a0 = 1, b0 = 1, g_init = 1),
    delta = 0,
    iter = 200L, warmup = 200L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  expect_true(!is.null(fit$posterior_mean_indicator))
  expect_equal(dim(fit$posterior_mean_indicator), c(p, p))
  # gg_diagnostics present + correct shape.
  expect_false(is.null(fit$gg_diagnostics))
  expect_equal(length(fit$gg_diagnostics$t), 1L)            # 1 chain
  expect_equal(length(fit$gg_diagnostics$t[[1]]), 200L)     # n_iter
  expect_true(all(fit$gg_diagnostics$g[[1]] > 0))
})


test_that("Zellner-Siow / hyper-g / hyper-g/n drive a varying g-trace", {
  # Each non-conjugate g-hyperprior should produce a g-trace with positive
  # variance (the MH-on-log(g) update is actually moving), and the
  # strong-signal edge should still be detected.
  set.seed(42)
  p = 4L; n = 200L
  Sigma = diag(p)
  Sigma[1, 2] = Sigma[2, 1] = 0.6
  X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  run_fit = function(prior_obj) {
    fit = bgm(
      x = X, variable_type = "continuous",
      interaction_prior = prior_obj,
      edge_selection = TRUE, delta = 0,
      iter = 600L, warmup = 400L, chains = 1L, cores = 1L,
      seed = 1L, display_progress = "none"
    )
    list(
      g_trace = fit$gg_diagnostics$g[[1]],
      bf      = inclusion_bayes_factor(fit)
    )
  }

  out_zs  = run_fit(graphical_g_prior(g_hyperprior = "zellner_siow",   b = 1))
  out_hg  = run_fit(graphical_g_prior(g_hyperprior = "hyper_g",        a = 3))
  out_hgn = run_fit(graphical_g_prior(g_hyperprior = "hyper_g_over_n", a = 3))

  # MH-on-log(g) is moving (positive variance).
  expect_gt(var(out_zs$g_trace),  0)
  expect_gt(var(out_hg$g_trace),  0)
  expect_gt(var(out_hgn$g_trace), 0)

  # Strong-edge detection holds across all three hyperpriors.
  expect_true(is.infinite(out_zs$bf$bf[1])  || out_zs$bf$bf[1]  > 5)
  expect_true(is.infinite(out_hg$bf$bf[1])  || out_hg$bf$bf[1]  > 5)
  expect_true(is.infinite(out_hgn$bf$bf[1]) || out_hgn$bf$bf[1] > 5)
})


test_that("hyper_g and hyper_g_over_n reject a <= 2 at the constructor", {
  expect_error(
    graphical_g_prior(g_hyperprior = "hyper_g", a = 2),
    "'a' must be > 2"
  )
  expect_error(
    graphical_g_prior(g_hyperprior = "hyper_g_over_n", a = 1.5),
    "'a' must be > 2"
  )
  # conjugate_gamma path doesn't care about `a`.
  expect_silent(graphical_g_prior(g_hyperprior = "conjugate_gamma", a = 1))
})


test_that("tCCH runs end-to-end and drives a varying g-trace", {
  # The truncated compound confluent hypergeometric branch shares the
  # MH-on-log(g) infrastructure with ZS / hyper-g / hyper-g/n; only the
  # log-prior shape differs.
  set.seed(42)
  p = 4L; n = 200L
  Sigma = diag(p)
  Sigma[1, 2] = Sigma[2, 1] = 0.6
  X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "tCCH"),
    edge_selection = TRUE, delta = 0,
    iter = 600L, warmup = 400L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  expect_gt(var(fit$gg_diagnostics$g[[1]]), 0)
  bf_strong = inclusion_bayes_factor(fit)$bf[1L]
  expect_true(is.infinite(bf_strong) || bf_strong > 5)
})


test_that("tCCH with a = 0.5, b = 1 (rest default) matches hyper_g(3)", {
  # tCCH(a, b, r=0, s=0, u=1) reduces to π(g) ∝ (1+g)^(-(a+b)).
  # Setting a + b = 1.5 reproduces hyper_g(α=3) up to the unnormalised
  # prior shape, so the two chains should agree on summary statistics
  # under the same data + seed.
  set.seed(42)
  p = 4L; n = 200L
  Sigma = diag(p)
  Sigma[1, 2] = Sigma[2, 1] = 0.6
  X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  run = function(prior_obj) {
    fit = bgm(
      x = X, variable_type = "continuous",
      interaction_prior = prior_obj,
      edge_selection = TRUE, delta = 0,
      iter = 600L, warmup = 400L, chains = 1L, cores = 1L,
      seed = 1L, display_progress = "none"
    )
    mean(fit$gg_diagnostics$g[[1]])
  }
  g_tcch = run(graphical_g_prior(g_hyperprior = "tCCH",
                                  tcch = list(a = 0.5, b = 1,
                                              r = 0, s = 0, u = 1)))
  g_hg   = run(graphical_g_prior(g_hyperprior = "hyper_g", a = 3))
  expect_equal(g_tcch, g_hg, tolerance = 0.05)
})


test_that("graphical_g_prior(g_hyperprior = 'tCCH') validates tcch list", {
  # Rejects negative r and s and non-positive u; rejects incomplete list.
  expect_error(
    graphical_g_prior(g_hyperprior = "tCCH",
                      tcch = list(a = 1, b = 1, r = -0.5, s = 0, u = 1)),
    "tcch\\$r"
  )
  expect_error(
    graphical_g_prior(g_hyperprior = "tCCH",
                      tcch = list(a = 1, b = 1, r = 0, s = 0, u = -1)),
    "tcch\\$u"
  )
  expect_error(
    graphical_g_prior(g_hyperprior = "tCCH", tcch = list(a = 1)),
    "components a, b, r, s, u"
  )
})


test_that("fixed g_hyperprior pins t = sqrt(g_fixed) across iterations", {
  set.seed(1)
  p = 4L; n = 60L
  X = matrix(stats::rnorm(n * p), n, p)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "fixed", g_fixed = 4),
    delta = 0,
    iter = 50L, warmup = 50L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  t_trace = fit$gg_diagnostics$t[[1]]
  g_trace = fit$gg_diagnostics$g[[1]]
  expect_equal(t_trace, rep(2.0, 50L))
  expect_equal(g_trace, rep(4.0, 50L))
})


test_that("inclusion_bayes_factor() returns a tidy table for a GG-prior fit", {
  set.seed(42)
  p = 4L; n = 200L
  Sigma = diag(p)
  Sigma[1, 2] = Sigma[2, 1] = 0.7        # strong edge between V1 and V2
  X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                          a0 = 1, b0 = 1, g_init = 1),
    edge_selection = TRUE, delta = 0,
    iter = 800L, warmup = 500L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  bf_tbl = inclusion_bayes_factor(fit)
  expect_s3_class(bf_tbl, "data.frame")
  expect_equal(nrow(bf_tbl), p * (p - 1L) / 2L)
  expect_equal(colnames(bf_tbl),
               c("edge", "posterior_inclusion", "prior_inclusion", "bf"))
  expect_true(all(bf_tbl$prior_inclusion > 0 &
                  bf_tbl$prior_inclusion < 1))

  # The strong-signal edge (V1-V2) is the first row under upper-triangle
  # row-major ordering.
  strong_idx = 1L
  expect_gt(bf_tbl$posterior_inclusion[strong_idx], 0.8)
  null_bfs   = bf_tbl$bf[-strong_idx]
  strong_bf  = bf_tbl$bf[strong_idx]
  median_null_bf = median(null_bfs)
  expect_true(is.finite(median_null_bf))
  expect_true(is.infinite(strong_bf) || strong_bf > 5 * median_null_bf)
})


test_that("inclusion_bayes_factor() works for non-GG continuous fits too", {
  # The prior-only chain machinery is prior-agnostic — it mutes the data
  # likelihood and samples the slab + diagonal + edge prior. Any
  # continuous interaction prior (cauchy, normal, ...) works.
  set.seed(7)
  p = 4L; n = 150L
  X = matrix(stats::rnorm(n * p), n, p)
  fit_cauchy = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = cauchy_prior(scale = 2.5),
    edge_selection = TRUE,
    iter = 400L, warmup = 200L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  bf_tbl = inclusion_bayes_factor(fit_cauchy)
  expect_s3_class(bf_tbl, "data.frame")
  expect_equal(nrow(bf_tbl), p * (p - 1L) / 2L)
  expect_true(all(is.finite(bf_tbl$prior_inclusion)))
  expect_true(all(bf_tbl$prior_inclusion > 0 &
                  bf_tbl$prior_inclusion < 1))
})


test_that("inclusion_bayes_factor() caches and respects recompute", {
  set.seed(7)
  p = 4L; n = 150L
  X = matrix(stats::rnorm(n * p), n, p)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(),
    edge_selection = TRUE, delta = 0,
    iter = 200L, warmup = 200L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  tbl_a = inclusion_bayes_factor(fit)
  cache = bgms:::get_fit_cache(fit)
  expect_false(is.null(cache$inclusion_bf_prior_pips))
  tbl_b = inclusion_bayes_factor(fit)               # cache hit
  expect_identical(tbl_a$prior_inclusion, tbl_b$prior_inclusion)
  # recompute = TRUE re-runs the prior-only chain (with the same seed,
  # so identical values — the point is that it didn't error).
  tbl_c = inclusion_bayes_factor(fit, recompute = TRUE)
  expect_equal(tbl_a$prior_inclusion, tbl_c$prior_inclusion)
})


test_that("inclusion_bayes_factor() errors clearly on misuse", {
  # No edge selection -> error.
  set.seed(1)
  p = 4L; n = 60L
  X = matrix(stats::rnorm(n * p), n, p)
  fit_noedge = bgm(
    x = X, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 100L, warmup = 100L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  expect_error(inclusion_bayes_factor(fit_noedge),
               "edge_selection = TRUE")
})


test_that("summary() no longer auto-computes prior_inclusion / inclusion_bf", {
  # The auto-trigger from summary() / coef() was removed. summary() now
  # returns just the standard $indicator columns; users opt in via
  # inclusion_bayes_factor(fit).
  set.seed(42)
  p = 4L; n = 100L
  X = matrix(stats::rnorm(n * p), n, p)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(),
    edge_selection = TRUE, delta = 0,
    iter = 200L, warmup = 200L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  s = summary(fit)
  expect_false("prior_inclusion" %in% colnames(s$indicator))
  expect_false("inclusion_bf"   %in% colnames(s$indicator))
})


# ---- Slow tests --------------------------------------------------------------

test_that("BF calibration under a null DGP: median BF over null edges ~ 1", {
  skip_unless_slow()

  # Generate independent multivariate normal data (true K = I, no edges)
  # and check that the inclusion BFs on the (all-null) edges do not
  # systematically favour inclusion. The joint-spec BF is a ratio of
  # PIPs estimated from chains with shared Z̃(Γ) weighting, so under
  # the null DGP the BFs should hover near 1 across edges.
  set.seed(2026)
  p = 5L; n = 200L
  X = matrix(stats::rnorm(n * p), n, p)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                          a0 = 1, b0 = 1, g_init = 1),
    edge_selection = TRUE, delta = 0,
    iter = 4000L, warmup = 1000L, chains = 1L, cores = 1L,
    seed = 2026L, display_progress = "none"
  )
  bfs = inclusion_bayes_factor(fit)$bf
  # The median BF on null edges should be on the order of 1
  # (lenient bounds because 1 chain at q=5 still has nontrivial MC error).
  expect_true(median(bfs) < 5,
    info = sprintf("median BF on null edges = %.2f (expected ~ 1)",
                   median(bfs)))
  expect_true(median(bfs) > 0.05,
    info = sprintf("median BF on null edges = %.2f (expected ~ 1)",
                   median(bfs)))
})


test_that("Cross-prior consistency: Fixed vs ConjugateGamma agree on strong signal", {
  skip_unless_slow()

  # Same data, two GG-prior configurations: Fixed g vs ConjugateGamma.
  # Both should identify the strong edge with a large BF; the two
  # configurations may differ on null edges but should agree on the
  # signal edge.
  set.seed(7)
  p = 4L; n = 250L
  Sigma = diag(p)
  Sigma[1, 2] = Sigma[2, 1] = 0.7
  X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

  fit_fixed = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "fixed",
                                          g_fixed = 1),
    edge_selection = TRUE, delta = 0,
    iter = 1500L, warmup = 500L, chains = 1L, cores = 1L,
    seed = 7L, display_progress = "none"
  )
  fit_cg = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                          a0 = 1, b0 = 1, g_init = 1),
    edge_selection = TRUE, delta = 0,
    iter = 1500L, warmup = 500L, chains = 1L, cores = 1L,
    seed = 7L, display_progress = "none"
  )

  bf_fixed = inclusion_bayes_factor(fit_fixed)$bf
  bf_cg    = inclusion_bayes_factor(fit_cg)$bf

  # Strong edge index = 1 (V1-V2).
  expect_true(is.infinite(bf_fixed[1L]) || bf_fixed[1L] > 10,
    info = sprintf("Fixed BF[strong] = %.3g", bf_fixed[1L]))
  expect_true(is.infinite(bf_cg[1L])    || bf_cg[1L]    > 10,
    info = sprintf("ConjugateGamma BF[strong] = %.3g", bf_cg[1L]))
})
