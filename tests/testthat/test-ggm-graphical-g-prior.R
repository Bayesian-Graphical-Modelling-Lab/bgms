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
      g_trace   = fit$gg_diagnostics$g[[1]],
      post_pip  = extract_posterior_inclusion_probabilities(fit),
      prior_pip = extract_prior_inclusion_probabilities(fit)
    )
  }

  out_zs  = run_fit(graphical_g_prior(g_hyperprior = "zellner_siow",   b = 1))
  out_hg  = run_fit(graphical_g_prior(g_hyperprior = "hyper_g",        a = 3))
  out_hgn = run_fit(graphical_g_prior(g_hyperprior = "hyper_g_over_n", a = 3))

  # MH-on-log(g) is moving (positive variance).
  expect_gt(var(out_zs$g_trace),  0)
  expect_gt(var(out_hg$g_trace),  0)
  expect_gt(var(out_hgn$g_trace), 0)

  # Posterior strong-edge PIP is high under all three hyperpriors.
  expect_gt(out_zs$post_pip[1, 2],  0.6)
  expect_gt(out_hg$post_pip[1, 2],  0.6)
  expect_gt(out_hgn$post_pip[1, 2], 0.6)
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
  # Posterior PIP on the strong edge should be high.
  post = extract_posterior_inclusion_probabilities(fit)
  expect_gt(post[1, 2], 0.6)
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


test_that("extract_prior_inclusion_probabilities() returns a symmetric p x p matrix for graphical_g_prior", {
  set.seed(42)
  p = 4L; n = 200L
  Sigma = diag(p)
  Sigma[1, 2] = Sigma[2, 1] = 0.7
  X = MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  fit = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = graphical_g_prior(g_hyperprior = "conjugate_gamma",
                                          a0 = 1, b0 = 1, g_init = 1),
    edge_selection = TRUE, delta = 0,
    iter = 800L, warmup = 500L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  prior_pip = extract_prior_inclusion_probabilities(fit)
  expect_true(is.matrix(prior_pip))
  expect_equal(dim(prior_pip), c(p, p))
  expect_equal(prior_pip, t(prior_pip))                       # symmetric
  expect_true(all(is.na(diag(prior_pip))))                    # NA diagonal
  off = prior_pip[upper.tri(prior_pip)]
  expect_true(all(off > 0 & off < 1))
  # GG-prior path varies per edge (simulated), not all identical.
  expect_gt(diff(range(off)), 0)
})


test_that("extract_prior_inclusion_probabilities() is analytic for non-GG priors", {
  # For any non-GG interaction prior the marginal is read off the edge
  # prior. Bernoulli(0.5) -> all off-diagonals are exactly 0.5.
  set.seed(7)
  p = 4L; n = 150L
  X = matrix(stats::rnorm(n * p), n, p)
  fit_cauchy = bgm(
    x = X, variable_type = "continuous",
    interaction_prior = cauchy_prior(scale = 2.5),
    edge_selection = TRUE,
    iter = 200L, warmup = 200L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  prior_pip = extract_prior_inclusion_probabilities(fit_cauchy)
  expect_equal(dim(prior_pip), c(p, p))
  off = prior_pip[upper.tri(prior_pip)]
  expect_true(all(off == 0.5))     # Bernoulli(0.5) is the default
})


test_that("extract_prior_inclusion_probabilities() caches the GG simulation", {
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
  m_a = extract_prior_inclusion_probabilities(fit)
  cache = bgms:::get_fit_cache(fit)
  expect_false(is.null(cache$prior_inclusion_pips))
  m_b = extract_prior_inclusion_probabilities(fit)            # cache hit
  expect_identical(m_a, m_b)
  m_c = extract_prior_inclusion_probabilities(fit, recompute = TRUE)
  expect_equal(m_a, m_c)
})


test_that("extract_prior_inclusion_probabilities() errors on misuse", {
  set.seed(1)
  p = 4L; n = 60L
  X = matrix(stats::rnorm(n * p), n, p)
  fit_noedge = bgm(
    x = X, variable_type = "continuous",
    edge_selection = FALSE,
    iter = 100L, warmup = 100L, chains = 1L, cores = 1L,
    seed = 1L, display_progress = "none"
  )
  expect_error(extract_prior_inclusion_probabilities(fit_noedge),
               "edge_selection = TRUE")
})


test_that("summary() does not auto-compute any BF-related columns", {
  # All BF / prior-PIP computation is now in the dedicated extractor.
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
  post = extract_posterior_inclusion_probabilities(fit)
  prior = extract_prior_inclusion_probabilities(fit)
  post_off  = post[upper.tri(post)]
  prior_off = prior[upper.tri(prior)]
  bfs = (post_off / (1 - post_off)) / (prior_off / (1 - prior_off))
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

  # Strong edge V1-V2: posterior PIP should be high under both priors.
  post_fixed = extract_posterior_inclusion_probabilities(fit_fixed)
  post_cg    = extract_posterior_inclusion_probabilities(fit_cg)
  expect_gt(post_fixed[1, 2], 0.8)
  expect_gt(post_cg[1, 2],    0.8)
})
