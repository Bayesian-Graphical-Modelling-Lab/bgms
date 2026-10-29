# Tests for prior_sensitivity_check().
#
# The check refits the model at a few multipliers of the chosen slab scale, at
# the original fit's run length, and classifies every edge from the fits' own
# statistics. Structural tests use small fits with short explicit run lengths:
# every refit is a full MCMC run, so the run length is the dominant cost. The
# classification itself is tested on synthetic per-fit statistics. The
# bgmCompare checks run nightly (T1, BGMS_RUN_SLOW_TESTS); the GGM, mixed and
# warm-start certifications run weekly (T2, BGMS_RUN_CERTIFICATION).

skip_unless_slow = function() {
  skip_if_not(
    identical(Sys.getenv("BGMS_RUN_SLOW_TESTS"), "true"),
    message = "Set BGMS_RUN_SLOW_TESTS=true to run prior-sensitivity certifications"
  )
}


# ------------------------------------------------------------------------------
# The classification, on synthetic per-fit statistics
# ------------------------------------------------------------------------------

# Per-fit statistics in the shape refit_edge_stats() returns, for n edges.
fake_stats = function(lbf, mcse, prior_q = 0.5, lthr = log(10),
                      unanimous = TRUE, zeroflip = FALSE) {
  n = length(lbf)
  prior_odds = rep(prior_q / (1 - prior_q), n)
  pip = plogis(lbf + log(prior_odds))
  list(
    edge = paste0("e", seq_len(n)), pip = pip, prior_odds = prior_odds,
    lbf = lbf, verdict = verdict_from_lbf(lbf, lthr), mcse_lbf = mcse,
    band_half = rep(0, n), unanimous = rep_len(unanimous, n),
    zeroflip = rep_len(zeroflip, n)
  )
}

# A bgms_prior_sensitivity object around a classification, enough for print().
fake_sensitivity = function(cls, multipliers, rhat_max = rep(1, length(multipliers)),
                            s_hat = NA_real_) {
  structure(
    list(
      edges = cls$edges,
      grid = data.frame(
        multiplier = multipliers, original_fit = multipliers == 1,
        iter = 1000L, warmup = 1000L, rhat_max = rhat_max
      ),
      scale_multipliers = multipliers, chosen_scale = 1,
      verdict_by_scale = cls$verdict_by_scale,
      preferred_scale = list(s_hat = s_hat, lo = s_hat / 1.3, hi = s_hat * 1.3),
      unit = list(
        noun = "edge", nouns = "edges",
        headline = "are the edge verdicts robust to the slab scale?"
      ),
      diagonal_eta = NA_real_, evidence_threshold = 10,
      refit_sampler = "nuts", warm = TRUE, runtime_seconds = 3
    ),
    class = "bgms_prior_sensitivity"
  )
}


test_that("the classification compares each verdict change with the two fits' Monte Carlo error", {
  lthr = log(10)
  multipliers = c(0.5, 1, 2)
  # Five edges at the chosen scale: presence, undecided, undecided, absence,
  # presence.
  s0 = fake_stats(c(3.0, 1.0, 2.0, -3.0, 2.6), mcse = c(0.2, 0.1, 0.2, 0.2, 0.3))
  # 0.5x: edge 2 becomes presence with a change of 1.5 against an error of
  # 2 * sqrt(0.1^2 + 0.1^2) = 0.28 (beyond); edge 5 stays presence.
  s_lo = fake_stats(c(3.2, 2.5, 2.1, -2.9, 2.8), mcse = c(0.2, 0.1, 0.2, 0.2, 0.3))
  # 2x: edge 3 becomes presence by a change of 0.35 against an error of
  # 2 * sqrt(0.2^2 + 0.2^2) = 0.57 (within); edge 5 becomes undecided by a
  # change of 0.4 against 2 * sqrt(0.3^2 + 0.3^2) = 0.85 (within).
  s_hi = fake_stats(c(2.9, 1.1, 2.35, -3.1, 2.2), mcse = c(0.2, 0.1, 0.2, 0.2, 0.3))

  cls = sensitivity_classification(list(s_lo, s0, s_hi), multipliers, 2L, lthr)
  e = cls$edges
  expect_equal(e$mover, c(
    "stable", "beyond-monte-carlo-error", "within-monte-carlo-error",
    "stable", "within-monte-carlo-error"
  ))
  expect_equal(e$moved_at, c(NA, 0.5, 2, NA, 2))
  # The bracket is the contiguous run of multipliers around 1x with the chosen
  # verdict.
  expect_equal(e$stability_lower, c(0.5, 1, 0.5, 0.5, 0.5))
  expect_equal(e$stability_upper, c(2, 2, 1, 2, 1))
  expect_equal(e$chosen_scale_verdict, s0$verdict)
  expect_equal(e$verdict_x0.5, s_lo$verdict)
  expect_equal(e$log_bf_x2, s_hi$lbf)
  expect_equal(e$mcse_x1, s0$mcse_lbf)
  expect_equal(dim(cls$verdict_by_scale), c(3L, 5L))
})


test_that("sufficiency at the chosen scale is the only not-certifiable criterion", {
  lthr = log(10)
  multipliers = c(0.5, 1, 2)
  # Edge 1 sits within 2 * MCSE of the presence threshold at 1x; edge 2 is
  # within error at 2x only; edge 3 has disagreeing chains at 1x.
  s0 = fake_stats(c(2.2, 1.0, 1.5), mcse = c(0.1, 0.1, 0.1), unanimous = c(TRUE, TRUE, FALSE))
  s_lo = fake_stats(c(2.2, 1.0, 1.5), mcse = c(0.1, 0.1, 0.1))
  s_hi = fake_stats(c(2.2, 2.5, 1.5), mcse = c(0.1, 0.1, 0.1))
  cls = sensitivity_classification(list(s_lo, s0, s_hi), multipliers, 2L, lthr)
  expect_equal(cls$edges$insufficient, c(TRUE, FALSE, TRUE))
  expect_equal(cls$edges$insufficient_noisy, c(TRUE, FALSE, FALSE))
  expect_equal(cls$edges$insufficient_disagree, c(FALSE, FALSE, TRUE))
  # Edge 2 is classified on its own change, not on its nearness to a threshold
  # at the other scale.
  expect_equal(cls$edges$mover[2], "beyond-monte-carlo-error")
})


test_that("saturated and verdict-less indicators are handled", {
  lthr = log(10)
  multipliers = c(0.5, 1, 2)
  # Edge 1 never flips in the original fit (infinite log BF, no MCSE) and
  # becomes finite presence at 2x; edge 2 has no verdict.
  s0 = fake_stats(c(Inf, NA), mcse = c(NA, NA), zeroflip = c(TRUE, FALSE))
  s0$verdict[2] = NA_character_
  s_lo = fake_stats(c(Inf, NA), mcse = c(NA, NA), zeroflip = c(TRUE, FALSE))
  s_hi = fake_stats(c(4, NA), mcse = c(0.2, NA))
  cls = sensitivity_classification(list(s_lo, s0, s_hi), multipliers, 2L, lthr)
  expect_equal(cls$edges$mover, c("stable", NA))
  expect_true(cls$edges$saturated[1])
  expect_false(cls$edges$insufficient[1])
  expect_true(is.na(cls$edges$stability_lower[2]))
  # A verdict change away from a saturated indicator is beyond error.
  s_hi2 = fake_stats(c(1, NA), mcse = c(0.2, NA))
  cls2 = sensitivity_classification(list(s_lo, s0, s_hi2), multipliers, 2L, lthr)
  expect_equal(cls2$edges$mover[1], "beyond-monte-carlo-error")
})


test_that("print() reports the counts, the movers and high R-hat fits", {
  lthr = log(10)
  multipliers = c(0.5, 1, 2)
  s0 = fake_stats(c(3.0, 1.0, 2.0, -3.0, 2.6), mcse = c(0.2, 0.1, 0.2, 0.2, 0.3))
  s_lo = fake_stats(c(3.2, 2.5, 2.1, -2.9, 2.8), mcse = c(0.2, 0.1, 0.2, 0.2, 0.3))
  s_hi = fake_stats(c(2.9, 1.1, 2.35, -3.1, 2.2), mcse = c(0.2, 0.1, 0.2, 0.2, 0.3))
  cls = sensitivity_classification(list(s_lo, s0, s_hi), multipliers, 2L, lthr)

  ps = fake_sensitivity(cls, multipliers, rhat_max = c(1.03, 1.005, 1.02), s_hat = 0.2)
  out = gsub("\\s+", " ", paste(utils::capture.output(print(ps)), collapse = " "))
  expect_match(out, "fits at 0.5x, 1x and 2x the chosen scale")
  # Edges 3 and 5 lie within 2 * MCSE of the presence threshold in the
  # original fit, so they are not certifiable rather than movers.
  expect_match(out, "verdict held at every scale 2")
  expect_match(out, "changed, within Monte Carlo error 0")
  expect_match(out, "changed, beyond Monte Carlo error 1")
  expect_match(out, "not certifiable at the chosen scale 2")
  expect_match(out, "2 edges are not certifiable at the chosen scale: e3, e5")
  expect_match(out, "within Monte Carlo error of a threshold, so a rerun")
  expect_match(out, "1 edge changed verdict beyond Monte Carlo error")
  expect_match(out, "e2 presence undecided undecided 0.5x")
  expect_match(out, "Largest split R-hat above 1.01: 0.5x fit 1.03 and 2x fit 1.02")
  expect_match(out, "the slab is 5 times that")
  expect_match(out, "toward absence")
  expect_match(out, "2 nuts refits")

  # Every fit at or below 1.01, and a scale within a factor two: neither
  # sentence appears, and the verdict table has every scale.
  ps2 = fake_sensitivity(cls, multipliers, s_hat = 0.8)
  out2 = gsub("\\s+", " ", paste(utils::capture.output(print(ps2)), collapse = " "))
  expect_false(grepl("Largest split R-hat", out2, fixed = TRUE))
  expect_match(out2, "the slab is 1.2 times that\\.")
  expect_false(grepl("toward", out2, fixed = TRUE))
  expect_match(out2, "0.5x 1x 2x")
})


# ------------------------------------------------------------------------------
# Engine pieces
# ------------------------------------------------------------------------------

test_that("diagonal_eta reads the standardized rate and refuses a raw one", {
  standardized = list(prior = list(
    pairwise_scale = 2, scale_rate = 0.5, scale_eta = 1
  ))
  raw_frame = list(prior = list(
    pairwise_scale = 2, scale_rate = 1, scale_eta = NA_real_
  ))
  discrete = list(prior = list(pairwise_scale = 1, scale_rate = NA_real_))

  expect_equal(diagonal_eta(standardized), 1)
  expect_true(is.na(diagonal_eta(discrete)))
  expect_error(diagonal_eta(raw_frame), "raw rate")
})


test_that("finite_reduce abstains on an all-NA diagnostic without warning", {
  expect_equal(finite_reduce(c(NA_real_, NA_real_), min, NA_real_), NA_real_)
  expect_equal(finite_reduce(numeric(0), max, 0), 0)
  expect_equal(finite_reduce(c(NA_real_, 2, 5), max, 0), 5)
  expect_silent(finite_reduce(rep(NA_real_, 3L), min, Inf))
  expect_equal(finite_reduce(c(-Inf, Inf), min, Inf), Inf)
})


test_that("a warm-start list with the wrong length errors", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:5],
    chains = 2, iter = 40, warmup = 40, seed = 3,
    display_progress = "none"
  )
  ws = extract_warm_state(fit)
  spec = get_fit_spec(fit)
  spec$sampler$warmup = 100L
  spec$sampler$iter = 100L
  spec$sampler$cores = 2L
  spec$sampler$display_progress = "none"
  spec$sampler$progress_type = 0L
  spec$sampler$progress_callback = NULL
  # three parameter vectors but only two chains
  spec$initial_state = list(parameters = c(ws$parameters, ws$parameters[1]))
  expect_error(run_sampler(spec), "one entry per chain")
})


# ------------------------------------------------------------------------------
# The check on fits
# ------------------------------------------------------------------------------

test_that("prior_sensitivity_check needs edge selection", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:5],
    edge_selection = FALSE,
    iter = 300, warmup = 300, chains = 1, seed = 7,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_error(prior_sensitivity_check(fit), "edge selection")
})


test_that("scale multipliers must be positive and include one other than 1", {
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:4],
    chains = 1, iter = 200, warmup = 200, seed = 3,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  expect_error(prior_sensitivity_check(fit, scale_multipliers = 1), "other than 1")
  expect_error(prior_sensitivity_check(fit, scale_multipliers = c(-1, 2)), "positive")
  expect_error(prior_sensitivity_check(fit, evidence_threshold = 1), "greater than 1")
})


test_that("a fit whose swept prior carries no scale errors before any refit", {
  # A beta-prime slab is scale-free, so bgm() records pairwise_scale = NA and
  # there is no scale to multiply.
  data("Wenchuan", package = "bgms")
  fit = suppressWarnings(bgm(Wenchuan[, 1:4],
    interaction_prior = beta_prime_prior(),
    iter = 200, warmup = 200, chains = 1, seed = 11,
    update_method = "adaptive-metropolis", display_progress = "none"
  ))
  expect_true(is.na(get_fit_spec(fit)$prior$pairwise_scale))
  t0 = Sys.time()
  expect_error(
    prior_sensitivity_check(fit),
    "no usable pairwise_scale"
  )
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 5)
})


test_that("prior_sensitivity_check builds the object from the fits' own statistics", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:6],
    chains = 2, iter = 400, warmup = 300, seed = 21,
    display_progress = "none"
  )
  ps = suppressMessages(prior_sensitivity_check(fit, seed = 21))

  expect_s3_class(ps, "bgms_prior_sensitivity")
  expect_equal(ps$scale_multipliers, c(0.5, 1, 2))
  # One grid row per fit (the 1x row is the original fit), every refit at the
  # original fit's run length.
  expect_equal(nrow(ps$grid), 3L)
  expect_equal(sum(ps$grid$original_fit), 1L)
  expect_equal(ps$grid$seconds[ps$grid$original_fit], 0)
  expect_true(all(ps$grid$iter == 400L))
  expect_true(all(ps$grid$warmup == 300L))
  expect_true(all(c("rhat_max", "ess_min", "inclusion_ess_min") %in% names(ps$grid)))
  expect_equal(nrow(ps$edges), 15L) # choose(6, 2)

  # The chosen-scale columns are the original fit's own Rao-Blackwellized
  # analysis: the reported PIP and its exact prior-odds transform.
  rb = extract_posterior_inclusion_probabilities(fit)
  rbv = unname(rb[upper.tri(rb)][order_upper_tri_rowmajor(6)])
  po = ps$edges$prior_inclusion_probability /
    (1 - ps$edges$prior_inclusion_probability)
  expect_equal(ps$edges$chosen_scale_pip, rbv, tolerance = 1e-10)
  pc = ps$edges$chosen_scale_pip
  expect_equal(
    ps$edges$chosen_scale_log_bf,
    log((pc / (1 - pc)) / po),
    tolerance = 1e-12
  )
  expect_equal(ps$edges$verdict_x1, ps$edges$chosen_scale_verdict)
  expect_equal(ps$edges$log_bf_x1, ps$edges$chosen_scale_log_bf)
  expect_true(all(ps$edges$mover %in%
    c("stable", "within-monte-carlo-error", "beyond-monte-carlo-error")))
  expect_true(all(is.na(ps$edges$moved_at[ps$edges$mover == "stable"])))
  expect_true(all(ps$edges$stability_lower %in% ps$scale_multipliers))
  expect_true(all(ps$edges$stability_upper %in% ps$scale_multipliers))
  expect_true(all(paste0(c("verdict_x", "log_bf_x", "mcse_x"), rep(ps$scale_multipliers, each = 3)) %in%
    names(ps$edges)))
  expect_equal(dim(ps$verdict_by_scale), c(3L, 15L))
  expect_true(is.finite(ps$preferred_scale$s_hat))
  expect_true(ps$warm)
  expect_output(print(ps), "are the edge verdicts robust")
})


test_that("the chosen-scale verdict uses the per-edge prior odds, not 1/2", {
  skip_on_cran()
  skip_unless_slow()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:5],
    edge_prior = bernoulli_prior(0.2),
    chains = 2, iter = 500, warmup = 400, seed = 5, display_progress = "none"
  )
  ps = suppressMessages(prior_sensitivity_check(fit, scale_multipliers = 2, seed = 5))
  expect_true(all(abs(ps$edges$prior_inclusion_probability - 0.2) < 1e-8))
  p = ps$edges$chosen_scale_pip
  expected = log((p / (1 - p)) / (0.2 / 0.8))
  expect_equal(ps$edges$chosen_scale_log_bf, expected, tolerance = 1e-8)
})


test_that("a warm-started refit at the chosen scale agrees with a cold one", {
  skip_on_cran()
  skip_unless_certification()
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:8],
    chains = 4, iter = 2000, warmup = 1000, seed = 11,
    display_progress = "none"
  )
  ws = extract_warm_state(fit)
  s0 = get_fit_spec(fit)$prior$pairwise_scale
  pipv = function(f) {
    m = extract_posterior_inclusion_probabilities(f)
    m[upper.tri(m)]
  }
  warm = suppressMessages(refit_at_scale(fit, s0, ws,
    warmup = 1000, iter = 2000,
    seed = 101, cores = 4L, sampler = "nuts"
  ))
  cold = suppressMessages(refit_at_scale(fit, s0,
    warm_state = NULL,
    warmup = 1000, iter = 2000, seed = 202, cores = 4L, sampler = "nuts"
  ))
  expect_lt(stats::median(abs(pipv(warm) - pipv(cold))), 0.03)
})


test_that("prior_sensitivity_check runs for GGM and mixed fits (cold refits)", {
  skip_on_cran()
  skip_unless_certification()
  set.seed(32)
  xg = matrix(rnorm(180 * 5), 180, 5)
  fg = bgm(xg,
    variable_type = "continuous",
    iter = 300, warmup = 250, chains = 2, seed = 8,
    precision_graph_prior = "joint",
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  psg = suppressMessages(prior_sensitivity_check(fg, scale_multipliers = 2, seed = 8))
  expect_s3_class(psg, "bgms_prior_sensitivity")
  expect_equal(psg$model_type, "ggm")
  expect_false(psg$warm)
  # A default GGM gives its diagonal prior as eta = 1; the refits hold it fixed
  # and the report says so.
  expect_equal(psg$diagonal_eta, 1)
  expect_match(
    paste(utils::capture.output(print(psg)), collapse = "\n"),
    "standardized rate eta"
  )
  expect_true(all(is.finite(psg$edges$log_bf_x2) | psg$edges$verdict_x2 %in% c("presence", "absence")))

  n = 180
  xm = data.frame(
    a = sample(0:2, n, TRUE), b = sample(0:2, n, TRUE),
    c = sample(0:1, n, TRUE), y1 = rnorm(n), y2 = rnorm(n)
  )
  vt = c("ordinal", "ordinal", "ordinal", "continuous", "continuous")
  fm = bgm(xm,
    variable_type = vt, iter = 300, warmup = 250, chains = 2, seed = 9,
    update_method = "adaptive-metropolis", display_progress = "none"
  )
  psm = suppressMessages(prior_sensitivity_check(fm, scale_multipliers = 2, seed = 9))
  expect_equal(psm$model_type, "mixed_mrf")
  expect_equal(psm$diagonal_eta, 1)
})


test_that("prior_sensitivity_check compares bgmCompare difference verdicts", {
  skip_on_cran()
  skip_unless_slow()
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[, 1:5], group_indicator = rep(1:2, length.out = nrow(Wenchuan)),
    iter = 800, warmup = 800, chains = 2, seed = 21,
    difference_selection = TRUE, display_progress = "none"
  )
  ps = suppressMessages(prior_sensitivity_check(fit, seed = 3))

  expect_s3_class(ps, "bgms_prior_sensitivity")
  nm = get_raw_samples(fit)$parameter_names$indicator
  # One row per difference indicator, main-effect differences included.
  expect_equal(ps$edges$edge, nm)
  # It sweeps the difference scale, not the interaction one.
  expect_equal(ps$unit$scale_field, "difference_scale")
  expect_equal(ps$chosen_scale, get_fit_spec(fit)$prior$difference_scale)
  # There is no data-preferred difference scale to compare against.
  expect_true(is.na(ps$preferred_scale$s_hat))

  out = paste(utils::capture.output(print(ps)), collapse = "\n")
  expect_match(out, "are the difference verdicts robust to the difference scale")
  # main_difference_selection is FALSE by default, so those indicators were
  # never updated; the report says so instead of counting them as undecided.
  expect_match(out, "never updated by the sampler and carry no verdict")
  expect_false(grepl("Note: the slab scale", out, fixed = TRUE))
})


test_that("prior_sensitivity_check needs difference selection", {
  skip_on_cran()
  data("Wenchuan", package = "bgms")
  fit = bgmCompare(
    x = Wenchuan[1:80, 1:4], group_indicator = rep(1:2, each = 40),
    iter = 150, warmup = 150, chains = 2, seed = 4,
    difference_selection = FALSE, display_progress = "none"
  )
  expect_error(prior_sensitivity_check(fit), "difference selection")
})


test_that("prior_sensitivity_check runs end to end on a single-indicator fit", {
  skip_on_cran()
  # A two-variable network has one edge, so apply()/vapply() over a one-column
  # matrix collapse to a vector; the shapes must survive that.
  data("Wenchuan", package = "bgms")
  fit = bgm(Wenchuan[, 1:2],
    chains = 2, iter = 300, warmup = 300, seed = 31,
    display_progress = "none"
  )
  ps = suppressMessages(prior_sensitivity_check(fit, seed = 31))
  expect_s3_class(ps, "bgms_prior_sensitivity")
  expect_equal(nrow(ps$edges), 1L)
  expect_equal(dim(ps$verdict_by_scale), c(3L, 1L))
  expect_true(all(paste0("verdict_x", ps$scale_multipliers) %in% names(ps$edges)))
  expect_output(print(ps), "are the edge verdicts robust")
})
