# --------------------------------------------------------------------------- #
# Cluster counts under the SBM prior: occupied clusters T and components B.
#
# Prior quantities are checked against independent derivations (enumeration
# over partitions, the Lah numbers, the closed form of P(T | B) at eps = 1, and
# the sampler's own partition coefficients); posterior quantities against
# arithmetic on known occupied-count posteriors; and the reconstruction of T
# from legacy summaries against round trips under each recognized convention.
# --------------------------------------------------------------------------- #

counts_enumerate_partitions = function(p) {
  rows = list(1L)
  for(k in seq_len(p - 1L)) {
    rows = unlist(lapply(rows, function(z) {
      lapply(seq_len(max(z) + 1L), function(j) c(z, j))
    }), recursive = FALSE)
  }
  rows
}

# sum over partitions with t blocks of prod_c eps^(|c|), by enumeration.
counts_partition_weights = function(p, eps) {
  w = numeric(p)
  for(z in counts_enumerate_partitions(p)) {
    w[max(z)] = w[max(z)] + exp(sum(lgamma(tabulate(z) + eps) - lgamma(eps)))
  }
  w
}

# P(T = t | B = b) = b_(t) / (b eps)^(p) * S_eps(p, t), by enumeration.
counts_t_given_b = function(b, p, eps) {
  S = counts_partition_weights(p, eps)
  vapply(seq_len(p), function(t) {
    if(b < t) {
      return(0)
    }
    exp(lgamma(b + 1) - lgamma(b - t + 1) - lgamma(b * eps + p) + lgamma(b * eps)) * S[t]
  }, numeric(1))
}

# Allocation draws whose numbers of occupied clusters follow `counts`.
allocations_with_counts = function(counts, p) {
  t = rep(seq_along(counts), counts)
  t(vapply(t, function(k) c(seq_len(k), rep(1L, p - k)), integer(p)))
}


test_that("occupied-cluster prior is normalized and matches enumeration", {
  for(cfg in list(c(1, 1), c(3, 1), c(1, 0.5), c(0.5, 2))) {
    for(p in 2:6) {
      k = bgms:::sbm_count_kernel(p, cfg[2], cfg[1])
      q = exp(k$log_q)
      expect_equal(sum(q), 1, tolerance = 1e-12)
      b = seq_len(400)
      q_enum = Reduce(`+`, lapply(b, function(bb) {
        dpois(bb - 1, cfg[1]) * counts_t_given_b(bb, p, cfg[2])
      }))
      expect_equal(q, q_enum, tolerance = 1e-10)
    }
  }
})

test_that("partition weights at eps = 1 are the Lah numbers", {
  p = 9
  lah = vapply(seq_len(p), function(t) {
    choose(p - 1, t - 1) * factorial(p) / factorial(t)
  }, numeric(1))
  expect_equal(exp(bgms:::sbm_log_partition_counts(p, 1)), lah, tolerance = 1e-12)
})

test_that("components given occupied clusters match the closed form at eps = 1", {
  p = 6
  lambda = 2
  k = bgms:::sbm_count_kernel(p, 1, lambda)
  b = seq_len(k$b_max)
  q = exp(k$log_q)
  closed = vapply(seq_len(p), function(t) {
    h = ifelse(b >= t, exp(lchoose(b, t) + lchoose(p - 1, t - 1) - lchoose(p + b - 1, p)), 0)
    h * dpois(b - 1, lambda) / q[t]
  }, numeric(length(b)))
  expect_equal(exp(k$log_b_given_t), closed, tolerance = 1e-12)
  expect_equal(colSums(exp(k$log_b_given_t)), rep(1, p), tolerance = 1e-12)
})

test_that("partition coefficients agree with the sampler's, including large lambda", {
  # The sampler used to stop at B = 501, which is too short past lambda ~ 200.
  for(cfg in list(c(1, 1, 6), c(3, 0.5, 10), c(400, 1, 8))) {
    k = bgms:::sbm_count_kernel(cfg[3], cfg[2], cfg[1])
    expect_equal(k$log_vp, as.numeric(bgms:::compute_Vn_mfm_sbm(cfg[3], cfg[2], cfg[3], cfg[1])),
      tolerance = 1e-12
    )
  }
})

test_that("partition coefficients agree with the sampler's for many clusters", {
  skip_on_cran()
  # B >= T, so the sum must reach past T as well as past lambda. A fixed
  # cut-off gave -Inf for T > 501 at p = 520, and was off by 1e-3 in log V(300)
  # at lambda = 250 with a small Dirichlet concentration.
  for(cfg in list(c(1, 1, 520), c(250, 0.01, 300))) {
    k = bgms:::sbm_count_kernel(cfg[3], cfg[2], cfg[1])
    expect_equal(k$log_vp, as.numeric(bgms:::compute_Vn_mfm_sbm(cfg[3], cfg[2], cfg[3], cfg[1])),
      tolerance = 1e-12
    )
  }
})

test_that("the component posterior keeps its mass beyond p", {
  # p = 3 and lambda = 5: most of the prior and posterior mass of B lies above
  # the number of variables, which the legacy summary cut off.
  p = 3
  z = allocations_with_counts(c(100, 300, 600), p)
  x = sbm_counts(allocations = z, lambda = 5, dirichlet_alpha = 1, version = "0.2.1.0")
  comp = x$components
  expect_gt(nrow(comp), p)
  expect_equal(comp$prior, dpois(comp$count - 1, 5), tolerance = 1e-14)
  expect_gt(sum(comp$posterior[comp$count > p]), 0.5)
  expect_lt(max(x$tail), 1e-10)
  expect_equal(sum(comp$posterior) + x$tail[["posterior"]], 1, tolerance = 1e-10)
  expect_equal(x$occupied$posterior, c(0.1, 0.3, 0.6))

  # Independently: P(B = b | X) = sum_t P(T = t | B = b) p_B(b) / q_t r_t.
  q = Reduce(`+`, lapply(seq_len(400), function(b) dpois(b - 1, 5) * counts_t_given_b(b, p, 1)))
  post = vapply(comp$count, function(b) {
    sum(counts_t_given_b(b, p, 1) * dpois(b - 1, 5) / q * c(0.1, 0.3, 0.6))
  }, numeric(1))
  expect_equal(comp$posterior, post, tolerance = 1e-10)
})

test_that("Bayes factors match independently derived odds", {
  p = 4
  lambda = 1.5
  eps = 0.7
  r = c(0.15, 0.45, 0.3, 0.1)
  z = allocations_with_counts(round(r * 1000), p)
  x = sbm_counts(allocations = z, lambda = lambda, dirichlet_alpha = eps, version = "0.2.1.0")

  q = Reduce(`+`, lapply(seq_len(400), function(b) dpois(b - 1, lambda) * counts_t_given_b(b, p, eps)))
  post_b = function(b) sum(counts_t_given_b(b, p, eps) * dpois(b - 1, lambda) / q * r)

  expect_equal(sbm_count_bf(x, "occupied")$bf,
    ((1 - r[1]) / r[1]) / ((1 - q[1]) / q[1]),
    tolerance = 1e-10
  )
  expect_equal(sbm_count_bf(x, "occupied", "point", 2, 3)$bf,
    (r[2] / r[3]) / (q[2] / q[3]),
    tolerance = 1e-10
  )
  expect_equal(sbm_count_bf(x, "components")$bf,
    ((1 - post_b(1)) / post_b(1)) / expm1(lambda),
    tolerance = 1e-10
  )
  expect_equal(sbm_count_bf(x, "components", "point", 2, 5)$bf,
    (post_b(2) / post_b(5)) / (dpois(1, lambda) / dpois(4, lambda)),
    tolerance = 1e-10
  )
})

test_that("the binary component factor equals a + (1 - a) times the occupied factor", {
  p = 5
  lambda = 1
  k = bgms:::sbm_count_kernel(p, 1, lambda)
  p_b = dpois(seq_len(k$b_max) - 1, lambda)
  t1_given_b = exp(k$log_b_given_t[, 1] + k$log_q[1]) / p_b
  a = sum(t1_given_b[-1] * p_b[-1]) / sum(p_b[-1])
  for(counts in list(c(10, 50, 30, 8, 2), c(70, 20, 8, 2, 0), c(1, 9, 60, 25, 5))) {
    x = sbm_counts(
      allocations = allocations_with_counts(counts, p),
      lambda = lambda, dirichlet_alpha = 1, version = "0.2.1.0"
    )
    expect_equal(sbm_count_bf(x, "components")$bf,
      a + (1 - a) * sbm_count_bf(x, "occupied")$bf,
      tolerance = 1e-10
    )
  }
})

test_that("printed tables and Bayes factors come from the same numbers", {
  fit = get_bgms_fit_sbm()
  x = extract_sbm_counts(fit)
  comp = x$components
  occ = x$occupied
  for(pair in list(c(1, 2), c(2, 3))) {
    expect_equal(
      sbm_count_bf(x, "components", "point", pair[1], pair[2])$bf,
      (comp$posterior[pair[1]] / comp$posterior[pair[2]]) /
        (comp$prior[pair[1]] / comp$prior[pair[2]])
    )
  }
  if(all(occ$posterior[1:2] > 0)) {
    expect_equal(
      sbm_count_bf(x, "occupied", "point", 1, 2)$bf,
      (occ$posterior[1] / occ$posterior[2]) / (occ$prior[1] / occ$prior[2])
    )
  }
  out = capture.output(print(x))
  expect_true(any(grepl("T > 1 vs T = 1", out, fixed = TRUE)))
  expect_true(any(grepl(bgms:::format_bf(sbm_count_bf(x, "components")), out, fixed = TRUE)))
})

test_that("unvisited counts give bounds, never infinite evidence", {
  p = 4
  z = allocations_with_counts(c(0, 600, 300, 100), p)
  x = sbm_counts(allocations = z, lambda = 1, dirichlet_alpha = 1, version = "0.2.1.0")
  for(target in c("occupied", "components")) {
    bf = sbm_count_bf(x, target)
    expect_identical(bf$status, "bound")
    expect_true(is.na(bf$bf))
    expect_identical(bf$bound_direction, ">")
    expect_true(is.finite(bf$bound) && bf$bound > 0)
    expect_match(bf$explanation, "run the chains longer")
  }
  expect_identical(sbm_count_bf(x, "occupied", "point", 3, 1)$bound_direction, ">")
  expect_identical(sbm_count_bf(x, "occupied", "point", 1, 3)$bound_direction, "<")
  y = sbm_counts(
    allocations = allocations_with_counts(c(0, 1000, 0, 0), p),
    lambda = 1, dirichlet_alpha = 1, version = "0.2.1.0"
  )
  expect_identical(sbm_count_bf(y, "occupied", "point", 3, 4)$status, "not estimable")
  expect_error(sbm_count_bf(x), "Choose target")
})

test_that("legacy summaries are inverted only under recognized conventions", {
  p = 6
  r = c(0.2, 0.4, 0.25, 0.1, 0.04, 0.01)
  k = bgms:::sbm_count_kernel(p, 1, 1)
  for(cv in list(c("shifted", "0.2.0.0"), c("zero_truncated", "0.1.6.3"))) {
    reported = as.numeric(bgms:::sbm_legacy_summary_matrix(k, cv[1]) %*% r)
    x = sbm_counts(
      num_blocks = data.frame(probability = reported),
      lambda = 1, dirichlet_alpha = 1, version = cv[2]
    )
    expect_equal(x$occupied$posterior, r, tolerance = 1e-12)
    expect_identical(x$source, cv[1])
  }
  expect_error(
    sbm_counts(num_blocks = reported, lambda = 1, dirichlet_alpha = 1, version = "0.1.5.0"),
    "different prior on the number of components"
  )
  expect_error(
    sbm_counts(num_blocks = reported, lambda = 1, dirichlet_alpha = 1),
    "cannot be established automatically"
  )
})

test_that("a fit's allocations and its stored summary give the same occupied counts", {
  fit = get_bgms_fit_sbm()
  from_draws = extract_sbm_counts(fit)
  from_summary = sbm_counts(
    num_blocks = fit$posterior_num_blocks,
    lambda = extract_arguments(fit)$lambda,
    dirichlet_alpha = extract_arguments(fit)$dirichlet_alpha,
    version = extract_arguments(fit)$version
  )
  expect_s3_class(from_draws, "bgms_sbm_counts")
  expect_identical(from_draws$source, "allocations")
  expect_equal(from_summary$occupied$posterior, from_draws$occupied$posterior,
    tolerance = 1e-10
  )
  expect_identical(from_draws$provenance$sampler, "fixed")
  expect_identical(
    sum(from_draws$chain_counts),
    sum(vapply(fit$raw_samples$allocations, nrow, integer(1)))
  )
  out = capture.output(print(summary(fit)))
  expect_true(any(grepl("T > 1 vs T = 1", out, fixed = TRUE)))
  expect_true(any(grepl("Estimated node allocations", out, fixed = TRUE)))
})

test_that("bgmCompare fits with an sbm difference prior give cluster counts", {
  skip_on_cran()
  set.seed(4)
  x = matrix(sample(0:2, 120 * 4, TRUE), 120, 4)
  fit = suppressWarnings(bgmCompare(x,
    group_indicator = rep(1:2, each = 60), variable_type = "ordinal",
    difference_prior = sbm_prior(), update_method = "adaptive-metropolis",
    iter = 100, warmup = 100, chains = 1, seed = 4,
    display_progress = "none", verbose = FALSE
  ))
  counts = extract_sbm_counts(fit)
  expect_s3_class(counts, "bgms_sbm_counts")
  expect_equal(sum(counts$occupied$posterior), 1)
})

test_that("provenance messages match the problem and appear once per session", {
  shown = bgms:::.sbm_messages_shown
  reset = function() rm(list = ls(shown), envir = shown)
  reset()
  on.exit(reset(), add = TRUE)

  old = list(version = package_version("0.2.0.0"), edge_prior = "Stochastic-Block")
  shown_text = capture.output(bgms:::sbm_provenance_message(old), type = "message")
  expect_match(paste(shown_text, collapse = " "), "resummarizing cannot repair the draws")
  expect_silent(bgms:::sbm_provenance_message(old))
  expect_identical(bgms:::sbm_fit_provenance(old)$sampler, "affected")
  expect_identical(bgms:::sbm_fit_provenance(old)$summary_convention, "shifted")
  expect_identical(
    bgms:::sbm_fit_provenance(list(version = "0.1.6.3"))$summary_convention,
    "zero_truncated"
  )

  expect_message(
    bgms:::sbm_provenance_message(list(edge_prior = "Stochastic-Block")),
    "cannot be established automatically"
  )

  current = list(version = package_version("0.2.1.0"), sbm_sampler_fix = TRUE)
  expect_silent(bgms:::sbm_provenance_message(current))
  # A caller can muffle one kind by class without using up its one showing.
  withCallingHandlers(
    bgms:::sbm_provenance_message(current, legacy_summary = TRUE),
    bgms_sbm_legacy_message = function(m) invokeRestart("muffleMessage")
  )
  expect_message(bgms:::sbm_provenance_message(current, legacy_summary = TRUE),
    "extract_sbm_counts",
    class = "bgms_sbm_legacy_message"
  )

  # knitr captures the message into the document and then muffles it; that
  # counts as shown, so a document repeats it at most once.
  reset()
  withr::with_options(list(knitr.in.progress = TRUE), {
    withCallingHandlers(
      bgms:::sbm_provenance_message(current, legacy_summary = TRUE),
      message = function(m) invokeRestart("muffleMessage")
    )
    expect_silent(bgms:::sbm_provenance_message(current, legacy_summary = TRUE))
  })

  reset()
  joint = list(version = package_version("0.2.0.0"), precision_graph_prior = "joint")
  expect_identical(bgms:::sbm_fit_provenance(joint)$sampler, "not assessed")
  expect_silent(bgms:::sbm_provenance_message(joint, legacy_summary = TRUE))
})

test_that("every fit records its SBM provenance", {
  a = extract_arguments(get_bgms_fit_sbm())
  expect_identical(a$sbm_prior_family, "shifted Poisson")
  expect_true(isTRUE(a$sbm_sampler_fix))
  expect_identical(a$sbm_summary_convention, "shifted")
  expect_true(is.na(a$build_commit) || grepl("^[0-9a-f]{40}$", a$build_commit))
})
