# ==============================================================================
# RB-proxy: indicator_accept_prob exposure
# ==============================================================================
#
# Verifies that the new $raw_samples$indicator_accept_prob field is populated
# when edge_selection = TRUE, has the expected per-chain matrix shape, lives
# in [0, 1], and is consistent with the indicator chain (an iteration in which
# an edge flipped must have alpha > 0).
# ==============================================================================

test_that("indicator_accept_prob is exposed and well-formed for OMRF", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[1:50, 1:4],
    iter = 40, warmup = 80, chains = 2,
    seed = 20260417,
    edge_selection = TRUE,
    display_progress = "none"
  )

  raw = fit$raw_samples
  expect_true(!is.null(raw$indicator))
  expect_true(!is.null(raw$indicator_accept_prob))
  expect_length(raw$indicator_accept_prob, length(raw$indicator))

  for(c in seq_along(raw$indicator_accept_prob)) {
    ind   = raw$indicator[[c]]
    alpha = raw$indicator_accept_prob[[c]]

    # Same per-chain shape as the indicator matrix.
    expect_equal(dim(alpha), dim(ind))

    # Valid probability range.
    expect_true(all(is.finite(alpha)))
    expect_true(all(alpha >= 0 - 1e-12))
    expect_true(all(alpha <= 1 + 1e-12))

    # Sanity: any iteration where an edge actually flipped must have had
    # a strictly positive acceptance probability.
    if(nrow(ind) >= 2) {
      flips = abs(diff(ind)) == 1L
      flipped_alpha = alpha[-1, , drop = FALSE][flips]
      if(length(flipped_alpha) > 0) {
        expect_true(all(flipped_alpha > 0))
      }
    }
  }
})


test_that("indicator_accept_prob aligns with empirical flip rate from state 1", {
  # For an interior edge with many 1 -> 0 transitions, the mean stored alpha
  # conditional on the prior state being 1 should approximate the empirical
  # flip rate from state 1. This is the RB-proxy alignment invariant.
  # Before the row-major indexing fix for GGM, this alignment was broken
  # because alpha values were written at column-major slots while the output
  # vector was read in row-major order, scrambling edge labels and collapsing
  # some slots under collisions.
  data("Wenchuan", package = "bgms")
  fit <- bgm(
    Wenchuan[, 1:5],
    iter = 2000, warmup = 500, chains = 2,
    seed = 20260417,
    edge_selection = TRUE,
    display_progress = "none"
  )

  raw <- fit$raw_samples
  n_chains <- length(raw$indicator)
  n_edges <- ncol(raw$indicator[[1]])

  flip_rate <- numeric(n_edges)
  mean_alpha <- numeric(n_edges)
  n_at_1 <- integer(n_edges)

  for(e in seq_len(n_edges)) {
    n1 <- 0L; flips <- 0L; alphas_at_1 <- numeric(0)
    for(c in seq_len(n_chains)) {
      ind_e   <- raw$indicator[[c]][, e]
      alpha_e <- raw$indicator_accept_prob[[c]][, e]
      if(length(ind_e) < 2) next
      prior_state <- ind_e[-length(ind_e)]
      next_state  <- ind_e[-1]
      alpha_next  <- alpha_e[-1]
      at_1 <- prior_state == 1
      n1 <- n1 + sum(at_1)
      flips <- flips + sum(at_1 & next_state == 0)
      alphas_at_1 <- c(alphas_at_1, alpha_next[at_1])
    }
    n_at_1[e] <- n1
    if(n1 > 0) {
      flip_rate[e]  <- flips / n1
      mean_alpha[e] <- mean(alphas_at_1)
    }
  }

  # Require at least a handful of edges with enough prior-state-1 samples to
  # compare (an "interior edge with many transitions" in the task statement).
  active <- n_at_1 >= 50
  expect_true(sum(active) >= 2,
              info = "Need at least two edges with >=50 prior-state-1 samples")

  # Tight per-edge agreement between flip_rate and mean_alpha.
  abs_err <- abs(flip_rate[active] - mean_alpha[active])
  expect_true(max(abs_err) < 0.05,
              info = sprintf("max |flip_rate - mean_alpha| = %.4f for active edges",
                             max(abs_err)))

  # Aggregate correlation across all edges with any prior-state-1 activity.
  eligible <- n_at_1 > 0
  if(sum(eligible) >= 3 && sd(flip_rate[eligible]) > 0 &&
     sd(mean_alpha[eligible]) > 0) {
    r <- cor(flip_rate[eligible], mean_alpha[eligible])
    expect_true(r > 0.95,
                info = sprintf("cor(flip_rate, mean_alpha) = %.4f, expected > 0.95", r))
  }
})


test_that("indicator_accept_prob is absent when edge_selection = FALSE", {
  data("Wenchuan", package = "bgms")
  fit = bgm(
    Wenchuan[1:50, 1:4],
    iter = 20, warmup = 40, chains = 1,
    seed = 20260417,
    edge_selection = FALSE,
    display_progress = "none"
  )

  expect_null(fit$raw_samples$indicator_accept_prob)
})
