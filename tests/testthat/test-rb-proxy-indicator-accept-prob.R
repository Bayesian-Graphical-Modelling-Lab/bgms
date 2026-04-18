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
