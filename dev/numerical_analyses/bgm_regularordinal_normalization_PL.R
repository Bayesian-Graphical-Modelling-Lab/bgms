################################################################################
# Reference: Numerical stability study for bounded vs. unbounded exponential sums
# Author: [Your Name]
# Date: [YYYY-MM-DD]
#
# Purpose:
#   Evaluate and compare four ways to compute the sum
#
#       S = 1 + Σ_{c=1..K} exp( m_c + (c+1)*r )
#
#   where r may vary widely. The goal is to identify numerically stable and
#   computationally efficient formulations for use in gradient calculations.
#
# Methods compared:
#   (1) direct        – naive computation using raw exp()
#   (2) bounded       – stabilized by subtracting a "bound" (i.e., scaled domain)
#   (3) preexp        – precomputes exp(m_c) and exp(r) to replace repeated calls
#   (4) preexp_bound  – preexp variant with the same "bound" scaling
#
# For each method, we compute both unscaled and scaled variants where relevant,
# and compare them against a high-precision MPFR reference.
#
# Key insight:
#   - For large negative r, preexp can lose precision (tiny multiplicative updates).
#   - For large positive r, bounded scaling avoids overflow.
#   - The combination (preexp + bound) gives the best general stability.
#
# Output:
#   - res: data frame with per-r results and relative errors
#   - Diagnostic plots and summary tables for numerical accuracy
################################################################################

library(Rmpfr)  # for arbitrary precision reference computations


################################################################################
# 1. Core comparison function
################################################################################
compare_all_methods <- function(K = 5,
                                r_vals = seq(-10, 10, length.out = 200),
                                m_vals = NULL,
                                mpfr_prec = 256) {
  # ---------------------------------------------------------------------------
  # Parameters:
  #   K          – number of categories (terms in the sum)
  #   r_vals     – vector of r values to evaluate over
  #   m_vals     – optional vector of m_c values; random if NULL
  #   mpfr_prec  – bits of precision for the high-precision reference
  #
  # Returns:
  #   A data.frame containing per-r computed values, reference values,
  #   relative errors, and failure flags.
  # ---------------------------------------------------------------------------

  if (is.null(m_vals)) m_vals <- runif(K, -1, 1)

  results <- data.frame(
    r = r_vals,
    direct = NA_real_,
    bounded = NA_real_,        # scaled-domain computation (exp(-bound) factor)
    preexp = NA_real_,
    preexp_bound = NA_real_,   # scaled-domain computation
    ref = NA_real_,            # unscaled MPFR reference
    ref_scaled = NA_real_,     # scaled reference
    err_direct = NA_real_,
    err_bounded = NA_real_,
    err_preexp = NA_real_,
    err_preexp_bound = NA_real_,
    ref_failed_unscaled = FALSE,
    ref_failed_scaled   = FALSE
  )

  # Loop over all r-values
  for (i in seq_along(r_vals)) {
    r <- r_vals[i]
    bound <- K * r  # can be unclipped; use max(0, K*r) for the clipped version

    # --- (0) High-precision MPFR reference -----------------------------------
    r_mp <- mpfr(r, precBits = mpfr_prec)
    m_mp <- mpfr(m_vals, precBits = mpfr_prec)
    b_mp <- mpfr(bound, precBits = mpfr_prec)

    ref_unscaled_mp <- 1 + sum(exp(m_mp + (1:K) * r_mp))
    ref_scaled_mp   <- exp(-b_mp) * ref_unscaled_mp

    # Convert to doubles for inspection
    ref_unscaled_num <- asNumeric(ref_unscaled_mp)
    ref_scaled_num   <- asNumeric(ref_scaled_mp)
    results$ref_failed_unscaled[i] <- !is.finite(ref_unscaled_num)
    results$ref_failed_scaled[i]   <- !is.finite(ref_scaled_num)
    results$ref[i]                 <- if (is.finite(ref_unscaled_num)) ref_unscaled_num else NA_real_
    results$ref_scaled[i]          <- if (is.finite(ref_scaled_num))   ref_scaled_num   else NA_real_

    # --- (1) Direct exponential sum (unscaled) -------------------------------
    results$direct[i] <- 1 + sum(exp(m_vals + (1:K) * r))

    # --- (2) Current bounded implementation (scaled) -------------------------
    eB <- exp(-bound)
    results$bounded[i] <- eB + sum(exp(m_vals + (1:K) * r - bound))

    # --- (3) Precomputed exp only (unscaled) ---------------------------------
    exp_r <- exp(r)
    exp_m <- exp(m_vals)
    powE  <- exp_r
    S_pre <- 1.0
    for (c in 1:K) {
      S_pre <- S_pre + exp_m[c] * powE
      powE  <- powE * exp_r
    }
    results$preexp[i] <- S_pre

    # --- (4) Precomputed exp + bound scaling (scaled) ------------------------
    exp_r <- exp(r)
    exp_m <- exp(m_vals)
    powE  <- exp_r
    S_preB <- eB
    for (c in 1:K) {
      S_preB <- S_preB + exp_m[c] * powE * eB
      powE   <- powE * exp_r
    }
    results$preexp_bound[i] <- S_preB

    # --- (5) Relative errors vs references -----------------------------------
    # Unscaled methods
    for (m in c("direct", "preexp")) {
      val <- results[[m]][i]
      if (is.finite(val)) {
        val_mp <- mpfr(val, precBits = mpfr_prec)
        err_mp <- abs((val_mp - ref_unscaled_mp) / ref_unscaled_mp)
        results[[paste0("err_", m)]][i] <- asNumeric(err_mp)
      }
    }

    # Scaled methods
    for (m in c("bounded", "preexp_bound")) {
      val <- results[[m]][i]
      if (is.finite(val)) {
        val_mp <- mpfr(val, precBits = mpfr_prec)
        err_mp <- abs((val_mp - ref_scaled_mp) / ref_scaled_mp)
        results[[paste0("err_", m)]][i] <- asNumeric(err_mp)
      }
    }
  }

  msg_a <- mean(results$ref_failed_unscaled)
  msg_b <- mean(results$ref_failed_scaled)
  message(sprintf("Ref (unscaled) non-finite in %.1f%%; Ref (scaled) non-finite in %.1f%% of r-values",
                  100 * msg_a, 100 * msg_b))
  results
}


################################################################################
# 2. Plotting: log-scale accuracy with failure marking
################################################################################
plot_errors <- function(res) {
  err_cols <- c("err_bounded", "err_direct", "err_preexp", "err_preexp_bound")
  cols <- c("gray", "black", "red", "blue")
  names(cols) <- err_cols

  # Compute a robust ylim (1st–99th percentile)
  finite_vals <- unlist(res[err_cols])
  finite_vals <- finite_vals[is.finite(finite_vals) & finite_vals > 0]
  if (length(finite_vals) > 0) {
    lower <- quantile(finite_vals, 0.01, na.rm = TRUE)
    upper <- quantile(finite_vals, 0.99, na.rm = TRUE)
    ylim <- c(lower / 10, upper * 10)
  } else {
    ylim <- c(1e-20, 1e-12)
  }

  # Baseline curve: bounded
  plot(res$r, res$err_bounded, type = "l", log = "y",
       col = cols["err_bounded"], lwd = 2,
       ylim = ylim,
       xlab = "r", ylab = "Relative error",
       main = "Accuracy and failure regions")

  # Add other methods
  for (e in setdiff(err_cols, "err_bounded"))
    lines(res$r, res[[e]], col = cols[e], lwd = 2)

  abline(h = .Machine$double.eps, col = "darkgray", lty = 2)

  legend("bottomright",
         legend = c("Current bounded", "Direct exp",
                    "Preexp only", "Preexp + bound"),
         col = cols, lwd = 2, bty = "n")

  # Mark numeric failures
  for (e in err_cols) {
    bad <- which(!is.finite(res[[e]]) | res[[e]] <= 0)
    if (length(bad) > 0)
      points(res$r[bad], rep(ylim[1], length(bad)),
             pch = 21, col = cols[e], bg = cols[e], cex = 0.6)
  }

  legend("bottomleft", legend = "dots = 0/Inf/NaN failures", bty = "n")
}


################################################################################
# 3. Summarize accuracy across r
################################################################################
summarize_accuracy <- function(res) {
  err_cols <- c("err_direct", "err_bounded", "err_preexp", "err_preexp_bound")

  summary <- data.frame(
    Method = c("Direct exp", "Current bounded",
               "Preexp only", "Preexp + bound"),
    Mean_error = NA_real_,
    Median_error = NA_real_,
    Max_error = NA_real_,
    Finite_fraction = NA_real_,
    Zero_or_Inf_fraction = NA_real_
  )

  for (j in seq_along(err_cols)) {
    e <- res[[err_cols[j]]]
    finite_mask <- is.finite(e) & e > 0
    summary$Mean_error[j]        <- mean(e[finite_mask], na.rm = TRUE)
    summary$Median_error[j]      <- median(e[finite_mask], na.rm = TRUE)
    summary$Max_error[j]         <- max(e[finite_mask], na.rm = TRUE)
    summary$Finite_fraction[j]   <- mean(finite_mask)
    summary$Zero_or_Inf_fraction[j] <- 1 - mean(finite_mask)
  }

  summary
}


################################################################################
# 4. Alternate jitter plot for fine-scale comparison
################################################################################
plot_errors_jitter <- function(res, offset_for_visibility = TRUE) {
  err_cols <- c("err_bounded", "err_direct", "err_preexp", "err_preexp_bound")
  cols <- c("gray", "black", "red", "blue")

  message("Plotting columns:")
  for (i in seq_along(err_cols))
    message(sprintf("  %-15s  ->  %s", err_cols[i], cols[i]))

  offset_factor <- if (offset_for_visibility) c(1, 5, 100, 1e4) else rep(1, 4)

  finite_vals <- unlist(res[err_cols])
  finite_vals <- finite_vals[is.finite(finite_vals) & finite_vals > 0]
  if (length(finite_vals) > 0) {
    lower <- quantile(finite_vals, 0.01, na.rm = TRUE)
    upper <- quantile(finite_vals, 0.99, na.rm = TRUE)
    ylim <- c(lower / 10, upper * 10)
  } else ylim <- c(1e-20, 1e-12)

  plot(res$r, res$err_bounded * offset_factor[1],
       type = "l", log = "y", lwd = 2, col = cols[1],
       ylim = ylim,
       xlab = "r", ylab = "Relative error",
       main = "Accuracy (offset for visibility)")

  for (j in 2:length(err_cols))
    lines(res$r, res[[err_cols[j]]] * offset_factor[j], col = cols[j], lwd = 2)

  abline(h = .Machine$double.eps, col = "darkgray", lty = 2)
  legend("bottomright",
         legend = c("Current bounded", "Direct exp", "Preexp only", "Preexp + bound"),
         col = cols, lwd = 2)
}


################################################################################
# 5. Example usage
################################################################################
# Run test for a moderate K and r-range.
# Expand range (e.g. seq(-100, 80, 1)) to probe overflow/underflow limits.
res <- compare_all_methods(K = 10, r_vals = seq(-71, 71, length.out = 1e4))

# Plot and summarize
plot_errors(res)
summary_table <- summarize_accuracy(res)
print(summary_table, digits = 3)
# plot_errors_jitter(res)   # optional visualization with offsets
################################################################################


################################################################################
# 6. Ratio stability check (direct vs preexp) × (bound vs clipped)
################################################################################
compare_prob_ratios <- function(K = 5,
                                r_vals = seq(-20, 20, length.out = 200),
                                m_vals = NULL,
                                mpfr_prec = 256) {

  if (!requireNamespace("Rmpfr", quietly = TRUE))
    stop("Please install Rmpfr: install.packages('Rmpfr')")

  if (is.null(m_vals)) m_vals <- runif(K, -1, 1)

  res <- data.frame(
    r = numeric(length(r_vals)),
    err_direct_bound = numeric(length(r_vals)),
    err_direct_clip  = numeric(length(r_vals)),
    err_preexp_bound = numeric(length(r_vals)),
    err_preexp_clip  = numeric(length(r_vals))
  )

  for (i in seq_along(r_vals)) {
    r <- r_vals[i]
    b_raw  <- K * r
    b_clip <- max(0, b_raw)

    # --- High-precision reference ---------------------------------------------
    r_mp <- Rmpfr::mpfr(r, precBits = mpfr_prec)
    m_mp <- Rmpfr::mpfr(m_vals, precBits = mpfr_prec)
    exp_terms_ref <- exp(m_mp + (1:K) * r_mp)
    denom_ref <- 1 + sum(exp_terms_ref)
    p_ref_num <- as.numeric(exp_terms_ref / denom_ref)

    # --- (1) Direct, un-clipped bound ----------------------------------------
    exp_terms_dB <- exp(m_vals + (1:K) * r - b_raw)
    denom_dB <- exp(-b_raw) + sum(exp_terms_dB)
    p_dB <- exp_terms_dB / denom_dB
    res$err_direct_bound[i] <- max(abs(p_dB - p_ref_num) / p_ref_num)

    # --- (2) Direct, clipped bound -------------------------------------------
    exp_terms_dC <- exp(m_vals + (1:K) * r - b_clip)
    denom_dC <- exp(-b_clip) + sum(exp_terms_dC)
    p_dC <- exp_terms_dC / denom_dC
    res$err_direct_clip[i] <- max(abs(p_dC - p_ref_num) / p_ref_num)

    # --- (3) Preexp, un-clipped bound ---------------------------------------
    eR <- exp(r)
    eM <- exp(m_vals)
    eB <- exp(-b_raw)
    powE <- eR
    S_preB <- eB
    terms_preB <- numeric(K)
    for (c in 1:K) {
      term <- eM[c] * powE * eB
      terms_preB[c] <- term
      S_preB <- S_preB + term
      powE <- powE * eR
    }
    p_preB <- terms_preB / S_preB
    res$err_preexp_bound[i] <- max(abs(p_preB - p_ref_num) / p_ref_num)

    # --- (4) Preexp, clipped bound ------------------------------------------
    eR <- exp(r)
    eM <- exp(m_vals)
    eB <- exp(-b_clip)
    powE <- eR
    S_preC <- eB
    terms_preC <- numeric(K)
    for (c in 1:K) {
      term <- eM[c] * powE * eB
      terms_preC[c] <- term
      S_preC <- S_preC + term
      powE <- powE * eR
    }
    p_preC <- terms_preC / S_preC
    res$err_preexp_clip[i] <- max(abs(p_preC - p_ref_num) / p_ref_num)

    res$r[i] <- r
  }

  return(res)
}






################################################################################
# 7. Example usage: compare probability ratio stability
################################################################################

K <- 10
r_vals <- seq(-75, 75, length.out = 1e4)
set.seed(123)
m_vals <- runif(K, -1, 1)

res_ratio <- compare_prob_ratios(K = K, r_vals = r_vals, m_vals = m_vals)

eps <- .Machine$double.eps
plot(res_ratio$r, pmax(res_ratio$err_direct_bound, eps),
     type = "l", log = "y", lwd = 2, col = "red",
     xlab = "r", ylab = "Relative error (vs MPFR reference)",
     main = "Numerical stability of p_c ratio computations — 4 variants")

lines(res_ratio$r, pmax(res_ratio$err_direct_clip, eps),  col = "blue",   lwd = 2)
lines(res_ratio$r, pmax(res_ratio$err_preexp_bound, eps), col = "orange", lwd = 2)
lines(res_ratio$r, pmax(res_ratio$err_preexp_clip, eps),  col = "purple", lwd = 2)

abline(h = .Machine$double.eps, col = "darkgray", lty = 2)
legend("top",
       legend = c("Direct + Bound", "Direct + Clipped Bound",
                  "Preexp + Bound", "Preexp + Clipped Bound"),
       col = c("red", "blue", "orange", "purple"),
       lwd = 2, bty = "n")

abline(v = -70)
abline(v = 70)

# Summarize numeric accuracy
summary_df <- data.frame(
  Method = c("Direct + Bound", "Direct + Clipped Bound",
             "Preexp + Bound", "Preexp + Clipped Bound"),
  Mean_error = c(mean(res_ratio$err_direct_bound, na.rm = TRUE),
                 mean(res_ratio$err_direct_clip,  na.rm = TRUE),
                 mean(res_ratio$err_preexp_bound, na.rm = TRUE),
                 mean(res_ratio$err_preexp_clip,  na.rm = TRUE)),
  Median_error = c(median(res_ratio$err_direct_bound, na.rm = TRUE),
                   median(res_ratio$err_direct_clip,  na.rm = TRUE),
                   median(res_ratio$err_preexp_bound, na.rm = TRUE),
                   median(res_ratio$err_preexp_clip,  na.rm = TRUE)),
  Max_error = c(max(res_ratio$err_direct_bound, na.rm = TRUE),
                max(res_ratio$err_direct_clip,  na.rm = TRUE),
                max(res_ratio$err_preexp_bound, na.rm = TRUE),
                max(res_ratio$err_preexp_clip,  na.rm = TRUE))
)
print(summary_df, digits = 3)


