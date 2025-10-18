# ==============================================================================
# Blume–Capel Numerical Stability Study (clean, documented)
# File: dev/numerical_analyses/BCvar_normalization_PL.r
#
# Goal
# ----
# Compare numerical stability of four ways to compute the Blume–Capel
# normalizing constant across a range of residual scores r:
#
#   Methods (exactly four):
#     1) Direct            — unbounded sum of exp(θ_lin*s + θ_quad*s^2 + s*r)
#     2) Preexp            — unbounded “power-chain” (precompute exp(θ-part), reuse exp(r))
#     3) Direct + bound    — SPLIT bound: b_theta + b_r, computing sum(exp(term - b_theta - b_r))
#     4) Preexp + bound    — SPLIT bound with power-chain on r-part only
#
#   Split bound (the only bound we use):
#     - b_theta = max_s (θ_lin*s + θ_quad*s^2)          # depends only on model, constant in r
#     - b_r     = C* * r, where C* = argmax_s |s|       # depends only on r and the score range
#
# References (for error calculation):
#   - ref_unscaled      = MPFR sum(exp(term))
#   - ref_split_scaled  = MPFR sum(exp( (θ_part - b_theta) + (s*r - b_r) ))
#
# Dependencies
# ------------
#   - Rmpfr
#
# Outputs
# -------
# - compare_bc_all_methods(...) returns a data.frame with:
#     r                    : grid of residual scores
#     direct               : numeric, Σ exp(term)
#     preexp               : numeric, Σ via power-chain (unbounded)
#     direct_bound         : numeric, Σ exp((θ-bθ) + (s*r - b_r))
#     preexp_bound         : numeric, Σ via power-chain with split bound
#     err_direct           : |(direct - ref_unscaled)/ref_unscaled|
#     err_preexp           : |(preexp - ref_unscaled)/ref_unscaled|
#     err_direct_bound     : |(direct_bound - ref_split_scaled)/ref_split_scaled|
#     err_preexp_bound     : |(preexp_bound - ref_split_scaled)/ref_split_scaled|
#     ref_unscaled         : numeric MPFR reference (unbounded)
#     ref_split_scaled     : numeric MPFR reference (split-scaled)
#
# Plotting helpers:
#   - plot_bc_four(res, ...)
#   - summarize_bc_four(res)
#
# ==============================================================================
library(Rmpfr)

# ------------------------------------------------------------------------------
# compare_bc_all_methods
# ------------------------------------------------------------------------------
# Compute all four methods and MPFR references over a vector of r-values.
#
# Args:
#   C          : integer, max category (scores are s = 0..C)
#   ref        : integer, baseline category index (scores centered by s <- 0:C - ref)
#   r_vals     : numeric vector of r values to scan
#   theta_lin  : numeric, linear θ parameter
#   theta_quad : numeric, quadratic θ parameter
#   mpfr_prec  : integer, MPFR precision (bits) for reference calculations
#
# Returns:
#   data.frame with columns described in the file header (see “Outputs”).
#
compare_bc_all_methods <- function(C = 10,
                                   ref = 3,
                                   r_vals = seq(-70, 70, length.out = 2000),
                                   theta_lin = 0.12,
                                   theta_quad = -0.02,
                                   mpfr_prec = 256) {

  # --- score grid and θ-part ---------------------------------------------------
  s_vals <- 0:C - ref
  smin   <- min(s_vals)

  # θ_part(s) = θ_lin*s + θ_quad*s^2
  theta_part <- theta_lin * s_vals + theta_quad * s_vals^2

  # --- split bound components (fixed model part, r-dependent rest part) -------
  b_theta <- max(theta_part)                       # constant over r
  Cstar   <- s_vals[ which.max(abs(s_vals)) ]      # farthest-from-zero score

  # Precomputed exponentials for the two chains we use
  exp_m                <- exp(theta_part)              # for unbounded preexp chain
  exp_m_theta_scaled   <- exp(theta_part - b_theta)    # for split-bounded variants

  # Output container
  res <- data.frame(
    r = r_vals,
    direct           = NA_real_,
    preexp           = NA_real_,
    direct_bound     = NA_real_,    # split bound
    preexp_bound     = NA_real_,    # split bound
    err_direct           = NA_real_,
    err_preexp           = NA_real_,
    err_direct_bound     = NA_real_,
    err_preexp_bound     = NA_real_,
    ref_unscaled     = NA_real_,
    ref_split_scaled = NA_real_
  )

  # --- MPFR constants independent of r ----------------------------------------
  tl    <- mpfr(theta_lin,  mpfr_prec)
  tq    <- mpfr(theta_quad, mpfr_prec)
  s_mp  <- mpfr(s_vals,     mpfr_prec)
  bth_mp<- mpfr(b_theta,    mpfr_prec)

  # --- Main loop over r --------------------------------------------------------
  for (i in seq_along(r_vals)) {
    r <- r_vals[i]
    term <- theta_part + s_vals * r              # double-precision exponents

    # r-dependent split bound
    b_r <- Cstar * r

    # ---------- MPFR references ----------
    r_mp <- mpfr(r, mpfr_prec)
    br_mp <- mpfr(b_r, mpfr_prec)

    term_mp <- tl*s_mp + tq*s_mp*s_mp + s_mp*r_mp
    ref_unscaled_mp     <- sum(exp(term_mp))                                 # Σ exp(term)
    ref_split_scaled_mp <- sum(exp((tl*s_mp + tq*s_mp*s_mp - bth_mp) +
                                     (s_mp*r_mp - br_mp)))                     # Σ exp((θ-bθ)+(s*r-b_r))

    # Store numeric references
    res$ref_unscaled[i]     <- asNumeric(ref_unscaled_mp)
    res$ref_split_scaled[i] <- asNumeric(ref_split_scaled_mp)

    # ---------- (1) Direct (unbounded) ----------
    v_direct <- sum(exp(term))
    res$direct[i] <- v_direct

    # ---------- (2) Preexp (unbounded) ----------
    # Power-chain on exp(r): start at s = smin
    eR  <- exp(r)
    pow <- exp(smin * r)
    S_pre <- 0.0
    for (j in seq_along(s_vals)) {
      S_pre <- S_pre + exp_m[j] * pow
      pow   <- pow * eR
    }
    res$preexp[i] <- S_pre

    # ---------- (3) Direct + bound (SPLIT) ----------
    # Σ exp( (θ - b_theta) + (s r - b_r) ), explicit loop for clarity
    sum_val <- 0.0
    for (j in seq_along(s_vals)) {
      sum_val <- sum_val + exp( (theta_part[j] - b_theta) + (s_vals[j]*r - b_r) )
    }
    res$direct_bound[i] <- sum_val

    # ---------- (4) Preexp + bound (SPLIT) ----------
    # Fold ONLY -b_r into the r-chain; θ-part is pre-scaled outside.
    pow_b <- exp(smin * r - b_r)
    S_pre_b <- 0.0
    for (j in seq_along(s_vals)) {
      S_pre_b <- S_pre_b + exp_m_theta_scaled[j] * pow_b
      pow_b   <- pow_b * eR
    }
    res$preexp_bound[i] <- S_pre_b

    # ---------- Errors (vs MPFR) ----------
    res$err_direct[i]       <- asNumeric(abs((mpfr(v_direct, mpfr_prec)      - ref_unscaled_mp)     / ref_unscaled_mp))
    res$err_preexp[i]       <- asNumeric(abs((mpfr(S_pre,   mpfr_prec)      - ref_unscaled_mp)     / ref_unscaled_mp))
    res$err_direct_bound[i] <- asNumeric(abs((mpfr(res$direct_bound[i], mpfr_prec) - ref_split_scaled_mp) / ref_split_scaled_mp))
    res$err_preexp_bound[i] <- asNumeric(abs((mpfr(res$preexp_bound[i], mpfr_prec) - ref_split_scaled_mp) / ref_split_scaled_mp))
  }

  res
}

# ------------------------------------------------------------------------------
# plot_bc_four
# ------------------------------------------------------------------------------
# Plot the four relative error curves on a log y-axis.
#
# Args:
#   res        : data.frame produced by compare_bc_all_methods()
#   draw_order : character vector with any ordering of:
#                c("err_direct","err_direct_bound","err_preexp_bound","err_preexp")
#   alpha      : named numeric vector (0..1) alphas for the same names
#   lwd        : line width
#
# Returns: (invisible) NULL. Draws a plot.
#
plot_bc_four <- function(res,
                         draw_order = c("err_direct","err_direct_bound","err_preexp_bound","err_preexp"),
                         alpha = c(err_direct=0.00, err_direct_bound=0.00, err_preexp_bound=0.40, err_preexp=0.40),
                         lwd = 2) {
  base_cols <- c(err_direct="#000000", err_preexp="#D62728",
                 err_direct_bound="#1F77B4", err_preexp_bound="#9467BD")
  to_rgba <- function(hex, a) rgb(t(col2rgb(hex))/255, alpha=a)
  cols <- mapply(to_rgba, base_cols[draw_order], alpha[draw_order], SIMPLIFY=TRUE, USE.NAMES=TRUE)

  vals <- unlist(res[draw_order]); vals <- vals[is.finite(vals)]
  ylim <- if (length(vals)) {
    q <- stats::quantile(vals, c(.01,.99), na.rm=TRUE); c(q[1]/10, q[2]*10)
  } else c(1e-20,1e-12)

  first <- draw_order[1]
  plot(res$r, res[[first]], type="l", log="y", col=cols[[1]], lwd=lwd, ylim=ylim,
       xlab="r", ylab="Relative error (vs MPFR)",
       main="Blume–Capel: Direct / Preexp / (Split) Bound")
  if (length(draw_order) > 1) {
    for (k in 2:length(draw_order)) lines(res$r, res[[draw_order[k]]], col=cols[[k]], lwd=lwd)
  }
  abline(h=.Machine$double.eps, col="gray70", lty=2)
  legend("top",
         legend=c("Direct","Direct + bound (split)","Preexp + bound (split)","Preexp")
         [match(draw_order, c("err_direct","err_direct_bound","err_preexp_bound","err_preexp"))],
         col=cols, lwd=lwd, bty="n")
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# summarize_bc_four
# ------------------------------------------------------------------------------
# Summarize accuracy per method.
#
# Args:
#   res : data.frame from compare_bc_all_methods()
#
# Returns:
#   data.frame with columns: Method, Mean, Median, Max, Finite
#
summarize_bc_four <- function(res) {
  cols <- c("err_direct","err_direct_bound","err_preexp_bound","err_preexp")
  labs <- c("Direct","Direct+Bound(split)","Preexp+Bound(split)","Preexp")
  mk <- function(v){
    f <- is.finite(v) & v > 0
    c(Mean=mean(v[f]), Median=median(v[f]), Max=max(v[f]), Finite=mean(f))
  }
  out <- t(sapply(cols, function(nm) mk(res[[nm]])))
  data.frame(Method=labs, out, row.names=NULL, check.names=FALSE)
}

# ==============================================================================
# Example usage (uncomment to run locally)
# ------------------------------------------------------------------------------
# res <- compare_bc_all_methods(
#   C = 10, ref = 0,
#   r_vals = seq(-80, -70, length.out = 2000),
#   theta_lin = 0.12,
#   theta_quad = 0.50,
#   mpfr_prec = 256
# )
# plot_bc_four(res,
#              draw_order = c("err_direct","err_direct_bound","err_preexp_bound","err_preexp"),
#              alpha = c(err_direct=0.00, err_direct_bound=0.50, err_preexp_bound=0.0, err_preexp=0.00),
#              lwd = 1)
# abline(v = 70); abline(v = -70)
# print(summarize_bc_four(res), digits = 3)
# ==============================================================================

# ==============================================================================
# Results Observations (curated runs)
# ------------------------------------------------------------------------------
# The blocks below capture observed behavior for representative settings.
# They are comments only—kept in-file so future readers see context quickly.
#
# ## ref = 5, C = 10
# res <- compare_bc_all_methods(C = 10, ref = 5,
#   r_vals = seq(-70, 70, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - Direct bounded a bit less variable
#  - Direct unbounded is significantly worse than the rest
# Conclusion: Don’t use direct unbounded.
#
# res <- compare_bc_all_methods(C = 10, ref = 5,
#   r_vals = seq(-100, -65, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - All methods continue to work decently even until r < -100.
#  - Direct unbounded is significantly worse than the rest
# Conclusion: Don’t use direct unbounded.
#
# res <- compare_bc_all_methods(C = 10, ref = 5,
#   r_vals = seq(65, 80, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - The bounded variants work until r < 71
#  - The unbounded variants continue to work decently even until r > 100.
#  - Direct unbounded is significantly worse than the rest
# Conclusion: Use Preexp bounded for r < 71; Use Preexp unbounded for r > 71.
#
# ## ref = 0, C = 10
# res <- compare_bc_all_methods(C = 10, ref = 0,
#   r_vals = seq(-70, 70, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - Direct bounded is less variable; Direct is worst overall
# Conclusion: Don’t use direct unbounded.
#
# res <- compare_bc_all_methods(C = 10, ref = 0,
#   r_vals = seq(-80, -70, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - The two unbounded variants fail when r < 0.
#  - Preexp bounded works until r < 71
#  - Direct bounded until r < -76.
# Conclusion: Use Direct bounded for r in [-75, -70]; Use Preexp bounded for r > -71.
#
# res <- compare_bc_all_methods(C = 10, ref = 0,
#   r_vals = seq(65, 80, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - The two unbounded variants fail when r > 65
#  - Direct bounded remains stable throughout.
#  - Preexp bounded blows up after r > 71
# Conclusion: Use Direct bounded for r > 71; Use Preexp bounded for r < 71.
#
# ## ref = 10, C = 10
# res <- compare_bc_all_methods(C = 10, ref = 10,
#   r_vals = seq(-70, 70, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - The unbounded variants become worse after r > 40
#  - Direct bounded is less variable
# Conclusion: Don’t use unbounded.
#
# res <- compare_bc_all_methods(C = 10, ref = 10,
#   r_vals = seq(-80, -65, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - The two unbounded variants fail when r < 66.
#  - The two bounded variants continue to work well even until r < -100.
# Conclusion: Don’t use unbounded.
#
# res <- compare_bc_all_methods(C = 10, ref = 10,
#   r_vals = seq(65, 80, length.out = 2000),
#   theta_lin = 0.12, theta_quad = 0.50, mpfr_prec = 256)
#  - The two unbounded variants fail when r > 75; direct unbounded already past r > 0.
#  - Direct bounded remains stable until r < 76
#  - Preexp bounded blows up after r > 71
# Conclusion: Use Direct bounded for r > 71; Use Preexp bounded for r < 71.
#
# ### Overall conclusions
# - Direct unbounded is never recommended.
# - For ref near center (e.g., ref = 5 when C = 10):
#     • Use Preexp bounded for moderate |r| (≈ |r| < 70)
#     • Use Preexp unbounded for large positive r (r > 71)
#     • Use Direct bounded for large negative r (r < -70)
# - For ref near edge (e.g., ref = 0 or ref = 10):
#     • Use Direct bounded for large |r| (|r| > 70)
#     • Use Preexp bounded for moderate |r|
# ==============================================================================

