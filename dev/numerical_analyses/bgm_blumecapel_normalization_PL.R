# ==============================================================================
# Blume–Capel Numerical Stability Study (reparametrized)
# File: dev/numerical_analyses/BCvar_normalization_PL.r
#
# Goal
# ----
# Compare numerical stability of four ways to compute the Blume–Capel
# normalizing constant across a range of residual scores r, using the
# reparametrized form
#
#   Z(r) = sum_{s=0}^C exp( θ_part(s) + s * r ),
#
# where
#
#   θ_part(s) = θ_lin * (s - ref) + θ_quad * (s - ref)^2.
#
# This corresponds to the reformulated denominator where:
#   - scores s are in {0, 1, ..., C},
#   - the quadratic/linear θ-part is in terms of the centered score (s - ref),
#   - the “residual” r enters only through s * r.
#
# Methods (exactly four):
#   1) Direct
#        Unbounded sum of exp(θ_part(s) + s * r).
#
#   2) Preexp
#        Unbounded “power-chain” over s, precomputing exp(θ_part(s)) and
#        reusing exp(r):
#          Z(r) = sum_s exp(θ_part(s)) * (exp(r))^s .
#
#   3) Direct + max-bound
#        Per-r max-term bound M(r) = max_s (θ_part(s) + s * r),
#        computing
#          Z(r) = exp(M(r)) * sum_s exp(θ_part(s) + s * r - M(r)),
#        but returning only the *scaled* sum:
#          sum_s exp(θ_part(s) + s * r - M(r)).
#
#   4) Preexp + max-bound
#        Same max-term bound M(r) as in (3), but using the power-chain:
#          sum_s exp(θ_part(s)) * exp(s * r - M(r)).
#
# References (for error calculation):
#   - ref_unscaled = MPFR sum_s exp(θ_part(s) + s * r)
#   - ref_scaled   = MPFR sum_s exp(θ_part(s) + s * r - M(r)),
#                    where M(r) = max_s (θ_part(s) + s * r) in MPFR.
#
# Dependencies
# ------------
#   - Rmpfr
#
# Outputs
# -------
# compare_bc_all_methods(...) returns a data.frame with:
#   r                : grid of residual scores
#   direct           : numeric, Σ_s exp(θ_part(s) + s * r)
#   preexp           : numeric, Σ_s via power-chain (unbounded)
#   direct_bound     : numeric, Σ_s exp(θ_part(s) + s * r - M(r))
#   preexp_bound     : numeric, Σ_s via power-chain with max-term bound
#   err_direct       : |(direct       - ref_unscaled)/ref_unscaled|
#   err_preexp       : |(preexp       - ref_unscaled)/ref_unscaled|
#   err_direct_bound : |(direct_bound - ref_scaled  )/ref_scaled  |
#   err_preexp_bound : |(preexp_bound - ref_scaled  )/ref_scaled  |
#   ref_unscaled     : numeric MPFR reference (unbounded)
#   ref_scaled       : numeric MPFR reference (max-term scaled)
#
# Plotting helpers (unchanged interface):
#   - plot_bc_four(res, ...)
#   - summarize_bc_four(res)
#
# ==============================================================================

library(Rmpfr)

# ------------------------------------------------------------------------------
# compare_bc_all_methods
# ------------------------------------------------------------------------------
# Compute all four methods and MPFR references over a vector of r-values
# for the reparametrized Blume–Capel normalizing constant
#
#   Z(r) = sum_{s=0}^C exp( θ_lin * (s - ref) + θ_quad * (s - ref)^2 + s * r ).
#
# Args:
#   max_cat    : integer, max category C (scores are s = 0..C)
#   ref        : integer, baseline category index for centering (s - ref)
#   r_vals     : numeric vector of r values to scan
#   theta_lin  : numeric, linear θ parameter
#   theta_quad : numeric, quadratic θ parameter
#   mpfr_prec  : integer, MPFR precision (bits) for reference calculations
#
# Returns:
#   data.frame with columns described in the file header (see “Outputs”).
# ------------------------------------------------------------------------------

compare_bc_all_methods = function(max_cat = 10,
                                  ref = 3,
                                  r_vals = seq(-70, 70, length.out = 2000),
                                  theta_lin = 0.12,
                                  theta_quad = -0.02,
                                  mpfr_prec = 256) {

  # --- score grid and θ-part ---------------------------------------------------
  scores   = 0:max_cat                 # s = 0..C
  centered = scores - ref              # (s - ref)

  # θ_part(s) = θ_lin*(s - ref) + θ_quad*(s - ref)^2
  theta_part = theta_lin * centered + theta_quad * centered^2

  # For the unbounded power-chain: exp(θ_part(s))
  exp_m = exp(theta_part)

  # Output container ------------------------------------------------------------
  res = data.frame(
    r = r_vals,
    direct = NA_real_,
    preexp = NA_real_,
    direct_bound = NA_real_,
    preexp_bound = NA_real_,
    err_direct = NA_real_,
    err_preexp = NA_real_,
    err_direct_bound = NA_real_,
    err_preexp_bound = NA_real_,
    ref_unscaled = NA_real_,
    ref_scaled = NA_real_,
    theta_lin = theta_lin,
    theta_quad = theta_quad,
    max_cat = max_cat,
    ref = ref
  )

  # --- MPFR constants independent of r ----------------------------------------
  tl_mpfr        = mpfr(theta_lin,  mpfr_prec)
  tq_mpfr        = mpfr(theta_quad, mpfr_prec)
  sc_center_mpfr = mpfr(centered,   mpfr_prec)   # (s - ref)
  sc_raw_mpfr    = mpfr(scores,     mpfr_prec)   # s

  # --- Main loop over r --------------------------------------------------------
  for (i in seq_along(r_vals)) {
    r = r_vals[i]

    # Standard double-precision exponents
    term = theta_part + scores * r

    # ---------- MPFR references ----------
    r_mpfr      = mpfr(r, mpfr_prec)
    term_mpfr   = tl_mpfr * sc_center_mpfr +
      tq_mpfr * sc_center_mpfr * sc_center_mpfr +
      sc_raw_mpfr * r_mpfr
    term_max_mpfr    = mpfr(max(asNumeric(term_mpfr)), mpfr_prec)
    ref_unscaled_mpfr = sum(exp(term_mpfr))
    ref_scaled_mpfr   = sum(exp(term_mpfr - term_max_mpfr))

    # Store numeric references
    res$ref_unscaled[i] = asNumeric(ref_unscaled_mpfr)
    res$ref_scaled[i]   = asNumeric(ref_scaled_mpfr)

    # ---------- (1) Direct (unbounded) ----------
    v_direct = sum(exp(term))
    res$direct[i] = v_direct

    # ---------- (2) Preexp (unbounded) ----------
    # Power-chain on exp(r): s = 0..max_cat, so start at s=0 with pow = 1
    eR  = exp(r)
    pow = 1.0
    S_pre = 0.0
    for (j in seq_along(scores)) {
      S_pre = S_pre + exp_m[j] * pow
      pow   = pow * eR
    }
    res$preexp[i] = S_pre

    # ---------- (3) Direct + max-bound ----------
    term_max = max(term)
    sum_val  = 0.0
    for (j in seq_along(scores)) {
      sum_val = sum_val + exp(theta_part[j] + scores[j] * r - term_max)
    }
    res$direct_bound[i] = sum_val

    # ---------- (4) Preexp + max-bound ----------
    pow_b   = exp(0 * r - term_max)   # = exp(-term_max), starting at s = 0
    S_pre_b = 0.0
    for (j in seq_along(scores)) {
      S_pre_b = S_pre_b + exp_m[j] * pow_b
      pow_b   = pow_b * eR
    }
    res$preexp_bound[i] = S_pre_b

    # ---------- Errors (vs MPFR) ----------
    res$err_direct[i] =
      asNumeric(abs((mpfr(v_direct, mpfr_prec) - ref_unscaled_mpfr) / ref_unscaled_mpfr))
    res$err_preexp[i] =
      asNumeric(abs((mpfr(S_pre,   mpfr_prec) - ref_unscaled_mpfr) / ref_unscaled_mpfr))
    res$err_direct_bound[i] =
      asNumeric(abs((mpfr(res$direct_bound[i], mpfr_prec) - ref_scaled_mpfr) / ref_scaled_mpfr))
    res$err_preexp_bound[i] =
      asNumeric(abs((mpfr(res$preexp_bound[i], mpfr_prec) - ref_scaled_mpfr) / ref_scaled_mpfr))
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
plot_bc_four = function(res,
                        draw_order = c("err_direct","err_direct_bound",
                                       "err_preexp_bound","err_preexp"),
                        alpha = c(err_direct       = 0.00,
                                  err_direct_bound = 0.00,
                                  err_preexp_bound = 0.40,
                                  err_preexp       = 0.40),
                        lwd = 2) {

  base_cols = c(err_direct       = "#000000",
                err_preexp       = "#D62728",
                err_direct_bound = "#1F77B4",
                err_preexp_bound = "#9467BD")

  to_rgba = function(hex, a) rgb(t(col2rgb(hex))/255, alpha = a)

  cols = mapply(to_rgba, base_cols[draw_order], alpha[draw_order],
                SIMPLIFY = TRUE, USE.NAMES = TRUE)

  vals = unlist(res[draw_order])
  vals = vals[is.finite(vals)]
  ylim = if (length(vals)) {
    q = stats::quantile(vals, c(.01, .99), na.rm = TRUE)
    c(q[1] / 10, q[2] * 10)
  } else c(1e-20, 1e-12)

  first = draw_order[1]
  plot(res$r, res[[first]], type = "l", log = "y",
       col = cols[[1]], lwd = lwd, ylim = ylim,
       xlab = "r", ylab = "Relative error (vs MPFR)",
       main = "Blume–Capel: Direct / Preexp / (Split) Bound")

  if (length(draw_order) > 1) {
    for (k in 2:length(draw_order)) {
      lines(res$r, res[[draw_order[k]]], col = cols[[k]], lwd = lwd)
    }
  }

  abline(h = .Machine$double.eps, col = "gray70", lty = 2)

  ## --- Theoretical overflow bound for the *sum* under the new parametrization
  ## Scores: s = 0..max_cat, centered: (s - ref)
  scores   = 0:res$max_cat[1]
  centered = scores - res$ref[1]

  # a_s = θ_lin*(s-ref) + θ_quad*(s-ref)^2
  a_s = res$theta_lin[1]  * centered +
    res$theta_quad[1] * centered * centered

  U  = 709
  N  = length(scores)
  Ue = U - log(N)   # effective bound for the *sum* (log N correction)

  pos = scores > 0
  if (any(pos)) {
    r_up  = min((Ue - a_s[pos]) / scores[pos])
  } else {
    r_up = Inf
  }
  r_low = -Inf  # no finite lower overflow bound with s >= 0

  if (is.finite(r_low)) abline(v = r_low, col = "darkblue",  lty = 2, lwd = 2)
  if (is.finite(r_up))  abline(v = r_up,  col = "darkgreen", lty = 2, lwd = 2)

  print(r_low)
  print(r_up)

  legend("top",
         legend = c("Direct",
                    "Direct + bound (split)",
                    "Preexp + bound (split)",
                    "Preexp")
         [match(draw_order,
                c("err_direct","err_direct_bound",
                  "err_preexp_bound","err_preexp"))],
         col = cols, lwd = lwd, bty = "n")

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
summarize_bc_four = function(res) {
  cols = c("err_direct","err_direct_bound","err_preexp_bound","err_preexp")
  labs = c("Direct","Direct+Bound(split)","Preexp+Bound(split)","Preexp")
  mk = function(v){
    f = is.finite(v) & v > 0
    c(Mean=mean(v[f]), Median=median(v[f]), Max=max(v[f]), Finite=mean(f))
  }
  out = t(sapply(cols, function(nm) mk(res[[nm]])))
  data.frame(Method=labs, out, row.names=NULL, check.names=FALSE)
}

# ==============================================================================
# Example usage (uncomment to run locally)
# ------------------------------------------------------------------------------
res = compare_bc_all_methods(
  max_cat = 10,
  ref = 0,
  r_vals = seq(60, 90, length.out = 1000),
  theta_lin = 0,
  theta_quad = 1.00,
  mpfr_prec = 256
)
plot_bc_four(res,
             draw_order = c("err_direct","err_direct_bound","err_preexp_bound","err_preexp"),
             alpha = c(err_direct = 0.00,
                       err_direct_bound = 1.00,
                       err_preexp_bound = 1.00,
                       err_preexp = 0.00),
             lwd = 1)
print(summarize_bc_four(res), digits = 3)
# ==============================================================================


