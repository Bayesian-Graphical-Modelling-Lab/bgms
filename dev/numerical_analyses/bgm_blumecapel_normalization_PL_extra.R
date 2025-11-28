############################################################
# Blume–Capel normalization analysis:
# Numerical comparison of FAST vs SAFE exponentiation methods
#
# Objective
# ---------
# This script provides a full numerical investigation of two methods
# to compute the *scaled* Blume–Capel partition sum:
#
#     Z_scaled(r) = sum_{s=0}^C exp( θ_part(s) + s*r - M(r) )
#
# where
#     θ_part(s) = θ_lin * (s - ref) + θ_quad * (s - ref)^2
# and
#     M(r) = max_s ( θ_part(s) + s*r ).
#
# We compare two computational approaches:
#
#   SAFE  = Direct computation     : sum_s exp(θ_part + s*r - M(r))
#   FAST  = Power-chain precompute : sum_s exp(θ_part(s)) * exp(s*r - M(r))
#
# MPFR (256-bit) is used as the ground-truth reference.
#
# The goals are:
#
#  1. Determine each method's numerical stability across a wide range
#     of (max_cat, ref, θ_lin, θ_quad, r).
#
#  2. Map all cases where FAST becomes inaccurate or produces NaN.
#
#  3. Identify the correct switching rule for the C++ implementation:
#        if (bound <= ~709) use FAST
#        else               use SAFE
#
#     where `bound = M(r)` is the maximum exponent before rescaling.
#
#  4. Produce plots and summary statistics to permanently document the
#     reasoning behind this rule.
#
# Key numerical fact
# ------------------
#   exp(x) in IEEE double precision overflows at x ≈ 709.782712893.
#   Therefore any exponent near ±709 is dangerous.
#
# Outcome summary
# ---------------
# - SAFE is stable across the entire tested range.
# - FAST is perfectly accurate **as long as bound ≤ ~709**
# - All FAST failures (NaN or large error) occur only when bound > ~709
# - No FAST failures were observed below this threshold.
#
# This provides strong empirical justification for the C++ switching rule.
############################################################

library(Rmpfr)
library(dplyr)
library(ggplot2)

############################################################
# 1. Simulation function
#
# Simulates FAST vs SAFE across:
#   - parameter grid (max_cat, ref, θ_lin, θ_quad)
#   - range of r values
#
# Returns one large data.frame containing:
#   - the computed bound M(r)
#   - FAST and SAFE values
#   - MPFR reference
#   - relative errors
#   - logical OK flags (err < tol)
############################################################

simulate_bc_fast_safe <- function(param_grid,
                                  r_vals = seq(-80, 80, length.out = 2000),
                                  mpfr_prec = 256,
                                  tol = 1e-12) {

  if (!all(c("max_cat", "ref", "theta_lin", "theta_quad") %in% names(param_grid))) {
    stop("param_grid must have columns: max_cat, ref, theta_lin, theta_quad")
  }

  out_list <- vector("list", nrow(param_grid))

  for (cfg_idx in seq_len(nrow(param_grid))) {
    cfg <- param_grid[cfg_idx, ]
    max_cat    <- as.integer(cfg$max_cat)
    ref        <- as.integer(cfg$ref)
    theta_lin  <- as.numeric(cfg$theta_lin)
    theta_quad <- as.numeric(cfg$theta_quad)

    # Score grid and θ(s)
    scores   <- 0:max_cat
    centered <- scores - ref
    theta_part <- theta_lin * centered + theta_quad * centered^2
    exp_m <- exp(theta_part)  # used by FAST

    # Build MPFR constants
    tl_mpfr        <- mpfr(theta_lin,  mpfr_prec)
    tq_mpfr        <- mpfr(theta_quad, mpfr_prec)
    sc_center_mpfr <- mpfr(centered,   mpfr_prec)
    sc_raw_mpfr    <- mpfr(scores,     mpfr_prec)

    n_r <- length(r_vals)

    res_cfg <- data.frame(
      config_id   = rep(cfg_idx, n_r),
      max_cat     = rep(max_cat, n_r),
      ref         = rep(ref, n_r),
      theta_lin   = rep(theta_lin, n_r),
      theta_quad  = rep(theta_quad, n_r),
      r           = r_vals,
      bound       = NA_real_,
      fast_val    = NA_real_,
      safe_val    = NA_real_,
      err_fast    = NA_real_,
      err_safe    = NA_real_,
      ok_fast     = NA,
      ok_safe     = NA,
      ref_scaled  = NA_real_,
      stringsAsFactors = FALSE
    )

    # Compute for all r
    for (i in seq_along(r_vals)) {
      r <- r_vals[i]

      term     <- theta_part + scores * r        # θ(s) + s*r
      term_max <- max(term)                      # numerical bound
      res_cfg$bound[i] <- term_max

      # MPFR scaled reference
      r_mpfr    <- mpfr(r, mpfr_prec)
      term_mpfr <- tl_mpfr * sc_center_mpfr +
        tq_mpfr * sc_center_mpfr * sc_center_mpfr +
        sc_raw_mpfr * r_mpfr
      term_max_mpfr   <- mpfr(max(asNumeric(term_mpfr)), mpfr_prec)
      ref_scaled_mpfr <- sum(exp(term_mpfr - term_max_mpfr))
      ref_scaled_num  <- asNumeric(ref_scaled_mpfr)
      res_cfg$ref_scaled[i] <- ref_scaled_num

      # SAFE method: direct evaluation
      safe_sum <- 0.0
      for (j in seq_along(scores)) {
        safe_sum <- safe_sum + exp(theta_part[j] + scores[j] * r - term_max)
      }
      res_cfg$safe_val[i] <- safe_sum

      # FAST method: preexp power-chain
      eR    <- exp(r)
      pow_b <- exp(-term_max)
      fast_sum <- 0.0
      for (j in seq_along(scores)) {
        fast_sum <- fast_sum + exp_m[j] * pow_b
        pow_b    <- pow_b * eR
      }
      res_cfg$fast_val[i] <- fast_sum

      # Relative errors
      if (is.finite(ref_scaled_num) && ref_scaled_num > 0) {
        res_cfg$err_safe[i] <- abs(safe_sum - ref_scaled_num) / ref_scaled_num
        res_cfg$err_fast[i] <- abs(fast_sum - ref_scaled_num) / ref_scaled_num
      }

      res_cfg$ok_safe[i] <- !is.na(res_cfg$err_safe[i]) &&
        is.finite(res_cfg$err_safe[i]) &&
        (res_cfg$err_safe[i] < tol)

      res_cfg$ok_fast[i] <- !is.na(res_cfg$err_fast[i]) &&
        is.finite(res_cfg$err_fast[i]) &&
        (res_cfg$err_fast[i] < tol)
    }

    out_list[[cfg_idx]] <- res_cfg
  }

  do.call(rbind, out_list)
}

############################################################
# 2. Parameter grid and simulation
############################################################

param_grid <- expand.grid(
  max_cat    = c(10),
  ref        = c(0, 5, 10),
  theta_lin  = c(-0.5, 0.0, 0.5),
  theta_quad = c(-0.2, 0.0, 0.2),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Very wide r-range so that bound covers deep negative and deep positive
r_vals <- seq(-100, 100, length.out = 5001)

tol <- 1e-12

sim_res <- simulate_bc_fast_safe(
  param_grid = param_grid,
  r_vals     = r_vals,
  mpfr_prec  = 256,
  tol        = tol
)

############################################################
# 3. Post-processing: classify regions, log-errors, abs(bound)
############################################################

df <- sim_res %>%
  mutate(
    err_fast_clipped = pmax(err_fast, 1e-300),
    err_safe_clipped = pmax(err_safe, 1e-300),

    log_err_fast = log10(err_fast_clipped),
    log_err_safe = log10(err_safe_clipped),

    abs_bound = abs(bound),

    region = case_when(
      ok_fast & ok_safe  ~ "both_ok",
      !ok_fast & ok_safe ~ "only_safe_ok",
      ok_fast & !ok_safe ~ "only_fast_ok",
      TRUE               ~ "neither_ok"
    )
  )

############################################################
# 4. NaN analysis for FAST
#
# We explicitly check:
#
#   Are there *any* NaN occurrences for FAST with |bound| < 709 ?
#
# This is essential: if NaN occurs for FAST even when |bound| is small,
# then the switching rule would fail.
############################################################

df_nan <- sim_res %>% filter(is.nan(err_fast))

nan_summary <- df_nan %>%
  summarise(
    n_nan     = n(),
    min_bound = min(bound, na.rm = TRUE),
    max_bound = max(bound, na.rm = TRUE)
  )

print(nan_summary)

df_nan_inside <- df_nan %>% filter(abs(bound) < 709)

cat("\nNumber of FAST NaN cases with |bound| < 709: ",
    nrow(df_nan_inside), "\n\n")

############################################################
# 5. FAST and SAFE plots vs bound
#
# We also explicitly count how many cases fail (ok_* == FALSE)
# while |bound| < 709. If the switching rule is correct, this
# number should be zero for FAST in the region where we intend
# to use it.
############################################################

# Count failures for FAST and SAFE when |bound| < 709
fast_fail_inside <- df %>%
  filter(abs(bound) < 709, !ok_fast) %>%
  nrow()

safe_fail_inside <- df %>%
  filter(abs(bound) < 709, !ok_safe) %>%
  nrow()

cat("\nFAST failures with |bound| < 709:", fast_fail_inside, "\n")
cat("SAFE failures with |bound| < 709:", safe_fail_inside, "\n\n")

# FAST
ggplot(df, aes(x = bound, y = log_err_fast, colour = region)) +
  geom_point(alpha = 0.3, size = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = log10(tol), linetype = 2) +
  geom_vline(xintercept = 709,  linetype = 2) +
  geom_vline(xintercept = -709, linetype = 2) +
  scale_color_manual(values = c(
    both_ok      = "darkgreen",
    only_safe_ok = "orange",
    only_fast_ok = "blue",
    neither_ok   = "red"
  )) +
  labs(
    x = "bound = max_s (theta_part(s) + s*r)",
    y = "log10(relative error) of FAST",
    colour = "region",
    subtitle = paste(
      "FAST failures with |bound| < 709:", fast_fail_inside
    )
  ) +
  ggtitle("FAST method vs bound") +
  theme_minimal()

# SAFE
ggplot(df, aes(x = bound, y = log_err_safe, colour = region)) +
  geom_point(alpha = 0.3, size = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = log10(tol), linetype = 2) +
  geom_vline(xintercept = 709,  linetype = 2) +
  geom_vline(xintercept = -709, linetype = 2) +
  scale_color_manual(values = c(
    both_ok      = "darkgreen",
    only_safe_ok = "orange",
    only_fast_ok = "blue",
    neither_ok   = "red"
  )) +
  labs(
    x = "bound = max_s (theta_part(s) + s*r)",
    y = "log10(relative error) of SAFE",
    colour = "region",
    subtitle = paste(
      "SAFE failures with |bound| < 709:", safe_fail_inside
    )
  ) +
  ggtitle("SAFE method vs bound") +
  theme_minimal()


############################################################
# 6. Fraction of configurations per |bound|-bin
############################################################

df_bins <- df %>%
  filter(is.finite(bound)) %>%
  mutate(
    abs_bound = abs(bound),
    bound_bin = cut(
      abs_bound,
      breaks = seq(0, max(abs_bound, na.rm = TRUE) + 10, by = 10),
      include_lowest = TRUE
    )
  ) %>%
  group_by(bound_bin) %>%
  summarise(
    mid_abs_bound = mean(abs_bound, na.rm = TRUE),
    frac_fast_ok  = mean(ok_fast, na.rm = TRUE),
    frac_safe_ok  = mean(ok_safe, na.rm = TRUE),
    n             = n(),
    .groups       = "drop"
  )

ggplot(df_bins, aes(x = mid_abs_bound)) +
  geom_line(aes(y = frac_fast_ok, colour = "FAST ok")) +
  geom_line(aes(y = frac_safe_ok, colour = "SAFE ok")) +
  geom_vline(xintercept = 709, linetype = 2) +
  scale_colour_manual(values = c("FAST ok" = "blue", "SAFE ok" = "darkgreen")) +
  labs(
    x = "|bound| bin center",
    y = "fraction of configurations with err < tol",
    colour = ""
  ) +
  ggtitle("FAST vs SAFE numerical stability by |bound|") +
  theme_minimal()

############################################################
# 7. Summary printed to console
############################################################

cat("\n================ SUMMARY =================\n")
print(nan_summary)

cat("\nFAST NaN cases with |bound| < 709: ",
    nrow(df_nan_inside), "\n\n")

cat("
Interpretation:
--------------
- The SAFE method (direct + bound) remains stable and accurate across the
  entire tested parameter and residual range.

- The FAST method (preexp + bound) is extremely accurate when the maximum
  exponent before rescaling, `bound = M(r)`, satisfies:

      |bound| ≤ ~709

- As soon as bound exceeds approximately +709, FAST becomes unstable:
    * large numerical error
    * or NaN (observed systematically)
    * No such failures appear below this threshold.

C++ Implementation Rule (recommended):
--------------------------------------
if (bound <= 709.0) {
    // FAST: preexp + bound (power-chain)
} else {
    // SAFE: direct + bound
}

This script constitutes the full reproducible analysis supporting the choice
of this switching threshold in the C++ Blume–Capel normalization code.
")

############################################################
# End of script
############################################################
