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

compare_bc_all_methods <- function(max_cat    = 10,
                                   ref        = 3,
                                   r_vals     = seq(-70, 70, length.out = 2000),
                                   theta_lin  = 0.12,
                                   theta_quad = -0.02,
                                   mpfr_prec  = 256) {

  # --- score grid and θ-part ---------------------------------------------------
  scores   <- 0:max_cat                 # s = 0..C
  centered <- scores - ref              # (s - ref)

  # θ_part(s) = θ_lin*(s - ref) + θ_quad*(s - ref)^2
  theta_part <- theta_lin * centered + theta_quad * centered^2

  # For the unbounded power-chain: exp(θ_part(s))
  exp_m <- exp(theta_part)

  # Output container ------------------------------------------------------------
  res <- data.frame(
    r                = r_vals,
    direct           = NA_real_,
    preexp           = NA_real_,
    direct_bound     = NA_real_,
    preexp_bound     = NA_real_,
    err_direct       = NA_real_,
    err_preexp       = NA_real_,
    err_direct_bound = NA_real_,
    err_preexp_bound = NA_real_,
    ref_unscaled     = NA_real_,
    ref_scaled       = NA_real_,
    bound            = NA_real_,   # term_max = M(r), puur ter inspectie
    theta_lin        = theta_lin,
    theta_quad       = theta_quad,
    max_cat          = max_cat,
    ref              = ref
  )

  # --- MPFR constants independent of r ----------------------------------------
  tl_mpfr        <- mpfr(theta_lin,  mpfr_prec)
  tq_mpfr        <- mpfr(theta_quad, mpfr_prec)
  sc_center_mpfr <- mpfr(centered,   mpfr_prec)   # (s - ref)
  sc_raw_mpfr    <- mpfr(scores,     mpfr_prec)   # s

  # --- Main loop over r --------------------------------------------------------
  for (i in seq_along(r_vals)) {
    r <- r_vals[i]

    # Standard double-precision exponents
    term <- theta_part + scores * r

    # ---------- MPFR references ----------
    r_mpfr    <- mpfr(r, mpfr_prec)
    term_mpfr <- tl_mpfr * sc_center_mpfr +
      tq_mpfr * sc_center_mpfr * sc_center_mpfr +
      sc_raw_mpfr * r_mpfr

    term_max_mpfr      <- mpfr(max(asNumeric(term_mpfr)), mpfr_prec)
    ref_unscaled_mpfr  <- sum(exp(term_mpfr))
    ref_scaled_mpfr    <- sum(exp(term_mpfr - term_max_mpfr))

    # Store numeric references
    res$ref_unscaled[i] <- asNumeric(ref_unscaled_mpfr)
    res$ref_scaled[i]   <- asNumeric(ref_scaled_mpfr)

    # ---------- (1) Direct (unbounded) ----------
    v_direct        <- sum(exp(term))
    res$direct[i]   <- v_direct

    # ---------- (2) Preexp (unbounded) ----------
    # Power-chain on exp(r): s = 0..max_cat, so start at s=0 with pow = 1
    eR    <- exp(r)
    pow   <- 1.0
    S_pre <- 0.0
    for (j in seq_along(scores)) {
      S_pre <- S_pre + exp_m[j] * pow
      pow   <- pow * eR
    }
    res$preexp[i] <- S_pre

    # ---------- (3) Direct + max-bound ----------
    term_max        <- max(term)     # M(r)
    res$bound[i]    <- term_max

    sum_direct_bound <- 0.0
    for (j in seq_along(scores)) {
      sum_direct_bound <- sum_direct_bound +
        exp(theta_part[j] + scores[j] * r - term_max)
    }
    res$direct_bound[i] <- sum_direct_bound

    # ---------- (4) Preexp + max-bound ----------
    pow_b   <- exp(-term_max)        # s = 0 → exp(0*r - term_max)
    S_pre_b <- 0.0
    for (j in seq_along(scores)) {
      S_pre_b <- S_pre_b + exp_m[j] * pow_b
      pow_b   <- pow_b * eR
    }
    res$preexp_bound[i] <- S_pre_b

    # ---------- Errors (vs MPFR) ----------
    res$err_direct[i] <-
      asNumeric(abs((mpfr(v_direct, mpfr_prec)      - ref_unscaled_mpfr) / ref_unscaled_mpfr))
    res$err_preexp[i] <-
      asNumeric(abs((mpfr(S_pre,   mpfr_prec)       - ref_unscaled_mpfr) / ref_unscaled_mpfr))
    res$err_direct_bound[i] <-
      asNumeric(abs((mpfr(sum_direct_bound, mpfr_prec) - ref_scaled_mpfr) / ref_scaled_mpfr))
    res$err_preexp_bound[i] <-
      asNumeric(abs((mpfr(S_pre_b, mpfr_prec)       - ref_scaled_mpfr) / ref_scaled_mpfr))
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

  ## --- Theoretical bound where max term hits exp(709)
  scores   <- 0:res$max_cat[1]
  centered <- scores - res$ref[1]

  # θ_part(s) = θ_lin*(s-ref) + θ_quad*(s-ref)^2
  theta_part <- res$theta_lin[1]  * centered +
    res$theta_quad[1] * centered * centered

  U <- 709
  pos <- scores > 0

  if (any(pos)) {
    r_up_vec <- (U - theta_part[pos]) / scores[pos]
    r_up <- min(r_up_vec)
  } else {
    r_up <- Inf
  }

  # Geen zinvolle beneden-grens voor overflow met s >= 0
  r_low <- -Inf

  if (is.finite(r_up)) {
    abline(v = r_up, col = "darkgreen", lty = 2, lwd = 2)
  }

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
# res = compare_bc_all_methods(
#   max_cat = 4,
#   ref = 0,
#   r_vals = seq(170, 175, length.out = 1000),
#   theta_lin = 0,
#   theta_quad = 1.00,
#   mpfr_prec = 256
# )
# plot_bc_four(res,
#              draw_order = c("err_direct","err_direct_bound","err_preexp_bound","err_preexp"),
#              alpha = c(err_direct = 0.00,
#                        err_direct_bound = 1.00,
#                        err_preexp_bound = 1.00,
#                        err_preexp = 0.00),
#              lwd = 1)
# print(summarize_bc_four(res), digits = 3)
# ==============================================================================

scan_bc_configs <- function(max_cat_vec    = c(4, 10),
                            ref_vec        = c(0, 2),
                            theta_lin_vec  = c(0.0, 0.12),
                            theta_quad_vec = c(-0.02, 0.0, 0.02),
                            r_vals         = seq(-80, 80, length.out = 2000),
                            mpfr_prec      = 256,
                            tol            = 1e-12) {

  cfg_grid <- expand.grid(
    max_cat    = max_cat_vec,
    ref        = ref_vec,
    theta_lin  = theta_lin_vec,
    theta_quad = theta_quad_vec,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  all_summaries <- vector("list", nrow(cfg_grid))

  for (i in seq_len(nrow(cfg_grid))) {
    cfg <- cfg_grid[i, ]
    cat("Config", i, "of", nrow(cfg_grid), ":",
        "max_cat =", cfg$max_cat,
        "ref =", cfg$ref,
        "theta_lin =", cfg$theta_lin,
        "theta_quad =", cfg$theta_quad, "\n")

    res_i <- compare_bc_all_methods(
      max_cat    = cfg$max_cat,
      ref        = cfg$ref,
      r_vals     = r_vals,
      theta_lin  = cfg$theta_lin,
      theta_quad = cfg$theta_quad,
      mpfr_prec  = mpfr_prec
    )

    summ_i <- summarize_bc_methods(res_i, tol = tol)
    all_summaries[[i]] <- summ_i
  }

  do.call(rbind, all_summaries)
}

classify_bc_bound_methods <- function(res, tol = 1e-12,
                                      eps_better = 1e-3) {
  # tol      : threshold for "good enough" relative error
  # eps_better : multiplicative margin to call one method "better" when both good

  r  <- res$r
  eD <- res$err_direct_bound
  eP <- res$err_preexp_bound

  finiteD <- is.finite(eD) & eD > 0
  finiteP <- is.finite(eP) & eP > 0

  goodD <- finiteD & (eD < tol)
  goodP <- finiteP & (eP < tol)

  state <- character(length(r))

  for (i in seq_along(r)) {
    if (!goodD[i] && !goodP[i]) {
      state[i] <- "neither_good"
    } else if (goodD[i] && !goodP[i]) {
      state[i] <- "only_direct_good"
    } else if (!goodD[i] && goodP[i]) {
      state[i] <- "only_preexp_good"
    } else {
      # both good: compare which is better
      # e.g. if preexp_bound error is at least eps_better times smaller than direct_bound
      if (eP[i] <= eD[i] * (1 - eps_better)) {
        state[i] <- "both_good_preexp_better"
      } else if (eD[i] <= eP[i] * (1 - eps_better)) {
        state[i] <- "both_good_direct_better"
      } else {
        # both good and within eps_better fraction: treat as "tie"
        state[i] <- "both_good_similar"
      }
    }
  }

  data.frame(
    r      = r,
    err_direct_bound  = eD,
    err_preexp_bound  = eP,
    state  = factor(state),
    bound  = res$bound,
    max_cat = res$max_cat[1],
    ref     = res$ref[1],
    theta_lin  = res$theta_lin[1],
    theta_quad = res$theta_quad[1],
    stringsAsFactors = FALSE
  )
}

summarize_bc_bound_classification <- function(class_df) {
  # class_df is the output of classify_bc_bound_methods()

  r     <- class_df$r
  state <- as.character(class_df$state)

  if (length(r) == 0) {
    return(class_df[FALSE, ])  # empty
  }

  # Identify run boundaries where state changes
  blocks <- list()
  start_idx <- 1
  current_state <- state[1]

  for (i in 2:length(r)) {
    if (state[i] != current_state) {
      # close previous block
      blocks[[length(blocks) + 1]] <- list(
        state = current_state,
        i_start = start_idx,
        i_end   = i - 1
      )
      # start new block
      start_idx <- i
      current_state <- state[i]
    }
  }
  # close last block
  blocks[[length(blocks) + 1]] <- list(
    state = current_state,
    i_start = start_idx,
    i_end   = length(r)
  )

  # Turn into a data.frame with r-intervals and some diagnostics
  out_list <- vector("list", length(blocks))
  for (k in seq_along(blocks)) {
    b <- blocks[[k]]
    idx <- b$i_start:b$i_end
    out_list[[k]] <- data.frame(
      state     = b$state,
      r_min     = min(r[idx]),
      r_max     = max(r[idx]),
      # a few handy diagnostics per block:
      max_err_direct_bound = max(class_df$err_direct_bound[idx], na.rm = TRUE),
      max_err_preexp_bound = max(class_df$err_preexp_bound[idx], na.rm = TRUE),
      min_bound = min(class_df$bound[idx], na.rm = TRUE),
      max_bound = max(class_df$bound[idx], na.rm = TRUE),
      n_points  = length(idx),
      max_cat   = class_df$max_cat[1],
      ref       = class_df$ref[1],
      theta_lin  = class_df$theta_lin[1],
      theta_quad = class_df$theta_quad[1],
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out_list)
}

# 1. Run the basic comparison
r_vals <- seq(0, 100, length.out = 2000)

res4 <- compare_bc_all_methods(
  max_cat    = 4,
  ref        = 0,
  r_vals     = r_vals,
  theta_lin  = 0.12,
  theta_quad = -0.02,
  mpfr_prec  = 256
)

# 2. Classify per-r which bound-method wins
class4 <- classify_bc_bound_methods(res4, tol = 1e-12, eps_better = 1e-3)

# 3. Compress into r-intervals
summary4 <- summarize_bc_bound_classification(class4)
print(summary4, digits = 3)




simulate_bc_fast_safe <- function(param_grid,
                                  r_vals = seq(-80, 80, length.out = 2000),
                                  mpfr_prec = 256,
                                  tol = 1e-12) {
  # param_grid: data.frame with columns
  #   max_cat, ref, theta_lin, theta_quad
  # r_vals    : vector of residual r values
  # tol       : tolerance for "ok" numerics (relative error)
  #
  # Returns one big data.frame with columns:
  #   config_id, max_cat, ref, theta_lin, theta_quad,
  #   r, bound, fast_val, safe_val,
  #   err_fast, err_safe, ok_fast, ok_safe,
  #   ref_scaled (MPFR reference)

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

    # --- score grid and θ-part for this config --------------------------------
    scores   <- 0:max_cat
    centered <- scores - ref

    theta_part <- theta_lin * centered + theta_quad * centered^2
    exp_m      <- exp(theta_part)  # for fast method

    # MPFR constants
    tl_mpfr        <- mpfr(theta_lin,  mpfr_prec)
    tq_mpfr        <- mpfr(theta_quad, mpfr_prec)
    sc_center_mpfr <- mpfr(centered,   mpfr_prec)
    sc_raw_mpfr    <- mpfr(scores,     mpfr_prec)

    # Storage for this config
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

    # --- main loop over r for this config -------------------------------------
    for (i in seq_along(r_vals)) {
      r <- r_vals[i]

      ## Double-precision exponents:
      term     <- theta_part + scores * r        # θ_part(s) + s*r
      term_max <- max(term)                      # M(r) = bound
      res_cfg$bound[i] <- term_max

      ## MPFR reference (scaled with max-term):
      r_mpfr    <- mpfr(r, mpfr_prec)
      term_mpfr <- tl_mpfr * sc_center_mpfr +
        tq_mpfr * sc_center_mpfr * sc_center_mpfr +
        sc_raw_mpfr * r_mpfr
      term_max_mpfr   <- mpfr(max(asNumeric(term_mpfr)), mpfr_prec)
      ref_scaled_mpfr <- sum(exp(term_mpfr - term_max_mpfr))
      ref_scaled_num  <- asNumeric(ref_scaled_mpfr)
      res_cfg$ref_scaled[i] <- ref_scaled_num

      # --- SAFE: Direct + max-bound ------------------------------------------
      # Z_safe = sum_s exp(θ_part(s) + s*r - term_max)
      safe_sum <- 0.0
      for (j in seq_along(scores)) {
        safe_sum <- safe_sum + exp(theta_part[j] + scores[j] * r - term_max)
      }
      res_cfg$safe_val[i] <- safe_sum

      # --- FAST: Preexp + max-bound (power-chain) ----------------------------
      # Z_fast = sum_s exp(θ_part(s)) * exp(s*r - term_max)
      eR    <- exp(r)
      pow_b <- exp(-term_max)    # s = 0 → exp(0*r - term_max)
      fast_sum <- 0.0
      for (j in seq_along(scores)) {
        fast_sum <- fast_sum + exp_m[j] * pow_b
        pow_b    <- pow_b * eR
      }
      res_cfg$fast_val[i] <- fast_sum

      # --- Relative errors vs MPFR (scaled) ----------------------------------
      if (is.finite(ref_scaled_num) && ref_scaled_num > 0) {
        res_cfg$err_safe[i] <- abs(safe_sum - ref_scaled_num) / ref_scaled_num
        res_cfg$err_fast[i] <- abs(fast_sum - ref_scaled_num) / ref_scaled_num
      } else {
        res_cfg$err_safe[i] <- NA_real_
        res_cfg$err_fast[i] <- NA_real_
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