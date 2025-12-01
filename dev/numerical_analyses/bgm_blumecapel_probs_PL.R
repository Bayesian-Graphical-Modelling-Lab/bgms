############################################################
# Blume–Capel probabilities:
# Numerical comparison of 4 methods vs MPFR reference
#
# Methods:
#   - direct_unscaled   : naive softmax
#   - direct_bound      : softmax with subtraction of M(r)
#   - preexp_unscaled   : preexp(theta_part) + power chain (no bound)
#   - preexp_bound      : preexp(theta_part) + power chain (with bound)
#
# Reference:
#   - MPFR softmax with scaling by M(r)
############################################################

library(Rmpfr)
library(dplyr)
library(ggplot2)

EXP_BOUND <- 709

############################################################
# 1. Compare 4 methods for one BC configuration
############################################################

compare_bc_prob_4methods_one <- function(max_cat,
                                         ref,
                                         theta_lin,
                                         theta_quad,
                                         r_vals,
                                         mpfr_prec = 256) {

  s_vals <- 0:max_cat
  c_vals <- s_vals - ref
  n_s    <- length(s_vals)
  n_r    <- length(r_vals)

  # theta_part(s)
  theta_part_num <- theta_lin * c_vals + theta_quad * c_vals^2

  # MPFR parameters
  tl_mp <- mpfr(theta_lin,  mpfr_prec)
  tq_mp <- mpfr(theta_quad, mpfr_prec)
  s_mp  <- mpfr(s_vals,     mpfr_prec)
  c_mp  <- mpfr(c_vals,     mpfr_prec)

  # Precompute for preexp methods
  exp_theta <- exp(theta_part_num)

  res <- data.frame(
    r               = r_vals,
    bound           = NA_real_,
    pow_bound       = NA_real_,
    err_direct      = NA_real_,
    err_bound       = NA_real_,
    err_preexp      = NA_real_,
    err_preexp_bound= NA_real_
  )

  for (i in seq_len(n_r)) {
    r    <- r_vals[i]
    r_mp <- mpfr(r, mpfr_prec)

    ## MPFR reference probabilities (softmax with scaling)
    term_mp <- tl_mp * c_mp +
      tq_mp * c_mp * c_mp +
      s_mp  * r_mp

    M_num <- max(asNumeric(term_mp))
    M_mp  <- mpfr(M_num, mpfr_prec)

    num_ref_mp <- exp(term_mp - M_mp)
    Z_ref_mp   <- sum(num_ref_mp)
    p_ref_mp   <- num_ref_mp / Z_ref_mp
    p_ref      <- asNumeric(p_ref_mp)

    ## Double: exponents
    term_num <- theta_part_num + s_vals * r
    M        <- max(term_num)
    res$bound[i] <- M
    res$pow_bound[i] <- max_cat * r - M

    ## (1) direct_unscaled
    num_dir  <- exp(term_num)
    den_dir  <- sum(num_dir)
    p_dir    <- num_dir / den_dir

    ## (2) direct_bound
    num_b    <- exp(term_num - M)
    den_b    <- sum(num_b)
    p_b      <- num_b / den_b

    ## (3) preexp_unscaled
    eR       <- exp(r)
    pow      <- eR
    num_pre  <- numeric(n_s)
    den_pre  <- 0.0

    # s = 0 term
    num_pre[1] <- exp_theta[1] * 1.0
    den_pre    <- den_pre + num_pre[1]

    if (max_cat >= 1) {
      for (s in 1:max_cat) {
        num_pre[s + 1] <- exp_theta[s + 1] * pow
        den_pre        <- den_pre + num_pre[s + 1]
        pow            <- pow * eR
      }
    }
    p_pre <- num_pre / den_pre

    ## (4) preexp_bound
    eR2       <- exp(r)
    pow_b     <- exp(-M)
    num_preB  <- numeric(n_s)
    den_preB  <- 0.0

    for (s in 0:max_cat) {
      idx <- s + 1
      num_preB[idx] <- exp_theta[idx] * pow_b
      den_preB      <- den_preB + num_preB[idx]
      pow_b         <- pow_b * eR2
    }
    p_preB <- num_preB / den_preB

    ## Relative errors vs MPFR reference on non-negligible support
    tau <- 1e-15  # <-- tweak this

    support_mask <- p_ref >= tau
    if (!any(support_mask)) {
      support_mask <- p_ref == max(p_ref)  # degenerate case: all tiny, pick the max
    }

    rel_direct <- abs(p_dir  - p_ref)[support_mask] / p_ref[support_mask]
    rel_bound  <- abs(p_b    - p_ref)[support_mask] / p_ref[support_mask]
    rel_preexp <- abs(p_pre  - p_ref)[support_mask] / p_ref[support_mask]
    rel_preB   <- abs(p_preB - p_ref)[support_mask] / p_ref[support_mask]

    res$err_direct[i]       <- max(rel_direct)
    res$err_bound[i]        <- max(rel_bound)
    res$err_preexp[i]       <- max(rel_preexp)
    res$err_preexp_bound[i] <- max(rel_preB)



  }

  res
}

############################################################
# 2. Sweep across param_grid × r_vals
############################################################

simulate_bc_prob_4methods <- function(param_grid,
                                      r_vals,
                                      mpfr_prec = 256,
                                      tol = 1e-12) {

  if (!all(c("max_cat", "ref", "theta_lin", "theta_quad") %in% names(param_grid))) {
    stop("param_grid must have columns: max_cat, ref, theta_lin, theta_quad")
  }

  out_list <- vector("list", nrow(param_grid))

  for (cfg_idx in seq_len(nrow(param_grid))) {
    cfg <- param_grid[cfg_idx, ]

    res_cfg <- compare_bc_prob_4methods_one(
      max_cat    = cfg$max_cat,
      ref        = cfg$ref,
      theta_lin  = cfg$theta_lin,
      theta_quad = cfg$theta_quad,
      r_vals     = r_vals,
      mpfr_prec  = mpfr_prec
    )

    res_cfg$config_id  <- cfg_idx
    res_cfg$max_cat    <- cfg$max_cat
    res_cfg$ref        <- cfg$ref
    res_cfg$theta_lin  <- cfg$theta_lin
    res_cfg$theta_quad <- cfg$theta_quad

    # simple ok flags
    res_cfg$ok_direct       <- is.finite(res_cfg$err_direct)       & (res_cfg$err_direct       < tol)
    res_cfg$ok_bound        <- is.finite(res_cfg$err_bound)        & (res_cfg$err_bound        < tol)
    res_cfg$ok_preexp       <- is.finite(res_cfg$err_preexp)       & (res_cfg$err_preexp       < tol)
    res_cfg$ok_preexp_bound <- is.finite(res_cfg$err_preexp_bound) & (res_cfg$err_preexp_bound < tol)

    out_list[[cfg_idx]] <- res_cfg
  }

  do.call(rbind, out_list)
}

############################################################
# 3. Example broad analysis (you can adjust this)
############################################################

param_grid <- expand.grid(
  max_cat    = c(4, 10),
  ref        = c(0, 2, 4, 5, 10),
  theta_lin  = c(-0.5, 0.0, 0.5),
  theta_quad = c(-0.2, 0.0, 0.2),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

r_vals <- seq(-80, 80, length.out = 2001)
tol    <- 1e-12

sim4 <- simulate_bc_prob_4methods(
  param_grid = param_grid,
  r_vals     = r_vals,
  mpfr_prec  = 256,
  tol        = tol
)

############################################################
# 4. Summaries: where each method fails, as a function of bound/pow_bound
############################################################

df4 <- sim4 %>%
  mutate(
    abs_bound = abs(bound),
    err_direct_cl       = pmax(err_direct,       1e-300),
    err_bound_cl        = pmax(err_bound,        1e-300),
    err_preexp_cl       = pmax(err_preexp,       1e-300),
    err_preexp_bound_cl = pmax(err_preexp_bound, 1e-300),
    log_err_direct      = log10(err_direct_cl),
    log_err_bound       = log10(err_bound_cl),
    log_err_preexp      = log10(err_preexp_cl),
    log_err_preexp_bound= log10(err_preexp_bound_cl)
  )

# Example: failures for each method inside |bound| <= 709 & pow_bound <= 709
inside <- df4 %>%
  filter(abs(bound) <= EXP_BOUND, pow_bound <= EXP_BOUND)

n_direct_fail       <- sum(!inside$ok_direct)
n_bound_fail        <- sum(!inside$ok_bound)
n_preexp_fail       <- sum(!inside$ok_preexp)
n_preexp_bound_fail <- sum(!inside$ok_preexp_bound)

cat("\nFailures inside fast region (|bound| <= 709 & pow_bound <= 709):\n")
cat("  direct_unscaled     :", n_direct_fail, "\n")
cat("  direct_bound        :", n_bound_fail, "\n")
cat("  preexp_unscaled     :", n_preexp_fail, "\n")
cat("  preexp_bound (FAST) :", n_preexp_bound_fail, "\n\n")
