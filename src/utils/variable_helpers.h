#ifndef BGMS_VARIABLE_HELPERS_H
#define BGMS_VARIABLE_HELPERS_H

#include <RcppArmadillo.h>
#include "math/explog_switch.h"


// -----------------------------------------------------------------------------
// Struct to hold both log-normalizer and probabilities from joint computation.
// Used by logp_and_gradient to avoid duplicate probability/denominator calculations.
// -----------------------------------------------------------------------------
struct LogZAndProbs {
  arma::vec log_Z;   // log-normalizer for each person
  arma::mat probs;   // (num_persons x num_cats+1) probability matrix
};


// -----------------------------------------------------------------------------
// Compute a numerically stable sum of the form:
//
//   denom = exp(-bound) + sum_{cat=0}^{K-1} exp(main_effect_param(cat)
//                 + (cat + 1) * residual_score - bound)
//
// but evaluated efficiently using precomputed exponentials:
//
//   exp_r = exp(residual_score)
//   exp_m = exp(main_effect_param)
//   denom = exp(-bound) * ( 1 + sum_c exp_m[c] * exp_r^(c+1) )
//
// If non-finite values arise (overflow, underflow, NaN), a safe fallback
// recomputes the naive version using direct exponentials.
// ----------------------------------------------------------------------------
arma::vec compute_denom_ordinal(
    const arma::vec& residual,
    const arma::vec& main_eff,
    const arma::vec& bound
);

// -----------------------------------------------------------------------------
// Compute denom = Σ_c exp( θ(c) + c*r - b ), with
//    θ(c) = lin_eff*(c-ref) + quad_eff*(c-ref)^2
//    b    = max_c( θ(c) + c*r )   (vectorized)
//
// Two modes:
//
// FAST (preexp + power-chain):
//    denom = Σ_c exp_theta[c] * exp(-b) * exp(r)^c
// Used only when all exponent terms are safe:
//    |b| ≤ EXP_BOUND,
//    underflow_bound ≥ -EXP_BOUND,
//    num_cats*r - b ≤ EXP_BOUND.
// This guarantees the recursive pow-chain stays finite.
//
// SAFE (direct evaluation):
//    denom = Σ_c exp(θ(c) + c*r - b)
// Used whenever any FAST-condition fails. Slower but always stable.
//
// FAST gives identical results when safe, otherwise SAFE is used.
// -----------------------------------------------------------------------------
arma::vec compute_denom_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
);

/**
 * Compute category probabilities in a numerically stable manner.
 *
 * Uses pre-exp or bounded formulations depending on the magnitude of `bound`.
 *  - If |bound| < 700: uses cheaper direct pre-exp computation
 *  - Else: clips bound at zero and applies stabilized scaling
 *
 * Empirical tests (see R/compare_prob_ratios.R) showed:
 *   - Clipping necessary for bound < -700
 *   - Bounds improve stability when large
 *
 * Returns:
 *   probs: num_persons × (num_cats + 1) matrix of probabilities (row-normalized)
 */
arma::mat compute_probs_ordinal(
    const arma::vec& main_param,
    const arma::vec& residual_score,
    const arma::vec& bound,
    int num_cats
);

// -----------------------------------------------------------------------------
// Blume–Capel probabilities, numerically stable via FAST/SAFE split.
//
// Model:
//   θ(c) = lin_eff * (c - ref) + quad_eff * (c - ref)^2,  c = 0..num_cats
//   exps_i(c) = θ(c) + c * r_i
//   b_i       = max_c exps_i(c)
//
// Probabilities:
//   p_i(c) ∝ exp( exps_i(c) - b_i )
//
// FAST (preexp + power-chain, same bounds as compute_denom_blume_capel):
//   used when |b_i| ≤ EXP_BOUND and pow_bound_i = num_cats * r_i - b_i ≤ EXP_BOUND
//
// SAFE (direct):
//   used otherwise: direct exp(θ(c) + (c-ref) * r_i - b_i)
//
// Under these conditions, denom is finite and > 0, so no one-hot fallback.
// -----------------------------------------------------------------------------
arma::mat compute_probs_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
);

// -----------------------------------------------------------------------------
// Joint computation of log-normalizer and probabilities for ordinal variables.
// Avoids redundant computation by computing both in a single pass.
// -----------------------------------------------------------------------------
LogZAndProbs compute_logZ_and_probs_ordinal(
    const arma::vec& main_param,
    const arma::vec& residual_score,
    const arma::vec& bound,
    int num_cats
);

// -----------------------------------------------------------------------------
// Joint computation of log-normalizer and probabilities for Blume-Capel variables.
// Avoids redundant computation by computing both in a single pass.
// -----------------------------------------------------------------------------
LogZAndProbs compute_logZ_and_probs_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
);

#endif // BGMS_VARIABLE_HELPERS_H