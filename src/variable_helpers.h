#include <RcppArmadillo.h>
#include "explog_switch.h"



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
inline arma::vec compute_denom_ordinal(const arma::vec& residual,
                                       const arma::vec& main_eff,
                                       const arma::vec& bound)
{
  constexpr double EXP_BOUND = 709.0;
  const int K = static_cast<int>(main_eff.n_elem);

  // --- Binary shortcut (K == 1) ---------------------------------------------
  if (K == 1) {
    return ARMA_MY_EXP(-bound) + ARMA_MY_EXP(main_eff[0] + residual - bound);
  }

  const arma::uword N = bound.n_elem;
  arma::vec denom(N, arma::fill::none);
  const arma::vec eM = ARMA_MY_EXP(main_eff);

  // Fast block: uses eB inside the loop (avoids intermediate overflow)
  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    arma::vec r = residual.rows(i0, i1);
    arma::vec b = bound.rows(i0, i1);
    arma::vec eR = ARMA_MY_EXP(r);
    arma::vec eB = ARMA_MY_EXP(-b);
    arma::vec pow = eR;

    arma::vec d = eB;
    for (int c = 0; c < K; ++c) {
      d += eM[c] * pow % eB;
      pow %= eR;
    }
    denom.rows(i0, i1) = d;
  };

  // Safe block: stabilized exponent; NO clamp here by design
  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    arma::vec r = residual.rows(i0, i1);
    arma::vec b = bound.rows(i0, i1);

    arma::vec d = ARMA_MY_EXP(-b);
    for (int c = 0; c < K; ++c) {
      arma::vec ex = main_eff[c] + (c + 1) * r - b;
      d += ARMA_MY_EXP(ex);
    }
    denom.rows(i0, i1) = d;
  };

  // Single linear scan over contiguous runs
  const double* bp = bound.memptr();
  arma::uword i = 0;
  while (i < N) {
    const bool fast = !(bp[i] < -EXP_BOUND || bp[i] > EXP_BOUND);
    arma::uword j = i + 1;
    while (j < N) {
      const bool fast_j = !(bp[j] < -EXP_BOUND || bp[j] > EXP_BOUND);
      if (fast_j != fast) break;
      ++j;
    }
    if (fast) do_fast_block(i, j - 1);
    else do_safe_block(i, j - 1);
    i = j;
  }

  return denom;
}

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
inline arma::vec compute_denom_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b          // update in place: per-person bound b[i]
) {

  constexpr double EXP_BOUND = 709.0;
  const arma::uword N = residual.n_elem;
  arma::vec denom(N);

  // ---- 1. Precompute theta_part[cat] and exp(theta_part) ----
  arma::vec cat = arma::regspace<arma::vec>(0, num_cats);
  arma::vec centered = cat - double(ref);
  arma::vec theta = lin_eff * centered + quad_eff * arma::square(centered);
  arma::vec exp_theta = ARMA_MY_EXP(theta);

  // ---- 2. Numerical bounds [b] ----
  b.set_size(N);
  b.fill(theta[0]);
  for (int c = 1; c <= num_cats; c++)
    b = arma::max(b, theta[c] + double(c) * residual);

  // ---- 3. Bounds for the FAST power chain: c*r - b ----
  // For fixed i, c*r[i] - b[i] ranges between -b[i] and num_cats*r[i] - b[i].
  // We need max_c (c*r[i] - b[i]) <= EXP_BOUND to avoid overflow in pow.
  arma::vec pow_bound = double(num_cats) * residual - b;

  // ---- 4. FAST BLOCK: Preexp + bounded power chain ----
  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);

    arma::vec eR = ARMA_MY_EXP(r);         // exp(r)
    arma::vec pow = ARMA_MY_EXP(-bb);       // start at cat=0 term: exp(0*r - b)
    arma::vec d = exp_theta[0] * pow;

    for (int c = 1; c <= num_cats; c++) {
      pow %= eR;                            // exp(c*r - b)
      d += exp_theta[c] * pow;
    }
    denom.rows(i0, i1) = d;
  };

  // ---- 5. SAFE BLOCK: direct exp(theta[c] + c*r - b) ----
  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);

    arma::vec d(bb.n_elem, arma::fill::zeros);
    for (int c = 0; c <= num_cats; c++) {
      arma::vec ex = theta[c] + double(c) * r - bb;
      d += ARMA_MY_EXP(ex);
    }



    denom.rows(i0, i1) = d;
  };

  // ---- 6. BLOCK SCAN: decide FAST vs SAFE per contiguous run ----
  const double* bp = b.memptr();
  const double* pp = pow_bound.memptr();

  arma::uword i = 0;
  while (i < N) {
    const bool fast_i = (std::abs(bp[i]) <= EXP_BOUND) && (std::abs(pp[i]) <= EXP_BOUND);

    arma::uword j = i + 1;
    while (j < N) {
      const bool fast_j = (std::abs(bp[j]) <= EXP_BOUND) && (std::abs(pp[j]) <= EXP_BOUND);
      if (fast_j != fast_i) break;
      ++j;
    }

    if (fast_i) do_fast_block(i, j - 1);
    else  do_safe_block(i, j - 1);

    i = j;
  }

  return denom;
}



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
inline arma::mat compute_probs_ordinal(const arma::vec& main_param,
                                       const arma::vec& residual_score,
                                       const arma::vec& bound,
                                       int num_cats)
{
  constexpr double EXP_BOUND = 709.0;
  const arma::uword N = bound.n_elem;

  if (num_cats == 1) {
    arma::vec b = arma::clamp(bound, 0.0, arma::datum::inf);
    arma::vec ex = main_param(0) + residual_score - b;
    arma::vec t = ARMA_MY_EXP(ex);
    arma::vec den = ARMA_MY_EXP(-b) + t;
    arma::mat probs(N, 2, arma::fill::none);
    probs.col(1) = t / den;
    probs.col(0) = 1.0 - probs.col(1);
    return probs;
  }

  arma::mat probs(N, num_cats + 1, arma::fill::none);
  const arma::vec eM = ARMA_MY_EXP(main_param);

  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    auto P = probs.rows(i0, i1).cols(1, num_cats);
    arma::vec r = residual_score.rows(i0, i1);
    arma::vec eR = ARMA_MY_EXP(r);
    arma::vec pow = eR;
    arma::vec den(P.n_rows, arma::fill::ones);
    for (int c = 0; c < num_cats; c++) {
      arma::vec term = eM[c] * pow;
      P.col(c) = term;
      den += term;
      pow %= eR;
    }
    P.each_col() /= den;
  };

  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    auto P = probs.rows(i0, i1).cols(1, num_cats);
    arma::vec r = residual_score.rows(i0, i1);
    arma::vec b = arma::clamp(bound.rows(i0, i1), 0.0, arma::datum::inf);
    arma::vec den = ARMA_MY_EXP(-b);
    for (int c = 0; c < num_cats; c++) {
      arma::vec ex = main_param(c) + (c + 1) * r - b;
      arma::vec t = ARMA_MY_EXP(ex);
      P.col(c) = t;
      den += t;
    }
    P.each_col() /= den;
  };

  // Single linear scan; no std::abs
  const double* bp = bound.memptr();
  arma::uword i = 0;
  while (i < N) {
    const bool fast = !(bp[i] < -EXP_BOUND || bp[i] > EXP_BOUND);
    arma::uword j = i + 1;
    while (j < N) {
      const bool fast_j = !(bp[j] < -EXP_BOUND || bp[j] > EXP_BOUND);
      if (fast_j != fast) break;
      j++;
    }
    if (fast) do_fast_block(i, j - 1);
    else do_safe_block(i, j - 1);
    i = j;
  }

  probs.col(0) = 1.0 - arma::sum(probs.cols(1, num_cats), 1);
  return probs;
}



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
//   used otherwise: direct exp(θ(c) + c * r_i - b_i)
//
// Under these conditions, denom is finite and > 0, so no one-hot fallback.
// -----------------------------------------------------------------------------
inline arma::mat compute_probs_blume_capel(const arma::vec& residual,
                                           const double lin_eff,
                                           const double quad_eff,
                                           const int ref,
                                           const int num_cats,
                                           arma::vec& b)   // updated in place
{
  constexpr double EXP_BOUND = 709.0;

  const arma::uword N = residual.n_elem;
  arma::mat probs(N, num_cats + 1, arma::fill::none);

  // 1. Precompute θ(c) and exp(θ(c))
  arma::vec cat = arma::regspace<arma::vec>(0, num_cats);
  arma::vec centered = cat - double(ref);
  arma::vec theta = lin_eff * centered + quad_eff * arma::square(centered);
  arma::vec exp_theta = ARMA_MY_EXP(theta);

  // 2. Compute bounds b[i] = max_c (θ(c) + c * r_i)
  b.set_size(N);
  b.fill(theta[0]);
  for (int c = 1; c <= num_cats; ++c) {
    b = arma::max(b, theta[c] + double(c) * residual);
  }

  // 3. Bound for the power chain: max_c (c * r_i - b_i) = num_cats * r_i - b_i
  arma::vec pow_bound = double(num_cats) * residual - b;

  // FAST block: preexp + bounded power chain
  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    auto P = probs.rows(i0, i1);
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);
    const arma::uword B = bb.n_elem;

    arma::vec eR = ARMA_MY_EXP(r);        // exp(r_i)
    arma::vec pow = ARMA_MY_EXP(-bb);      // exp(0 * r_i - b_i)
    arma::vec denom(B, arma::fill::zeros);

    // c = 0
    arma::vec col0 = exp_theta[0] * pow;
    P.col(0) = col0;
    denom += col0;

    // c = 1..num_cats
    for (int c = 1; c <= num_cats; ++c) {
      pow %= eR;                           // exp(c * r_i - b_i)
      arma::vec col = exp_theta[c] * pow;
      P.col(c) = col;
      denom += col;
    }

    P.each_col() /= denom;
  };

  // SAFE block: direct exp(θ(c) + c * r_i - b_i)
  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    auto P = probs.rows(i0, i1);
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);
    const arma::uword B = bb.n_elem;
    arma::vec denom(B, arma::fill::zeros);

    for (int c = 0; c <= num_cats; ++c) {
      arma::vec ex = theta[c] + double(c) * r - bb;
      arma::vec col = ARMA_MY_EXP(ex);
      P.col(c) = col;
      denom += col;
    }
    P.each_col() /= denom;
  };

  // 4. Single linear scan over contiguous FAST/SAFE runs (same as denom)
  const double* bp = b.memptr();
  const double* pp = pow_bound.memptr();
  arma::uword i = 0;
  while (i < N) {
    const bool fast_i =
      (std::abs(bp[i]) <= EXP_BOUND) && (std::abs(pp[i]) <= EXP_BOUND);

    arma::uword j = i + 1;
    while (j < N) {
      const bool fast_j =
        (std::abs(bp[j]) <= EXP_BOUND) && (std::abs(pp[j]) <= EXP_BOUND);
      if (fast_j != fast_i) break;
      j++;
    }

    if (fast_i) do_fast_block(i, j - 1);
    else do_safe_block(i, j - 1);

    i = j;
  }

  return probs;
}