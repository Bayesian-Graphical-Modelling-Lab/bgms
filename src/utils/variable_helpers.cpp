#include "utils/variable_helpers.h"


// =============================================================================
// compute_denom_ordinal
// =============================================================================
arma::vec compute_denom_ordinal(
    const arma::vec& residual,
    const arma::vec& main_eff,
    const arma::vec& bound
) {
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


// =============================================================================
// compute_denom_blume_capel
// =============================================================================
arma::vec compute_denom_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
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
  b = theta[0] + double(centered[0]) * residual;
  for (int c = 1; c <= num_cats; c++)
    b = arma::max(b, theta[c] + double(centered[c]) * residual);

  // ---- 3. Bounds for the FAST power chain ----
  arma::vec pow_bound_low = double(-ref) * residual - b;
  arma::vec pow_bound_high = double(num_cats - ref) * residual - b;
  arma::vec pow_bound = arma::max(arma::abs(pow_bound_low), arma::abs(pow_bound_high));

  // ---- 4. FAST BLOCK ----
  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);

    arma::vec eR = ARMA_MY_EXP(r);
    arma::vec pow = ARMA_MY_EXP(double(-ref) * r - bb);
    arma::vec d = exp_theta[0] * pow;

    for (int c = 1; c <= num_cats; c++) {
      pow %= eR;
      d += exp_theta[c] * pow;
    }
    denom.rows(i0, i1) = d;
  };

  // ---- 5. SAFE BLOCK ----
  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);

    arma::vec d(bb.n_elem, arma::fill::zeros);
    for (int c = 0; c <= num_cats; c++) {
      arma::vec ex = theta[c] + double(centered[c]) * r - bb;
      d += ARMA_MY_EXP(ex);
    }

    denom.rows(i0, i1) = d;
  };

  // ---- 6. BLOCK SCAN ----
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


// =============================================================================
// compute_probs_ordinal
// =============================================================================
arma::mat compute_probs_ordinal(
    const arma::vec& main_param,
    const arma::vec& residual_score,
    const arma::vec& bound,
    int num_cats
) {
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


// =============================================================================
// compute_probs_blume_capel
// =============================================================================
arma::mat compute_probs_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
) {
  constexpr double EXP_BOUND = 709.0;

  const arma::uword N = residual.n_elem;
  arma::mat probs(N, num_cats + 1, arma::fill::none);

  arma::vec cat = arma::regspace<arma::vec>(0, num_cats);
  arma::vec centered = cat - double(ref);
  arma::vec theta = lin_eff * centered + quad_eff * arma::square(centered);
  arma::vec exp_theta = ARMA_MY_EXP(theta);

  b.set_size(N);
  b = theta[0] + double(centered[0]) * residual;
  for (int c = 1; c <= num_cats; ++c) {
    b = arma::max(b, theta[c] + double(centered[c]) * residual);
  }

  arma::vec pow_bound_low = double(-ref) * residual - b;
  arma::vec pow_bound_high = double(num_cats - ref) * residual - b;
  arma::vec pow_bound = arma::max(arma::abs(pow_bound_low), arma::abs(pow_bound_high));

  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    auto P = probs.rows(i0, i1);
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);
    const arma::uword B = bb.n_elem;

    arma::vec eR = ARMA_MY_EXP(r);
    arma::vec pow = ARMA_MY_EXP(double(-ref) * r - bb);
    arma::vec denom(B, arma::fill::zeros);

    arma::vec col0 = exp_theta[0] * pow;
    P.col(0) = col0;
    denom += col0;

    for (int c = 1; c <= num_cats; ++c) {
      pow %= eR;
      arma::vec col = exp_theta[c] * pow;
      P.col(c) = col;
      denom += col;
    }

    P.each_col() /= denom;
  };

  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    auto P = probs.rows(i0, i1);
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);
    const arma::uword B = bb.n_elem;
    arma::vec denom(B, arma::fill::zeros);

    for (int c = 0; c <= num_cats; ++c) {
      arma::vec ex = theta[c] + double(centered[c]) * r - bb;
      arma::vec col = ARMA_MY_EXP(ex);
      P.col(c) = col;
      denom += col;
    }
    P.each_col() /= denom;
  };

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


// =============================================================================
// compute_logZ_and_probs_ordinal
// =============================================================================
LogZAndProbs compute_logZ_and_probs_ordinal(
    const arma::vec& main_param,
    const arma::vec& residual_score,
    const arma::vec& bound,
    int num_cats
) {
  constexpr double EXP_BOUND = 709.0;
  const arma::uword N = bound.n_elem;

  LogZAndProbs result;
  result.probs.set_size(N, num_cats + 1);
  result.log_Z.set_size(N);

  if (num_cats == 1) {
    arma::vec b = arma::clamp(bound, 0.0, arma::datum::inf);
    arma::vec ex = main_param(0) + residual_score - b;
    arma::vec t = ARMA_MY_EXP(ex);
    arma::vec eB = ARMA_MY_EXP(-b);
    arma::vec den = eB + t;
    result.probs.col(1) = t / den;
    result.probs.col(0) = 1.0 - result.probs.col(1);
    result.log_Z = b + ARMA_MY_LOG(den);
    return result;
  }

  const arma::vec eM = ARMA_MY_EXP(main_param);

  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    auto P = result.probs.rows(i0, i1).cols(1, num_cats);
    arma::vec r = residual_score.rows(i0, i1);
    arma::vec b = bound.rows(i0, i1);
    arma::vec eR = ARMA_MY_EXP(r);
    arma::vec eB = ARMA_MY_EXP(-b);
    arma::vec pow = eR;
    arma::vec den = eB;
    for (int c = 0; c < num_cats; c++) {
      arma::vec term = eM[c] * pow % eB;
      P.col(c) = term;
      den += term;
      pow %= eR;
    }
    P.each_col() /= den;
    result.log_Z.rows(i0, i1) = b + ARMA_MY_LOG(den);
  };

  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    auto P = result.probs.rows(i0, i1).cols(1, num_cats);
    arma::vec r = residual_score.rows(i0, i1);
    arma::vec b = arma::clamp(bound.rows(i0, i1), 0.0, arma::datum::inf);
    arma::vec eB = ARMA_MY_EXP(-b);
    arma::vec den = eB;
    for (int c = 0; c < num_cats; c++) {
      arma::vec ex = main_param(c) + (c + 1) * r - b;
      arma::vec t = ARMA_MY_EXP(ex);
      P.col(c) = t;
      den += t;
    }
    P.each_col() /= den;
    result.log_Z.rows(i0, i1) = b + ARMA_MY_LOG(den);
  };

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

  result.probs.col(0) = 1.0 - arma::sum(result.probs.cols(1, num_cats), 1);
  return result;
}


// =============================================================================
// compute_logZ_and_probs_blume_capel
// =============================================================================
LogZAndProbs compute_logZ_and_probs_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b
) {
  constexpr double EXP_BOUND = 709.0;
  const arma::uword N = residual.n_elem;

  LogZAndProbs result;
  result.probs.set_size(N, num_cats + 1);
  result.log_Z.set_size(N);

  arma::vec cat = arma::regspace<arma::vec>(0, num_cats);
  arma::vec centered = cat - double(ref);
  arma::vec theta = lin_eff * centered + quad_eff * arma::square(centered);
  arma::vec exp_theta = ARMA_MY_EXP(theta);

  b.set_size(N);
  b = theta[0] + double(centered[0]) * residual;
  for (int c = 1; c <= num_cats; ++c) {
    b = arma::max(b, theta[c] + double(centered[c]) * residual);
  }

  arma::vec pow_bound_low = double(-ref) * residual - b;
  arma::vec pow_bound_high = double(num_cats - ref) * residual - b;
  arma::vec pow_bound = arma::max(arma::abs(pow_bound_low), arma::abs(pow_bound_high));

  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    auto P = result.probs.rows(i0, i1);
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);
    const arma::uword B = bb.n_elem;

    arma::vec eR = ARMA_MY_EXP(r);
    arma::vec pow = ARMA_MY_EXP(double(-ref) * r - bb);
    arma::vec denom(B, arma::fill::zeros);

    arma::vec col0 = exp_theta[0] * pow;
    P.col(0) = col0;
    denom += col0;

    for (int c = 1; c <= num_cats; ++c) {
      pow %= eR;
      arma::vec col = exp_theta[c] * pow;
      P.col(c) = col;
      denom += col;
    }

    P.each_col() /= denom;
    result.log_Z.rows(i0, i1) = bb + ARMA_MY_LOG(denom);
  };

  auto do_safe_block = [&](arma::uword i0, arma::uword i1) {
    auto P = result.probs.rows(i0, i1);
    arma::vec r = residual.rows(i0, i1);
    arma::vec bb = b.rows(i0, i1);
    const arma::uword B = bb.n_elem;
    arma::vec denom(B, arma::fill::zeros);

    for (int c = 0; c <= num_cats; ++c) {
      arma::vec ex = theta[c] + double(centered[c]) * r - bb;
      arma::vec col = ARMA_MY_EXP(ex);
      P.col(c) = col;
      denom += col;
    }
    P.each_col() /= denom;
    result.log_Z.rows(i0, i1) = bb + ARMA_MY_LOG(denom);
  };

  const double* bp = b.memptr();
  const double* pp = pow_bound.memptr();
  arma::uword i = 0;
  while (i < N) {
    const bool fast_i = (std::abs(bp[i]) <= EXP_BOUND) && (std::abs(pp[i]) <= EXP_BOUND);
    arma::uword j = i + 1;
    while (j < N) {
      const bool fast_j = (std::abs(bp[j]) <= EXP_BOUND) && (std::abs(pp[j]) <= EXP_BOUND);
      if (fast_j != fast_i) break;
      j++;
    }
    if (fast_i) do_fast_block(i, j - 1);
    else do_safe_block(i, j - 1);
    i = j;
  }

  return result;
}
