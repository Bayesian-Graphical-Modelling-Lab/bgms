#include <RcppArmadillo.h>
#include "bgm_helper.h"
#include "bgm_logp_and_grad.h"
#include "common_helpers.h"
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
  constexpr double EXP_BOUND = 709;
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
// Compute a numerically stable sum of the form
//
//    denom = sum_{cat=0}^{num_cats} exp( theta_part(cat) + cat * r - b )
//
// where:
//    theta_part(cat) = lin_eff * (cat - ref) + quad_eff * (cat - ref) ^ 2
//    b = max_cat (theta_part(cat) + cat * r)     // vectorized
//
// Two evaluation modes:
//   * FAST (preexp + power chain): used when residual <= r_cut
//          denom = sum_c exp_theta[c] * exp( -b ) * exp(r)^c
//
//   * SAFE (direct + bound): fallback when exp(r)^c becomes unstable
//          denom = sum_c exp( theta_part(c) + c * r - b )
// -----------------------------------------------------------------------------
inline arma::vec compute_denom_blume_capel(
    const arma::vec& residual,
    const double lin_eff,
    const double quad_eff,
    const int ref,
    const int num_cats,
    arma::vec& b          // update: per-person bound b[i]
) {

  const arma::uword N = residual.n_elem;
  arma::vec denom(N);

  // ---- Precompute theta_part[cat] and exp_theta[cat] ----
  arma::vec cat = arma::regspace<arma::vec>(0, num_cats);
  arma::vec centered = cat - double(ref);
  arma::vec theta = lin_eff * centered + quad_eff * arma::square(centered);
  arma::vec exp_theta = ARMA_MY_EXP(theta);

  // ---- Bound: b[i] = max_cat( theta_part(cat) + cat*r[i] ) ----
  b.set_size(residual.n_elem);
  b.fill(theta[0]);
  for (int c = 1; c <= num_cats; c++)
    b = arma::max(b, theta[c] + double(c) * residual);

  // ---- Cutoff r_cut ----
  const double Ueff = 709.0 - MY_LOG(double(num_cats + 1));
  double r_cut = std::numeric_limits<double>::infinity();
  for (int c = 1; c <= num_cats; c++) {
    const double cand = (Ueff - theta[c]) / double(c);
    if (cand < r_cut) r_cut = cand;
  }

  // ---- FAST BLOCK ----
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

  // ---- SAFE BLOCK ----
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

  // ---- BLOCK SCAN ----
  const double* rp = residual.memptr();
  arma::uword i = 0;
  while (i < N) {
    const bool fast = (rp[i] <= r_cut);
    arma::uword j = i + 1;
    while (j < N) {
      const bool fast_j = (rp[j] <= r_cut);
      if (fast_j != fast) break;
      ++j;
    }
    if (fast) do_fast_block(i, j - 1);
    else do_safe_block(i, j - 1);
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
 * Empirical tests (see dev/numerical_analyses/bgm_regularordinal_normalization.r) showed:
 *   - Clipping necessary for bound < -700
 *   - Bounds improve stability when large
 *
 * Returns:
 *   probs: num_persons × num_cats matrix of probabilities (row-normalized)
 */
inline arma::mat compute_probs_ordinal(const arma::vec& main_param,
                                       const arma::vec& residual_score,
                                       const arma::vec& bound,
                                       int num_cats)
{
  constexpr double EXP_BOUND = 709;
  const arma::uword N = bound.n_elem;

  if (num_cats == 1) {
    arma::vec b = arma::clamp(bound, 0.0, arma::datum::inf);
    arma::vec ex = main_param(0) + residual_score - b;
    arma::vec t = ARMA_MY_EXP(ex);
    arma::vec den = ARMA_MY_EXP(-b) + t;
    arma::mat probs(N, 1, arma::fill::none);
    probs.col(0) = t / den;
    return probs;
  }

  arma::mat probs(N, num_cats, arma::fill::none);
  const arma::vec eM = ARMA_MY_EXP(main_param);

  auto do_fast_block = [&](arma::uword i0, arma::uword i1) {
    auto P = probs.rows(i0, i1);
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
    auto P = probs.rows(i0, i1);
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

  return probs;
}



/**
 * Computes the log-pseudoposterior contribution for a single main-effect parameter (bgm model).
 *
 * For the specified variable, this function evaluates the log-pseudoposterior of either:
 *  - an ordinal threshold parameter (category-specific), or
 *  - a Blume–Capel main-effect parameter (linear or quadratic).
 *
 * The log-pseudoposterior combines:
 *  - Prior contribution: Beta prior on the logistic scale.
 *  - Sufficient statistic contribution: from category counts or Blume–Capel statistics.
 *  - Likelihood contribution: vectorized across all persons using the residual matrix.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (used for ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference category for Blume–Capel centering.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Prior hyperparameters for the Beta prior.
 *  - variable: Index of the variable under consideration.
 *  - category: Category index (ordinal variables only).
 *  - parameter: Parameter index (Blume–Capel only: 0 = linear, 1 = quadratic).
 *
 * Returns:
 *  - The log-pseudoposterior value for the specified parameter.
 *
 * Notes:
 *  - Exactly one of `category` or `parameter` is relevant depending on variable type.
 *  - Uses a numerically stable denominator with exponential bounding.
 *  - This function is used within Metropolis and gradient-based updates.
 */
double log_pseudoposterior_main_effects_component (
    const arma::mat& main_effects,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const int variable,
    const int category,
    const int parameter
) {
  const int num_persons = residual_matrix.n_rows;
  double log_posterior = 0.0;

  auto log_beta_prior = [&](double main_effect_param) {
    return main_effect_param * main_alpha - std::log1p (MY_EXP (main_effect_param)) * (main_alpha + main_beta);
  };

  const int num_cats = num_categories(variable);

  if (is_ordinal_variable(variable)) {
    // Prior contribution + sufficient statistic
    const double value = main_effects(variable, category);
    log_posterior += value * counts_per_category(category + 1, variable);
    log_posterior += log_beta_prior (value);

    arma::vec residual_score = residual_matrix.col (variable);                  // rest scores for all persons
    arma::vec bound = num_cats * residual_score;                                // numerical bound vector
    arma::vec main_effect_param = main_effects.row (variable).cols (0, num_cats - 1).t ();   // main_effect parameters

    arma::vec denom = compute_denom_ordinal(
      residual_score, main_effect_param, bound
    );

    // We then compute the total log-likelihood contribution as:
    //   log_posterior -= bound + log (denom), summed over all persons
    log_posterior -= arma::accu (bound + ARMA_MY_LOG (denom));                    // total contribution
  } else {
    const double value = main_effects(variable, parameter);
    const double linear_main_effect = main_effects(variable, 0);
    const double quadratic_main_effect = main_effects(variable, 1);
    const int ref = baseline_category(variable);

    // Prior contribution + sufficient statistic
    log_posterior += value * blume_capel_stats(parameter, variable);
    log_posterior += log_beta_prior(value);

    if(false) {
      // Vectorized likelihood contribution
      // For each person, we compute the unnormalized log-likelihood denominator:
      //   denom = sum_c exp (θ_lin * c + θ_quad * (c - ref)^2 + c * residual_score - bound)
      // Where:
      //   - θ_lin, θ_quad are linear and quadratic main_effects
      //   - ref is the reference category (used for centering)
      //   - bound = num_cats * residual_score (stabilizes exponentials)
      arma::vec residual_score = residual_matrix.col(variable);                     // rest scores for all persons
      arma::vec bound = num_cats * residual_score;                                // numerical bound vector
      arma::vec denom(num_persons, arma::fill::zeros);                          // initialize denominator
      for (int cat = 0; cat <= num_cats; cat++) {
        int score = cat - ref;                                               // centered category
        double lin_term = linear_main_effect * score;                                      // precompute linear term
        double quad_term = quadratic_main_effect * score * score;                    // precompute quadratic term

        arma::vec exponent = lin_term + quad_term + score * residual_score - bound;
        denom += ARMA_MY_EXP (exponent);                                           // accumulate over categories
      }

      // The final log-likelihood contribution is then:
      //   log_posterior -= bound + log (denom), summed over all persons
      log_posterior -= arma::accu (bound + ARMA_MY_LOG (denom));                    // total contribution

    } else {
      arma::vec residual_score = residual_matrix.col(variable);                     // rest scores for all persons
      arma::vec bound(residual_score.n_elem);



      arma::vec denom = compute_denom_blume_capel(
          residual_score, linear_main_effect, quadratic_main_effect, ref,
          num_cats, bound
      );



      // residual_score is r_i per person
      log_posterior -= arma::accu(bound + ARMA_MY_LOG(denom));
      log_posterior += ref * arma::accu(residual_score);

    }
  }

  return log_posterior;
}



/**
 * Computes the log-pseudoposterior contribution for a single pairwise interaction (bgm model).
 *
 * The contribution consists of:
 *  - Sufficient statistic term: interaction × pairwise count.
 *  - Likelihood term: summed over all observations, using either
 *    * ordinal thresholds, or
 *    * Blume–Capel quadratic/linear main effects.
 *  - Prior term: Cauchy prior on the interaction coefficient (if active).
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters.
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Sufficient statistics for pairwise counts.
 *  - var1, var2: Indices of the variable pair being updated.
 *
 * Returns:
 *  - The log-pseudoposterior value for the specified interaction parameter.
 *
 * Notes:
 *  - Bounds are applied for numerical stability in exponential terms.
 *  - The function assumes that `pairwise_effects` is symmetric.
 *  - Used within Metropolis and gradient-based updates of pairwise effects.
 */
double log_pseudoposterior_interactions_component (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const arma::imat& pairwise_stats,
    const int var1,
    const int var2
) {
  const int num_observations = observations.n_rows;

  double log_pseudo_posterior = 2.0 * pairwise_effects(var1, var2) * pairwise_stats(var1, var2);

  for (int var : {var1, var2}) {
    int num_cats = num_categories (var);

    // Compute rest score: contribution from other variables
    arma::vec residual_score = observations * pairwise_effects.col (var);
    arma::vec denominator = arma::zeros (num_observations);

    if (is_ordinal_variable (var)) {
      arma::vec bound = num_cats * residual_score;                                // numerical bound vector
      arma::vec main_effect_param = main_effects.row (var).cols (0, num_cats - 1).t ();   // main_effect parameters
      denominator += compute_denom_ordinal(
        residual_score, main_effect_param, bound
      );
      // Subtract log partition function and bounds adjustment
      log_pseudo_posterior -= arma::accu (ARMA_MY_LOG (denominator));
      log_pseudo_posterior -= arma::accu (bound);
    } else {
      const int ref = baseline_category (var);
      if (false) {
        arma::vec bound = num_cats * residual_score;                                // numerical bound vector
        // Binary/categorical variable: quadratic + linear term
        for (int category = 0; category <= num_cats; category++) {
          int score = category - ref;
          double lin_term = main_effects (var, 0) * score;
          double quad_term = main_effects (var, 1) * score * score;
          arma::vec exponent = lin_term + quad_term + score * residual_score - bound;
          denominator += ARMA_MY_EXP (exponent);
        }
        // Subtract log partition function and bounds adjustment
        log_pseudo_posterior -= arma::accu (ARMA_MY_LOG (denominator));
        log_pseudo_posterior -= arma::accu (bound);
      } else {
        arma::vec bound(residual_score.n_elem);



        denominator = compute_denom_blume_capel(
          residual_score, main_effects (var, 0), main_effects (var, 1), ref,
          num_cats, bound
        );



        // residual_score is r_i per person
        log_pseudo_posterior -= arma::accu(bound + ARMA_MY_LOG(denominator));
        log_pseudo_posterior += ref * arma::accu(residual_score);
      }
    }
  }

  // Add Cauchy prior terms for included pairwise effects
  if (inclusion_indicator (var1, var2) == 1) {
    log_pseudo_posterior += R::dcauchy (pairwise_effects (var1, var2), 0.0, pairwise_scale, true);
  }

  return log_pseudo_posterior;
}



/**
 * Computes the full log-pseudoposterior for the bgm model.
 *
 * The log-pseudoposterior combines:
 *  - Main-effect contributions:
 *    * Ordinal variables: one parameter per category with Beta prior.
 *    * Blume–Capel variables: linear and quadratic parameters with Beta priors.
 *  - Pairwise-effect contributions:
 *    * Included interactions (per inclusion_indicator) with Cauchy prior.
 *  - Likelihood contributions:
 *    * Vectorized over all persons, with numerically stabilized denominators.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (for ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Hyperparameters for the Beta priors.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Pairwise sufficient statistics.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *
 * Returns:
 *  - The scalar log-pseudoposterior value for the full model.
 *
 * Notes:
 *  - Exponential terms are bounded with nonnegative `bound` values for stability.
 *  - Pairwise effects are included only when marked in `inclusion_indicator`.
 *  - This is the top-level function combining both main and interaction components.
 */
double log_pseudoposterior (
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const arma::imat& pairwise_stats,
    const arma::mat& residual_matrix
) {

  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;

  double log_pseudoposterior = 0.0;

  // Calculate the contribution from the data and the prior
  auto log_beta_prior = [&](double main_effect_param) {
    return main_effect_param * main_alpha - std::log1p (MY_EXP (main_effect_param)) * (main_alpha + main_beta);
  };

  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        double value = main_effects(variable, cat);
        log_pseudoposterior += counts_per_category(cat + 1, variable) * value;
        log_pseudoposterior += log_beta_prior(value);
      }
    } else {
      double value = main_effects(variable, 0);
      log_pseudoposterior += log_beta_prior(value);
      log_pseudoposterior += blume_capel_stats(0, variable) * value;

      value = main_effects(variable, 1);
      log_pseudoposterior += log_beta_prior(value);
      log_pseudoposterior += blume_capel_stats(1, variable) * value;
    }
  }
  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      if (inclusion_indicator(var1, var2) == 0) continue;

      double value = pairwise_effects(var1, var2);
      log_pseudoposterior += 2.0 * pairwise_stats(var1, var2) * value;
      log_pseudoposterior += R::dcauchy(value, 0.0, pairwise_scale, true); // Cauchy prior
    }
  }

  // Calculate the log denominators
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);
    arma::vec residual_score = residual_matrix.col (variable);                  // rest scores for all persons

    arma::vec denom(num_persons, arma::fill::zeros);
    if (is_ordinal_variable(variable)) {
      arma::vec bound = num_cats * residual_score;                              // numerical bound vector
      arma::vec main_effect_param = main_effects.row (variable).cols (0, num_cats - 1).t ();   // main_effect parameters for variable
      denom += compute_denom_ordinal(
        residual_score, main_effect_param, bound
      );
    } else {
      const int ref = baseline_category(variable);
      const double lin_effect = main_effects(variable, 0);
      const double quad_effect = main_effects(variable, 1);

      if(false) {
        // ----
        arma::vec bound = num_cats * residual_score;                              // numerical bound vector

        const int score_min = -ref;
        const int score_max = num_cats - ref;
        const int max_diff = (std::abs(score_min) >= std::abs(score_max)) ? score_min : score_max;
        bound = max_diff * residual_score;                              // numerical bound vector

        double main_bound = lin_effect * score_min + quad_effect * score_min * score_min;
        for (int cat = 1; cat <= num_cats; cat++) {
          const int score = cat - ref;
          const double tmp = lin_effect * score + quad_effect * score * score;
          if (std::abs(tmp) > std::abs(main_bound)) main_bound = tmp;
        }
        bound += main_bound;  // final bound adjustment
        // ----

        for (int cat = 0; cat <= num_cats; cat++) {
          int score = cat - ref;                                                  // centered category
          double lin = lin_effect * score;                                        // precompute linear term
          double quad = quad_effect * score * score;                              // precompute quadratic term
          arma::vec exponent = lin + quad + score * residual_score - bound;
          denom += ARMA_MY_EXP (exponent);                                        // accumulate over categories
        }

        log_pseudoposterior -= arma::accu (bound + ARMA_MY_LOG (denom));            // total contribution

      } else {
        arma::vec bound(residual_score.n_elem);



        denom = compute_denom_blume_capel(
          residual_score, lin_effect, quad_effect, ref, num_cats, bound
        );



        // residual_score is r_i per person
        log_pseudoposterior -= arma::accu(bound + ARMA_MY_LOG(denom));
        log_pseudoposterior += ref * arma::accu(residual_score);
      }
    }
  }

  return log_pseudoposterior;
}



std::pair<arma::vec, arma::imat> gradient_observed_active(
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const arma::imat& pairwise_stats
) {
  const int num_variables = observations.n_cols;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::imat index_matrix(num_variables, num_variables, arma::fill::zeros);

  // Count active pairwise effects + Index map for pairwise parameters
  int num_active = 0;
  for (int i = 0; i < num_variables - 1; i++) {
    for (int j = i + 1; j < num_variables; j++) {
      if (inclusion_indicator(i, j) == 1){
        index_matrix(i, j) = num_main + num_active++;
        index_matrix(j, i) = index_matrix(i, j);
      }
    }
  }

  // Allocate gradient vector (main + active pairwise only)
  arma::vec gradient(num_main + num_active, arma::fill::zeros);

  // ---- STEP 1: Observed statistics ----
  int offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) = counts_per_category(cat + 1, variable);
      }
      offset += num_cats;
    } else {
      gradient(offset) = blume_capel_stats(0, variable);
      gradient(offset + 1) = blume_capel_stats(1, variable);
      offset += 2;
    }
  }
  for (int i = 0; i < num_variables - 1; i++) {
    for (int j = i + 1; j < num_variables; j++) {
      if (inclusion_indicator(i, j) == 0) continue;
      int location = index_matrix(i, j);
      gradient(location) = 2.0 * pairwise_stats(i, j);
    }
  }

  return {gradient, index_matrix};
}



/**
 * Computes the gradient of the log-pseudoposterior for main and active pairwise parameters.
 *
 * Gradient components:
 *  - Observed sufficient statistics (from counts_per_category, blume_capel_stats, pairwise_stats).
 *  - Minus expected sufficient statistics (computed via probabilities over categories).
 *  - Plus gradient contributions from priors:
 *    * Beta priors on main effects.
 *    * Cauchy priors on active pairwise effects.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (for ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Hyperparameters for Beta priors.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *
 * Returns:
 *  - A vector containing the gradient of the log-pseudoposterior with respect to
 *    all main and active pairwise parameters, in the same order as
 *    `vectorize_model_parameters_bgm()`.
 */
arma::vec gradient_log_pseudoposterior(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const arma::mat& residual_matrix,
    const arma::imat index_matrix,
    const arma::vec grad_obs
) {
  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;

  // Allocate gradient vector (main + active pairwise only)
  arma::vec gradient = grad_obs;

  // ---- STEP 2: Expected statistics ----
  int offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_cats = num_categories(variable);
    arma::vec residual_score = residual_matrix.col(variable);

    if (is_ordinal_variable(variable)) {
      arma::vec bound = num_cats * residual_score;
      arma::vec main_param = main_effects.row(variable).cols(0, num_cats - 1).t();
      arma::mat probs = compute_probs_ordinal(
        main_param, residual_score, bound, num_cats
      );

      // main effects
      for (int cat = 0; cat < num_cats; cat++) {
        gradient(offset + cat) -= arma::accu(probs.col(cat));
      }

      // pairwise effects
      for (int j = 0; j < num_variables; j++) {
        if (inclusion_indicator(variable, j) == 0 || variable == j) continue;
        arma::vec expected_value = arma::zeros(num_persons);
        for (int cat = 0; cat < num_cats; cat++) {
          expected_value += (cat + 1) * probs.col(cat) % observations.col(j);
        }
        int location = (variable < j) ? index_matrix(variable, j) : index_matrix(j, variable);
        gradient(location) -= arma::accu(expected_value);
      }
      offset += num_cats;
    } else {

      const int ref = baseline_category(variable);
      const double lin_eff = main_effects(variable, 0);
      const double quad_eff = main_effects(variable, 1);

      // Compute bounds
      // const int score_min = -ref;
      // const int score_max = num_cats - ref;
      // const int max_diff  = (std::abs(score_min) >= std::abs(score_max)) ? score_min : score_max;
      // arma::vec bound = max_diff * residual_score;
      //
      // double main_bound = lin_eff * score_min + quad_eff * score_min * score_min;
      // for (int cat = 1; cat <= num_cats; cat++) {
      //   const int score = cat - ref;
      //   const double tmp = lin_eff * score + quad_eff * score * score;
      //   if (std::abs(tmp) > std::abs(main_bound)) main_bound = tmp;
      // }
      // bound += main_bound;

      // arma::mat exponents(num_persons, num_cats + 1);
      // for (int cat = 0; cat <= num_cats; cat++) {
      //   int score = cat - ref;
      //   double lin  = lin_eff * score;
      //   double quad = quad_eff * score * score;
      //   exponents.col(cat) = lin + quad + score * residual_score - bound;
      // }
      // arma::mat probs = ARMA_MY_EXP(exponents);
      // arma::vec denom = arma::sum(probs, 1);
      // probs.each_col() /= denom;

      // arma::ivec lin_score = arma::regspace<arma::ivec>(0 - ref, num_cats - ref);
      // arma::ivec quad_score = arma::square(lin_score);

      // Compute exponents
      arma::vec scores = arma::regspace<arma::vec>(0 - ref, num_cats - ref);
      arma::rowvec offsets = lin_eff * scores.t() + quad_eff * arma::square(scores.t());
      arma::mat exponents = residual_score * scores.t();
      exponents.each_row() += offsets;
      //exponents.each_col() -= bound;
      arma::vec row_max = arma::max(exponents, /*dim=*/1);
      exponents.each_col() -= row_max;

      // Compute probabilities
      arma::mat probs = ARMA_MY_EXP(exponents);
      arma::vec denom = arma::sum(probs, 1);
      probs.each_col() /= denom;

      arma::ivec lin_score = arma::conv_to<arma::ivec>::from(scores);
      arma::ivec quad_score = arma::square(lin_score);

      // main effects
      gradient(offset) -= arma::accu(probs * lin_score);
      gradient(offset + 1) -= arma::accu(probs * quad_score);

      // pairwise effects
      for (int j = 0; j < num_variables; j++) {
        if (inclusion_indicator(variable, j) == 0 || variable == j) continue;
        arma::vec expected_value = arma::zeros(num_persons);
        for (int cat = 0; cat <= num_cats; cat++) {
          int score = cat - ref;
          expected_value += score * probs.col(cat) % observations.col(j);
        }
        int location = (variable < j) ? index_matrix(variable, j) : index_matrix(j, variable);
        gradient(location) -= arma::accu(expected_value);
      }
      offset += 2;
    }
  }

  // ---- STEP 3: Priors ----
  offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    if (is_ordinal_variable(variable)) {
      const int num_cats = num_categories(variable);
      for (int cat = 0; cat < num_cats; cat++) {
        const double p = 1.0 / (1.0 + MY_EXP(-main_effects(variable, cat)));
        gradient(offset + cat) += main_alpha - (main_alpha + main_beta) * p;
      }
      offset += num_cats;
    } else {
      for (int k = 0; k < 2; k++) {
        const double param = main_effects(variable, k);
        const double p = 1.0 / (1.0 + MY_EXP(-param));
        gradient(offset + k) += main_alpha - (main_alpha + main_beta) * p;
      }
      offset += 2;
    }
  }
  for (int i = 0; i < num_variables - 1; i++) {
    for (int j = i + 1; j < num_variables; j++) {
      if (inclusion_indicator(i, j) == 0) continue;
      int location = index_matrix(i, j);
      const double effect = pairwise_effects(i, j);
      gradient(location) -= 2.0 * effect / (effect * effect + pairwise_scale * pairwise_scale);
    }
  }
  return gradient;
}



/**
 * Computes the log-likelihood ratio for updating a single variable’s parameter (bgm model).
 *
 * The ratio compares the likelihood of the current parameter state versus a proposed state,
 * given the observed data, main effects, and residual contributions from other variables.
 *
 * Calculation:
 *  - Removes the current interaction contribution from the residual scores.
 *  - Recomputes denominators of the softmax likelihood under both current and proposed states.
 *  - Returns the accumulated log difference across all persons.
 *
 * Inputs:
 *  - variable: Index of the variable being updated.
 *  - interacting_score: Integer vector of interaction scores with other variables.
 *  - proposed_state: Candidate value for the parameter being updated.
 *  - current_state: Current value of the parameter.
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - num_categories: Number of categories per variable.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *
 * Returns:
 *  - The log-likelihood ratio (current − proposed) for the specified variable.
 *
 * Notes:
 *  - For ordinal variables, the likelihood includes thresholds and category-specific scores.
 *  - For Blume–Capel variables, linear and quadratic terms are applied with centered categories.
 *  - Bounds are used to stabilize exponentials in the softmax denominator.
 */
double compute_log_likelihood_ratio_for_variable (
    int variable,
    const arma::ivec& interacting_score,
    double proposed_state,
    double current_state,
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::mat& residual_matrix,
    const arma::imat& observations,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category
) {
  // Convert interaction score vector to double precision
  arma::vec interaction = arma::conv_to<arma::vec>::from (interacting_score);

  const int num_persons = residual_matrix.n_rows;
  const int num_cats = num_categories (variable);

  // Compute adjusted linear predictors without the current interaction
  arma::vec residual_score = residual_matrix.col (variable) - interaction * current_state;
  arma::vec bounds = residual_score * num_cats;

  arma::vec denom_current = arma::zeros (num_persons);
  arma::vec denom_proposed = arma::zeros (num_persons);

  if (is_ordinal_variable (variable)) {
    arma::vec main_param = main_effects.row(variable).cols(0, num_cats - 1).t();

    // ---- main change: use safe helper ----
    denom_current += compute_denom_ordinal(
      residual_score + interaction * current_state, main_param, bounds
    );
    denom_proposed += compute_denom_ordinal(
      residual_score + interaction * proposed_state, main_param, bounds
    );

  } else {
    // Binary or categorical variable: linear + quadratic score
    const int ref = baseline_category (variable);

    if(false) {
      for (int category = 0; category <= num_cats; category++) {
        int score = category - ref;
        double lin_term = main_effects (variable, 0) * score;
        double quad_term = main_effects (variable, 1) * score * score;
        arma::vec exponent = lin_term + quad_term + score * residual_score - bounds;

        denom_current += ARMA_MY_EXP (exponent + score * interaction * current_state);
        denom_proposed += ARMA_MY_EXP (exponent + score * interaction * proposed_state);
      }
    } else {
      arma::vec residual_current = residual_score + interaction * current_state;
      arma::vec residual_proposed = residual_score + interaction * proposed_state;

      arma::vec b_cur(residual_score.n_elem);
      arma::vec b_prp(residual_score.n_elem);

      arma::vec denom_cur  = compute_denom_blume_capel(
        residual_current, main_effects(variable, 0), main_effects(variable, 1),
        ref, num_cats, b_cur
      );

      arma::vec denom_prp = compute_denom_blume_capel(
        residual_proposed,  main_effects(variable, 0), main_effects(variable, 1),
        ref, num_cats, b_prp
      );

      arma::vec logZ_cur = b_cur + ARMA_MY_LOG(denom_cur) -
        double(ref) * residual_current;
      arma::vec logZ_prp = b_prp + ARMA_MY_LOG(denom_prp) -
        double(ref) * residual_proposed;

      return arma::accu(logZ_cur - logZ_prp);
    }
  }

  // Accumulated log-likelihood difference across persons
  return arma::accu (ARMA_MY_LOG (denom_current) - ARMA_MY_LOG (denom_proposed));
}



/**
 * Computes the log-pseudolikelihood ratio for updating a single pairwise interaction (bgm model).
 *
 * The ratio compares the pseudo-likelihood under a proposed value versus the current value
 * of an interaction parameter between two variables.
 *
 * Calculation:
 *  1. Direct contribution from the interaction term:
 *       Δβ × ∑(score_var1 × score_var2) over all persons,
 *     where Δβ = proposed_state − current_state.
 *  2. Change in pseudo-likelihood for variable1, accounting for its updated interaction with variable2.
 *  3. Symmetric change in pseudo-likelihood for variable2, accounting for its updated interaction with variable1.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of pairwise interaction parameters.
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - num_persons: Number of persons (rows of observations).
 *  - variable1, variable2: Indices of the variable pair being updated.
 *  - proposed_state: Candidate value for the interaction parameter.
 *  - current_state: Current value of the interaction parameter.
 *  - residual_matrix: Matrix of residual scores (persons × variables).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - pairwise_stats: Sufficient statistics for pairwise counts.
 *
 * Returns:
 *  - The log-pseudolikelihood ratio (current − proposed) for the specified interaction.
 *
 * Notes:
 *  - Calls `compute_log_likelihood_ratio_for_variable()` for both variables involved.
 *  - The factor of 2.0 in the direct term reflects symmetry of the pairwise sufficient statistic.
 *  - Used in Metropolis updates for pairwise interaction parameters.
 */
double log_pseudolikelihood_ratio_interaction (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const int num_persons,
    const int variable1,
    const int variable2,
    const double proposed_state,
    const double current_state,
    const arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& pairwise_stats
) {
  double log_ratio = 0.0;
  const double delta = proposed_state - current_state;

  // Extract score vectors for both variables across all persons
  arma::ivec score1 = observations.col(variable1);
  arma::ivec score2 = observations.col(variable2);

  // (1) Direct interaction contribution to the linear predictor:
  //     Δβ × ∑(score1_i × score2_i) for all persons i
  log_ratio += 2.0 * pairwise_stats(variable1, variable2) * delta;

  // (2) Change in pseudo-likelihood for variable1 due to the update in its interaction with variable2
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable1, score2, proposed_state, current_state, main_effects,
    num_categories, residual_matrix, observations, is_ordinal_variable,
    baseline_category
  );

  // (3) Symmetric change for variable2 due to its interaction with variable1
  log_ratio += compute_log_likelihood_ratio_for_variable (
    variable2, score1, proposed_state, current_state, main_effects,
    num_categories, residual_matrix, observations, is_ordinal_variable,
    baseline_category
  );

  return log_ratio;
}