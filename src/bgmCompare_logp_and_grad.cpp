#include <RcppArmadillo.h>
#include "bgmCompare_helper.h"
#include "bgmCompare_logp_and_grad.h"
#include <cmath>
#include "explog_switch.h"

using namespace Rcpp;



/**
 * Computes the log pseudoposterior for the bgmCompare model.
 *
 * The pseudoposterior combines:
 *  - Data contributions (main effects, pairwise effects, normalizing constants),
 *    evaluated separately for each group.
 *  - Prior contributions on main effects, pairwise effects, and group differences.
 *
 * Procedure:
 *  - For each group:
 *    * Construct group-specific main and pairwise effects using
 *      `compute_group_main_effects()` and `compute_group_pairwise_effects()`.
 *    * Add linear contributions from sufficient statistics.
 *    * Add quadratic contributions from pairwise sufficient statistics.
 *    * Subtract log normalizing constants computed from residual scores.
 *  - Add prior contributions:
 *    * Logistic–Beta prior for all main-effect baselines.
 *    * Cauchy priors for main-effect group differences (if active).
 *    * Cauchy priors for pairwise effects (baseline + group differences if active).
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
 *  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable in main_effects.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - observations: Observation matrix (persons × variables).
 *  - group_indices: Matrix of row ranges [start,end] for each group in observations.
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category_group: Per-group category counts (for ordinal variables).
 *  - blume_capel_stats_group: Per-group sufficient statistics (for Blume–Capel variables).
 *  - pairwise_stats_group: Per-group pairwise sufficient statistics.
 *  - num_groups: Number of groups.
 *  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - interaction_scale: Scale parameter for Cauchy priors on baseline pairwise effects.
 *  - difference_scale: Scale parameter for Cauchy priors on group differences.
 *
 * Returns:
 *  - The scalar value of the log pseudoposterior.
 *
 * Notes:
 *  - This function generalizes the bgm pseudoposterior to multiple groups.
 *  - Group differences are only included if marked in inclusion_indicator.
 *  - Residual scores are recomputed separately for each group.
 */
double log_pseudoposterior(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const double difference_scale
) {
  const int num_variables = observations.n_cols;
  const int max_num_categories = num_categories.max();
  double log_pp = 0.0;

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::imat counts_per_category = counts_per_category_group[group];
    const arma::imat blume_capel_stats = blume_capel_stats_group[group];

    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);

    const arma::vec proj_g = projection.row(group).t(); // length = num_groups-1

    // ---- build group-specific main & pairwise effects ----
    for (int v = 0; v < num_variables; ++v) {
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );

      // store into row v (padded with zeros if variable has < max_num_categories params)
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();

      // upper triangle incl. base value; mirror to keep symmetry
      for (int u = v; u < num_variables; ++u) { // Combines with loop over v
        if(u == v) continue;
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }

    // ---- data contribution pseudolikelihood (linear terms) ----
      const int num_cats = num_categories(v);
      if (is_ordinal_variable(v)) {
        // use group-specific main_effects
        for (int c = 0; c < num_cats; ++c) {
          const double val = main_group(v, c);
          log_pp += static_cast<double>(counts_per_category(c, v)) * val;
        }
      } else {
        // two sufficient stats for binary-ish coding? keep original shape
        log_pp += static_cast<double>(blume_capel_stats(0, v)) * main_group(v, 0);
        log_pp += static_cast<double>(blume_capel_stats(1, v)) * main_group(v, 1);
      }
    }

    // ---- data contribution pseudolikelihood (quadratic terms) ----
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    const arma::mat pairwise_stats = pairwise_stats_group[group];

    log_pp += arma::accu(pairwise_group % pairwise_stats); // trace(X' * W * X) = sum(W %*% (X'X))

    // ---- pseudolikelihood normalizing constants (per variable) ----
    const arma::mat residual_matrix = obs * pairwise_group;
    for (int v = 0; v < num_variables; ++v) {
      const int num_cats = num_categories(v);
      const arma::vec rest_score = residual_matrix.col(v);

      // bound to stabilize exp; use group-specific params consistently
      arma::vec bound = num_cats * rest_score;
      bound = arma::clamp(bound, 0.0, arma::datum::inf);

      arma::vec denom(rest_score.n_elem, arma::fill::zeros);

      if (is_ordinal_variable(v)) {
        // base term exp(-bound)
        denom = ARMA_MY_EXP(-bound);
        // main_effects from main_group
        for (int c = 0; c < num_cats; ++c) {
          const double th = main_group(v, c);
          const arma::vec exponent = th + (c + 1) * rest_score - bound;
          denom += ARMA_MY_EXP(exponent);
        }
      } else {
        // linear/quadratic main effects from main_group
        const double lin_effect  = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        const int ref = baseline_category(v);
        for (int c = 0; c <= num_cats; ++c) {
          const int centered = c - ref;
          const double quad = quad_effect * centered * centered;
          const double lin  = lin_effect * c;
          const arma::vec exponent = lin + quad + c * rest_score - bound;
          denom += ARMA_MY_EXP(exponent);
        }
      }
      // - sum_i [ bound_i + log denom_i ]
      log_pp -= arma::accu(bound + ARMA_MY_LOG(denom));
    }
  }

  // ---- priors ----
  auto log_beta_prior = [&](double x) {
    return x * main_alpha - std::log1p(MY_EXP(x)) * (main_alpha + main_beta);
  };

  // Main effects prior
  for (int v = 0; v < num_variables; ++v) {
    const int row0 = main_effect_indices(v, 0);
    const int row1 = main_effect_indices(v, 1);
    for (int r = row0; r <= row1; ++r) {
      log_pp += log_beta_prior(main_effects(r, 0));

      if (inclusion_indicator(v, v) == 0) continue;
      for (int eff = 1; eff < num_groups; ++eff) {
        log_pp += R::dcauchy(main_effects(r, eff), 0.0, difference_scale, true);
      }
    }
  }

  // Pairwise effects prior
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      const int idx = pairwise_effect_indices(v1, v2);
      log_pp += R::dcauchy(pairwise_effects(idx, 0), 0.0, interaction_scale, true);

      if (inclusion_indicator(v1, v2) == 0) continue;
      for (int eff = 1; eff < num_groups; ++eff) {
        log_pp += R::dcauchy(pairwise_effects(idx, eff), 0.0, difference_scale, true);
      }
    }
  }

  return log_pp;
}



/**
 * Computes the gradient of the log pseudoposterior for the bgmCompare model.
 *
 * The gradient combines contributions from:
 *  - Observed sufficient statistics (counts, Blume–Capel stats, pairwise stats).
 *  - Expected sufficient statistics under the model (via softmax probabilities).
 *  - Priors on main effects, pairwise effects, and group differences.
 *
 * Procedure:
 *  - Loop over groups:
 *    * Add observed sufficient statistics to gradient.
 *    * Build group-specific main and pairwise effects
 *      (`compute_group_main_effects()`, `compute_group_pairwise_effects()`).
 *    * Compute expected sufficient statistics from residual scores and subtract them.
 *  - Add prior contributions:
 *    * Logistic–Beta prior gradient for main-effect baselines.
 *    * Cauchy prior gradient for main-effect differences and pairwise effects.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
 *  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - observations: Observation matrix (persons × variables).
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category_group: Per-group category counts (for ordinal variables).
 *  - blume_capel_stats_group: Per-group sufficient statistics (for Blume–Capel variables).
 *  - pairwise_stats_group: Per-group pairwise sufficient statistics.
 *  - num_groups: Number of groups.
 *  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - interaction_scale: Scale parameter for Cauchy prior on baseline pairwise effects.
 *  - difference_scale: Scale parameter for Cauchy prior on group differences.
 *  - main_index: Index map for main-effect parameters (from build_index_maps()).
 *  - pair_index: Index map for pairwise-effect parameters (from build_index_maps()).
 *
 * Returns:
 *  - A vector containing the gradient of the log pseudoposterior
 *    with respect to all active parameters (layout matches vectorization).
 *
 * Notes:
 *  - Must be consistent with `vectorize_model_parameters_bgmcompare()` and
 *    `unvectorize_model_parameters_bgmcompare()`.
 *  - Uses logistic–Beta gradient for main baselines, and Cauchy gradients for differences.
 *  - Index maps (`main_index`, `pair_index`) are required to map parameters
 *    back into the flat gradient vector.
 */
arma::vec gradient(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double main_alpha,
    const double main_beta,
    const double interaction_scale,
    const double difference_scale,
    const arma::imat& main_index,
    const arma::imat& pair_index
) {
  const int num_variables      = observations.n_cols;
  const int max_num_categories = num_categories.max();
  const int n_main_rows        = main_effects.n_rows;
  const int n_pair_rows        = pairwise_effects.n_rows;

  // total length (must match vectorize_model_parameters_bgmcompare)
  long long total_len = 0;
  total_len += n_main_rows;      // main col 0
  total_len += n_pair_rows;      // pair col 0
  for (int v = 0; v < num_variables; ++v) {
    if (inclusion_indicator(v, v) == 1) {
      const int r0 = main_effect_indices(v, 0);
      const int r1 = main_effect_indices(v, 1);
      total_len += static_cast<long long>(r1 - r0 + 1) * (num_groups - 1);
    }
  }
  for (int v2 = 0; v2 < num_variables - 1; ++v2) {
    for (int v1 = v2 + 1; v1 < num_variables; ++v1) {
      total_len += (inclusion_indicator(v1, v2) == 1) * (num_groups - 1);
    }
  }

  arma::vec grad(total_len, arma::fill::zeros);
  int off;

  // -------------------------------------------------
  // Allocate temporaries ONCE (reused inside loops)
  // -------------------------------------------------
  arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);
  arma::mat exponents; // resized per variable
  arma::vec lin_score, quad_score;

  // -------------------------------
  // Observed sufficient statistics
  // -------------------------------
  for (int g = 0; g < num_groups; ++g) {
    const arma::imat counts_per_category = counts_per_category_group[g];
    const arma::imat blume_capel_stats   = blume_capel_stats_group[g];
    const arma::mat pairwise_stats       = pairwise_stats_group[g];

    // Main effects
    for (int v = 0; v < num_variables; ++v) {
      const int base     = main_effect_indices(v, 0);
      const int num_cats = num_categories(v);

      if (is_ordinal_variable(v)) {
        for (int c = 0; c < num_cats; ++c) {
          off = main_index(base + c, 0);
          grad(off) += counts_per_category(c, v);

          if (inclusion_indicator(v, v) != 0) {
            const int cnt = counts_per_category(c, v);
            for (int k = 1; k < num_groups; ++k) {
              off = main_index(base + c, k);
              grad(off) += cnt * projection(g, k - 1);
            }
          }
        }
      } else {
        // overall (2 stats)
        off = main_index(base, 0);
        grad(off) += blume_capel_stats(0, v);

        off = main_index(base + 1, 0);
        grad(off) += blume_capel_stats(1, v);

        if (inclusion_indicator(v, v) != 0) {
          for (int k = 1; k < num_groups; ++k) {
            off = main_index(base, k);
            grad(off) += blume_capel_stats(0, v) * projection(g, k - 1);

            off = main_index(base + 1, k);
            grad(off) += blume_capel_stats(1, v) * projection(g, k - 1);
          }
        }
      }
    }

    // Pairwise (observed)
    for (int v1 = 0; v1 < num_variables - 1; ++v1) {
      for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
        const int row = pairwise_effect_indices(v1, v2);

        off = pair_index(row, 0);
        grad(off) += 2.0 * pairwise_stats(v1, v2); // upper tri counted once

        if (inclusion_indicator(v1, v2) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = pair_index(row, k);
          grad(off) += 2.0 * pairwise_stats(v1, v2) * projection(g, k - 1);
        }
      }
    }
  }

  // --------------------------------
  // Expected sufficient statistics
  // --------------------------------
  for (int g = 0; g < num_groups; ++g) {
    const int r0 = group_indices(g, 0);
    const int r1 = group_indices(g, 1);

    main_group.zeros();
    pairwise_group.zeros();
    arma::vec proj_g = projection.row(g).t(); // length = num_groups-1

    // build group-specific params
    for (int v = 0; v < num_variables; ++v) {
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();

      for (int u = v + 1; u < num_variables; ++u) {
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }
    }

    // group slice and rest matrix
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    const arma::mat residual_matrix = obs * pairwise_group;
    const int num_group_obs = obs.n_rows;

    for (int v = 0; v < num_variables; ++v) {
      const int K   = num_categories(v);
      const int ref = baseline_category(v);

      arma::vec rest_score = residual_matrix.col(v);
      arma::vec bound      = K * rest_score;
      bound.clamp(0.0, arma::datum::inf);

      exponents.set_size(num_group_obs, K + 1);
      exponents.zeros();

      if (is_ordinal_variable(v)) {
        exponents.col(0) -= bound;
        for (int j = 0; j < K; ++j) {
          exponents.col(j + 1) = main_group(v, j) + (j + 1) * rest_score - bound;
        }
      } else {
        const double lin_effect  = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        for (int s = 0; s <= K; ++s) {
          const int centered = s - ref;
          const double lin  = lin_effect * s;
          const double quad = quad_effect * centered * centered;
          exponents.col(s) = lin + quad + s * rest_score - bound;
        }
      }

      // exponentiate + normalize row-wise (in-place)
      for (arma::uword i = 0; i < exponents.n_rows; ++i) {
        double denom = 0.0;
        for (arma::uword j = 0; j < exponents.n_cols; ++j) {
          double val = MY_EXP(exponents(i, j));
          exponents(i, j) = val;
          denom += val;
        }
        double inv = 1.0 / denom;
        for (arma::uword j = 0; j < exponents.n_cols; ++j)
          exponents(i, j) *= inv;
      }

      // ---- MAIN expected ----
      const int base = main_effect_indices(v, 0);
      if (is_ordinal_variable(v)) {
        for (int s = 1; s <= K; ++s) {
          const int j = s - 1;
          double sum_col_s = arma::accu(exponents.col(s));

          off = main_index(base + j, 0);
          grad(off) -= sum_col_s;

          if (inclusion_indicator(v, v) == 0) continue;
          for (int k = 1; k < num_groups; ++k) {
            off = main_index(base + j, k);
            grad(off) -= projection(g, k - 1) * sum_col_s;
          }
        }
      } else {
        if (lin_score.n_elem != K + 1) {
          lin_score  = arma::regspace<arma::vec>(0, K);
          quad_score = arma::square(lin_score - ref);
        }

        double sum_lin  = arma::accu(exponents * lin_score);
        double sum_quad = arma::accu(exponents * quad_score);

        off = main_index(base, 0);
        grad(off) -= sum_lin;
        off = main_index(base + 1, 0);
        grad(off) -= sum_quad;

        if (inclusion_indicator(v, v) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = main_index(base, k);
          grad(off) -= projection(g, k - 1) * sum_lin;

          off = main_index(base + 1, k);
          grad(off) -= projection(g, k - 1) * sum_quad;
        }
      }

      // ---- PAIRWISE expected ----
      for (int v2 = 0; v2 < num_variables; ++v2) {
        if (v == v2) continue;

        const int row = (v < v2) ? pairwise_effect_indices(v, v2)
                                 : pairwise_effect_indices(v2, v);

        double total = 0.0;
        for (int s = 1; s <= K; ++s)
          total += s * arma::dot(exponents.col(s), obs.col(v2));

        off = pair_index(row, 0);
        grad(off) -= total;

        if (inclusion_indicator(v, v2)) {
          for (int k = 1; k < num_groups; ++k) {
            off = pair_index(row, k);
            grad(off) -= projection(g, k - 1) * total;
          }
        }
      }
    }
  }

  // -------------------------------
  // Priors
  // -------------------------------
  // Main
  for (int v = 0; v < num_variables; ++v) {
    const int base     = main_effect_indices(v, 0);
    const int num_cats = num_categories(v);

    if (is_ordinal_variable(v)) {
      for (int c = 0; c < num_cats; ++c) {
        off = main_index(base + c, 0);
        double value = main_effects(base + c, 0);
        const double p = 1.0 / (1.0 + MY_EXP(-value));
        grad(off) += main_alpha - (main_alpha + main_beta) * p;

        if (inclusion_indicator(v, v) == 0) continue;
        for (int k = 1; k < num_groups; ++k) {
          off = main_index(base + c, k);
          double value = main_effects(base + c, k);
          grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);
        }

      }
    } else {
      off = main_index(base, 0);
      double value = main_effects(base, 0);
      double p = 1.0 / (1.0 + MY_EXP(-value));
      grad(off) += main_alpha - (main_alpha + main_beta) * p;

      off = main_index(base + 1, 0);
      value = main_effects(base + 1, 0);
      p = 1.0 / (1.0 + MY_EXP(-value));
      grad(off) += main_alpha - (main_alpha + main_beta) * p;


      if (inclusion_indicator(v, v) == 0) continue;
      for (int k = 1; k < num_groups; ++k) {
        off = main_index(base, k);
        double value = main_effects(base, k);
        grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);

        off = main_index(base + 1, k);
        value = main_effects(base + 1, k);
        grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);
      }
    }
  }

  // Pairwise (Cauchy)
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      const int row = pairwise_effect_indices(v1, v2);

      // overall uses interaction_scale
      off = pair_index(row, 0);
      double value = pairwise_effects(row, 0);
      grad(off) -= 2.0 * value / (value * value + interaction_scale * interaction_scale);


      if (inclusion_indicator(v1, v2) == 0) continue;
      for (int k = 1; k < num_groups; ++k) {
        off = pair_index(row, k);
        double value = pairwise_effects(row, k);
        grad(off) -= 2.0 * value / (value * value + difference_scale * difference_scale);
      }
    }
  }

  return grad;
}



/**
 * Computes the log pseudoposterior contribution of a single main-effect parameter (bgmCompare model).
 *
 * This function isolates the contribution of one main-effect parameter,
 * either the overall (baseline) effect or one of its group-specific differences.
 *
 * Procedure:
 *  - For each group:
 *    * Construct group-specific main effects for the selected variable
 *      with `compute_group_main_effects()`.
 *    * Construct group-specific pairwise effects for the variable.
 *    * Add linear contributions from sufficient statistics.
 *    * Subtract log normalizing constants from the group-specific likelihood.
 *  - Add prior contribution:
 *    * Logistic–Beta prior for baseline (h == 0).
 *    * Cauchy prior for group differences (h > 0), if included.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
 *  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - observations: Observation matrix (persons × variables).
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category_group: Per-group category counts (for ordinal variables).
 *  - blume_capel_stats_group: Per-group sufficient statistics (for Blume–Capel variables).
 *  - num_groups: Number of groups.
 *  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - difference_scale: Scale parameter for Cauchy priors on group differences.
 *  - variable: Index of the variable of interest.
 *  - category: Category index (only used if variable is ordinal).
 *  - par: Parameter index (0 = linear, 1 = quadratic; used for Blume–Capel).
 *  - h: Column index (0 = overall baseline, >0 = group difference).
 *
 * Returns:
 *  - The scalar log pseudoposterior contribution of the selected parameter.
 *
 * Notes:
 *  - If h > 0 but inclusion_indicator(variable, variable) == 0,
 *    the function returns 0.0 (no contribution).
 *  - This component function is used in parameter-wise Metropolis updates.
 *  - Consistent with the full `log_pseudoposterior()` for bgmCompare.
 */
double log_pseudoposterior_main_component(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double main_alpha,
    const double main_beta,
    const double difference_scale,
    int variable,
    int category, // for ordinal variables only
    int par, // for Blume-Capel variables only
    int h // Overall = 0, differences are 1, ....
) {
  if(h > 0 && inclusion_indicator(variable, variable) == 0) {
    return 0.0; // No contribution if differences not included
  }

  const int num_variables = observations.n_cols;
  const int max_num_categories = num_categories.max();
  double log_pp = 0.0;

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::imat counts_per_category = counts_per_category_group[group];
    const arma::imat blume_capel_stats = blume_capel_stats_group[group];

    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);

    const arma::vec proj_g = projection.row(group).t(); // length = num_groups-1

    // ---- build group-specific main & pairwise effects ----
    arma::vec me = compute_group_main_effects(
      variable, num_groups, main_effects, main_effect_indices, proj_g
    );

    // store into row v (padded with zeros if variable has < max_num_categories params)
    main_group(variable, arma::span(0, me.n_elem - 1)) = me.t();

    // upper triangle incl. base value; mirror to keep symmetry
    for (int u = 0; u < num_variables; u++) {
      if(u == variable) continue;
      double w = compute_group_pairwise_effects(
        variable, u, num_groups, pairwise_effects, pairwise_effect_indices,
        inclusion_indicator, proj_g
      );
      pairwise_group(variable, u) = w;
      pairwise_group(u, variable) = w;
    }

    // ---- data contribution pseudolikelihood (linear terms) ----
    if (is_ordinal_variable(variable)) {
      const double val = main_group(variable, category);
      log_pp += static_cast<double>(counts_per_category(category, variable)) *
        val;
    } else {
      log_pp += static_cast<double>(blume_capel_stats(par, variable)) *
        main_group(variable, par);
    }

    // ---- data contribution pseudolikelihood (quadratic terms) ----
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));

    // ---- pseudolikelihood normalizing constants (per variable) ----
    const arma::vec rest_score = obs * pairwise_group.col(variable);
    const int num_cats = num_categories(variable);

    // bound to stabilize exp; use group-specific params consistently
    arma::vec bound = num_cats * rest_score;
    bound = arma::clamp(bound, 0.0, arma::datum::inf);

    arma::vec denom(rest_score.n_elem, arma::fill::zeros);
    if (is_ordinal_variable(variable)) {
      // base term exp(-bound)
      denom = ARMA_MY_EXP(-bound);
      // main_effects from main_group
      for (int cat = 0; cat < num_cats; cat++) {
        const double th = main_group(variable, cat);
        const arma::vec exponent = th + (cat + 1) * rest_score - bound;
        denom += ARMA_MY_EXP(exponent);
      }
    } else {
      // linear/quadratic main effects from main_group
      const double lin_effect  = main_group(variable, 0);
      const double quad_effect = main_group(variable, 1);
      const int ref = baseline_category(variable);
      for (int cat = 0; cat <= num_cats; cat++) {
        const int centered = cat - ref;
        const double quad = quad_effect * centered * centered;
        const double lin  = lin_effect * cat;
        const arma::vec exponent = lin + quad + cat * rest_score - bound;
        denom += ARMA_MY_EXP(exponent);
      }
    }
    // - sum_i [ bound_i + log denom_i ]
    log_pp -= arma::accu(bound + ARMA_MY_LOG(denom));
  }

  // ---- priors ----
  if (h == 0) {
    // overall
    auto log_beta_prior = [&](double x) {
      return x * main_alpha - std::log1p(MY_EXP(x)) * (main_alpha + main_beta);
    };

    // Main effects prior
    if(is_ordinal_variable(variable)) {
      int r = main_effect_indices(variable, 0) + category;
      log_pp += log_beta_prior(main_effects(r, 0));
    } else {
      int r = main_effect_indices(variable, 0) + par;
      log_pp += log_beta_prior(main_effects(r, 0));
    }
  } else {
    if(is_ordinal_variable(variable)) {
      int r = main_effect_indices(variable, 0) + category;
      log_pp += R::dcauchy(main_effects(r, h), 0.0, difference_scale, true);
    } else {
      int r = main_effect_indices(variable, 0) + par;
      log_pp += R::dcauchy(main_effects(r, h), 0.0, difference_scale, true);
    }
  }

  return log_pp;
}



/**
 * Computes the log pseudoposterior contribution of a single pairwise-effect parameter (bgmCompare model).
 *
 * This function isolates the contribution of one pairwise-effect parameter
 * between two variables, either the baseline effect (h == 0) or a group-specific
 * difference (h > 0).
 *
 * Procedure:
 *  - For each group:
 *    * Construct group-specific main effects for the two variables.
 *    * Construct group-specific pairwise effects for all pairs.
 *    * Add linear contributions from the pairwise sufficient statistic.
 *      - Baseline (h == 0): contribution = 2 * suff_pair * effect.
 *      - Difference (h > 0): scaled by projection value proj_g(h-1).
 *    * Subtract log normalizing constants from both variables’ likelihoods.
 *  - Add prior contribution:
 *    * Cauchy prior for baseline (scale = interaction_scale).
 *    * Cauchy prior for group differences (scale = difference_scale).
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (rows = categories, cols = groups).
 *  - pairwise_effects: Matrix of pairwise-effect parameters (rows = pairs, cols = groups).
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - observations: Observation matrix (persons × variables).
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - num_categories: Number of categories per variable.
 *  - pairwise_stats_group: Per-group pairwise sufficient statistics.
 *  - num_groups: Number of groups.
 *  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - interaction_scale: Scale parameter for Cauchy prior on baseline pairwise effects.
 *  - difference_scale: Scale parameter for Cauchy prior on group differences.
 *  - variable1, variable2: Indices of the variable pair.
 *  - h: Column index (0 = baseline, >0 = group difference).
 *
 * Returns:
 *  - The scalar log pseudoposterior contribution of the selected pairwise effect.
 *
 * Notes:
 *  - If h > 0 but inclusion_indicator(variable1, variable2) == 0,
 *    the function returns 0.0 (no contribution).
 *  - Symmetry is enforced: pairwise effects are mirrored across (var1,var2).
 *  - This component is used in parameter-wise Metropolis updates for pairwise terms.
 */
double log_pseudoposterior_pair_component(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::mat>&  pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double interaction_scale,
    const double difference_scale,
    int variable1,
    int variable2,
    int h // Overall = 0, differences are 1, ....
) {
  if(h > 0 && inclusion_indicator(variable1, variable2) == 0) {
    return 0.0; // No contribution if differences not included
  }

  const int num_variables = observations.n_cols;
  const int max_num_categories = num_categories.max();
  double log_pp = 0.0;
  int idx = pairwise_effect_indices(variable1, variable2);

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    arma::mat main_group(num_variables, max_num_categories, arma::fill::zeros);
    arma::mat pairwise_group(num_variables, num_variables, arma::fill::zeros);

    const arma::vec proj_g = projection.row(group).t(); // length = num_groups-1

    // ---- build group-specific main & pairwise effects ----
    for (int v : {variable1, variable2}) { // Only populate two rows
      arma::vec me = compute_group_main_effects(
        v, num_groups, main_effects, main_effect_indices, proj_g
      );
      main_group(v, arma::span(0, me.n_elem - 1)) = me.t();
    }
    for (int v = 0; v < num_variables; v++) {
      for (int u = v; u < num_variables; ++u) { // Combines with loop over v
        if(u == v) continue;
        double w = compute_group_pairwise_effects(
          v, u, num_groups, pairwise_effects, pairwise_effect_indices,
          inclusion_indicator, proj_g
        );
        pairwise_group(v, u) = w;
        pairwise_group(u, v) = w;
      }
    }

    // ---- data contribution pseudolikelihood ----
    const arma::mat pairwise_stats = pairwise_stats_group[group];
    const double suff_pair = pairwise_stats(variable1, variable2);

    if(h == 0) {
      log_pp += 2.0 * suff_pair * pairwise_effects(idx, h);
    } else {
      log_pp += 2.0 * suff_pair * proj_g(h-1) * pairwise_effects(idx, h);
    }

    // ---- pseudolikelihood normalizing constants (per variable) ----
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));
    const arma::mat residual_matrix = obs * pairwise_group;

    for (int v : {variable1, variable2}) {
      const int num_cats = num_categories(v);
      const arma::vec rest_score = residual_matrix.col(v);

      // bound to stabilize exp; use group-specific params consistently
      arma::vec bound = num_cats * rest_score;
      bound = arma::clamp(bound, 0.0, arma::datum::inf);

      arma::vec denom(rest_score.n_elem, arma::fill::zeros);

      if (is_ordinal_variable(v)) {
        // base term exp(-bound)
        denom = ARMA_MY_EXP(-bound);
        // main_effects from main_group
        for (int c = 0; c < num_cats; ++c) {
          const double th = main_group(v, c);
          const arma::vec exponent = th + (c + 1) * rest_score - bound;
          denom += ARMA_MY_EXP(exponent);
        }
      } else {
        // linear/quadratic main effects from main_group
        const double lin_effect  = main_group(v, 0);
        const double quad_effect = main_group(v, 1);
        const int ref = baseline_category(v);
        for (int c = 0; c <= num_cats; ++c) {
          const int centered = c - ref;
          const double quad = quad_effect * centered * centered;
          const double lin  = lin_effect * c;
          const arma::vec exponent = lin + quad + c * rest_score - bound;
          denom += ARMA_MY_EXP(exponent);
        }
      }
      // - sum_i [ bound_i + log denom_i ]
      log_pp -= arma::accu(bound + ARMA_MY_LOG(denom));
    }
  }

  // ---- priors ----
  if (h == 0) {
    log_pp += R::dcauchy(pairwise_effects(idx, 0), 0.0, interaction_scale, true);
  } else {
    log_pp += R::dcauchy(pairwise_effects(idx, h), 0.0, difference_scale, true);
  }
  return log_pp;
}



/**
 * Computes the log-ratio of pseudolikelihood normalizing constants
 * for a single variable under current vs. proposed parameters (bgmCompare model).
 *
 * This function is used in Metropolis–Hastings updates for main-effect parameters.
 * It evaluates how the normalizing constant (denominator of the pseudolikelihood)
 * changes when switching from the current to the proposed parameter values.
 *
 * Procedure:
 *  - For each group:
 *    * Construct group-specific main effects (current vs. proposed).
 *    * Construct group-specific pairwise weights for the variable.
 *    * Compute residual scores for observations under both models.
 *    * Calculate denominators with stability bounds (ordinal vs. Blume–Capel cases).
 *    * Accumulate the log-ratio contribution across all observations.
 *
 * Inputs:
 *  - current_main_effects, proposed_main_effects: Matrices of main-effect parameters
 *    (rows = categories, cols = groups).
 *  - current_pairwise_effects, proposed_pairwise_effects: Matrices of pairwise-effect parameters
 *    (rows = pairs, cols = groups).
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - observations: Observation matrix (persons × variables).
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - num_categories: Number of categories per variable.
 *  - num_groups: Number of groups.
 *  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - variable: Index of the variable being updated.
 *
 * Returns:
 *  - The scalar log-ratio of pseudolikelihood constants
 *    (current model vs. proposed model).
 *
 * Notes:
 *  - For ordinal variables, denominators include exp(-bound) and category terms.
 *  - For Blume–Capel variables, denominators use linear/quadratic scores
 *    with baseline centering.
 *  - Stability bounds (`bound_current`, `bound_proposed`) are applied to avoid overflow.
 */
double log_ratio_pseudolikelihood_constant_variable(
    const arma::mat& current_main_effects,
    const arma::mat& current_pairwise_effects,
    const arma::mat& proposed_main_effects,
    const arma::mat& proposed_pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int variable
) {
  const int num_cats = num_categories(variable);
  const int num_variables = observations.n_cols;

  double log_ratio = 0.0;

  // --- per group ---
  for (int group = 0; group < num_groups; ++group) {
    const arma::vec proj_g = projection.row(group).t();

    // --- group-specific main effects (current/proposed) ---
    const arma::vec main_current = compute_group_main_effects(
      variable, num_groups, current_main_effects, main_effect_indices, proj_g
    );
    const arma::vec main_proposed = compute_group_main_effects(
      variable, num_groups, proposed_main_effects, main_effect_indices, proj_g
    );

    // --- group-specific pairwise effects for this variable (column) ---
    arma::vec weights_current(num_variables, arma::fill::zeros);
    arma::vec weights_proposed(num_variables, arma::fill::zeros);
    for (int u = 0; u < num_variables; ++u) {
      if (u == variable) continue;
      weights_current(u) = compute_group_pairwise_effects(
        variable, u, num_groups, current_pairwise_effects,
        pairwise_effect_indices, inclusion_indicator, proj_g
      );
      weights_proposed(u) = compute_group_pairwise_effects(
        variable, u, num_groups, proposed_pairwise_effects,
        pairwise_effect_indices, inclusion_indicator, proj_g
      );
    }

    // --- group observations and rest scores ---
    const int r0 = group_indices(group, 0);
    const int r1 = group_indices(group, 1);
    const arma::mat obs = arma::conv_to<arma::mat>::from(observations.rows(r0, r1));

    const arma::vec rest_current = obs * weights_current;
    const arma::vec rest_proposed = obs * weights_proposed;

    // --- denominators with stability bounds ---
    arma::vec bound_current;
    arma::vec bound_proposed;
    arma::vec denom_current(rest_current.n_elem, arma::fill::zeros);
    arma::vec denom_proposed(rest_proposed.n_elem, arma::fill::zeros);

    if (is_ordinal_variable(variable)) {
      // regular ordinal/binary
      bound_current = num_cats * arma::clamp(rest_current, 0.0, arma::datum::inf);
      bound_proposed = num_cats * arma::clamp(rest_proposed, 0.0, arma::datum::inf);

      denom_current = ARMA_MY_EXP(-bound_current);
      denom_proposed = ARMA_MY_EXP(-bound_proposed);

      for (int c = 0; c < num_cats; ++c) {
        denom_current += ARMA_MY_EXP(main_current(c) + (c + 1) * rest_current - bound_current);
        denom_proposed += ARMA_MY_EXP(main_proposed(c) + (c + 1) * rest_proposed - bound_proposed);
      }
    } else {
      // Blume-Capel: linear + quadratic
      const int ref = baseline_category(variable);

      arma::vec const_current(num_cats + 1, arma::fill::zeros);
      arma::vec const_proposed(num_cats + 1, arma::fill::zeros);
      for (int s = 0; s <= num_cats; ++s) {
        const int centered = s - ref;
        const_current(s) = main_current(0) * s + main_current(1) * centered * centered;
        const_proposed(s) = main_proposed(0) * s + main_proposed(1) * centered * centered;
      }

      double lbound = std::max(const_current.max(), const_proposed.max());
      if (lbound < 0.0) lbound = 0.0;

      bound_current = lbound + num_cats * arma::clamp(rest_current, 0.0, arma::datum::inf);
      bound_proposed = lbound + num_cats * arma::clamp(rest_proposed, 0.0, arma::datum::inf);

      for (int s = 0; s <= num_cats; ++s) {
        denom_current += ARMA_MY_EXP(const_current(s) + s * rest_current - bound_current);
        denom_proposed += ARMA_MY_EXP(const_proposed(s) + s * rest_proposed - bound_proposed);
      }
    }

    // --- accumulate contribution ---
    log_ratio += arma::accu((bound_current - bound_proposed) +
      ARMA_MY_LOG(denom_current) - ARMA_MY_LOG(denom_proposed));
  }

  return log_ratio;
}



/**
 * Computes the log pseudolikelihood ratio for updating a single main-effect parameter (bgmCompare model).
 *
 * This function is used in Metropolis–Hastings updates for main effects.
 * It compares the likelihood of the data under the current vs. proposed
 * value of a single variable’s main-effect parameter, while keeping
 * all other parameters fixed.
 *
 * Procedure:
 *  - For each group:
 *    * Compute group-specific main effects for the variable (current vs. proposed).
 *    * Add contributions from observed sufficient statistics
 *      (category counts or Blume–Capel stats).
 *  - Add the ratio of pseudolikelihood normalizing constants by calling
 *    `log_ratio_pseudolikelihood_constant_variable()`.
 *
 * Inputs:
 *  - current_main_effects: Matrix of main-effect parameters (current state).
 *  - proposed_main_effects: Matrix of main-effect parameters (candidate state).
 *  - current_pairwise_effects: Matrix of pairwise-effect parameters (fixed at current state).
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - observations: Observation matrix (persons × variables).
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category_group: Per-group category counts (for ordinal variables).
 *  - blume_capel_stats_group: Per-group sufficient statistics (for Blume–Capel variables).
 *  - num_groups: Number of groups.
 *  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - variable: Index of the variable being updated.
 *
 * Returns:
 *  - The scalar log pseudolikelihood ratio (proposed vs. current).
 *
 * Notes:
 *  - A temporary copy of `inclusion_indicator` is made to ensure the
 *    variable’s self-term (diagonal entry) is included.
 *  - Only the variable under update changes between current and proposed states;
 *    all other variables and pairwise effects remain fixed.
 *  - This function does not add prior contributions — only pseudolikelihood terms.
 */
double log_pseudolikelihood_ratio_main(
    const arma::mat& current_main_effects,
    const arma::mat& proposed_main_effects,
    const arma::mat& current_pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat&  projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::imat>& counts_per_category_group,
    const std::vector<arma::imat>& blume_capel_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int variable
) {
  double lr = 0.0;
  arma::imat tmp_ind = inclusion_indicator;
  tmp_ind(variable, variable) = 1; // Ensure self-interaction is included

  // Add data contribution (group-specific parameters via projection)
  for (int g = 0; g < num_groups; ++g) {
    const arma::vec proj_g = projection.row(g).t();

    const arma::vec main_cur  = compute_group_main_effects(
      variable, num_groups, current_main_effects,  main_effect_indices, proj_g
    );
    const arma::vec main_prop = compute_group_main_effects(
      variable, num_groups, proposed_main_effects, main_effect_indices, proj_g
    );

    if (is_ordinal_variable(variable)) {
      const arma::imat& num_obs = counts_per_category_group[g];
      const int num_cats = num_categories(variable);
      for (int c = 0; c < num_cats; ++c) {
        lr += (main_prop(c) - main_cur(c)) * static_cast<double>(num_obs(c, variable));
      }
    } else {
      const arma::imat& suff = blume_capel_stats_group[g];
      lr += (main_prop(0) - main_cur(0)) * static_cast<double>(suff(0, variable));
      lr += (main_prop(1) - main_cur(1)) * static_cast<double>(suff(1, variable));
    }
  }

  // Add ratio of normalizing constants
  lr += log_ratio_pseudolikelihood_constant_variable(
    current_main_effects, current_pairwise_effects, proposed_main_effects,
    /* same */ current_pairwise_effects, main_effect_indices,
    pairwise_effect_indices, projection, observations, group_indices,
    num_categories, num_groups, tmp_ind, is_ordinal_variable,
    baseline_category, variable
  );

  return lr;
}


/**
 * Computes the log pseudolikelihood ratio for updating a single pairwise-effect parameter (bgmCompare model).
 *
 * This function is used in Metropolis–Hastings updates for pairwise effects.
 * It compares the likelihood of the data under the current vs. proposed
 * value of a single interaction (var1,var2), while keeping all other
 * parameters fixed.
 *
 * Procedure:
 *  - Ensure the interaction is included in a temporary copy of inclusion_indicator.
 *  - For each group:
 *    * Compute group-specific pairwise effect for (var1,var2), current vs. proposed.
 *    * Add linear contribution from the pairwise sufficient statistic.
 *  - Add the ratio of pseudolikelihood normalizing constants for both variables:
 *    * Call `log_ratio_pseudolikelihood_constant_variable()` separately for var1 and var2,
 *      comparing current vs. proposed pairwise weights.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (fixed).
 *  - current_pairwise_effects: Matrix of pairwise-effect parameters (current state).
 *  - proposed_pairwise_effects: Matrix of pairwise-effect parameters (candidate state).
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - observations: Observation matrix (persons × variables).
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - num_categories: Number of categories per variable.
 *  - pairwise_stats_group: Per-group pairwise sufficient statistics.
 *  - num_groups: Number of groups.
 *  - inclusion_indicator: Symmetric binary matrix of active variables (diag) and pairs (off-diag).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - var1, var2: Indices of the variable pair being updated.
 *
 * Returns:
 *  - The scalar log pseudolikelihood ratio (proposed vs. current).
 *
 * Notes:
 *  - A temporary copy of `inclusion_indicator` is used to force the edge (var1,var2) as active.
 *  - Only the selected pair changes between current and proposed states;
 *    all other effects remain fixed.
 *  - This function does not add prior contributions — only pseudolikelihood terms.
 */
double log_pseudolikelihood_ratio_pairwise(
    const arma::mat& main_effects,
    const arma::mat& current_pairwise_effects,
    const arma::mat& proposed_pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::imat& observations,
    const arma::imat& group_indices,
    const arma::ivec& num_categories,
    const std::vector<arma::mat>& pairwise_stats_group,
    const int num_groups,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int var1,
    const int var2
) {
  double lr = 0.0;
  // Ensure interaction is included
  arma::imat tmp_ind = inclusion_indicator;
  tmp_ind(var1, var2) = 1;
  tmp_ind(var2, var1) = 1;

  // Add data contribution
  for (int g = 0; g < num_groups; ++g) {
    const arma::vec proj_g = projection.row(g).t();
    const arma::mat& suff  = pairwise_stats_group[g];

    const double w_cur = compute_group_pairwise_effects(
      var1, var2, num_groups, current_pairwise_effects,
      pairwise_effect_indices, tmp_ind, proj_g
    );
    const double w_prop = compute_group_pairwise_effects(
      var1, var2, num_groups, proposed_pairwise_effects,
      pairwise_effect_indices, tmp_ind, proj_g
    );

    lr += 2.0 * (w_prop - w_cur) * suff(var1, var2);
  }

  // Add ratio of normalizing constant for `var1`
  lr += log_ratio_pseudolikelihood_constant_variable(
    main_effects, current_pairwise_effects, /* same */ main_effects,
    proposed_pairwise_effects, main_effect_indices, pairwise_effect_indices,
    projection, observations, group_indices, num_categories, num_groups,
    tmp_ind, is_ordinal_variable, baseline_category, var1
  );

  // Add ratio of normalizing constant for `var2`
  lr += log_ratio_pseudolikelihood_constant_variable(
    main_effects, current_pairwise_effects, /* same */ main_effects,
    proposed_pairwise_effects, main_effect_indices, pairwise_effect_indices,
    projection, observations, group_indices, num_categories, num_groups,
    tmp_ind, is_ordinal_variable, baseline_category, var2
  );

  return lr;
}