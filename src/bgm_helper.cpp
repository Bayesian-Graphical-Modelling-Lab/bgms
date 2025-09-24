#include <RcppArmadillo.h>
#include "bgm_helper.h"
#include "common_helpers.h"

using namespace Rcpp;



/**
 * Flattens the main-effect parameter matrix into a vector (bgm model).
 *
 * For each variable:
 *  - Ordinal variables contribute one parameter per category.
 *  - Blume–Capel variables contribute two parameters (linear and quadratic).
 *
 * The parameters are written into the output vector in variable order.
 * This function is the inverse of `unvectorize_main_effects_bgm()`.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - num_categories: Number of categories per variable.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *
 * Returns:
 *  - A vector containing all main-effect parameters in row-major order.
 */
arma::vec vectorize_main_effects_bgm (
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_parameters = count_num_main_effects (num_categories, is_ordinal_variable);
  arma::vec main_effect_vector (num_parameters);
  int offset = 0;

  for (int variable = 0; variable < static_cast<int>(main_effects.n_rows); variable++) {
    const int num_pars = is_ordinal_variable (variable) ? num_categories(variable) : 2;
    main_effect_vector.subvec (offset, offset + num_pars - 1) =
      main_effects.row (variable).cols (0, num_pars - 1).t ();
    offset += num_pars;
  }

  return main_effect_vector;
}



/**
 * Reconstructs a main-effect parameter matrix from a flat vector (bgm model).
 *
 * For each variable:
 *  - Ordinal variables read one parameter per category.
 *  - Blume–Capel variables read two parameters (linear and quadratic).
 *
 * The parameters are written into the output matrix in variable order.
 * This function is the inverse of `vectorize_main_effects_bgm()`.
 *
 * Inputs:
 *  - main_effect_vector: Flat vector of main-effect parameters.
 *  - num_categories: Number of categories per variable.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *
 * Returns:
 *  - A matrix of main-effect parameters (variables × max_categories).
 *    Unused entries are filled with zeros.
 */
arma::mat unvectorize_main_effects_bgm (
    const arma::vec& main_effect_vector,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = num_categories.n_elem;
  const int max_categories = num_categories.max ();

  arma::mat matrix (num_variables, max_categories, arma::fill::zeros);

  int offset = 0;
  for (int variable = 0; variable < num_variables; variable++) {
    const int num_pars = is_ordinal_variable[variable] ? num_categories[variable] : 2;
    matrix.row (variable).cols (0, num_pars - 1) =
      main_effect_vector.subvec (offset, offset + num_pars - 1).t ();
    offset += num_pars;
  }

  return matrix;
}



/**
 * Flattens all main-effect and pairwise-effect parameters into a single vector (bgm model).
 *
 * The vector is constructed in two parts:
 *  1. Main effects: each variable contributes either
 *     - one parameter per category (if ordinal), or
 *     - two parameters (linear and quadratic, if Blume–Capel).
 *  2. Pairwise effects: upper-triangle interactions (v1 < v2) are added
 *     only if the pair is active in the inclusion indicator.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - pairwise_effects: Symmetric matrix of interaction strengths.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - num_categories: Number of categories per variable.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *
 * Returns:
 *  - A flat vector containing all main-effect and active pairwise-effect parameters.
 *
 * Note: This function is the inverse of `unvectorize_model_parameters_bgm()`.
 */
arma::vec vectorize_model_parameters_bgm(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = main_effects.n_rows;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);

  // Compute total number of active interactions
  int num_active = 0;
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) {
        num_active++;
      }
    }
  }

  arma::vec param_vec(num_main + num_active, arma::fill::zeros);
  int offset = 0;

  // main_effect parameters
  for (int v = 0; v < num_variables; ++v) {
    if (is_ordinal_variable(v)) {
      int num_cats = num_categories(v);
      for (int c = 0; c < num_cats; ++c) {
        param_vec(offset++) = main_effects(v, c);
      }
    } else {
      // Blume-Capel: always two params
      param_vec(offset++) = main_effects(v, 0);  // linear
      param_vec(offset++) = main_effects(v, 1);  // quadratic
    }
  }

  // Interaction parameters (upper triangle)
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) {
        param_vec(offset++) = pairwise_effects(v1, v2);
      }
    }
  }

  return param_vec;
}



/**
 * Reconstructs main-effect and pairwise-effect matrices from a parameter vector (bgm model).
 *
 * The input vector must have the layout created by `vectorize_model_parameters_bgm()`:
 *  - Main effects first, in variable order
 *    (one per category if ordinal, or two if Blume–Capel).
 *  - Then active pairwise effects, in upper-triangle order (v1 < v2).
 *
 * Inputs:
 *  - param_vec: Flattened parameter vector.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - num_categories: Number of categories per variable.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *
 * Outputs:
 *  - main_effects_out: Matrix of main-effect parameters (variables × max_categories).
 *  - pairwise_effects_out: Symmetric matrix of pairwise effects.
 *    Only active interactions are filled; inactive entries remain zero.
 *
 * Note: This function is the inverse of `vectorize_model_parameters_bgm()`.
 */
void unvectorize_model_parameters_bgm(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,
    arma::mat& pairwise_effects_out,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
) {
  const int num_variables = num_categories.n_elem;
  const int max_num_categories = num_categories.max();

  main_effects_out.set_size(num_variables, max_num_categories);
  pairwise_effects_out.set_size(num_variables, num_variables);
  main_effects_out.fill(0.0);
  pairwise_effects_out.fill(0.0);

  int offset = 0;

  // --- Reconstruct main_effects
  for (int v = 0; v < num_variables; ++v) {
    if (is_ordinal_variable(v)) {
      int num_cats = num_categories(v);
      for (int c = 0; c < num_cats; ++c) {
        main_effects_out(v, c) = param_vec(offset++);
      }
    } else {
      main_effects_out(v, 0) = param_vec(offset++);
      main_effects_out(v, 1) = param_vec(offset++);
    }
  }

  // --- Reconstruct only active interactions
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) {
        double val = param_vec(offset++);
        pairwise_effects_out(v1, v2) = val;
        pairwise_effects_out(v2, v1) = val;
      }
    }
  }
}



/**
 * Extracts the inverse mass matrix entries for active parameters only (bgm model).
 *
 * When selection is disabled, the full diagonal vector is returned unchanged.
 * When selection is enabled, the output vector is restricted to:
 *  - all main-effect parameters, and
 *  - only those pairwise effects marked as active in the inclusion indicator.
 *
 * Inputs:
 *  - inv_diag: Full inverse mass matrix diagonal (main + all pairwise effects).
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - num_categories: Number of categories per variable.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - selection: If false, return inv_diag unchanged; if true, restrict to active parameters.
 *
 * Returns:
 *  - Vector of inverse mass entries for main effects and active pairwise effects.
 *
 * Note: The layout of parameters must match `vectorize_model_parameters_bgm()`.
 */
arma::vec inv_mass_active(
    const arma::vec& inv_diag,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const bool& selection
) {
  if(selection == false)
    return inv_diag;

  const int num_variables = inclusion_indicator.n_rows;
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);

  // Compute total number of active interactions
  int num_active = 0;
  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) {
        num_active++;
      }
    }
  }

  arma::vec active_inv_diag(num_main + num_active, arma::fill::zeros);
  active_inv_diag.head(num_main) = inv_diag.head(num_main);

  int offset_full = num_main;
  int offset_active = num_main;

  for (int v1 = 0; v1 < num_variables - 1; ++v1) {
    for (int v2 = v1 + 1; v2 < num_variables; ++v2) {
      if (inclusion_indicator(v1, v2) == 1) {
        active_inv_diag(offset_active) = inv_diag(offset_full);
        offset_active++;
      }
      offset_full++;
    }
  }
  return active_inv_diag;
}