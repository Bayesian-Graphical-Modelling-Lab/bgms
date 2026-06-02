#pragma once

/**
 * @file bgmCompare_helper.h
 * @brief Utility functions for the bgmCompare multi-group comparison model.
 *
 * Functions for computing group-specific effects via contrast projection,
 * vectorizing/unvectorizing the parameter space (restricting group-difference
 * entries to active edges), building index maps, and initializing the graph.
 */

#include <RcppArmadillo.h>
#include "rng/rng_utils.h"


/**
 * Compute the group-specific main-effect vector for a single variable.
 *
 * Combines the baseline column with group-deviation columns weighted
 * by the contrast vector for one group.
 *
 * @param variable            Variable index
 * @param num_groups           Number of groups (G)
 * @param main_effects         Main-effect parameter matrix (n_main_rows x G)
 * @param main_effect_indices  Start/end row indices per variable (V x 2)
 * @param proj_group           Contrast vector for the target group (length G-1)
 * @return Group-specific main-effect vector (length = categories for variable)
 */
arma::vec compute_group_main_effects(
    const int variable,
    const int num_groups,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::vec& proj_group
);

/**
 * Compute the group-specific pairwise interaction for a variable pair.
 *
 * Adds the baseline effect and (if the edge is active) the group deviation
 * projected through the contrast vector. Returns zero for inactive edges.
 *
 * @param var1                     First variable index
 * @param var2                     Second variable index
 * @param num_groups               Number of groups (G)
 * @param pairwise_effects         Pairwise-effect parameter matrix (n_pair_rows x G)
 * @param pairwise_effect_indices  Row index per variable pair (V x V)
 * @param inclusion_indicator      Edge inclusion indicators (V x V)
 * @param proj_group               Contrast vector for the target group (length G-1)
 * @return Group-specific pairwise interaction (scalar)
 */
double compute_group_pairwise_effects(
    const int var1,
    const int var2,
    const int num_groups,
    const arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::vec& proj_group
);

/**
 * Flatten the main-effect and pairwise-effect matrices into a single vector.
 *
 * Layout: overall mains, overall pairs, active main differences,
 * active pair differences. Group-difference entries are included only
 * for variables/pairs whose inclusion indicator is 1.
 *
 * @param main_effects             Main-effect matrix (n_main_rows x G)
 * @param pairwise_effects         Pairwise-effect matrix (n_pair_rows x G)
 * @param inclusion_indicator      Edge inclusion indicators (V x V)
 * @param main_effect_indices      Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices  Row index per variable pair (V x V)
 * @param num_categories           Number of categories per variable
 * @param is_ordinal_variable      1 = ordinal, 0 = Blume-Capel
 * @return Flat parameter vector
 */
arma::vec vectorize_model_parameters_bgmcompare(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

/**
 * Reconstruct parameter matrices from a flat parameter vector.
 *
 * Inverse of vectorize_model_parameters_bgmcompare(). Writes group-difference
 * columns only for active variables/pairs; inactive columns remain zero.
 *
 * @param param_vec                Flat parameter vector
 * @param[out] main_effects_out    Reconstructed main-effect matrix (n_main_rows x G)
 * @param[out] pairwise_effects_out Reconstructed pairwise-effect matrix (n_pair_rows x G)
 * @param inclusion_indicator      Edge inclusion indicators (V x V)
 * @param main_effect_indices      Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices  Row index per variable pair (V x V)
 * @param num_groups               Number of groups (G)
 * @param num_categories           Number of categories per variable
 * @param is_ordinal_variable      1 = ordinal, 0 = Blume-Capel
 */
void unvectorize_model_parameters_bgmcompare(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,
    arma::mat& pairwise_effects_out,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const int num_groups,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

/**
 * Build index maps from parameter-matrix entries to flat-vector positions.
 *
 * Returns two integer matrices mapping each (row, group-column) entry of the
 * main-effect and pairwise-effect matrices to its index in the flat parameter
 * vector, or -1 if the entry is inactive.
 *
 * @param main_effects             Main-effect matrix (n_main_rows x G)
 * @param pairwise_effects         Pairwise-effect matrix (n_pair_rows x G)
 * @param inclusion_indicator      Edge inclusion indicators (V x V)
 * @param main_effect_indices      Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices  Row index per variable pair (V x V)
 * @param num_categories           Number of categories per variable
 * @param is_ordinal_variable      1 = ordinal, 0 = Blume-Capel
 * @return Pair of (main_index, pair_index) matrices
 */
std::pair<arma::imat, arma::imat> build_index_maps(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

/**
 * Extract inverse-mass diagonal entries for active parameters only.
 *
 * When selection is enabled, restricts to entries whose inclusion indicator
 * is 1. When selection is disabled, returns the full diagonal.
 *
 * @param inv_diag                 Full inverse-mass diagonal
 * @param inclusion_indicator      Edge inclusion indicators (V x V)
 * @param num_groups               Number of groups (G)
 * @param num_categories           Number of categories per variable
 * @param is_ordinal_variable      1 = ordinal, 0 = Blume-Capel
 * @param main_index               Main-effect index map from build_index_maps()
 * @param pair_index               Pairwise index map from build_index_maps()
 * @param main_effect_indices      Start/end row indices per variable (V x 2)
 * @param pairwise_effect_indices  Row index per variable pair (V x V)
 * @param selection                true restricts to active entries
 * @return Inverse-mass vector for active parameters
 */
arma::vec inv_mass_active(
    const arma::vec& inv_diag,
    const arma::imat& inclusion_indicator,
    const int num_groups,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const bool& selection
);

/**
 * Initialize the graph structure for a bgmCompare model.
 *
 * Draws each pairwise edge on/off from its prior inclusion probability.
 * For excluded edges, zeros out the group-difference columns. Optionally
 * draws main-effect inclusion indicators when main_difference_selection
 * is enabled.
 *
 * @param[in,out] indicator                Edge/main inclusion indicators (V x V)
 * @param[in,out] main                     Main-effect matrix (modified in place)
 * @param[in,out] pairwise                 Pairwise-effect matrix (modified in place)
 * @param main_indices                     Start/end row indices per variable (V x 2)
 * @param pairwise_indices                 Row index per variable pair (V x V)
 * @param incl_prob                        Prior inclusion probabilities (V x V)
 * @param main_difference_selection        Enable main-effect difference selection
 * @param rng                              Random number generator
 */

inline void initialise_graph_bgmcompare(
    arma::imat& indicator,
    arma::mat& main,
    arma::mat& pairwise,
    const arma::imat& main_indices,
    const arma::imat& pairwise_indices,
    const arma::mat& incl_prob,
    const bool main_difference_selection,
    SafeRNG& rng
) {
  int V = indicator.n_rows;
  int G = main.n_cols;
  // Initialize pairwise indicators
  for (int i = 0; i < V-1; ++i) {
    for (int j = i+1; j < V; ++j) {
      double p = incl_prob(i,j);
      int draw = (runif(rng) < p) ? 1 : 0;
      indicator(i,j) = indicator(j,i) = draw;
      if (!draw) {
        int row = pairwise_indices(i, j);
        pairwise(row, arma::span(1, G-1)).zeros();
      }
    }
  }
  // Initialize main effect indicators (only if main_difference_selection is enabled)
  for(int i = 0; i < V; i++) {
    if (main_difference_selection) {
      double p = incl_prob(i,i);
      int draw = (runif(rng) < p) ? 1 : 0;
      indicator(i,i) = draw;
      if(!draw) {
        int start = main_indices(i,0);
        int end = main_indices(i,1);
        main(arma::span(start, end), arma::span(1, G - 1)).zeros();
      }
    } else {
      // Keep main effect indicators at 1 (all differences included)
      indicator(i,i) = 1;
    }
  }
};