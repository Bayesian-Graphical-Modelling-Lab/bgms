#pragma once

#include <RcppArmadillo.h>
#include "rng_utils.h"



arma::vec compute_group_main_effects(
    const int variable,
    const int num_groups,
    const arma::mat& main_effects,
    const arma::imat& main_effect_indices,
    const arma::vec& proj_group
);

double compute_group_pairwise_effects(
    const int var1,
    const int var2,
    const int num_groups,
    const arma::mat& pairwise_effects,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::vec& proj_group
);

arma::vec vectorize_model_parameters_bgmcompare(
    const arma::mat& main_effects,                 // [n_main_rows × G]
    const arma::mat& pairwise_effects,             // [n_pair_rows × G]
    const arma::imat& inclusion_indicator,         // [V × V]
    const arma::imat& main_effect_indices,         // [V × 2], inclusive [start,end]
    const arma::imat& pairwise_effect_indices,      // [V × V]
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

void unvectorize_model_parameters_bgmcompare(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,                 // [n_main_rows × G]
    arma::mat& pairwise_effects_out,             // [n_pair_rows × G]
    const arma::imat& inclusion_indicator,       // [V × V]
    const arma::imat& main_effect_indices,       // [V × 2], inclusive [start,end]
    const arma::imat& pairwise_effect_indices,   // [V × V]
    const int num_groups,                          // G
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

std::pair<arma::imat, arma::imat> build_index_maps(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

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

inline void initialise_graph_bgmcompare(
    arma::imat& indicator,
    arma::mat& main,
    arma::mat& pairwise,
    const arma::imat& main_indices,
    const arma::imat& pairwise_indices,
    const arma::mat& incl_prob,
    SafeRNG& rng
) {
  int V = indicator.n_rows;
  int G = main.n_cols;
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
  for(int i = 0; i < V; i++) {
    double p = incl_prob(i,i);
    int draw = (runif(rng) < p) ? 1 : 0;
    indicator(i,i) = draw;
    if(!draw) {
      int start = main_indices(i,0);
      int end = main_indices(i,1);
      main(arma::span(start, end), arma::span(1, G - 1)).zeros();
    }
  }
};
