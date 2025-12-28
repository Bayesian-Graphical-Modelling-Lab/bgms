#pragma once

#include <RcppArmadillo.h>
#include "rng_utils.h"

// Vectorize main_effect matrix
arma::vec vectorize_main_effects_bgm(
    const arma::mat& main_effects,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

// Unvectorize main_effects into matrix
arma::mat unvectorize_main_effects_bgm(
    const arma::vec& main_effect_vector,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

// Vectorize model parameters (main + interaction effects)
arma::vec vectorize_model_parameters_bgm(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

// Unvectorize model parameters back into matrices
void unvectorize_model_parameters_bgm(
    const arma::vec& param_vec,
    arma::mat& main_effects_out,
    arma::mat& pairwise_effects_out,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable
);

arma::vec inv_mass_active(
    const arma::vec& inv_diag,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const bool& selection
);

inline void initialise_graph_bgm(
    arma::imat& indicator,
    arma::mat&  pairwise,
    const arma::mat& incl_prob,
    arma::mat&  rest,
    const arma::imat& X,
    SafeRNG& rng
) {
  int V = indicator.n_rows;
  for (int i = 0; i < V-1; ++i) {
    for (int j = i+1; j < V; ++j) {
      double p = incl_prob(i,j);
      int draw = (runif(rng) < p) ? 1 : 0;
      indicator(i,j) = indicator(j,i) = draw;
      if (!draw)
        pairwise(i,j) = pairwise(j,i) = 0.0;
    }
  }
  rest = X * pairwise;
}
