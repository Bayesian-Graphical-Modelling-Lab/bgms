#pragma once

#include <RcppArmadillo.h>

// Log posterior for a single component of main effect
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
);

// Log posterior for a single component of interactions
double log_pseudoposterior_interactions_component (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    const arma::imat& pairwise_stats,
    const int var1,
    const int var2
);

// Log posterior for a single component of interactions (optimized with residual matrix)
double log_pseudoposterior_interactions_component (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::mat& residual_matrix,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    const arma::imat& pairwise_stats,
    const int var1,
    const int var2,
    const double delta
);

std::pair<arma::vec, arma::imat>  gradient_observed_active(
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const arma::imat& pairwise_stats
);

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
    const arma::mat& pairwise_scaling_factors,
    const arma::mat& residual_matrix,
    const arma::imat index_matrix,
    const arma::vec grad_obs
);

// Joint log-pseudoposterior and gradient (for efficient NUTS/HMC)
std::pair<double, arma::vec> logp_and_gradient(
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
    const arma::mat& pairwise_scaling_factors,
    const arma::imat& pairwise_stats,
    const arma::mat& residual_matrix,
    const arma::imat index_matrix,
    const arma::vec grad_obs
);

// Pseudolikelihood ratio for a single variable
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
);

// Pseudolikelihood ratio for an interaction
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
);