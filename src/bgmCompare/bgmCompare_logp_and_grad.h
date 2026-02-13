#pragma once

#include <RcppArmadillo.h>



arma::uword total_length(
    const int num_variables,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::imat& inclusion_indicator,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const int num_groups
);

arma::vec gradient_observed_active(
  const arma::imat& main_effect_indices,
  const arma::imat& pairwise_effect_indices,
  const arma::mat& projection,
  const arma::imat& observations,
  const arma::imat& group_indices,
  const arma::ivec& num_categories,
  const arma::imat& inclusion_indicator,
  const std::vector<arma::imat>& counts_per_category_group,
  const std::vector<arma::imat>& blume_capel_stats_group,
  const std::vector<arma::mat>&  pairwise_stats_group,
  const int num_groups,
  const arma::uvec& is_ordinal_variable,
  const arma::ivec& baseline_category,
  const arma::imat main_index,
  const arma::imat pair_index
);

double log_pseudoposterior(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::mat& observations_double,
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
    const arma::mat& pairwise_scaling_factors,
    const double difference_scale
);

arma::vec gradient(
    const arma::mat& main_effects,
    const arma::mat& pairwise_effects,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const arma::mat& projection,
    const arma::mat& observations_double,
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
    const arma::mat& pairwise_scaling_factors,
    const double difference_scale,
    const arma::imat& main_index,
    const arma::imat& pair_index,
    const arma::vec& grad_obs
);

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
);

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
    const arma::mat& pairwise_scaling_factors,
    const double difference_scale,
    int variable1,
    int variable2,
    int h // Overall = 0, differences are 1, ....
);


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
);

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
);