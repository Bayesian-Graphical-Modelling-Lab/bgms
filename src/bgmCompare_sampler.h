#pragma once

#include <RcppArmadillo.h>
#include <string>

struct SamplerOutput;
struct SafeRNG;
class ProgressManager;

SamplerOutput run_gibbs_sampler_bgmCompare(
    int chain_id,
    arma::imat observations,
    const int num_groups,
    std::vector<arma::imat>& counts_per_category,
    std::vector<arma::imat>& blume_capel_stats,
    std::vector<arma::mat>& pairwise_stats,
    const arma::ivec& num_categories,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const double difference_scale,
    const double difference_selection_alpha,
    const double difference_selection_beta,
    const std::string& difference_prior,
    const int iter,
    const int burnin,
    const bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const bool difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    const double target_accept,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    arma::mat inclusion_probability,
    SafeRNG& rng,
    const std::string& update_method,
    const int hmc_num_leapfrogs,
    ProgressManager& pm
);