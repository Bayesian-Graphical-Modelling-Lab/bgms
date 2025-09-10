#pragma once
#include <RcppArmadillo.h>

struct SafeRNG;

Rcpp::List run_gibbs_sampler_for_bgmCompare(
    int chain_id,
    arma::imat observations,
    const int num_groups,
    Rcpp::List num_obs_categories,
    Rcpp::List sufficient_blume_capel,
    Rcpp::List sufficient_pairwise,
    const arma::ivec& num_categories,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const double difference_scale,//new
    const double difference_selection_alpha,//new
    const double difference_selection_beta,//new
    const std::string difference_prior,//new
    const int iter,
    const int burnin,
    const bool na_impute,
    const arma::imat& missing_data_indices,//updated
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const bool difference_selection,//new
    const arma::imat main_effect_indices,
    const arma::imat pairwise_effect_indices,
    const double target_accept,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const arma::mat projection,//new
    const arma::ivec& group_membership,//new
    const arma::imat& group_indices,//new
    const arma::imat& interaction_index_matrix,//new
    arma::mat inclusion_probability,
    SafeRNG& rng
);