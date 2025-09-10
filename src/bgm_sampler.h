#pragma once
#include <RcppArmadillo.h>

// forward declaration
struct SafeRNG;

Rcpp::List run_gibbs_sampler_for_bgm(
    int chain_id,
    arma::imat observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    const std::string& edge_prior,
    arma::mat inclusion_probability,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& interaction_index_matrix,
    const int iter,
    const int burnin,
    arma::imat num_obs_categories,
    arma::imat sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat pairwise_effect_indices,
    const double target_accept,
    arma::imat sufficient_pairwise,
    const int hmc_num_leapfrogs,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    SafeRNG& rng
);