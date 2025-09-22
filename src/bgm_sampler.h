#pragma once
#include <RcppArmadillo.h>
// forward declaration
struct SafeRNG;
class ProgressManager;

Rcpp::List run_gibbs_sampler_bgm(
    int chain_id,
    arma::imat observations,
    const arma::ivec& num_categories,
    const double pairwise_scale,
    const std::string& edge_prior,
    arma::mat inclusion_probability,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& interaction_index_matrix,
    const int iter,
    const int burnin,
    arma::imat counts_per_category,
    arma::imat blume_capel_stats,
    const double main_alpha,
    const double main_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat pairwise_effect_indices,
    const double target_accept,
    arma::imat pairwise_stats,
    const int hmc_num_leapfrogs,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    SafeRNG& rng,
    ProgressManager& pm
);