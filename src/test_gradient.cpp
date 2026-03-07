// test_gradient.cpp - Rcpp-exported helper for numerical gradient checking
//
// Constructs a MixedMRFModel from R inputs and evaluates logp_and_gradient
// at a given parameter vector. Used by tests to compare analytical gradients
// against finite-difference approximations.
#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"

// Evaluate log-posterior and gradient at a given NUTS parameter vector.
//
// @param inputFromR  List with model specification (same as sample_mixed_mrf)
// @param theta       NUTS parameter vector (excludes Kyy)
// @param edge_selection Whether to do edge selection
//
// @return List with "logp" (double) and "gradient" (numeric vector)
// [[Rcpp::export]]
Rcpp::List test_mixed_gradient(
    const Rcpp::List& inputFromR,
    const arma::vec& theta,
    const bool edge_selection = false
) {
    arma::imat discrete_obs = Rcpp::as<arma::imat>(inputFromR["discrete_observations"]);
    arma::mat continuous_obs = Rcpp::as<arma::mat>(inputFromR["continuous_observations"]);
    arma::ivec num_categories = Rcpp::as<arma::ivec>(inputFromR["num_categories"]);
    arma::uvec is_ordinal = Rcpp::as<arma::uvec>(inputFromR["is_ordinal_variable"]);
    arma::ivec baseline_cat = Rcpp::as<arma::ivec>(inputFromR["baseline_category"]);
    double main_alpha = Rcpp::as<double>(inputFromR["main_alpha"]);
    double main_beta = Rcpp::as<double>(inputFromR["main_beta"]);
    double pairwise_scale = Rcpp::as<double>(inputFromR["pairwise_scale"]);
    std::string pseudolikelihood = Rcpp::as<std::string>(inputFromR["pseudolikelihood"]);

    int p = discrete_obs.n_cols;
    int q = continuous_obs.n_cols;
    int total = p + q;

    arma::mat inclusion_prob = arma::ones<arma::mat>(total, total);
    arma::imat edge_indicators = arma::ones<arma::imat>(total, total);

    MixedMRFModel model(
        discrete_obs, continuous_obs,
        num_categories, is_ordinal, baseline_cat,
        inclusion_prob, edge_indicators,
        edge_selection, pseudolikelihood,
        main_alpha, main_beta, pairwise_scale,
        42
    );

    // Set the NUTS block to the provided parameters
    model.set_vectorized_parameters(theta);

    auto result = model.logp_and_gradient(theta);

    return Rcpp::List::create(
        Rcpp::Named("logp") = result.first,
        Rcpp::Named("gradient") = Rcpp::NumericVector(result.second.begin(), result.second.end())
    );
}
