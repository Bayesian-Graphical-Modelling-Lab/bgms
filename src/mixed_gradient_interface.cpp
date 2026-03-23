// Test interface for the mixed MRF gradient engine.
//
// Exposes logp_and_gradient to R for finite-difference validation.

#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"

// [[Rcpp::export]]
Rcpp::List mixed_test_logp_and_gradient(
    const arma::vec& params,
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& edge_indicators,
    const std::string& pseudolikelihood,
    double pairwise_scale)
{
    size_t p = discrete_observations.n_cols;
    size_t q = continuous_observations.n_cols;
    size_t total = p + q;

    arma::mat inc_prob(total, total, arma::fill::value(0.5));
    bool edge_selection = false;

    MixedMRFModel model(
        discrete_observations, continuous_observations,
        num_categories, is_ordinal_variable, baseline_category,
        inc_prob, edge_indicators, edge_selection,
        pseudolikelihood,
        1.0, 1.0,   // main_alpha, main_beta
        pairwise_scale,
        42           // seed
    );

    auto result = model.logp_and_gradient(params);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}


// [[Rcpp::export]]
Rcpp::List mixed_test_logp_and_gradient_full(
    const arma::vec& params,
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& edge_indicators,
    const std::string& pseudolikelihood,
    double pairwise_scale)
{
    size_t p = discrete_observations.n_cols;
    size_t q = continuous_observations.n_cols;
    size_t total = p + q;

    arma::mat inc_prob(total, total, arma::fill::value(0.5));
    bool edge_selection = true;

    MixedMRFModel model(
        discrete_observations, continuous_observations,
        num_categories, is_ordinal_variable, baseline_category,
        inc_prob, edge_indicators, edge_selection,
        pseudolikelihood,
        1.0, 1.0,
        pairwise_scale,
        42
    );

    auto result = model.logp_and_gradient_full(params);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}


// [[Rcpp::export]]
Rcpp::List mixed_test_project_position(
    const arma::vec& x,
    const arma::vec& inv_mass,
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& edge_indicators,
    const std::string& pseudolikelihood,
    double pairwise_scale)
{
    size_t p = discrete_observations.n_cols;
    size_t q = continuous_observations.n_cols;
    size_t total = p + q;

    arma::mat inc_prob(total, total, arma::fill::value(0.5));
    bool edge_selection = true;

    MixedMRFModel model(
        discrete_observations, continuous_observations,
        num_categories, is_ordinal_variable, baseline_category,
        inc_prob, edge_indicators, edge_selection,
        pseudolikelihood,
        1.0, 1.0,
        pairwise_scale,
        42
    );

    arma::vec x_proj = x;
    model.project_position(x_proj, inv_mass);

    return Rcpp::List::create(
        Rcpp::Named("projected") = Rcpp::wrap(x_proj)
    );
}


// [[Rcpp::export]]
Rcpp::List mixed_test_project_momentum(
    const arma::vec& r,
    const arma::vec& x,
    const arma::vec& inv_mass,
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& edge_indicators,
    const std::string& pseudolikelihood,
    double pairwise_scale)
{
    size_t p = discrete_observations.n_cols;
    size_t q = continuous_observations.n_cols;
    size_t total = p + q;

    arma::mat inc_prob(total, total, arma::fill::value(0.5));
    bool edge_selection = true;

    MixedMRFModel model(
        discrete_observations, continuous_observations,
        num_categories, is_ordinal_variable, baseline_category,
        inc_prob, edge_indicators, edge_selection,
        pseudolikelihood,
        1.0, 1.0,
        pairwise_scale,
        42
    );

    arma::vec r_proj = r;
    model.project_momentum(r_proj, x, inv_mass);

    return Rcpp::List::create(
        Rcpp::Named("projected") = Rcpp::wrap(r_proj)
    );
}
