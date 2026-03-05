// test_mixed_mrf.cpp — Rcpp test helpers for MixedMRFModel skeleton (Phase A.3)
//
// These are lightweight exported functions that construct a MixedMRFModel from
// R inputs and return diagnostic information for testthat assertions.

#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"

// [[Rcpp::export]]
Rcpp::List test_mixed_mrf_skeleton(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection,
    const std::string& pseudolikelihood,
    int seed
) {
    MixedMRFModel model(
        discrete_observations,
        continuous_observations,
        num_categories,
        is_ordinal_variable,
        baseline_category,
        inclusion_probability,
        initial_edge_indicators,
        edge_selection,
        pseudolikelihood,
        1.0, 1.0, 2.5,
        seed
    );

    // --- Dimension checks ---
    size_t param_dim = model.parameter_dimension();
    size_t full_param_dim = model.full_parameter_dimension();
    int num_variables = model.get_num_variables();
    int num_pairwise = model.get_num_pairwise();

    // --- Vectorized parameter round-trip ---
    arma::vec params_before = model.get_vectorized_parameters();
    model.set_vectorized_parameters(params_before);
    arma::vec params_after = model.get_vectorized_parameters();
    double roundtrip_max_diff = arma::max(arma::abs(params_before - params_after));

    // --- Full-parameter vector ---
    arma::vec full_params = model.get_full_vectorized_parameters();

    // --- Edge indicators ---
    arma::ivec indicator_vec = model.get_vectorized_indicator_parameters();
    const arma::imat& edge_mat = model.get_edge_indicators();

    // --- Non-trivial parameter round-trip ---
    // Set some non-zero parameters and verify round-trip
    arma::vec test_params(full_param_dim);
    for(size_t i = 0; i < full_param_dim; ++i) {
        test_params(i) = static_cast<double>(i + 1) * 0.01;
    }
    // Kyy must be SPD: overwrite the Kyy block with a diagonal-dominant SPD matrix.
    // Kyy occupies indices: num_main + num_pairwise_xx + q ... num_main + num_pairwise_xx + q + q(q+1)/2 - 1
    // For simplicity, just make it diagonal = large positive
    int p = discrete_observations.n_cols;
    int q = continuous_observations.n_cols;
    size_t num_main = 0;
    for(int s = 0; s < p; ++s) {
        if(is_ordinal_variable(s)) {
            num_main += num_categories(s);
        } else {
            num_main += 2;
        }
    }
    size_t kyy_start = num_main + p * (p - 1) / 2 + q;
    // Write identity-like SPD into Kyy block
    size_t kyy_idx = kyy_start;
    for(int i = 0; i < q; ++i) {
        for(int j = i; j < q; ++j) {
            if(i == j) {
                test_params(kyy_idx) = 2.0 + i * 0.1;  // positive diagonal
            } else {
                test_params(kyy_idx) = 0.01;  // small off-diagonal
            }
            kyy_idx++;
        }
    }

    model.set_vectorized_parameters(test_params);
    arma::vec recovered = model.get_vectorized_parameters();
    double nontrivial_max_diff = arma::max(arma::abs(test_params - recovered));

    // --- Clone round-trip ---
    auto cloned = model.clone();
    arma::vec cloned_params = cloned->get_vectorized_parameters();
    double clone_max_diff = arma::max(arma::abs(recovered - cloned_params));

    return Rcpp::List::create(
        Rcpp::Named("parameter_dimension") = param_dim,
        Rcpp::Named("full_parameter_dimension") = full_param_dim,
        Rcpp::Named("num_variables") = num_variables,
        Rcpp::Named("num_pairwise") = num_pairwise,
        Rcpp::Named("params_length") = params_before.n_elem,
        Rcpp::Named("roundtrip_max_diff") = roundtrip_max_diff,
        Rcpp::Named("nontrivial_roundtrip_max_diff") = nontrivial_max_diff,
        Rcpp::Named("clone_max_diff") = clone_max_diff,
        Rcpp::Named("indicator_length") = indicator_vec.n_elem,
        Rcpp::Named("edge_indicators_rows") = edge_mat.n_rows,
        Rcpp::Named("edge_indicators_cols") = edge_mat.n_cols,
        Rcpp::Named("has_edge_selection") = model.has_edge_selection(),
        Rcpp::Named("has_adaptive_metropolis") = model.has_adaptive_metropolis(),
        Rcpp::Named("has_gradient") = model.has_gradient(),
        Rcpp::Named("has_missing_data") = model.has_missing_data()
    );
}


// [[Rcpp::export]]
Rcpp::List test_mixed_mrf_likelihoods(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection,
    const std::string& pseudolikelihood,
    const arma::vec& params,
    int seed
) {
    MixedMRFModel model(
        discrete_observations,
        continuous_observations,
        num_categories,
        is_ordinal_variable,
        baseline_category,
        inclusion_probability,
        initial_edge_indicators,
        edge_selection,
        pseudolikelihood,
        1.0, 1.0, 2.5,
        seed
    );

    // Set parameters (caller must ensure Kyy block is SPD)
    model.set_vectorized_parameters(params);

    int p = discrete_observations.n_cols;

    // Evaluate log_conditional_omrf for each discrete variable
    Rcpp::NumericVector omrf_ll(p);
    for(int s = 0; s < p; ++s) {
        omrf_ll[s] = model.log_conditional_omrf(s);
    }

    // Evaluate log_conditional_ggm
    double ggm_ll = model.log_conditional_ggm();

    return Rcpp::List::create(
        Rcpp::Named("omrf_ll") = omrf_ll,
        Rcpp::Named("ggm_ll") = ggm_ll
    );
}
