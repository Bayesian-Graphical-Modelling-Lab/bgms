// Test interface for the mixed MRF gradient engine.
//
// Exposes logp_and_gradient, project_position, project_momentum,
// and constrained leapfrog to R for validation.

#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "mcmc/algorithms/leapfrog.h"

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


// [[Rcpp::export]]
Rcpp::List mixed_test_leapfrog_constrained(
    const arma::vec& x0,
    const arma::vec& r0,
    double step_size,
    int n_steps,
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& edge_indicators,
    const std::string& pseudolikelihood,
    double pairwise_scale,
    Rcpp::Nullable<Rcpp::NumericVector> inv_mass_in = R_NilValue)
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

    Memoizer::JointFn joint = [&model](const arma::vec& x)
        -> std::pair<double, arma::vec> {
        return model.logp_and_gradient_full(x);
    };
    Memoizer memo(joint);

    arma::vec inv_mass;
    if(inv_mass_in.isNotNull()) {
        inv_mass = Rcpp::as<arma::vec>(inv_mass_in);
    } else {
        inv_mass = arma::ones<arma::vec>(x0.n_elem);
    }

    ProjectPositionFn proj_pos = [&model, &inv_mass](arma::vec& x) {
        model.project_position(x, inv_mass);
    };
    ProjectMomentumFn proj_mom = [&model, &inv_mass](arma::vec& r, const arma::vec& x) {
        model.project_momentum(r, x, inv_mass);
    };

    arma::vec x = x0;
    arma::vec r = r0;
    double logp0 = memo.cached_log_post(x);

    for(int s = 0; s < n_steps; ++s) {
        std::tie(x, r) = leapfrog_constrained(
            x, r, step_size, memo, inv_mass, proj_pos, proj_mom
        );
    }

    double logp_final = memo.cached_log_post(x);
    double kin0 = 0.5 * arma::dot(r0, inv_mass % r0);
    double kin_final = 0.5 * arma::dot(r, inv_mass % r);
    double H0 = -logp0 + kin0;
    double H_final = -logp_final + kin_final;

    return Rcpp::List::create(
        Rcpp::Named("x") = Rcpp::wrap(x),
        Rcpp::Named("r") = Rcpp::wrap(r),
        Rcpp::Named("logp0") = logp0,
        Rcpp::Named("logp_final") = logp_final,
        Rcpp::Named("H0") = H0,
        Rcpp::Named("H_final") = H_final,
        Rcpp::Named("dH") = H_final - H0
    );
}
