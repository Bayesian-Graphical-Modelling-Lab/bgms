// Test interface for the GGM gradient engine and RATTLE projection.
//
// Exposes logp_and_gradient, forward_map, project_position,
// project_momentum, and constrained leapfrog to R for validation.

#include <RcppArmadillo.h>
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "models/ggm/ggm_model.h"
#include "mcmc/algorithms/leapfrog.h"

// [[Rcpp::export]]
Rcpp::List ggm_test_logp_and_gradient(
    const arma::vec& theta,
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    double pairwise_scale)
{
    GraphConstraintStructure cs;
    cs.build(edge_indicators);

    GGMGradientEngine engine;
    engine.rebuild(cs, static_cast<size_t>(n), suf_stat, pairwise_scale);

    auto result = engine.logp_and_gradient(theta);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}

// [[Rcpp::export]]
Rcpp::List ggm_test_forward_map(
    const arma::vec& theta,
    const arma::imat& edge_indicators)
{
    GraphConstraintStructure cs;
    cs.build(edge_indicators);

    // Minimal engine just for forward map (no suf_stat needed)
    arma::mat dummy_S(edge_indicators.n_rows, edge_indicators.n_rows, arma::fill::zeros);
    GGMGradientEngine engine;
    engine.rebuild(cs, 100, dummy_S, 1.0);

    ForwardMapResult fm = engine.forward_map(theta);

    return Rcpp::List::create(
        Rcpp::Named("Phi") = Rcpp::wrap(fm.Phi),
        Rcpp::Named("K") = Rcpp::wrap(fm.K),
        Rcpp::Named("log_det_jacobian") = fm.log_det_jacobian,
        Rcpp::Named("psi") = Rcpp::wrap(fm.psi)
    );
}

// [[Rcpp::export]]
Rcpp::List ggm_test_project_position(
    const arma::vec& x,
    const arma::imat& edge_indicators)
{
    size_t p = edge_indicators.n_rows;

    // Build a minimal GGMModel from sufficient statistics
    arma::mat suf_stat = arma::eye(p, p);
    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true, 2.5);

    // Unpack x into the model's Cholesky factor
    model.set_full_position(x);

    // Project
    arma::vec x_proj = x;
    model.project_position(x_proj);

    // Unpack projected x to get Phi and K
    model.set_full_position(x_proj);

    // Compute K from projected Phi to verify constraints
    arma::vec full_pos = model.get_full_position();

    // Reconstruct Phi for output
    arma::mat Phi(p, p, arma::fill::zeros);
    GraphConstraintStructure cs;
    cs.build(edge_indicators);
    for (size_t q = 0; q < p; ++q) {
        size_t offset = cs.full_theta_offsets[q];
        for (size_t i = 0; i < q; ++i) {
            Phi(i, q) = x_proj(offset + i);
        }
        Phi(q, q) = std::exp(x_proj(offset + q));
    }
    arma::mat K = Phi.t() * Phi;

    return Rcpp::List::create(
        Rcpp::Named("x_projected") = Rcpp::wrap(x_proj),
        Rcpp::Named("Phi") = Rcpp::wrap(Phi),
        Rcpp::Named("K") = Rcpp::wrap(K)
    );
}

// [[Rcpp::export]]
arma::vec ggm_test_get_full_position(
    const arma::mat& Phi,
    const arma::imat& edge_indicators)
{
    size_t p = Phi.n_rows;

    arma::mat suf_stat = arma::eye(p, p);
    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true, 2.5);

    // Set the model's Cholesky factor directly, then get full position
    // Use a full-edge graph to set the Cholesky (no constraints bite)
    arma::vec x(p * (p + 1) / 2);
    GraphConstraintStructure cs;
    cs.build(edge_indicators);
    for (size_t q = 0; q < p; ++q) {
        size_t offset = cs.full_theta_offsets[q];
        for (size_t i = 0; i < q; ++i) {
            x(offset + i) = Phi(i, q);
        }
        x(offset + q) = std::log(Phi(q, q));
    }

    return x;
}

// [[Rcpp::export]]
Rcpp::List ggm_test_logp_and_gradient_full(
    const arma::vec& x,
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    double pairwise_scale)
{
    GraphConstraintStructure cs;
    cs.build(edge_indicators);

    GGMGradientEngine engine;
    engine.rebuild(cs, static_cast<size_t>(n), suf_stat, pairwise_scale);

    auto result = engine.logp_and_gradient_full(x);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}

// [[Rcpp::export]]
arma::vec ggm_test_project_momentum(
    const arma::vec& r,
    const arma::vec& x,
    const arma::imat& edge_indicators)
{
    size_t p = edge_indicators.n_rows;

    arma::mat suf_stat = arma::eye(p, p);
    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true, 2.5);

    arma::vec r_proj = r;
    model.project_momentum(r_proj, x);

    return r_proj;
}

// [[Rcpp::export]]
Rcpp::List ggm_test_leapfrog_constrained(
    const arma::vec& x0,
    const arma::vec& r0,
    double step_size,
    int n_steps,
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    double pairwise_scale)
{
    size_t p = edge_indicators.n_rows;

    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true, pairwise_scale);

    // Joint function: logp_and_gradient_full
    Memoizer::JointFn joint = [&model](const arma::vec& x)
        -> std::pair<double, arma::vec> {
        return model.logp_and_gradient_full(x);
    };
    Memoizer memo(joint);

    // Projection callback
    ProjectFn project = [&model](arma::vec& x, arma::vec& r) {
        model.project_position(x);
        model.project_momentum(r, x);
    };

    // Identity mass
    arma::vec inv_mass = arma::ones<arma::vec>(x0.n_elem);

    // Run n_steps constrained leapfrog steps
    arma::vec x = x0;
    arma::vec r = r0;
    double logp0 = memo.cached_log_post(x);

    for (int s = 0; s < n_steps; ++s) {
        std::tie(x, r) = leapfrog_constrained(
            x, r, step_size, memo, inv_mass, project
        );
    }

    double logp_final = memo.cached_log_post(x);
    double kin0 = 0.5 * arma::dot(r0, r0);
    double kin_final = 0.5 * arma::dot(r, r);
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
