// Test interface for the GGM gradient engine.
//
// Exposes logp_and_gradient and forward_map to R for validation
// against finite differences.

#include <RcppArmadillo.h>
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"

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
