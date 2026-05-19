// Test interface for the closed-form log Z(G) approximations used by the
// hierarchical-spec MH (Stage 3 of dev/plans/backlog/hierarchical-ggm-degord-rr.md).
// Exposes the C++ entry points to R for parity and incremental checks.

#include <RcppArmadillo.h>
#include "models/ggm/log_z_nlo.h"


// [[Rcpp::export]]
double log_Z_NLO_gamma_cpp(
    const arma::imat& G,
    double alpha, double beta, double sigma,
    bool include_F = false,
    double delta = 0.0
) {
    return log_Z_NLO_gamma(G, alpha, beta, sigma, include_F, delta);
}


// [[Rcpp::export]]
double log_Z_NLO_gamma_degord_cpp(
    const arma::imat& G, int i, int j,
    double alpha, double beta, double sigma,
    bool include_F = false,
    double delta = 0.0
) {
    return log_Z_NLO_gamma_degord(G, i, j, alpha, beta, sigma, include_F, delta);
}


// [[Rcpp::export]]
double log_Z_NLO_gamma_delta_incr_alpha1_cpp(
    const arma::imat& G_before, int i, int j,
    double beta, double sigma, double delta,
    bool include_F = false
) {
    return log_Z_NLO_gamma_delta_incr_alpha1(
        G_before, i, j, beta, sigma, delta, include_F);
}


// [[Rcpp::export]]
double log_Z_NLO_gamma_delta_incr_alphaN_cpp(
    const arma::imat& G_before, int i, int j,
    double alpha, double beta, double sigma, double delta,
    bool include_F = false
) {
    return log_Z_NLO_gamma_delta_incr_alphaN(
        G_before, i, j, alpha, beta, sigma, delta, include_F);
}
