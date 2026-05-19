// Test interface for the closed-form log Z(G) approximations used by the
// hierarchical-spec MH (Stage 3 of dev/plans/backlog/hierarchical-ggm-degord-rr.md).
// Exposes the C++ entry points to R for parity and incremental checks.

#include <RcppArmadillo.h>
#include "models/ggm/log_z_nlo.h"
#include "models/ggm/degord_sampler.h"
#include "models/ggm/z_ratio_estimator.h"
#include "rng/rng_utils.h"


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


// ---- DEGORD sampler test interface ------------------------------------

// [[Rcpp::export]]
Rcpp::List degord_chain_aux_cpp(
    int q, double alpha, double beta, double sigma, double delta
) {
    auto c = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    return Rcpp::List::create(
        Rcpp::Named("q")                        = c.q,
        Rcpp::Named("alpha")                    = c.alpha,
        Rcpp::Named("beta")                     = c.beta,
        Rcpp::Named("sigma")                    = c.sigma,
        Rcpp::Named("delta")                    = c.delta,
        Rcpp::Named("sigma_diag")               = c.sigma_diag,
        Rcpp::Named("nu_chi_df")                = Rcpp::wrap(c.nu_chi_df),
        Rcpp::Named("nu_mu_l")                  = Rcpp::wrap(c.nu_mu_l),
        Rcpp::Named("nu_H_e_saddle")            = Rcpp::wrap(c.nu_H_e_saddle),
        Rcpp::Named("nu_lgamma_half_k")         = Rcpp::wrap(c.nu_lgamma_half_k),
        Rcpp::Named("nu_diag_const")            = Rcpp::wrap(c.nu_diag_const),
        Rcpp::Named("nu_slab_const_saddle")     = Rcpp::wrap(c.nu_slab_const_saddle),
        Rcpp::Named("nu_per_vertex")            = Rcpp::wrap(c.nu_per_vertex)
    );
}


// [[Rcpp::export]]
Rcpp::List degord_pi_aux_cpp(
    const arma::imat& G_pi,
    double alpha, double beta, double sigma, double delta
) {
    int q = G_pi.n_rows;
    auto c = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    auto a = degord::make_pi_aux(G_pi, c);
    return Rcpp::List::create(
        Rcpp::Named("q")       = a.q,
        Rcpp::Named("nu_pi")   = Rcpp::wrap(a.nu_pi),
        Rcpp::Named("E_count") = a.E_count,
        Rcpp::Named("log_C0")  = a.log_C0
    );
}


// [[Rcpp::export]]
arma::imat degord_permute_graph_cpp(const arma::imat& G, int i, int j) {
    int q = G.n_rows;
    auto pi = degord::degord_permutation(q, i, j);
    return degord::permute_graph(G, pi);
}


// [[Rcpp::export]]
double degord_log_Zhat_pi_from_pool_cpp(
    const arma::mat& noise_pool_t,
    const arma::imat& G_pi,
    double alpha, double beta, double sigma, double delta,
    int slab_tilt_mode = 0
) {
    int q = G_pi.n_rows;
    auto c = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    c.slab_tilt_mode = slab_tilt_mode;
    auto a = degord::make_pi_aux(G_pi, c);
    return degord::log_Zhat_pi_from_pool(noise_pool_t, a, c);
}


// [[Rcpp::export]]
double degord_delta_log_Zhat_pi_toggle_cpp(
    const arma::mat& noise_pool,
    const arma::mat& noise_pool_t,
    const arma::imat& G_curr,
    int i, int j,
    double alpha, double beta, double sigma, double delta,
    int slab_tilt_mode = 0
) {
    int q = G_curr.n_rows;
    auto c = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    c.slab_tilt_mode = slab_tilt_mode;
    return degord::delta_log_Zhat_pi_toggle(
        noise_pool, noise_pool_t, G_curr, i, j, c);
}


// [[Rcpp::export]]
arma::mat degord_draw_bartlett_pool_cpp(int q, int M_inner, int seed) {
    SafeRNG rng(seed);
    return degord::draw_bartlett_pool(rng, q, M_inner);
}


// ---- Phase 3: V(Γ, U) Russian-Roulette estimator test interface --------

// [[Rcpp::export]]
double degord_V_at_Gamma_pi_cpp(
    int K_depth,
    const Rcpp::List& pools_t,
    const arma::imat& G_pi,
    double alpha, double beta, double sigma, double delta,
    double c_val, double rho,
    int slab_tilt_mode = 0
) {
    int q = G_pi.n_rows;
    auto chain_aux = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    chain_aux.slab_tilt_mode = slab_tilt_mode;
    std::vector<arma::mat> pools_t_cpp;
    pools_t_cpp.reserve(static_cast<size_t>(K_depth));
    for (int n = 0; n < K_depth; ++n)
        pools_t_cpp.push_back(Rcpp::as<arma::mat>(pools_t[n]));
    return degord::V_at_Gamma_pi_degord(
        K_depth, pools_t_cpp, G_pi, chain_aux, c_val, rho);
}


// [[Rcpp::export]]
Rcpp::List degord_draw_U_rr_cpp(int M_inner, int q, double rho, int seed) {
    SafeRNG rng(seed);
    int K_depth = 0;
    std::vector<arma::mat> pools_t;
    degord::draw_U_degord_rr(rng, K_depth, pools_t, M_inner, q, rho);
    Rcpp::List pools_R(K_depth);
    for (int n = 0; n < K_depth; ++n) pools_R[n] = pools_t[n];
    return Rcpp::List::create(
        Rcpp::Named("K_depth") = K_depth,
        Rcpp::Named("pools_t") = pools_R
    );
}
