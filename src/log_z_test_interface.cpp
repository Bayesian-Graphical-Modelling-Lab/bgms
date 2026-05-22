// Test interface for the closed-form log Z(G) approximations used by the
// hierarchical-spec MH (Stage 3 of dev/plans/backlog/hierarchical-ggm-degord-rr.md).
// Exposes the C++ entry points to R for parity and incremental checks.

#include <RcppArmadillo.h>
#include "models/ggm/log_z_nlo.h"
#include "models/ggm/manuscript_nlo.h"
#include "models/ggm/degord_sampler.h"
#include "models/ggm/sd_density_at_zero.h"
#include "models/ggm/z_ratio_estimator.h"
#include "models/ggm/ggm_model.h"
#include "rng/rng_utils.h"
#include "math/cholesky_helpers.h"


// [[Rcpp::export]]
double log_Z_NLO_gamma_cpp(
    const arma::imat& G,
    double alpha, double beta, double sigma,
    bool include_F = false,
    double delta = 0.0
) {
    return log_Z_NLO_gamma(G, alpha, beta, sigma, include_F, delta);
}


// Manuscript App C NLO at alpha = 1 (companion-AI delivery 2026-05-21).
// Tracks ~/SV/Z/R/src/manuscript_NLO.h::log_Z_manuscript_NLO_alpha1.
//
// [[Rcpp::export]]
double log_Z_manuscript_NLO_alpha1_cpp(
    const arma::imat& G,
    double beta, double sigma, double delta
) {
    return ggm_nlo::log_Z_manuscript_NLO_alpha1(G, beta, sigma, delta);
}


// Manuscript NLO under DEGORD reordering (relabel toggle (i, j) -> (0, 1)).
//
// [[Rcpp::export]]
double log_Z_manuscript_NLO_alpha1_degord_cpp(
    const arma::imat& G, int i, int j,
    double beta, double sigma, double delta
) {
    return ggm_nlo::log_Z_manuscript_NLO_alpha1_degord(
        G, i, j, beta, sigma, delta);
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


// ---- Savage-Dickey 1D conditional density at K_ij = 0 -----------------
//
// SD spec sourced from ~/SV/Z/notes/2026-05-22_message-to-bgms-companion-SD-pivot.md
// and ~/SV/Z/R/src/sd_density_at_zero.cpp. (i, j) are 1-indexed at the R
// boundary to match Z's exposed function signature for the bit-for-bit test.

// [[Rcpp::export]]
Rcpp::List sd_log_density_at_zero_cpp(
    const arma::mat& K, int i, int j,
    const arma::mat& S, int n_obs,
    double delta, double sigma,
    bool nlo = true,
    bool apply_pd_truncation = false,
    int newton_max_iter = 50,
    double newton_tol = 1e-10
) {
    ggm_sd::SDResult r = ggm_sd::density_at_zero_one(
        K, i - 1, j - 1, S, n_obs, delta, sigma,
        nlo, apply_pd_truncation, newton_max_iter, newton_tol);
    return Rcpp::List::create(
        Rcpp::Named("log_density") = r.log_density,
        Rcpp::Named("x_mode")      = r.x_mode,
        Rcpp::Named("curvature")   = r.curvature,
        Rcpp::Named("x_minus")     = r.x_minus,
        Rcpp::Named("x_plus")      = r.x_plus,
        Rcpp::Named("log_Z_trunc") = r.log_Z_trunc,
        Rcpp::Named("status")      = r.status);
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


// Log-space V test interface. Returns (log_abs, sign) for log|V| and sign(V).
// sign = 0 with log_abs = NaN signals an auto-reject sentinel (S = 0 or any
// non-finite intermediate).
//
// [[Rcpp::export]]
Rcpp::List degord_V_log_at_Gamma_pi_cpp(
    int K_depth,
    const Rcpp::List& pools_t,
    const arma::imat& G_pi,
    double alpha, double beta, double sigma, double delta,
    double log_c, double rho,
    int slab_tilt_mode = 0
) {
    int q = G_pi.n_rows;
    auto chain_aux = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    chain_aux.slab_tilt_mode = slab_tilt_mode;
    std::vector<arma::mat> pools_t_cpp;
    pools_t_cpp.reserve(static_cast<size_t>(K_depth));
    for (int n = 0; n < K_depth; ++n)
        pools_t_cpp.push_back(Rcpp::as<arma::mat>(pools_t[n]));
    auto res = degord::V_log_at_Gamma_pi_degord(
        K_depth, pools_t_cpp, G_pi, chain_aux, log_c, rho);
    return Rcpp::List::create(
        Rcpp::Named("log_abs") = res.first,
        Rcpp::Named("sign")    = res.second
    );
}


// Paired log-space V at (Γ_curr, Γ_star) with within-pool cache reuse.
// Returns a list with `curr` and `star` sub-lists, each carrying
// (log_abs, sign).
//
// [[Rcpp::export]]
Rcpp::List degord_V_log_pair_at_Gamma_curr_star_cpp(
    int K_depth,
    const Rcpp::List& pools_t,
    const arma::imat& G_pi_curr,
    const arma::imat& G_pi_star,
    double alpha, double beta, double sigma, double delta,
    double log_c_curr, double log_c_star, double rho,
    int slab_tilt_mode = 0
) {
    int q = G_pi_curr.n_rows;
    auto chain_aux = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    chain_aux.slab_tilt_mode = slab_tilt_mode;
    std::vector<arma::mat> pools_t_cpp;
    pools_t_cpp.reserve(static_cast<size_t>(K_depth));
    for (int n = 0; n < K_depth; ++n)
        pools_t_cpp.push_back(Rcpp::as<arma::mat>(pools_t[n]));
    auto res = degord::V_log_pair_at_Gamma_curr_star_degord(
        K_depth, pools_t_cpp, G_pi_curr, G_pi_star, chain_aux,
        log_c_curr, log_c_star, rho);
    return Rcpp::List::create(
        Rcpp::Named("curr") = Rcpp::List::create(
            Rcpp::Named("log_abs") = res.curr.first,
            Rcpp::Named("sign")    = res.curr.second),
        Rcpp::Named("star") = Rcpp::List::create(
            Rcpp::Named("log_abs") = res.star.first,
            Rcpp::Named("sign")    = res.star.second)
    );
}


// log_Zhat under G_pi_star using cached state built under G_pi_curr.
// Used by tests to validate that the cache adapter matches a fresh
// log_Zhat_pi_from_pool at G_pi_star to FP-reordering tolerance.
//
// [[Rcpp::export]]
double degord_log_Zhat_star_from_cache_cpp(
    const arma::mat& noise_pool_t,
    const arma::imat& G_pi_curr,
    const arma::imat& G_pi_star,
    double alpha, double beta, double sigma, double delta,
    int slab_tilt_mode = 0
) {
    int q = G_pi_curr.n_rows;
    auto chain_aux = degord::make_chain_aux(q, alpha, beta, sigma, delta);
    chain_aux.slab_tilt_mode = slab_tilt_mode;
    auto a_curr = degord::make_pi_aux(G_pi_curr, chain_aux);
    auto a_star = degord::make_pi_aux(G_pi_star, chain_aux);
    degord::PoolCache cache_curr;
    (void) degord::log_Zhat_pi_from_pool_cache(
        noise_pool_t, a_curr, chain_aux, cache_curr);
    return degord::log_Zhat_star_from_cache(
        noise_pool_t, a_star, chain_aux, cache_curr);
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


// ---- Phase 4a: hierarchical-spec smoke test ----------------------------

// Constructs a small GGMModel with Normal slab + Gamma diagonal, switches to
// hierarchical-spec, runs n_sweeps of MH, and returns the final edge
// indicators and a few summary statistics. Crashes here are a regression
// in the Phase 4a wiring.

// [[Rcpp::export]]
Rcpp::List ggm_hierarchical_smoke_cpp(
    const arma::mat& observations,
    double inclusion_prob,
    double interaction_scale,    // sigma for Normal slab
    double diagonal_shape,       // alpha for Gamma diag
    double diagonal_rate,        // beta for Gamma diag
    double delta,                // determinant tilt
    int    M_inner,
    double kappa,
    double rho,
    int    n_sweeps,
    int    seed,
    bool   use_manuscript_nlo = false
) {
    int p = observations.n_cols;
    arma::mat inclusion_probability(p, p, arma::fill::value(inclusion_prob));
    arma::imat initial_edges(p, p, arma::fill::zeros);
    // Start at empty graph: edge_indicators are 0 off-diagonal, 1 on the
    // diagonal (by GGMModel convention).
    for (int i = 0; i < p; ++i) initial_edges(i, i) = 1;

    auto slab = std::make_unique<NormalPrior>(interaction_scale);
    auto diag = std::make_unique<GammaScalePrior>(diagonal_shape, diagonal_rate);

    GGMModel model(observations,
                   inclusion_probability,
                   initial_edges,
                   /*edge_selection=*/true,
                   std::move(slab),
                   std::move(diag),
                   /*na_impute=*/false);
    model.set_seed(seed);
    model.set_determinant_tilt(delta);
    model.set_z_ratio_tuning(M_inner, kappa, rho);
    model.set_use_manuscript_nlo(use_manuscript_nlo);
    model.set_graph_prior_spec(GraphPriorSpec::Hierarchical);

    arma::ivec n_edges(n_sweeps, arma::fill::zeros);
    for (int s = 0; s < n_sweeps; ++s) {
        model.prepare_iteration();
        model.update_edge_indicators();
        const arma::imat& E = model.get_edge_indicators();
        n_edges[s] = arma::accu(E) / 2;  // off-diagonal edges (E is symmetric)
        // Discount the diagonal-1 convention if it applies. Standard ggm
        // counts edges as accu(upper-tri) which is accu(E)/2 here.
    }

    return Rcpp::List::create(
        Rcpp::Named("final_edges")  = model.get_edge_indicators(),
        Rcpp::Named("n_edges_path") = n_edges
    );
}


// ---- Savage-Dickey between-Γ smoke test (handoff 2026-05-22) -----------
//
// Constructs a small GGMModel with Normal slab + Gamma diagonal, switches
// to the SD between-step variant, runs n_sweeps of (within-K + between-Γ)
// MH, and returns the per-edge inclusion-frequency matrix plus the
// n_edges path. Prior-only mode disables the likelihood contribution in
// every MH ratio; under prior_only the inclusion frequencies should
// converge to inclusion_prob.
//
// observations: n × p data matrix (Y). When prior_only = true, the entries
// are ignored except for n_ and p_; pass any conformable matrix.
//
// Plug-in DEGORD smoke test — apples-to-apples counterpart of ggm_sd_smoke_cpp.
// Runs the hierarchical-prior chain with plug-in mNLO ratio (closed-form,
// V/RR machinery bypassed), prior-only optional, PIP accumulator over the
// post-warmup window. Used to compare β-Bernoulli edge marginals against the
// SD chain under identical priors.
//
// [[Rcpp::export]]
Rcpp::List ggm_plug_in_smoke_cpp(
    const arma::mat& observations,
    double inclusion_prob,
    double interaction_scale,    // sigma for Normal slab
    double diagonal_shape,       // alpha for Gamma diag
    double diagonal_rate,        // beta for Gamma diag
    double delta,                // determinant tilt
    int    M_inner,              // V/RR tuning (unused under plug-in; kept for parity)
    double kappa,
    double rho,
    int    n_warmup,
    int    n_sweeps,
    int    seed,
    bool   prior_only       = false,
    bool   include_within_k = true,
    bool   use_manuscript_nlo = false
) {
    int p = observations.n_cols;
    arma::mat inclusion_probability(p, p, arma::fill::value(inclusion_prob));
    arma::imat initial_edges(p, p, arma::fill::zeros);
    for (int i = 0; i < p; ++i) initial_edges(i, i) = 1;

    auto slab = std::make_unique<NormalPrior>(interaction_scale);
    auto diag = std::make_unique<GammaScalePrior>(diagonal_shape, diagonal_rate);

    GGMModel model(observations,
                   inclusion_probability,
                   initial_edges,
                   /*edge_selection=*/true,
                   std::move(slab),
                   std::move(diag),
                   /*na_impute=*/false);
    model.set_seed(seed);
    model.set_determinant_tilt(delta);
    model.set_z_ratio_tuning(M_inner, kappa, rho);
    model.set_use_manuscript_nlo(use_manuscript_nlo);
    model.set_graph_prior_spec(GraphPriorSpec::Hierarchical);
    model.set_plug_in_nlo(true);
    model.set_prior_only(prior_only);

    arma::mat pip_counts(p, p, arma::fill::zeros);
    arma::ivec n_edges_path(n_sweeps, arma::fill::zeros);

    int n_total = n_warmup + n_sweeps;
    for (int s = 0; s < n_total; ++s) {
        model.prepare_iteration();
        if (include_within_k) {
            model.do_one_metropolis_step(s);
        }
        model.update_edge_indicators();
        if (s >= n_warmup) {
            const arma::imat& E = model.get_edge_indicators();
            for (int ii = 0; ii < p; ++ii) {
                for (int jj = 0; jj < p; ++jj) {
                    if (E(ii, jj) == 1 && ii != jj) {
                        pip_counts(ii, jj) += 1.0;
                    }
                }
            }
            n_edges_path[s - n_warmup] = arma::accu(E) / 2;
        }
    }

    arma::mat pip = pip_counts / static_cast<double>(n_sweeps);
    return Rcpp::List::create(
        Rcpp::Named("pip")          = pip,
        Rcpp::Named("final_edges")  = model.get_edge_indicators(),
        Rcpp::Named("n_edges_path") = n_edges_path
    );
}


// Test interface for cholesky_helpers::perm_to_trailing_2x2.
// Accepts a precision matrix K and 1-based edge indices (i, j); internally
// computes U = chol(K, "upper") and delegates. Returns the trailing 2×2 of
// the permuted lower-triangular factor: (l_ii, l_ji, l_jj) and an ok flag.
//
// [[Rcpp::export]]
Rcpp::List chol_perm_trailing_2x2_cpp(
    const arma::mat& K, int i_1based, int j_1based
) {
    arma::mat U;
    if (!arma::chol(U, K, "upper")) {
        return Rcpp::List::create(
            Rcpp::Named("l_ii") = NA_REAL,
            Rcpp::Named("l_ji") = NA_REAL,
            Rcpp::Named("l_jj") = NA_REAL,
            Rcpp::Named("ok")   = false
        );
    }
    const std::size_t i0 = static_cast<std::size_t>(i_1based - 1);
    const std::size_t j0 = static_cast<std::size_t>(j_1based - 1);
    auto r = cholesky_helpers::perm_to_trailing_2x2(U, i0, j0);
    return Rcpp::List::create(
        Rcpp::Named("l_ii") = r.l_ii,
        Rcpp::Named("l_ji") = r.l_ji,
        Rcpp::Named("l_jj") = r.l_jj,
        Rcpp::Named("ok")   = r.ok
    );
}


// [[Rcpp::export]]
Rcpp::List ggm_sd_smoke_cpp(
    const arma::mat& observations,
    double inclusion_prob,
    double interaction_scale,    // sigma for Normal slab
    double diagonal_shape,       // alpha for Gamma diag
    double diagonal_rate,        // beta for Gamma diag
    double delta,                // determinant tilt
    int    n_warmup,
    int    n_sweeps,
    int    seed,
    bool   prior_only = false,
    bool   include_within_k = true,
    bool   use_lspace = false
) {
    int p = observations.n_cols;
    arma::mat inclusion_probability(p, p, arma::fill::value(inclusion_prob));
    arma::imat initial_edges(p, p, arma::fill::zeros);
    for (int i = 0; i < p; ++i) initial_edges(i, i) = 1;

    auto slab = std::make_unique<NormalPrior>(interaction_scale);
    auto diag = std::make_unique<GammaScalePrior>(diagonal_shape, diagonal_rate);

    GGMModel model(observations,
                   inclusion_probability,
                   initial_edges,
                   /*edge_selection=*/true,
                   std::move(slab),
                   std::move(diag),
                   /*na_impute=*/false);
    model.set_seed(seed);
    model.set_determinant_tilt(delta);
    model.set_use_sd_between_step(true);
    model.set_use_sd_lspace(use_lspace);
    model.set_prior_only(prior_only);
    model.reset_pd_revert_count();

    // PIP accumulator: count visits with γ_ij = 1 across post-warmup sweeps.
    arma::mat pip_counts(p, p, arma::fill::zeros);
    arma::ivec n_edges_path(n_sweeps, arma::fill::zeros);

    int n_total = n_warmup + n_sweeps;
    for (int s = 0; s < n_total; ++s) {
        model.prepare_iteration();
        if (include_within_k) {
            model.do_one_metropolis_step(s);
        }
        model.update_edge_indicators();
        if (s >= n_warmup) {
            const arma::imat& E = model.get_edge_indicators();
            for (int ii = 0; ii < p; ++ii) {
                for (int jj = 0; jj < p; ++jj) {
                    if (E(ii, jj) == 1 && ii != jj) {
                        pip_counts(ii, jj) += 1.0;
                    }
                }
            }
            n_edges_path[s - n_warmup] = arma::accu(E) / 2;
        }
    }

    arma::mat pip = pip_counts / static_cast<double>(n_sweeps);

    const long n_pd_reverts = model.n_pd_reverts();
    if (n_pd_reverts > 0) {
        Rcpp::warning(
            "GGM SD L-space chain: PD-revert defense fired %ld time(s) during "
            "this run. Canary for K landing numerically off the PD cone near "
            "the boundary; the chain reverted the offending toggle.",
            n_pd_reverts);
    }
    return Rcpp::List::create(
        Rcpp::Named("pip")           = pip,
        Rcpp::Named("final_edges")   = model.get_edge_indicators(),
        Rcpp::Named("n_edges_path")  = n_edges_path,
        Rcpp::Named("n_pd_reverts")  = n_pd_reverts
    );
}
