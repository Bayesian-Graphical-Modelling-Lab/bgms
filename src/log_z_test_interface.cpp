// Test fixtures for the L-space Savage-Dickey primitives and the GGM SD
// chain. The bgm() production path routes hierarchical-spec to the same
// SD machinery; this file exists so the SD primitives and a thin chain
// driver are reachable from R for SBC / PD diagnostics that don't want to
// go through the full bgm() wrapper.

#include <RcppArmadillo.h>
#include "models/ggm/sd_density_l_space.h"
#include "models/ggm/sd_density_l_space_quad.h"
#include "models/ggm/ggm_model.h"
#include "rng/rng_utils.h"
#include "math/cholesky_helpers.h"
#include "mcmc/algorithms/nuts.h"
#include "mcmc/algorithms/hmc.h"
#include "mcmc/samplers/nuts_adaptation.h"


// Test interface for ggm_sd::density_at_l_ji_one (Laplace + optional NLO
// 1D conditional density of l_ji at x_eval).
//
// [[Rcpp::export]]
Rcpp::List sd_log_density_at_l_ji_cpp(
    double x_eval, double A, double B, double s_jj, double alpha,
    bool nlo = true, int newton_max_iter = 50, double newton_tol = 1e-10
) {
    auto r = ggm_sd::density_at_l_ji_one(x_eval, A, B, s_jj, alpha,
                                          nlo, newton_max_iter, newton_tol);
    return Rcpp::List::create(
        Rcpp::Named("log_density") = r.log_density,
        Rcpp::Named("x_mode")      = r.x_mode,
        Rcpp::Named("curvature")   = r.curvature,
        Rcpp::Named("status")      = r.status
    );
}


// Test interface for ggm_sd::density_at_l_ji_gh (Gauss-Hermite quadrature
// variant; ~64 evaluations per call but more reliable than Laplace+NLO
// across all chain configurations).
//
// [[Rcpp::export]]
Rcpp::List sd_log_density_at_l_ji_gh_cpp(
    double x_eval, double A, double B, double s_jj, double alpha,
    int num_nodes = 64
) {
    auto r = ggm_sd::density_at_l_ji_gh(x_eval, A, B, s_jj, alpha, num_nodes);
    return Rcpp::List::create(
        Rcpp::Named("log_density") = r.log_density,
        Rcpp::Named("log_Z")       = r.log_Z,
        Rcpp::Named("status")      = r.status
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


// Thin GGM SD chain driver, used by SBC and PD diagnostic scripts that
// want PIP, K samples, and PD-revert counts without going through the
// bgm() wrapper. The production bgm() path with
// prior_factorization = "hierarchical" hits the same code (the L-space SD
// between-step + the configurable within-K sampler).
//
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
    int    sample_thin = 0,      // 0 = don't save K samples; k>0 = save every k-th
    bool   edge_selection = true, // false = skip SD between-step, fixed graph
    std::string within_k_mode = "am",  // "am" = Roverato+RW, "nuts" = unconstrained NUTS on active theta
    int    nuts_max_depth = 10,
    double nuts_target_accept = 0.8,
    Rcpp::Nullable<Rcpp::IntegerMatrix> user_edges = R_NilValue
) {
    int p = observations.n_cols;
    arma::mat inclusion_probability(p, p, arma::fill::value(inclusion_prob));
    arma::imat initial_edges(p, p, arma::fill::zeros);
    for (int i = 0; i < p; ++i) initial_edges(i, i) = 1;
    // For fixed-graph (no edge-selection) mode the caller can override the
    // default via `user_edges` (a p×p {0,1} matrix). If omitted, default to
    // full graph when edge_selection = false, or identity when true.
    if (user_edges.isNotNull()) {
        Rcpp::IntegerMatrix ue(user_edges);
        if (ue.nrow() != p || ue.ncol() != p) {
            Rcpp::stop("ggm_sd_smoke_cpp: user_edges must be %d × %d", p, p);
        }
        for (int i = 0; i < p; ++i)
            for (int j = 0; j < p; ++j)
                initial_edges(i, j) = ue(i, j);
    } else if (!edge_selection) {
        initial_edges.ones();
    }

    auto slab = std::make_unique<NormalPrior>(interaction_scale);
    auto diag = std::make_unique<GammaScalePrior>(diagonal_shape, diagonal_rate);

    // Under prior_only the chain must target π(K, Γ) with no likelihood. The
    // AM within-K path checks `prior_only_` at the per-step MH ratio, but the
    // NUTS gradient engine bakes `n_` and `suf_stat_` into the log-posterior
    // without that flag. Build the model with n = 0, S = 0 so the prior is
    // the engine's actual target (mirrors the trick sample_ggm_prior uses).
    std::unique_ptr<GGMModel> model_ptr;
    if (prior_only) {
        arma::mat suf_stat_zero(p, p, arma::fill::zeros);
        model_ptr = std::make_unique<GGMModel>(
            /*n=*/0, suf_stat_zero, inclusion_probability,
            initial_edges, edge_selection,
            std::move(slab), std::move(diag));
    } else {
        model_ptr = std::make_unique<GGMModel>(
            observations, inclusion_probability, initial_edges,
            edge_selection, std::move(slab), std::move(diag),
            /*na_impute=*/false);
    }
    GGMModel& model = *model_ptr;
    model.set_seed(seed);
    model.set_determinant_tilt(delta);
    model.set_graph_prior_spec(GraphPriorSpec::Hierarchical);
    model.set_prior_only(prior_only);
    model.reset_pd_revert_count();

    // NUTS-within setup (only when within_k_mode == "nuts").
    const bool use_nuts = (within_k_mode == "nuts");
    if (!use_nuts && within_k_mode != "am") {
        Rcpp::stop("ggm_sd_smoke_cpp: within_k_mode must be \"am\" or \"nuts\" "
                   "(got: %s)", within_k_mode.c_str());
    }
    double nuts_step_size = 0.1;
    std::unique_ptr<DualAveraging> nuts_da;
    if (use_nuts) {
        SafeRNG& rng = model.get_rng();
        arma::vec theta0 = model.get_vectorized_parameters();
        auto joint_fn0 = [&model](const arma::vec& th)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(th);
        };
        auto grad_fn0 = [&model](const arma::vec& th) -> arma::vec {
            return model.logp_and_gradient(th).second;
        };
        nuts_step_size = heuristic_initial_step_size(
            theta0, grad_fn0, joint_fn0, rng, nuts_target_accept);
        nuts_da = std::make_unique<DualAveraging>(nuts_step_size);
    }

    // PIP accumulator + edge-count trajectory.
    arma::mat pip_counts(p, p, arma::fill::zeros);
    arma::ivec n_edges_path(n_sweeps, arma::fill::zeros);

    // Optional K-sample accumulators (for SBC). Off-diagonals stored
    // upper-triangle row-major to match the rest of bgms.
    const int n_pairs   = p * (p - 1) / 2;
    const int n_samples = (sample_thin > 0) ? (n_sweeps / sample_thin) : 0;
    arma::mat K_offdiag_samples(n_samples, n_pairs, arma::fill::zeros);
    arma::mat K_diag_samples   (n_samples, p,       arma::fill::zeros);

    int n_total = n_warmup + n_sweeps;
    int sample_idx = 0;
    for (int s = 0; s < n_total; ++s) {
        model.prepare_iteration();
        if (include_within_k) {
            if (use_nuts) {
                arma::vec theta  = model.get_vectorized_parameters();
                arma::vec inv_m  = arma::ones<arma::vec>(theta.n_elem);
                auto joint_fn = [&model](const arma::vec& th)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient(th);
                };
                StepResult res = nuts_step(
                    theta, nuts_step_size, joint_fn, inv_m,
                    model.get_rng(), nuts_max_depth);
                model.set_vectorized_parameters(res.state);
                if (s < n_warmup) {
                    nuts_da->update(res.accept_prob, nuts_target_accept);
                    nuts_step_size = nuts_da->current();
                } else if (s == n_warmup) {
                    nuts_step_size = nuts_da->averaged();
                }
            } else {
                model.do_one_metropolis_step(s);
            }
        }
        if (edge_selection) {
            model.update_edge_indicators();
        }
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
            const int post_idx = s - n_warmup;
            if (sample_thin > 0 && (post_idx + 1) % sample_thin == 0 &&
                sample_idx < n_samples) {
                const arma::mat& K = model.get_precision_matrix();
                int off_idx = 0;
                for (int ii = 0; ii < p; ++ii) {
                    K_diag_samples(sample_idx, ii) = K(ii, ii);
                    for (int jj = ii + 1; jj < p; ++jj) {
                        // Store K_ij * γ_ij so inactive edges are exact 0 (not
                        // ~1e-18 floating-point noise from the Roverato DEL
                        // cancellation), giving a clean spike-slab mixture
                        // sample for downstream SBC analysis.
                        K_offdiag_samples(sample_idx, off_idx) =
                            (E(ii, jj) == 1) ? K(ii, jj) : 0.0;
                        ++off_idx;
                    }
                }
                ++sample_idx;
            }
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
        Rcpp::Named("pip")               = pip,
        Rcpp::Named("final_edges")       = model.get_edge_indicators(),
        Rcpp::Named("final_K")           = model.get_precision_matrix(),
        Rcpp::Named("n_edges_path")      = n_edges_path,
        Rcpp::Named("n_pd_reverts")      = n_pd_reverts,
        Rcpp::Named("K_offdiag_samples") = K_offdiag_samples,
        Rcpp::Named("K_diag_samples")    = K_diag_samples,
        Rcpp::Named("within_k_mode")     = within_k_mode,
        Rcpp::Named("nuts_step_size")    = nuts_step_size
    );
}
