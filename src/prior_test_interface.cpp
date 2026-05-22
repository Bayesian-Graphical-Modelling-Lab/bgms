// Test interface for the polymorphic parameter prior classes.
//
// Exposes logp, grad, and scaled variants to R for unit testing.

#include <RcppArmadillo.h>
#include "priors/parameter_prior.h"
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "models/ggm/ggm_model.h"


// [[Rcpp::export]]
Rcpp::List test_parameter_prior(
    const std::string& type,
    double x,
    double scale = 1.0,
    double alpha = 0.5,
    double beta = 0.5,
    double scale_factor = 1.0
) {
    auto prior = create_parameter_prior(type, scale, alpha, beta);

    return Rcpp::List::create(
        Rcpp::Named("logp") = prior->logp(x),
        Rcpp::Named("grad") = prior->grad(x),
        Rcpp::Named("logp_scaled") = prior->logp(x, scale_factor),
        Rcpp::Named("grad_scaled") = prior->grad(x, scale_factor)
    );
}


// Test entry for the edgewise GraphicalGPrior. Builds a small prior over
// the supplied V_ij matrix and t, then evaluates logp/grad at (x, i, j).
// Used to confirm bit-equality with NormalPrior(scale = t * sqrt(V_ij(i,j))).
//
// [[Rcpp::export]]
Rcpp::List test_graphical_g_prior(
    double x,
    const arma::mat& V_ij,
    double t,
    int i,
    int j
) {
    GraphicalGPrior prior(&V_ij, &t);
    return Rcpp::List::create(
        Rcpp::Named("logp_edge")     = prior.logp(x, i, j),
        Rcpp::Named("grad_edge")     = prior.grad(x, i, j),
        Rcpp::Named("logp_no_edge")  = prior.logp(x),
        Rcpp::Named("grad_no_edge")  = prior.grad(x),
        // For comparison: equivalent NormalPrior at the resolved scale.
        Rcpp::Named("logp_ref")      =
            R::dnorm(x, 0.0, t * std::sqrt(V_ij(i, j)), true)
    );
}


// Test entry for the edgewise GraphicalGDiag (Gamma diagonal with rate 1/t).
//
// [[Rcpp::export]]
Rcpp::List test_graphical_g_diag(double x, double t) {
    GraphicalGDiag prior(&t);
    return Rcpp::List::create(
        Rcpp::Named("logp") = prior.logp(x),
        Rcpp::Named("grad") = prior.grad(x),
        // Reference Gamma(shape = 1, rate = 1/t)  ≡  Exp(1/t) density.
        Rcpp::Named("logp_ref") = R::dgamma(x, 1.0, t, true)
    );
}


// [[Rcpp::export]]
Rcpp::List test_scale_prior(
    const std::string& type,
    double x,
    double shape = 1.0,
    double rate = 1.0
) {
    auto prior = create_scale_prior(type, shape, rate);

    return Rcpp::List::create(
        Rcpp::Named("logp") = prior->logp(x),
        Rcpp::Named("grad") = prior->grad(x)
    );
}


// [[Rcpp::export]]
Rcpp::List ggm_test_logp_and_gradient_prior(
    const arma::vec& theta,
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    const std::string& interaction_prior_type = "cauchy",
    double interaction_scale = 1.0,
    double interaction_alpha = 0.5,
    double interaction_beta = 0.5,
    const std::string& diagonal_prior_type = "gamma",
    double diagonal_shape = 1.0,
    double diagonal_rate = 1.0
) {
    GraphConstraintStructure cs;
    cs.build(edge_indicators);

    auto ip = create_parameter_prior(
        interaction_prior_type, interaction_scale,
        interaction_alpha, interaction_beta);
    auto dp = create_scale_prior(diagonal_prior_type, diagonal_shape, diagonal_rate);

    GGMGradientEngine engine;
    engine.rebuild(cs, static_cast<size_t>(n), suf_stat, *ip, *dp);

    auto result = engine.logp_and_gradient(theta);

    return Rcpp::List::create(
        Rcpp::Named("value") = result.first,
        Rcpp::Named("gradient") = Rcpp::wrap(result.second)
    );
}


// End-to-end Graphical G-prior smoke chain.
//
// Builds a GGMModel from observations, enables the GG-prior with the chosen
// hyperprior, runs `n_iter` Metropolis sweeps, and returns the K trajectory,
// edge-indicator trajectory, the V_ij table, and the t (= √g) trajectory.
//
// Used for bit-equality + dynamics validation before the full R API is built.
// Once the R API is in (bgm() → graphical_g_prior()) this entry point can be
// dropped or kept under .StatsVault. Not part of the public API.
//
// [[Rcpp::export]]
Rcpp::List ggm_gg_prior_smoke_cpp(
    const arma::mat& observations,
    double inclusion_prob,
    int    n_iter,
    int    n_warmup,
    int    seed,
    int    g_hyperprior,
    double g_init      = 1.0,
    double tcch_a      = 1.0,
    double tcch_b      = 1.0,
    double delta       = 0.0
) {
    // g_hyperprior: 0 = Fixed, 1 = ConjugateGamma
    const int p = observations.n_cols;
    arma::mat inclusion_probability(p, p, arma::fill::value(inclusion_prob));
    arma::imat initial_edge_indicators(p, p, arma::fill::zeros);
    for (int i = 0; i < p; ++i) initial_edge_indicators(i, i) = 1;

    // Placeholder priors — enable_gg_prior() will replace them with the
    // GraphicalG-bound versions.
    auto slab = std::make_unique<NormalPrior>(1.0);
    auto diag = std::make_unique<GammaScalePrior>(1.0, 1.0);

    GGMModel model(observations,
                   inclusion_probability,
                   initial_edge_indicators,
                   /*edge_selection=*/true,
                   std::move(slab),
                   std::move(diag),
                   /*na_impute=*/false);
    model.set_seed(seed);
    model.set_determinant_tilt(delta);
    model.enable_gg_prior(
        static_cast<GGMModel::GGHyperprior>(g_hyperprior),
        g_init, tcch_a, tcch_b);

    const int n_edges_off = p * (p - 1) / 2;

    arma::cube K_samples(p, p, n_iter, arma::fill::zeros);
    arma::imat ind_samples(n_edges_off, n_iter, arma::fill::zeros);
    arma::vec t_trace(n_iter, arma::fill::zeros);
    arma::ivec n_edges_per_iter(n_iter, arma::fill::zeros);

    // Manual Metropolis loop. No NUTS path here — we keep this entry
    // minimal so it directly exercises the MH-path slab/diag substitution.
    for (int it = 0; it < n_warmup; ++it) {
        model.prepare_iteration();
        model.set_edge_selection_active(true);
        model.update_edge_indicators();
        model.do_one_metropolis_step(it);
    }
    for (int it = 0; it < n_iter; ++it) {
        model.prepare_iteration();
        model.set_edge_selection_active(true);
        model.update_edge_indicators();
        model.do_one_metropolis_step(it);

        // Snapshot.
        K_samples.slice(it) = model.get_precision_matrix();
        const arma::imat& E = model.get_edge_indicators();
        int e = 0;
        for (int i = 0; i < p - 1; ++i) {
            for (int j = i + 1; j < p; ++j) {
                ind_samples(e, it) = E(i, j);
                ++e;
            }
        }
        n_edges_per_iter(it) = arma::accu(E) / 2;
        t_trace(it) = model.gg_current_t();
    }

    return Rcpp::List::create(
        Rcpp::Named("K_samples")        = Rcpp::wrap(K_samples),
        Rcpp::Named("indicator_samples")= Rcpp::wrap(ind_samples),
        Rcpp::Named("t_trace")          = Rcpp::wrap(t_trace),
        Rcpp::Named("n_edges_per_iter") = Rcpp::wrap(n_edges_per_iter),
        Rcpp::Named("V_ij")             = Rcpp::wrap(model.gg_V_ij()),
        Rcpp::Named("g_final")          = model.gg_current_g(),
        Rcpp::Named("t_final")          = model.gg_current_t()
    );
}

