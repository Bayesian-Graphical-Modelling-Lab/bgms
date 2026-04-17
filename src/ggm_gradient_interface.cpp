// Test interface for the GGM gradient engine and RATTLE projection.
//
// Exposes logp_and_gradient, forward_map, project_position,
// project_momentum, and constrained leapfrog to R for validation.
// Also exposes sample_ggm_prior() for sampling precision matrices
// from the GGM prior using NUTS.

#include <RcppArmadillo.h>
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "models/ggm/ggm_model.h"
#include "mcmc/algorithms/leapfrog.h"
#include "mcmc/algorithms/nuts.h"
#include "rng/rng_utils.h"
#include "priors/parameter_prior.h"
#include "priors/edge_prior.h"
#include "utils/common_helpers.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/execution/chain_runner.h"
#include "mcmc/execution/chain_result.h"
#include "utils/progress_manager.h"

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

    CauchyPrior ip(pairwise_scale);
    GammaScalePrior dp(1.0, 1.0);
    GGMGradientEngine engine;
    engine.rebuild(cs, static_cast<size_t>(n), suf_stat, ip, dp);

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
    CauchyPrior ip(1.0);
    GammaScalePrior dp(1.0, 1.0);
    GGMGradientEngine engine;
    engine.rebuild(cs, 100, dummy_S, ip, dp);

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
    const arma::imat& edge_indicators,
    Rcpp::Nullable<Rcpp::NumericVector> inv_mass_in = R_NilValue)
{
    size_t p = edge_indicators.n_rows;

    // Build a minimal GGMModel from sufficient statistics
    arma::mat suf_stat = arma::eye(p, p);
    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true,
                   std::make_unique<CauchyPrior>(2.5),
                   std::make_unique<GammaScalePrior>(1.0, 1.0));

    // Unpack x into the model's Cholesky factor
    model.set_full_position(x);

    // Project
    arma::vec x_proj = x;
    arma::vec inv_mass;
    if(inv_mass_in.isNotNull()) {
        inv_mass = Rcpp::as<arma::vec>(inv_mass_in);
    } else {
        inv_mass = arma::ones<arma::vec>(x.n_elem);
    }
    model.project_position(x_proj, inv_mass);

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
                   edge_indicators, true,
                   std::make_unique<CauchyPrior>(2.5),
                   std::make_unique<GammaScalePrior>(1.0, 1.0));

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

    CauchyPrior ip(pairwise_scale);
    GammaScalePrior dp(1.0, 1.0);
    GGMGradientEngine engine;
    engine.rebuild(cs, static_cast<size_t>(n), suf_stat, ip, dp);

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
    const arma::imat& edge_indicators,
    Rcpp::Nullable<Rcpp::NumericVector> inv_mass_in = R_NilValue)
{
    size_t p = edge_indicators.n_rows;

    arma::mat suf_stat = arma::eye(p, p);
    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true,
                   std::make_unique<CauchyPrior>(2.5),
                   std::make_unique<GammaScalePrior>(1.0, 1.0));

    arma::vec inv_mass;
    if(inv_mass_in.isNotNull()) {
        inv_mass = Rcpp::as<arma::vec>(inv_mass_in);
    } else {
        inv_mass = arma::ones<arma::vec>(r.n_elem);
    }

    arma::vec r_proj = r;
    model.project_momentum(r_proj, x, inv_mass);

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
    double pairwise_scale,
    Rcpp::Nullable<Rcpp::NumericVector> inv_mass_in = R_NilValue)
{
    size_t p = edge_indicators.n_rows;

    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true,
                   std::make_unique<CauchyPrior>(pairwise_scale),
                   std::make_unique<GammaScalePrior>(1.0, 1.0));

    // Joint function: logp_and_gradient_full
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

    // Projection callbacks using mass-weighted overloads
    ProjectPositionFn proj_pos = [&model, &inv_mass](arma::vec& x) {
        model.project_position(x, inv_mass);
    };
    ProjectMomentumFn proj_mom = [&model, &inv_mass](arma::vec& r, const arma::vec& x) {
        model.project_momentum(r, x, inv_mass);
    };

    // Run n_steps constrained leapfrog steps
    arma::vec x = x0;
    arma::vec r = r0;
    double logp0 = memo.cached_log_post(x);

    for (int s = 0; s < n_steps; ++s) {
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


// [[Rcpp::export]]
Rcpp::List ggm_test_leapfrog_constrained_checked(
    const arma::vec& x0,
    const arma::vec& r0,
    double step_size,
    int n_steps,
    const arma::mat& suf_stat,
    int n,
    const arma::imat& edge_indicators,
    double pairwise_scale,
    double reverse_check_tol = 0.5,
    Rcpp::Nullable<Rcpp::NumericVector> inv_mass_in = R_NilValue)
{
    size_t p = edge_indicators.n_rows;

    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(static_cast<int>(p), suf_stat, inc_prob,
                   edge_indicators, true,
                   std::make_unique<CauchyPrior>(pairwise_scale),
                   std::make_unique<GammaScalePrior>(1.0, 1.0));

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
    int non_reversible_count = 0;

    for (int s = 0; s < n_steps; ++s) {
        auto result = leapfrog_constrained_checked(
            x, r, step_size, memo, inv_mass, proj_pos, proj_mom,
            reverse_check_tol
        );
        x = std::move(result.theta);
        r = std::move(result.r);
        if (!result.reversible) {
            non_reversible_count++;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("x") = Rcpp::wrap(x),
        Rcpp::Named("r") = Rcpp::wrap(r),
        Rcpp::Named("non_reversible_count") = non_reversible_count
    );
}


// -----------------------------------------------------------------------------
// sample_ggm_prior: Sample precision matrices from the GGM prior using NUTS
// -----------------------------------------------------------------------------
// Uses the Cholesky parameterization with NUTS. By setting n=0 and S=0,
// the likelihood vanishes and the sampler targets the prior:
//   K_ij | graph ~ Cauchy(0, scale) or Normal(0, scale)  (included edges)
//   K_ij = 0                                             (excluded edges)
//   K_ii ~ Gamma(1, 1)                                   (diagonal)
//
// edge_indicators: p x p integer matrix with 1 = edge included, 0 = excluded.
//   Defaults to all-ones (full graph). For edge selection SBC, pass the
//   graph drawn from the edge prior so K is sampled conditional on that graph.
//
// The prior is the product of the element-wise priors plus the
// Jacobian from theta -> K induced by the Cholesky parameterization.
//
// Uses the same NUTS infrastructure (windowed warmup, dual-averaging step-size
// adaptation, Welford diagonal mass-matrix adaptation) as sample_ggm and bgm.
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List sample_ggm_prior(
    int p,
    int n_samples,
    int n_warmup = 1000,
    double pairwise_scale = 2.5,
    const std::string& interaction_prior_type = "cauchy",
    double step_size = 0.1,
    int max_depth = 10,
    int seed = 1,
    bool verbose = true,
    Rcpp::Nullable<Rcpp::IntegerMatrix> edge_indicators_nullable = R_NilValue)
{
    // Build edge indicators (default: full graph, no constraints)
    arma::imat edge_indicators;
    if(edge_indicators_nullable.isNotNull()) {
        edge_indicators = Rcpp::as<arma::imat>(
            Rcpp::IntegerMatrix(edge_indicators_nullable.get()));
    } else {
        edge_indicators.ones(p, p);
    }

    auto ip = create_parameter_prior(interaction_prior_type, pairwise_scale);
    auto dp = std::make_unique<GammaScalePrior>(1.0, 1.0);

    // Create model with n=0, S=0 so likelihood is flat (prior-only).
    // edge_selection=false ensures the graph stays fixed throughout:
    // no update_edge_indicators() calls, required for RATTLE correctness
    // when the graph is sparse.
    arma::mat suf_stat(p, p, arma::fill::zeros);
    arma::mat inc_prob(p, p, arma::fill::value(0.5));
    GGMModel model(0, suf_stat, inc_prob, edge_indicators,
                   false, std::move(ip), std::move(dp));

    // Configure NUTS with the standard windowed warmup (dual averaging +
    // diagonal Welford mass-matrix adaptation). This is the same path used
    // by sample_ggm / bgm, so prior draws match the sampler used for the
    // full posterior.
    SamplerConfig config;
    config.sampler_type = "nuts";
    config.no_iter = n_samples;
    config.no_warmup = n_warmup;
    config.edge_selection = false;
    config.seed = seed;
    config.target_acceptance = 0.8;
    config.max_tree_depth = max_depth;
    config.initial_step_size = step_size;
    config.na_impute = false;

    // Edge prior is required by the API but is never consulted because
    // edge_selection=false. Use a cheap Bernoulli placeholder.
    auto edge_prior_obj = create_edge_prior(
        Bernoulli, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

    // Silent progress manager (no R console spam). progress_type=0 => none.
    ProgressManager pm(1, n_samples, n_warmup, 0, 0, false, R_NilValue);

    // Run a single chain.
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, /*no_chains=*/1, /*no_threads=*/1, pm);
    pm.finish();

    if(results.empty() || results[0].error) {
        std::string msg = results.empty() ? "empty result" : results[0].error_msg;
        Rcpp::stop("sample_ggm_prior: chain failed (%s)", msg.c_str());
    }

    // Samples are stored as the upper triangle of K (p(p+1)/2 elements per
    // iteration, column-wise iteration index): see GGMModel::extract_upper_triangle.
    // Layout: for i in 0..p-1, for j in i..p-1, entry K(i,j).
    const arma::mat& samples = results[0].samples;  // dim x no_iter

    int n_edges = p * (p - 1) / 2;
    arma::mat K_offdiag_samples(n_samples, n_edges);
    arma::mat K_diag_samples(n_samples, p);

    for(int s = 0; s < n_samples; ++s) {
        const arma::vec& upper = samples.col(s);
        size_t e = 0;
        int off_idx = 0;
        for(int i = 0; i < p; ++i) {
            for(int j = i; j < p; ++j) {
                if(i == j) {
                    K_diag_samples(s, i) = upper(e);
                } else {
                    K_offdiag_samples(s, off_idx++) = upper(e);
                }
                ++e;
            }
        }
    }

    // Build column names
    Rcpp::CharacterVector offdiag_names(n_edges);
    int idx = 0;
    for(int i = 0; i < p - 1; ++i) {
        for(int j = i + 1; j < p; ++j) {
            offdiag_names[idx++] = "K_" + std::to_string(i + 1) + "_" +
                                   std::to_string(j + 1);
        }
    }

    Rcpp::CharacterVector diag_names(p);
    for(int i = 0; i < p; ++i) {
        diag_names[i] = "K_" + std::to_string(i + 1) + "_" +
                        std::to_string(i + 1);
    }

    return Rcpp::List::create(
        Rcpp::Named("K_offdiag") = K_offdiag_samples,
        Rcpp::Named("K_diag") = K_diag_samples,
        Rcpp::Named("offdiag_names") = offdiag_names,
        Rcpp::Named("diag_names") = diag_names,
        Rcpp::Named("step_size") = config.initial_step_size,
        Rcpp::Named("edge_indicators") = Rcpp::wrap(edge_indicators)
    );
}
