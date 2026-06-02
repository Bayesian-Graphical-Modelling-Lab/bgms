// omrf_residual_test_interface.cpp - test-only interface
//
// Exposes a minimal harness to verify the OMRF residual-matrix invariant
//   residual_matrix_ == 2 * observations_double_ * pairwise_effects_
// is preserved across a warmup run of tune_proposal_sd(). Used by
// tests/testthat to characterize the stage-3b proposal-SD tuner.
#include <RcppArmadillo.h>

#include "models/omrf/omrf_model.h"
#include "priors/parameter_prior.h"
#include "mcmc/execution/warmup_schedule.h"

// Build an all-ordinal OMRF model on the given data, set the pairwise
// effects (which rebuilds residual_matrix_ exactly), then run `warmup`
// iterations of tune_proposal_sd(). Returns the residual matrix the model
// holds afterwards alongside the value the invariant requires, so the
// caller can assert they match.
//
// [[Rcpp::export]]
Rcpp::List test_omrf_residual_invariant(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& pairwise,
    const int warmup,
    const int seed,
    const double target_accept = 0.44,
    const bool enable_selection = true,
    const bool learn_sd = true
) {
    const int p = static_cast<int>(observations.n_cols);

    // All-ordinal so observations_double_ == observations (no BC centering),
    // making the invariant residual == 2 * X * pairwise hold exactly.
    arma::uvec is_ordinal = arma::ones<arma::uvec>(p);
    arma::ivec baseline   = arma::zeros<arma::ivec>(p);
    arma::mat  incl_prob  = 0.5 * arma::ones<arma::mat>(p, p);
    arma::imat edges      = arma::ones<arma::imat>(p, p);
    edges.diag().zeros();

    auto interaction_prior = create_parameter_prior("cauchy", 2.5, NA_REAL, NA_REAL);
    auto threshold_prior   = create_parameter_prior("beta-prime", 1.0, 0.5, 0.5);

    OMRFModel model(
        observations, num_categories, incl_prob, edges,
        is_ordinal, baseline,
        std::move(interaction_prior), std::move(threshold_prior),
        /*edge_selection=*/false);

    model.set_seed(seed);
    model.set_metropolis_target_accept(target_accept);
    model.set_pairwise_effects(pairwise);  // residual now exactly 2*X*pairwise

    WarmupSchedule schedule(warmup, enable_selection, learn_sd);
    model.init_metropolis_adaptation(schedule);

    for (int iter = 0; iter < warmup; ++iter) {
        model.tune_proposal_sd(iter, schedule);
    }

    const arma::mat residual = model.get_residual_matrix();
    const arma::mat pw       = model.get_pairwise_effects();
    const arma::mat obs_d    = arma::conv_to<arma::mat>::from(observations);
    const arma::mat expected = 2.0 * obs_d * pw;

    return Rcpp::List::create(
        Rcpp::Named("residual")            = residual,
        Rcpp::Named("expected")            = expected,
        Rcpp::Named("pairwise_after")      = pw,
        Rcpp::Named("pairwise_before")     = pairwise,
        Rcpp::Named("max_abs_discrepancy") = arma::abs(residual - expected).max(),
        Rcpp::Named("any_move_accepted")   = (arma::abs(pw - pairwise).max() > 0.0)
    );
}
