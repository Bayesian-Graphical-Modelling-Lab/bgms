#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>

#include "models/ggm/ggm_model.h"
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/edge_prior.h"
#include "mcmc/chain_result.h"
#include "mcmc/mcmc_runner.h"
#include "mcmc/sampler_config.h"

// [[Rcpp::export]]
Rcpp::List sample_ggm(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const int no_iter,
    const int no_warmup,
    const int no_chains,
    const bool edge_selection,
    const int seed,
    const int no_threads,
    const int progress_type,
    const std::string& edge_prior = "Bernoulli",
    const double beta_bernoulli_alpha = 1.0,
    const double beta_bernoulli_beta = 1.0,
    const double beta_bernoulli_alpha_between = 1.0,
    const double beta_bernoulli_beta_between = 1.0,
    const double dirichlet_alpha = 1.0,
    const double lambda = 1.0
) {

    // Create model from R input
    GGMModel model = createGGMModelFromR(
        inputFromR, prior_inclusion_prob, initial_edge_indicators, edge_selection);

    // Configure sampler - GGM only supports MH
    SamplerConfig config;
    config.sampler_type = "mh";
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.seed = seed;

    // Set up progress manager
    ProgressManager pm(no_chains, no_iter, no_warmup, 50, progress_type);

    // Create edge prior
    EdgePrior edge_prior_enum = edge_prior_from_string(edge_prior);
    auto edge_prior_obj = create_edge_prior(
        edge_prior_enum,
        beta_bernoulli_alpha, beta_bernoulli_beta,
        beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
        dirichlet_alpha, lambda
    );

    // Run MCMC using unified infrastructure
    std::vector<ChainResultNew> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}