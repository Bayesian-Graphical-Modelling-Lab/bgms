#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>

#include "ggm_model.h"
#include "utils/progress_manager.h"
#include "chainResultNew.h"
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
    const int progress_type
) {

    // Create model from R input
    GaussianVariables model = createGaussianVariablesFromR(
        inputFromR, prior_inclusion_prob, initial_edge_indicators, edge_selection);

    // Configure sampler - GGM only supports MH
    SamplerConfig config;
    config.sampler_type = "mh";
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.seed = seed;
    // Edge selection starts at no_warmup/2 by default (handled by get_edge_selection_start())

    // Set up progress manager
    ProgressManager pm(no_chains, no_iter, no_warmup, 50, progress_type);

    // Run MCMC using unified infrastructure
    std::vector<ChainResultNew> results = run_mcmc_sampler(
        model, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}