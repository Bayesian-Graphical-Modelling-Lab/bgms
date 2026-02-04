/**
 * sample_omrf.cpp - R interface for OMRF model sampling
 *
 * Uses the unified MCMC runner infrastructure to sample from OMRF models.
 * Supports MH, NUTS, and HMC samplers with optional edge selection.
 */
#include <vector>
#include <memory>
#include <RcppArmadillo.h>

#include "omrf_model.h"
#include "utils/progress_manager.h"
#include "chainResultNew.h"
#include "mcmc/mcmc_runner.h"
#include "mcmc/sampler_config.h"

/**
 * R-exported function to sample from an OMRF model
 *
 * @param inputFromR          List with model specification
 * @param prior_inclusion_prob Prior inclusion probabilities (p × p matrix)
 * @param initial_edge_indicators Initial edge indicators (p × p integer matrix)
 * @param no_iter             Number of post-warmup iterations
 * @param no_warmup           Number of warmup iterations
 * @param no_chains           Number of parallel chains
 * @param edge_selection      Whether to do edge selection (spike-and-slab)
 * @param sampler_type        "mh", "nuts", or "hmc"
 * @param seed                Random seed
 * @param no_threads          Number of threads for parallel execution
 * @param progress_type       Progress bar type
 * @param target_acceptance   Target acceptance rate for NUTS/HMC (default: 0.8)
 * @param max_tree_depth      Maximum tree depth for NUTS (default: 10)
 * @param num_leapfrogs       Number of leapfrog steps for HMC (default: 10)
 * @param edge_selection_start Iteration to start edge selection (-1 = no_warmup/2)
 *
 * @return List with per-chain results including samples and diagnostics
 */
// [[Rcpp::export]]
Rcpp::List sample_omrf(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const int no_iter,
    const int no_warmup,
    const int no_chains,
    const bool edge_selection,
    const std::string& sampler_type,
    const int seed,
    const int no_threads,
    const int progress_type,
    const double target_acceptance = 0.8,
    const int max_tree_depth = 10,
    const int num_leapfrogs = 10,
    const int edge_selection_start = -1
) {
    // Create model from R input
    OMRFModel model = createOMRFModelFromR(
        inputFromR, prior_inclusion_prob, initial_edge_indicators, edge_selection);

    // Configure sampler
    SamplerConfig config;
    config.sampler_type = sampler_type;
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.edge_selection_start = edge_selection_start;  // -1 means use default (no_warmup/2)
    config.seed = seed;
    config.target_acceptance = target_acceptance;
    config.max_tree_depth = max_tree_depth;
    config.num_leapfrogs = num_leapfrogs;

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
