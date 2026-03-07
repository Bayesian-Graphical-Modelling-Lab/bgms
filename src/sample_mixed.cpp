// sample_mixed.cpp - R interface for Mixed MRF model sampling
//
// Uses the unified MCMC runner infrastructure to sample from models with
// both discrete (ordinal / Blume-Capel) and continuous variables.
// Supports MH and hybrid-nuts (NUTS for unconstrained block + MH for Kyy)
// samplers, with optional edge selection.
#include <vector>
#include <memory>
#include <RcppArmadillo.h>

#include "models/mixed/mixed_mrf_model.h"
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/edge_prior.h"
#include "mcmc/execution/chain_result.h"
#include "mcmc/execution/chain_runner.h"
#include "mcmc/execution/sampler_config.h"

// R-exported function to sample from a Mixed MRF model.
//
// @param inputFromR              List with model specification:
//                                  discrete_observations (integer matrix n x p),
//                                  continuous_observations (numeric matrix n x q),
//                                  num_categories (integer vector, length p),
//                                  is_ordinal_variable (integer vector, length p),
//                                  baseline_category (integer vector, length p),
//                                  main_alpha, main_beta, pairwise_scale (doubles),
//                                  pseudolikelihood (string: "conditional" or "marginal")
// @param prior_inclusion_prob    Prior inclusion probabilities ((p+q) x (p+q) matrix)
// @param initial_edge_indicators Initial edge indicators ((p+q) x (p+q) integer matrix)
// @param no_iter                 Number of post-warmup iterations
// @param no_warmup               Number of warmup iterations
// @param no_chains               Number of parallel chains
// @param edge_selection          Whether to do edge selection (spike-and-slab)
// @param seed                    Random seed
// @param no_threads              Number of threads for parallel execution
// @param progress_type           Progress bar type
// @param edge_prior              Edge prior type
// @param beta_bernoulli_alpha         Beta-Bernoulli alpha hyperparameter
// @param beta_bernoulli_beta          Beta-Bernoulli beta hyperparameter
// @param beta_bernoulli_alpha_between SBM between-cluster alpha
// @param beta_bernoulli_beta_between  SBM between-cluster beta
// @param dirichlet_alpha         Dirichlet alpha for SBM
// @param lambda                  Lambda for SBM
// @param sampler_type            Sampler type string ("mh", "hybrid-nuts", etc.)
// @param target_acceptance       Target acceptance rate for gradient-based samplers
// @param max_tree_depth          Maximum tree depth for NUTS
// @param num_leapfrogs           Number of leapfrog steps for HMC
//
// @return List with per-chain results including samples and diagnostics
// [[Rcpp::export]]
Rcpp::List sample_mixed_mrf(
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
    const double lambda = 1.0,
    const std::string& sampler_type = "mh",
    const double target_acceptance = 0.80,
    const int max_tree_depth = 10,
    const int num_leapfrogs = 100
) {
    // Extract model inputs from R list
    arma::imat discrete_obs = Rcpp::as<arma::imat>(inputFromR["discrete_observations"]);
    arma::mat continuous_obs = Rcpp::as<arma::mat>(inputFromR["continuous_observations"]);
    arma::ivec num_categories = Rcpp::as<arma::ivec>(inputFromR["num_categories"]);
    arma::uvec is_ordinal = Rcpp::as<arma::uvec>(inputFromR["is_ordinal_variable"]);
    arma::ivec baseline_cat = Rcpp::as<arma::ivec>(inputFromR["baseline_category"]);
    double main_alpha = Rcpp::as<double>(inputFromR["main_alpha"]);
    double main_beta = Rcpp::as<double>(inputFromR["main_beta"]);
    double pairwise_scale = Rcpp::as<double>(inputFromR["pairwise_scale"]);
    std::string pseudolikelihood = Rcpp::as<std::string>(inputFromR["pseudolikelihood"]);

    // Create model
    MixedMRFModel model(
        discrete_obs, continuous_obs,
        num_categories, is_ordinal, baseline_cat,
        prior_inclusion_prob, initial_edge_indicators,
        edge_selection, pseudolikelihood,
        main_alpha, main_beta, pairwise_scale,
        seed
    );

    // Create edge prior
    EdgePrior edge_prior_enum = edge_prior_from_string(edge_prior);
    auto edge_prior_obj = create_edge_prior(
        edge_prior_enum,
        beta_bernoulli_alpha, beta_bernoulli_beta,
        beta_bernoulli_alpha_between, beta_bernoulli_beta_between,
        dirichlet_alpha, lambda
    );

    // Configure sampler
    SamplerConfig config;
    config.sampler_type = sampler_type;
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.seed = seed;
    config.na_impute = false;
    config.target_acceptance = target_acceptance;
    config.max_tree_depth = max_tree_depth;
    config.num_leapfrogs = num_leapfrogs;

    // Set up progress manager
    ProgressManager pm(no_chains, no_iter, no_warmup, 50, progress_type);

    // Run MCMC using unified infrastructure
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}
