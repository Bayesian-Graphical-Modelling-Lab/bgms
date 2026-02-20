/**
 * sample_omrf.cpp - R interface for OMRF model sampling
 *
 * Uses the unified MCMC runner infrastructure to sample from OMRF models.
 * Supports MH, NUTS, and HMC samplers with optional edge selection.
 */
#include <vector>
#include <memory>
#include <RcppArmadillo.h>

#include "models/omrf/omrf_model.h"
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/edge_prior.h"
#include "mcmc/chain_result.h"
#include "mcmc/mcmc_runner.h"
#include "mcmc/sampler_config.h"

/**
 * R-exported function to sample from an OMRF model
 *
 * @param inputFromR          List with model specification
 * @param prior_inclusion_prob Prior inclusion probabilities (p x p matrix)
 * @param initial_edge_indicators Initial edge indicators (p x p integer matrix)
 * @param no_iter             Number of post-warmup iterations
 * @param no_warmup           Number of warmup iterations
 * @param no_chains           Number of parallel chains
 * @param edge_selection      Whether to do edge selection (spike-and-slab)
 * @param sampler_type        "mh", "nuts", or "hmc"
 * @param seed                Random seed
 * @param no_threads          Number of threads for parallel execution
 * @param progress_type       Progress bar type
 * @param edge_prior          Edge prior type: "Bernoulli", "Beta-Bernoulli", "Stochastic-Block"
 * @param na_impute           Whether to impute missing data
 * @param missing_index       Matrix of missing data indices (n_missing x 2, 0-based)
 * @param beta_bernoulli_alpha     Beta-Bernoulli alpha hyperparameter
 * @param beta_bernoulli_beta      Beta-Bernoulli beta hyperparameter
 * @param beta_bernoulli_alpha_between SBM between-cluster alpha
 * @param beta_bernoulli_beta_between  SBM between-cluster beta
 * @param dirichlet_alpha     Dirichlet alpha for SBM
 * @param lambda              Lambda for SBM
 * @param target_acceptance   Target acceptance rate for NUTS/HMC (default: 0.8)
 * @param max_tree_depth      Maximum tree depth for NUTS (default: 10)
 * @param num_leapfrogs       Number of leapfrog steps for HMC (default: 10)
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
    const std::string& edge_prior = "Bernoulli",
    const bool na_impute = false,
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_nullable = R_NilValue,
    const double beta_bernoulli_alpha = 1.0,
    const double beta_bernoulli_beta = 1.0,
    const double beta_bernoulli_alpha_between = 1.0,
    const double beta_bernoulli_beta_between = 1.0,
    const double dirichlet_alpha = 1.0,
    const double lambda = 1.0,
    const double target_acceptance = 0.8,
    const int max_tree_depth = 10,
    const int num_leapfrogs = 10,
    const Rcpp::Nullable<Rcpp::NumericMatrix> pairwise_scaling_factors_nullable = R_NilValue
) {
    // Create model from R input
    OMRFModel model = createOMRFModelFromR(
        inputFromR, prior_inclusion_prob, initial_edge_indicators, edge_selection);

    // Set pairwise scaling factors (if provided)
    if (pairwise_scaling_factors_nullable.isNotNull()) {
        arma::mat sf = Rcpp::as<arma::mat>(
            Rcpp::NumericMatrix(pairwise_scaling_factors_nullable.get()));
        model.set_pairwise_scaling_factors(sf);
    }

    // Set up missing data imputation
    if (na_impute && missing_index_nullable.isNotNull()) {
        arma::imat missing_index = Rcpp::as<arma::imat>(
            Rcpp::IntegerMatrix(missing_index_nullable.get()));
        model.set_missing_data(missing_index);
    }

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
    config.target_acceptance = target_acceptance;
    config.max_tree_depth = max_tree_depth;
    config.num_leapfrogs = num_leapfrogs;
    config.na_impute = na_impute;

    // Set up progress manager
    ProgressManager pm(no_chains, no_iter, no_warmup, 50, progress_type);

    // Run MCMC using unified infrastructure
    std::vector<ChainResultNew> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}
