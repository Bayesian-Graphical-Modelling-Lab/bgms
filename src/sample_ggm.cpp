#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>

#include "models/ggm/ggm_model.h"
#include "utils/progress_manager.h"
#include "utils/common_helpers.h"
#include "priors/edge_prior.h"
#include "mcmc/execution/chain_result.h"
#include "mcmc/execution/chain_runner.h"
#include "mcmc/execution/sampler_config.h"

// [[Rcpp::export]]
Rcpp::List sample_ggm(
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
    const double beta_bernoulli_alpha = 1.0,
    const double beta_bernoulli_beta = 1.0,
    const double beta_bernoulli_alpha_between = 1.0,
    const double beta_bernoulli_beta_between = 1.0,
    const double dirichlet_alpha = 1.0,
    const double lambda = 1.0,
    const double target_acceptance = 0.8,
    const int max_tree_depth = 10,
    const bool na_impute = false,
    const Rcpp::Nullable<Rcpp::IntegerMatrix> missing_index_nullable = R_NilValue
) {

    // Create model from R input
    double pairwise_scale = Rcpp::as<double>(inputFromR["pairwise_scale"]);

    GGMModel model = createGGMModelFromR(
        inputFromR, prior_inclusion_prob, initial_edge_indicators,
        edge_selection, pairwise_scale, na_impute);

    // Set up missing data imputation (same pattern as OMRF)
    if (na_impute && missing_index_nullable.isNotNull()) {
        arma::imat missing_index = Rcpp::as<arma::imat>(
            Rcpp::IntegerMatrix(missing_index_nullable.get()));
        model.set_missing_data(missing_index);
    }

    // Configure sampler
    SamplerConfig config;
    config.sampler_type = sampler_type;
    config.no_iter = no_iter;
    config.no_warmup = no_warmup;
    config.edge_selection = edge_selection;
    config.seed = seed;
    config.target_acceptance = target_acceptance;
    config.max_tree_depth = max_tree_depth;
    config.na_impute = na_impute;

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
    std::vector<ChainResult> results = run_mcmc_sampler(
        model, *edge_prior_obj, config, no_chains, no_threads, pm);

    // Convert to R list format
    Rcpp::List output = convert_results_to_list(results);

    pm.finish();

    return output;
}