#pragma once

#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>

#include "models/base_model.h"
#include "mcmc/execution/chain_result.h"
#include "priors/edge_prior.h"
#include "utils/progress_manager.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/samplers/sampler_base.h"
#include "mcmc/execution/warmup_schedule.h"


/**
 * Create a sampler matching config.sampler_type
 *
 * @param config   Sampler configuration (type, step size, tree depth, etc.)
 * @param schedule Shared warmup schedule for adaptation staging
 * @return Owning pointer to a concrete SamplerBase subclass
 */
std::unique_ptr<SamplerBase> create_sampler(const SamplerConfig& config, WarmupSchedule& schedule);

/**
 * Run a single MCMC chain (warmup + sampling)
 *
 * @param chain_result  Output container for samples and diagnostics
 * @param model         Model to sample from (state is mutated in place)
 * @param edge_prior    Edge prior for edge-selection updates
 * @param config        Sampler configuration
 * @param chain_id      Chain identifier (0-based)
 * @param pm            Progress manager for user-facing status updates
 */
void run_mcmc_chain(
    ChainResult& chain_result,
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    int chain_id,
    ProgressManager& pm
);


/**
 * MCMCChainRunner - TBB worker for parallel chain execution
 *
 * Each chain gets its own model clone and edge prior. The worker dispatches
 * chains to threads via RcppParallel::parallelFor.
 */
struct MCMCChainRunner : public RcppParallel::Worker {
    std::vector<ChainResult>& results_;
    std::vector<std::unique_ptr<BaseModel>>& models_;
    std::vector<std::unique_ptr<BaseEdgePrior>>& edge_priors_;
    const SamplerConfig& config_;
    ProgressManager& pm_;

    MCMCChainRunner(
        std::vector<ChainResult>& results,
        std::vector<std::unique_ptr<BaseModel>>& models,
        std::vector<std::unique_ptr<BaseEdgePrior>>& edge_priors,
        const SamplerConfig& config,
        ProgressManager& pm
    ) :
        results_(results),
        models_(models),
        edge_priors_(edge_priors),
        config_(config),
        pm_(pm)
    {}

    void operator()(std::size_t begin, std::size_t end);
};


/**
 * Run multi-chain MCMC (parallel or sequential based on thread count)
 *
 * @param model       Prototype model (cloned per chain)
 * @param edge_prior  Prototype edge prior (cloned per chain)
 * @param config      Sampler configuration
 * @param no_chains   Number of chains to run
 * @param no_threads  Number of threads (1 = sequential)
 * @param pm          Progress manager for user-facing status updates
 * @return Vector of ChainResult, one per chain
 */
std::vector<ChainResult> run_mcmc_sampler(
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    int no_chains,
    int no_threads,
    ProgressManager& pm
);

/**
 * Convert chain results to an Rcpp::List for return to R
 *
 * @param results  Vector of completed chain results
 * @return Named list with samples, diagnostics, and metadata
 */
Rcpp::List convert_results_to_list(const std::vector<ChainResult>& results);
