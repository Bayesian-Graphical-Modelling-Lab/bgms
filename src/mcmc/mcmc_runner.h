#pragma once

#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>

#include "models/base_model.h"
#include "mcmc/chain_result.h"
#include "priors/edge_prior.h"
#include "utils/progress_manager.h"
#include "mcmc/sampler_config.h"
#include "mcmc/base_sampler.h"
#include "mcmc/mcmc_utils.h"
#include "mcmc/mcmc_adaptation.h"


/// Create a sampler matching config.sampler_type.
std::unique_ptr<BaseSampler> create_sampler(const SamplerConfig& config, WarmupSchedule& schedule);

/// Run a single MCMC chain (warmup + sampling) writing into chain_result.
void run_mcmc_chain(
    ChainResultNew& chain_result,
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    int chain_id,
    ProgressManager& pm
);


/// Worker struct for TBB parallel chain execution.
struct MCMCChainRunner : public RcppParallel::Worker {
    std::vector<ChainResultNew>& results_;
    std::vector<std::unique_ptr<BaseModel>>& models_;
    std::vector<std::unique_ptr<BaseEdgePrior>>& edge_priors_;
    const SamplerConfig& config_;
    ProgressManager& pm_;

    MCMCChainRunner(
        std::vector<ChainResultNew>& results,
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


/// Run multi-chain MCMC (parallel or sequential based on no_threads).
std::vector<ChainResultNew> run_mcmc_sampler(
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    int no_chains,
    int no_threads,
    ProgressManager& pm
);

/// Convert chain results to Rcpp::List for return to R.
Rcpp::List convert_results_to_list(const std::vector<ChainResultNew>& results);
