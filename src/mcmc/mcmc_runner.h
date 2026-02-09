#pragma once

#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>

#include "../base_model.h"
#include "../chainResultNew.h"
#include "../priors/edge_prior.h"
#include "../utils/progress_manager.h"
#include "sampler_config.h"
#include "base_sampler.h"
#include "nuts_sampler.h"
#include "hmc_sampler.h"
#include "mh_sampler.h"
#include "mcmc_utils.h"


/**
 * Create a sampler based on configuration
 *
 * Factory function that returns the appropriate sampler type.
 *
 * @param config  Sampler configuration
 * @return Unique pointer to the created sampler
 */
inline std::unique_ptr<BaseSampler> create_sampler(const SamplerConfig& config) {
    if (config.sampler_type == "nuts") {
        return std::make_unique<NUTSSampler>(config, config.no_warmup);
    } else if (config.sampler_type == "hmc" || config.sampler_type == "hamiltonian-mc") {
        return std::make_unique<HMCSampler>(config);
    } else if (config.sampler_type == "mh" || config.sampler_type == "adaptive-metropolis") {
        return std::make_unique<MHSampler>(config);
    } else {
        Rcpp::stop("Unknown sampler_type: '%s'", config.sampler_type.c_str());
    }
}


/**
 * Run MCMC sampling for a single chain
 *
 * Supports MH, NUTS, and HMC samplers with optional edge selection.
 * Handles warmup adaptation and diagnostic collection.
 *
 * @param chain_result  Output storage for this chain
 * @param model         The model to sample from
 * @param config        Sampler configuration
 * @param chain_id      Chain identifier (0-based)
 * @param pm            Progress manager for user feedback
 */
inline void run_mcmc_chain(
    ChainResultNew& chain_result,
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int chain_id,
    ProgressManager& pm
) {
    chain_result.chain_id = chain_id + 1;

    const int edge_start = config.get_edge_selection_start();

    // Create sampler for this chain
    auto sampler = create_sampler(config);

    // =========================================================================
    // Warmup phase
    // =========================================================================
    for (int iter = 0; iter < config.no_warmup; ++iter) {

        // Impute missing data if applicable
        if (config.na_impute && model.has_missing_data()) {
            model.impute_missing();
        }

        // Edge selection starts after edge_start iterations
        if (config.edge_selection && iter >= edge_start && model.has_edge_selection()) {
            model.update_edge_indicators();
        }

        // Sampler step (unified interface)
        sampler->warmup_step(model);

        // Update edge prior parameters (Beta-Bernoulli, SBM, etc.)
        if (config.edge_selection && iter >= edge_start && model.has_edge_selection()) {
            edge_prior.update(
                model.get_edge_indicators(),
                model.get_inclusion_probability(),
                model.get_num_variables(),
                model.get_num_pairwise(),
                model.get_rng()
            );
        }

        // Progress and interrupt check
        pm.update(chain_id);
        if (pm.shouldExit()) {
            chain_result.userInterrupt = true;
            return;
        }
    }

    // Finalize warmup (samplers fix their adapted parameters)
    sampler->finalize_warmup();

    // =========================================================================
    // Activate edge selection mode (if enabled)
    // =========================================================================
    if (config.edge_selection && model.has_edge_selection()) {
        model.set_edge_selection_active(true);
        model.initialize_graph();  // Randomly initialize graph structure
    }

    // =========================================================================
    // Sampling phase
    // =========================================================================
    for (int iter = 0; iter < config.no_iter; ++iter) {

        // Impute missing data if applicable
        if (config.na_impute && model.has_missing_data()) {
            model.impute_missing();
        }

        // Edge selection continues during sampling
        if (config.edge_selection && model.has_edge_selection()) {
            model.update_edge_indicators();
        }

        // Sampler step (unified interface)
        SamplerResult result = sampler->sample_step(model);

        // Update edge prior parameters (Beta-Bernoulli, SBM, etc.)
        if (config.edge_selection && model.has_edge_selection()) {
            edge_prior.update(
                model.get_edge_indicators(),
                model.get_inclusion_probability(),
                model.get_num_variables(),
                model.get_num_pairwise(),
                model.get_rng()
            );
        }

        // Store NUTS diagnostics if available
        if (chain_result.has_nuts_diagnostics && sampler->has_nuts_diagnostics()) {
            auto* diag = dynamic_cast<NUTSDiagnostics*>(result.diagnostics.get());
            if (diag) {
                chain_result.store_nuts_diagnostics(iter, diag->tree_depth, diag->divergent, diag->energy);
            }
        }

        // Store samples
        chain_result.store_sample(iter, model.get_full_vectorized_parameters());

        // Store edge indicators if applicable
        if (chain_result.has_indicators) {
            chain_result.store_indicators(iter, model.get_vectorized_indicator_parameters());
        }

        // Progress and interrupt check
        pm.update(chain_id);
        if (pm.shouldExit()) {
            chain_result.userInterrupt = true;
            return;
        }
    }
}


/**
 * Worker struct for parallel chain execution
 */
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

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            ChainResultNew& chain_result = results_[i];
            BaseModel& model = *models_[i];
            BaseEdgePrior& edge_prior = *edge_priors_[i];
            model.set_seed(config_.seed + static_cast<int>(i));

            try {
                run_mcmc_chain(chain_result, model, edge_prior, config_, static_cast<int>(i), pm_);
            } catch (std::exception& e) {
                chain_result.error = true;
                chain_result.error_msg = e.what();
            } catch (...) {
                chain_result.error = true;
                chain_result.error_msg = "Unknown error";
            }
        }
    }
};


/**
 * Run MCMC sampling with parallel chains
 *
 * Main entry point for multi-chain MCMC. Handles:
 * - Chain allocation and model cloning
 * - Parallel or sequential execution based on no_threads
 * - Result collection
 *
 * @param model       Template model (will be cloned for each chain)
 * @param config      Sampler configuration
 * @param no_chains   Number of chains to run
 * @param no_threads  Number of threads (1 = sequential)
 * @param pm          Progress manager
 * @return Vector of chain results
 */
inline std::vector<ChainResultNew> run_mcmc_sampler(
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int no_chains,
    const int no_threads,
    ProgressManager& pm
) {
    const bool has_nuts_diag = (config.sampler_type == "nuts");

    // Allocate result storage
    std::vector<ChainResultNew> results(no_chains);
    for (int c = 0; c < no_chains; ++c) {
        results[c].reserve(model.full_parameter_dimension(), config.no_iter);

        if (config.edge_selection) {
            size_t n_edges = model.get_vectorized_indicator_parameters().n_elem;
            results[c].reserve_indicators(n_edges, config.no_iter);
        }

        if (has_nuts_diag) {
            results[c].reserve_nuts_diagnostics(config.no_iter);
        }
    }

    if (no_threads > 1) {
        // Multi-threaded execution
        std::vector<std::unique_ptr<BaseModel>> models;
        std::vector<std::unique_ptr<BaseEdgePrior>> edge_priors;
        models.reserve(no_chains);
        edge_priors.reserve(no_chains);
        for (int c = 0; c < no_chains; ++c) {
            models.push_back(model.clone());
            models[c]->set_seed(config.seed + c);
            edge_priors.push_back(edge_prior.clone());
        }

        MCMCChainRunner runner(results, models, edge_priors, config, pm);
        tbb::global_control control(tbb::global_control::max_allowed_parallelism, no_threads);
        RcppParallel::parallelFor(0, static_cast<size_t>(no_chains), runner);

    } else {
        // Single-threaded execution
        model.set_seed(config.seed);
        for (int c = 0; c < no_chains; ++c) {
            auto chain_model = model.clone();
            chain_model->set_seed(config.seed + c);
            auto chain_edge_prior = edge_prior.clone();
            run_mcmc_chain(results[c], *chain_model, *chain_edge_prior, config, c, pm);
        }
    }

    return results;
}


/**
 * Convert chain results to Rcpp::List format
 *
 * Creates a standardized output format for both GGM and OMRF models.
 * Each chain is a list with:
 *   - chain_id: Chain identifier
 *   - samples: Parameter samples matrix (param_dim × n_iter)
 *   - indicator_samples: Edge indicators (if edge_selection)
 *   - treedepth / divergent / energy: NUTS diagnostics (if NUTS/HMC)
 *   - error / error_msg: Error information (if error occurred)
 *
 * @param results  Vector of chain results
 * @return Rcpp::List with per-chain output
 */
inline Rcpp::List convert_results_to_list(const std::vector<ChainResultNew>& results) {
    Rcpp::List output(results.size());

    for (size_t i = 0; i < results.size(); ++i) {
        const ChainResultNew& chain = results[i];
        Rcpp::List chain_list;

        chain_list["chain_id"] = chain.chain_id;

        if (chain.error) {
            chain_list["error"] = true;
            chain_list["error_msg"] = chain.error_msg;
        } else {
            chain_list["error"] = false;
            chain_list["samples"] = chain.samples;
            chain_list["userInterrupt"] = chain.userInterrupt;

            if (chain.has_indicators) {
                chain_list["indicator_samples"] = chain.indicator_samples;
            }

            if (chain.has_nuts_diagnostics) {
                chain_list["treedepth"] = chain.treedepth_samples;
                chain_list["divergent"] = chain.divergent_samples;
                chain_list["energy"] = chain.energy_samples;
            }
        }

        output[i] = chain_list;
    }

    return output;
}
