#include "mcmc/mcmc_runner.h"

#include "mcmc/nuts_sampler.h"
#include "mcmc/hmc_sampler.h"
#include "mcmc/mh_sampler.h"


std::unique_ptr<BaseSampler> create_sampler(const SamplerConfig& config) {
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


void run_mcmc_chain(
    ChainResultNew& chain_result,
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int chain_id,
    ProgressManager& pm
) {
    chain_result.chain_id = chain_id + 1;

    const int edge_start = config.get_edge_selection_start();

    auto sampler = create_sampler(config);

    // Warmup phase
    for (int iter = 0; iter < config.no_warmup; ++iter) {

        if (config.na_impute && model.has_missing_data()) {
            model.impute_missing();
        }

        if (config.edge_selection && iter >= edge_start && model.has_edge_selection()) {
            model.update_edge_indicators();
        }

        sampler->warmup_step(model);

        if (config.edge_selection && iter >= edge_start && model.has_edge_selection()) {
            edge_prior.update(
                model.get_edge_indicators(),
                model.get_inclusion_probability(),
                model.get_num_variables(),
                model.get_num_pairwise(),
                model.get_rng()
            );
        }

        pm.update(chain_id);
        if (pm.shouldExit()) {
            chain_result.userInterrupt = true;
            return;
        }
    }

    sampler->finalize_warmup();

    // Activate edge selection mode
    if (config.edge_selection && model.has_edge_selection()) {
        model.set_edge_selection_active(true);
        model.initialize_graph();
    }

    // Sampling phase
    for (int iter = 0; iter < config.no_iter; ++iter) {

        if (config.na_impute && model.has_missing_data()) {
            model.impute_missing();
        }

        if (config.edge_selection && model.has_edge_selection()) {
            model.update_edge_indicators();
        }

        SamplerResult result = sampler->sample_step(model);

        if (config.edge_selection && model.has_edge_selection()) {
            edge_prior.update(
                model.get_edge_indicators(),
                model.get_inclusion_probability(),
                model.get_num_variables(),
                model.get_num_pairwise(),
                model.get_rng()
            );
        }

        if (chain_result.has_nuts_diagnostics && sampler->has_nuts_diagnostics()) {
            auto* diag = dynamic_cast<NUTSDiagnostics*>(result.diagnostics.get());
            if (diag) {
                chain_result.store_nuts_diagnostics(iter, diag->tree_depth, diag->divergent, diag->energy);
            }
        }

        chain_result.store_sample(iter, model.get_full_vectorized_parameters());

        if (chain_result.has_indicators) {
            chain_result.store_indicators(iter, model.get_vectorized_indicator_parameters());
        }

        if (chain_result.has_allocations && edge_prior.has_allocations()) {
            chain_result.store_allocations(iter, edge_prior.get_allocations());
        }

        pm.update(chain_id);
        if (pm.shouldExit()) {
            chain_result.userInterrupt = true;
            return;
        }
    }
}


void MCMCChainRunner::operator()(std::size_t begin, std::size_t end) {
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


std::vector<ChainResultNew> run_mcmc_sampler(
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int no_chains,
    const int no_threads,
    ProgressManager& pm
) {
    const bool has_nuts_diag = (config.sampler_type == "nuts");
    const bool has_sbm_alloc = edge_prior.has_allocations() ||
        (config.edge_selection && dynamic_cast<StochasticBlockEdgePrior*>(&edge_prior) != nullptr);

    std::vector<ChainResultNew> results(no_chains);
    for (int c = 0; c < no_chains; ++c) {
        results[c].reserve(model.full_parameter_dimension(), config.no_iter);

        if (config.edge_selection) {
            size_t n_edges = model.get_vectorized_indicator_parameters().n_elem;
            results[c].reserve_indicators(n_edges, config.no_iter);
        }

        if (has_sbm_alloc) {
            results[c].reserve_allocations(model.get_num_variables(), config.no_iter);
        }

        if (has_nuts_diag) {
            results[c].reserve_nuts_diagnostics(config.no_iter);
        }
    }

    if (no_threads > 1) {
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


Rcpp::List convert_results_to_list(const std::vector<ChainResultNew>& results) {
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

            if (chain.has_allocations) {
                chain_list["allocation_samples"] = chain.allocation_samples;
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
