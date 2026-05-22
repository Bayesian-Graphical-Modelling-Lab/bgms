#include "mcmc/execution/chain_runner.h"

#include <chrono>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <string>
#include <unistd.h>
#include <tbb/global_control.h>
#include "mcmc/samplers/nuts_sampler.h"
#include "mcmc/samplers/metropolis_sampler.h"


std::unique_ptr<SamplerBase> create_sampler(const SamplerConfig& config, WarmupSchedule& schedule) {
    if (config.sampler_type == "nuts") {
        return std::make_unique<NUTSSampler>(config, schedule);
    } else if (config.sampler_type == "mh" || config.sampler_type == "adaptive-metropolis") {
        return std::make_unique<MetropolisSampler>(config, schedule);
    } else {
        Rcpp::stop("Unknown sampler_type: '%s'", config.sampler_type.c_str());
    }
}


void run_mcmc_chain(
    ChainResult& chain_result,
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int chain_id,
    ProgressManager& pm
) {
    chain_result.chain_id = chain_id + 1;

    // Construct warmup schedule (shared by runner and sampler)
    const bool learn_sd = (config.sampler_type == "nuts");
    WarmupSchedule schedule(config.no_warmup, config.edge_selection, learn_sd);

    auto sampler = create_sampler(config, schedule);

    // Initialize sampler (step-size heuristic) before the main loop
    sampler->initialize(model);

    const int total_iter = config.no_warmup + config.no_iter;

    // Live-diagnostic envelope. Activated when BGMS_LIVE_DIAG=<path> and the
    // model carries V-ratio diagnostics (hierarchical-spec). Each line is
    // appended to <path>.<pid>.<chain> with explicit flush, so `tail -f`
    // shows real-time per-iter state. Format CSV.
    const char* env_live  = std::getenv("BGMS_LIVE_DIAG");
    const bool  live_diag = (env_live != nullptr)
                            && model.has_v_ratio_diagnostics();
    const char* env_every = std::getenv("BGMS_LIVE_DIAG_EVERY");
    const int   live_every = (env_every != nullptr)
                             ? std::max(1, std::atoi(env_every)) : 25;
    auto chain_t0 = std::chrono::steady_clock::now();
    std::ofstream diag_file;
    // Snapshot of mh_U counters at the LAST print, so the live diag can
    // emit per-window accept rate / per-window auto-reject rate (rather
    // than just cumulative totals).
    long long last_mh_U_attempts = 0;
    long long last_mh_U_accepts  = 0;
    if (live_diag) {
        std::string path = std::string(env_live)
                         + "." + std::to_string(::getpid())
                         + "." + std::to_string(chain_id + 1);
        diag_file.open(path, std::ios::out | std::ios::trunc);
        if (diag_file.is_open()) {
            diag_file << "iter,phase,K_depth,sign_V,wall_per_iter_s,"
                      << "wall_total_s,"
                      << "mh_U_attempts_window,mh_U_accepts_window,"
                      << "mh_U_attempts_total,mh_U_accepts_total\n";
            diag_file.flush();
        }
    }

    // ---- Main MCMC loop (warmup + sampling) ----
    for (int iter = 0; iter < total_iter; ++iter) {
        // Iteration wall-clock timer (only used when V-ratio diagnostics are
        // tracked; cheap regardless).
        auto iter_t0 = std::chrono::steady_clock::now();

        // Per-iteration preparation (e.g., shuffle edge order)
        model.prepare_iteration();

        // Auxiliary-U refresh (hierarchical-spec V/RR machinery). Gated by
        // the schedule so the U pool doesn't drift via PMMH dynamics during
        // early warmup when no Γ moves consume it.
        if (schedule.u_refresh_enabled(iter)) {
            model.refresh_auxiliary_u();
        }

        // Optional missing-data imputation
        if (config.na_impute && model.has_missing_data()) {
            model.impute_missing();
        }

        // Edge selection
        if (schedule.selection_enabled(iter) && model.has_edge_selection()) {
            if (iter == schedule.stage3c_start) {
                model.set_edge_selection_active(true);
            }
            model.update_edge_indicators();
        }

        // Main parameter update — adaptation is internal to sampler
        StepResult result = sampler->step(model, iter);

        // Stage 3b: proposal-SD tuning
        model.tune_proposal_sd(iter, schedule);

        // Edge prior update
        if (schedule.selection_enabled(iter) && model.has_edge_selection()) {
            edge_prior.update(
                model.get_edge_indicators(),
                model.get_inclusion_probability(),
                model.get_num_variables(),
                model.get_num_pairwise(),
                model.get_rng()
            );
        }

        // Store samples (only during sampling phase)
        if (schedule.sampling(iter)) {
            int sample_index = iter - config.no_warmup;

            if (chain_result.has_nuts_diagnostics && sampler->has_nuts_diagnostics()) {
                auto* diag = dynamic_cast<NUTSDiagnostics*>(result.diagnostics.get());
                if (diag) {
                    chain_result.store_nuts_diagnostics(sample_index, diag->tree_depth, diag->divergent, diag->non_reversible, diag->energy, diag->accept_prob);
                }
            }

            if (chain_result.has_am_diagnostics) {
                chain_result.store_am_diagnostics(sample_index, result.accept_prob);
            }

            if (chain_result.has_v_ratio_diagnostics) {
                chain_result.store_v_ratio_diagnostics(
                    sample_index,
                    model.current_sign_V(),
                    model.current_log_abs_V());
            }

            chain_result.store_sample(sample_index, model.get_storage_vectorized_parameters());

            if (chain_result.has_indicators) {
                chain_result.store_indicators(sample_index, model.get_vectorized_indicator_parameters());
            }

            if (chain_result.has_allocations && edge_prior.has_allocations()) {
                chain_result.store_allocations(sample_index, edge_prior.get_allocations());
            }
        }

        // Capture per-iter wall delta (always cheap; only used by the
        // K_depth_samples store below and the live-diag print).
        auto iter_t1 = std::chrono::steady_clock::now();
        double wall_iter = std::chrono::duration<double>(iter_t1 - iter_t0).count();

        // Per-iter K_depth + wall in the sampling phase, for post-hoc analysis.
        if (schedule.sampling(iter) && chain_result.has_v_ratio_diagnostics) {
            int sample_index = iter - config.no_warmup;
            chain_result.K_depth_samples(sample_index)   = model.current_K_depth();
            chain_result.iter_wall_samples(sample_index) = wall_iter;
        }

        // Live diagnostic: append one CSV row every live_every iters,
        // regardless of phase, so we can watch warmup vs sampling cost
        // transition in real time. Explicit flush per line.
        if (live_diag && diag_file.is_open()
            && (iter % live_every == 0 || iter + 1 == total_iter)) {
            double wall_total =
                std::chrono::duration<double>(iter_t1 - chain_t0).count();
            const char* phase = schedule.sampling(iter)
                ? "sample"
                : (schedule.selection_enabled(iter) ? "stage3c" : "warm");

            // Pull cumulative mh_U counters from the model's diag summary
            // (cheap: it's just a struct read). Compute per-window delta.
            Rcpp::List ds = model.get_diagnostics_summary();
            long long mh_att_tot = 0, mh_acc_tot = 0;
            if (ds.containsElementNamed("mh_U_attempts")) {
                mh_att_tot = static_cast<long long>(
                    Rcpp::as<double>(ds["mh_U_attempts"]));
            }
            if (ds.containsElementNamed("mh_U_accepts")) {
                mh_acc_tot = static_cast<long long>(
                    Rcpp::as<double>(ds["mh_U_accepts"]));
            }
            long long mh_att_win = mh_att_tot - last_mh_U_attempts;
            long long mh_acc_win = mh_acc_tot - last_mh_U_accepts;
            last_mh_U_attempts = mh_att_tot;
            last_mh_U_accepts  = mh_acc_tot;

            diag_file << iter
                      << "," << phase
                      << "," << model.current_K_depth()
                      << "," << model.current_sign_V()
                      << "," << wall_iter
                      << "," << wall_total
                      << "," << mh_att_win
                      << "," << mh_acc_win
                      << "," << mh_att_tot
                      << "," << mh_acc_tot << "\n";
            diag_file.flush();
        }

        pm.update(chain_id);
        if (pm.shouldExit()) {
            chain_result.userInterrupt = true;
            return;
        }
    }

    // Capture end-of-chain diagnostic snapshot from the model. For GGMModel
    // this surfaces the hierarchical auto-reject counters; for other models
    // the override returns an empty list and we just store that.
    chain_result.diagnostics_summary = model.get_diagnostics_summary();
}


void MCMCChainRunner::operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
        ChainResult& chain_result = results_[i];
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


std::vector<ChainResult> run_mcmc_sampler(
    BaseModel& model,
    BaseEdgePrior& edge_prior,
    const SamplerConfig& config,
    const int no_chains,
    const int no_threads,
    ProgressManager& pm
) {
    const bool has_nuts_diag = (config.sampler_type == "nuts");
    const bool has_am_diag = (config.sampler_type == "adaptive-metropolis" ||
                              config.sampler_type == "mh");
    const bool has_sbm_alloc = edge_prior.has_allocations() ||
        (config.edge_selection && dynamic_cast<StochasticBlockEdgePrior*>(&edge_prior) != nullptr);

    std::vector<ChainResult> results(no_chains);
    for (int c = 0; c < no_chains; ++c) {
        results[c].reserve(model.storage_dimension(), config.no_iter);

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

        if (has_am_diag) {
            results[c].reserve_am_diagnostics(config.no_iter);
        }

        if (model.has_v_ratio_diagnostics()) {
            results[c].reserve_v_ratio_diagnostics(config.no_iter);
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


Rcpp::List convert_results_to_list(const std::vector<ChainResult>& results) {
    Rcpp::List output(results.size());

    for (size_t i = 0; i < results.size(); ++i) {
        const ChainResult& chain = results[i];
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
                chain_list["non_reversible"] = chain.non_reversible_samples;
                chain_list["energy"] = chain.energy_samples;
                chain_list["accept_prob"] = chain.accept_prob_samples;
            }

            if (chain.has_am_diagnostics) {
                chain_list["am_accept_prob"] = chain.am_accept_prob_samples;
            }

            if (chain.has_v_ratio_diagnostics) {
                chain_list["v_sign"]    = chain.v_sign_samples;
                chain_list["v_log_abs"] = chain.v_log_abs_samples;
                chain_list["K_depth"]   = chain.K_depth_samples;
                chain_list["iter_wall"] = chain.iter_wall_samples;
            }

            if (chain.diagnostics_summary.size() > 0) {
                chain_list["diagnostics_summary"] = chain.diagnostics_summary;
            }
        }

        output[i] = chain_list;
    }

    return output;
}
