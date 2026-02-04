#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <vector>
#include <memory>
#include "base_sampler.h"
#include "mcmc_utils.h"
#include "mcmc_nuts.h"
#include "mcmc_adaptation.h"
#include "sampler_config.h"
#include "../base_model.h"

/**
 * NUTSSampler - No-U-Turn Sampler implementation
 *
 * Provides a clean interface to the NUTS algorithm for any BaseModel
 * with gradient support. Handles:
 * - Step size adaptation via dual averaging during warmup
 * - Mass matrix adaptation using Welford's algorithm during warmup
 * - Trajectory simulation with the no-U-turn criterion
 * - Diagnostics collection (tree depth, divergence, energy)
 *
 * The sampler fully owns all sampling logic. The model only provides:
 * - logp_and_gradient(theta): Compute log posterior and gradient
 * - get_vectorized_parameters(): Get current state as vector
 * - set_vectorized_parameters(theta): Update model state from vector
 * - get_active_inv_mass(): Get inverse mass diagonal (used for initialization)
 * - get_rng(): Get random number generator
 *
 * Warmup Schedule (Stan-style):
 * - Stage 1 (7.5%): Initial adaptation, step size only
 * - Stage 2 (82.5%): Mass matrix learning in doubling windows
 * - Stage 3 (10%): Final step size tuning with fixed mass
 *
 * Usage:
 *   NUTSSampler nuts(config, n_warmup);
 *   for (iter in warmup) {
 *       auto result = nuts.warmup_step(model);
 *   }
 *   nuts.finalize_warmup();
 *   for (iter in sampling) {
 *       auto result = nuts.sample_step(model);
 *   }
 */
class NUTSSampler : public BaseSampler {
public:
    /**
     * Construct NUTS sampler with configuration
     * @param config   Sampler configuration (step size, target acceptance, etc.)
     * @param n_warmup Number of warmup iterations (for scheduling mass matrix adaptation)
     */
    explicit NUTSSampler(const SamplerConfig& config, int n_warmup = 1000)
        : step_size_(config.initial_step_size),
          target_acceptance_(config.target_acceptance),
          max_tree_depth_(config.max_tree_depth),
          no_warmup_(config.no_warmup),
          n_warmup_(n_warmup),
          warmup_iteration_(0),
          initialized_(false),
          step_adapter_(config.initial_step_size)
    {
        build_warmup_schedule(n_warmup);
    }

    /**
     * Perform one NUTS step during warmup (with step size and mass matrix adaptation)
     * @param model  The model to sample from
     * @return SamplerResult with state and diagnostics
     */
    SamplerResult warmup_step(BaseModel& model) override {
        // Initialize on first warmup iteration
        if (!initialized_) {
            initialize(model);
            initialized_ = true;
        }

        SamplerResult result = do_nuts_step(model);

        // Adapt step size during all warmup phases
        step_adapter_.update(result.accept_prob, target_acceptance_);
        step_size_ = step_adapter_.current();

        // During Stage 2, accumulate samples for mass matrix estimation
        if (in_stage2()) {
            mass_accumulator_->update(result.state);

            // Check if we're at the end of a window
            if (at_window_end()) {
                // Update mass matrix from accumulated samples
                inv_mass_ = mass_accumulator_->inverse_mass();
                mass_accumulator_->reset();

                // Restart step size adaptation with new mass matrix
                step_adapter_.restart(step_size_);
            }
        }

        warmup_iteration_++;
        return result;
    }

    /**
     * Finalize warmup phase (fix step size to averaged value)
     */
    void finalize_warmup() override {
        step_size_ = step_adapter_.averaged();
    }

    /**
     * Perform one NUTS step during sampling (fixed step size and mass matrix)
     * @param model  The model to sample from
     * @return SamplerResult with state and diagnostics
     */
    SamplerResult sample_step(BaseModel& model) override {
        return do_nuts_step(model);
    }

    /**
     * NUTS produces tree depth, divergence, and energy diagnostics
     */
    bool has_nuts_diagnostics() const override { return true; }

    /**
     * Get the current (or final) step size
     */
    double get_step_size() const { return step_size_; }

    /**
     * Get the averaged step size (for reporting after warmup)
     */
    double get_averaged_step_size() const {
        return step_adapter_.averaged();
    }

    /**
     * Get the current inverse mass matrix diagonal
     */
    const arma::vec& get_inv_mass() const { return inv_mass_; }

private:
    /**
     * Build Stan-style warmup schedule with doubling windows
     */
    void build_warmup_schedule(int n_warmup) {
        // Stage 1: 7.5% of warmup
        stage1_end_ = static_cast<int>(0.075 * n_warmup);

        // Stage 3 starts at 90% of warmup
        stage3_start_ = n_warmup - static_cast<int>(0.10 * n_warmup);

        // Stage 2: build doubling windows between stage1_end and stage3_start
        window_ends_.clear();
        int cur = stage1_end_;
        int wsize = 25;  // Initial window size

        while (cur < stage3_start_) {
            int win = std::min(wsize, stage3_start_ - cur);
            window_ends_.push_back(cur + win);
            cur += win;
            wsize = std::min(wsize * 2, stage3_start_ - cur);
        }
    }

    /**
     * Check if we're in Stage 2 (mass matrix learning phase)
     */
    bool in_stage2() const {
        return warmup_iteration_ >= stage1_end_ && warmup_iteration_ < stage3_start_;
    }

    /**
     * Check if we're at the end of a Stage 2 window
     */
    bool at_window_end() const {
        for (int end : window_ends_) {
            if (warmup_iteration_ + 1 == end) {
                return true;
            }
        }
        return false;
    }

    /**
     * Initialize step size and mass matrix on first iteration
     */
    void initialize(BaseModel& model) {
        arma::vec theta = model.get_vectorized_parameters();
        SafeRNG& rng = model.get_rng();

        // Initialize inverse mass to identity (or from model)
        inv_mass_ = model.get_active_inv_mass();

        // Initialize mass matrix accumulator
        mass_accumulator_ = std::make_unique<DiagMassMatrixAccumulator>(
            static_cast<int>(theta.n_elem));

        // Create log posterior and gradient functions
        auto log_post = [&model](const arma::vec& params) -> double {
            return model.logp_and_gradient(params).first;
        };
        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };

        // Use heuristic to find good initial step size
        step_size_ = heuristic_initial_step_size(
            theta, log_post, grad_fn, rng, target_acceptance_);

        // Restart dual averaging with the heuristic step size
        step_adapter_.restart(step_size_);
    }

    /**
     * Execute one NUTS step using the sampler's learned mass matrix
     */
    SamplerResult do_nuts_step(BaseModel& model) {
        // Get current state
        arma::vec theta = model.get_vectorized_parameters();
        SafeRNG& rng = model.get_rng();

        // Create log posterior and gradient functions that call the model
        auto log_post = [&model](const arma::vec& params) -> double {
            return model.logp_and_gradient(params).first;
        };
        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };

        // Call the NUTS free function with our learned inverse mass
        SamplerResult result = nuts_sampler(
            theta,
            step_size_,
            log_post,
            grad_fn,
            inv_mass_,
            rng,
            max_tree_depth_
        );

        // Update model state with new parameters
        model.set_vectorized_parameters(result.state);

        return result;
    }

    // Configuration
    double step_size_;
    double target_acceptance_;
    int max_tree_depth_;
    int no_warmup_;
    int n_warmup_;

    // State tracking
    int warmup_iteration_;
    bool initialized_;

    // Step size adaptation
    DualAveraging step_adapter_;

    // Mass matrix adaptation
    arma::vec inv_mass_;
    std::unique_ptr<DiagMassMatrixAccumulator> mass_accumulator_;

    // Warmup schedule
    int stage1_end_;
    int stage3_start_;
    std::vector<int> window_ends_;
};
