#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "base_sampler.h"
#include "mcmc_utils.h"
#include "mcmc_hmc.h"
#include "sampler_config.h"
#include "../base_model.h"

/**
 * HMCSampler - Hamiltonian Monte Carlo sampler
 *
 * Uses fixed-length leapfrog integration with optional step size
 * adaptation during warmup via dual averaging.
 *
 * The sampler fully owns all sampling logic. The model only provides:
 * - logp_and_gradient(theta): Compute log posterior and gradient
 * - get_vectorized_parameters(): Get current state as vector
 * - set_vectorized_parameters(theta): Update model state from vector
 * - get_active_inv_mass(): Get inverse mass diagonal
 * - get_rng(): Get random number generator
 */
class HMCSampler : public BaseSampler {
public:
    /**
     * Construct HMC sampler with configuration
     * @param config  Sampler configuration
     */
    explicit HMCSampler(const SamplerConfig& config)
        : step_size_(config.initial_step_size),
          target_acceptance_(config.target_acceptance),
          num_leapfrogs_(config.num_leapfrogs),
          no_warmup_(config.no_warmup),
          warmup_iteration_(0)
    {
        // Initialize dual averaging state
        dual_avg_state_.set_size(3);
        dual_avg_state_(0) = std::log(step_size_);
        dual_avg_state_(1) = std::log(step_size_);
        dual_avg_state_(2) = 0.0;
    }

    /**
     * Perform one HMC step during warmup (with step size adaptation)
     */
    SamplerResult warmup_step(BaseModel& model) override {
        SamplerResult result = do_hmc_step(model);

        // Update step size via dual averaging
        warmup_iteration_++;
        update_step_size_with_dual_averaging(
            step_size_,
            result.accept_prob,
            warmup_iteration_,
            dual_avg_state_,
            target_acceptance_
        );
        step_size_ = std::exp(dual_avg_state_(0));

        return result;
    }

    /**
     * Finalize warmup phase (fix step size to averaged value)
     */
    void finalize_warmup() override {
        step_size_ = std::exp(dual_avg_state_(1));
    }

    /**
     * Perform one HMC step during sampling (fixed step size)
     */
    SamplerResult sample_step(BaseModel& model) override {
        return do_hmc_step(model);
    }

    double get_step_size() const { return step_size_; }
    double get_averaged_step_size() const { return std::exp(dual_avg_state_(1)); }

private:
    /**
     * Execute one HMC step using the model's interface
     */
    SamplerResult do_hmc_step(BaseModel& model) {
        // Get current state
        arma::vec theta = model.get_vectorized_parameters();
        arma::vec inv_mass = model.get_active_inv_mass();
        SafeRNG& rng = model.get_rng();

        // Create log posterior and gradient functions that call the model
        auto log_post = [&model](const arma::vec& params) -> double {
            return model.logp_and_gradient(params).first;
        };
        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };

        // Call the HMC free function
        SamplerResult result = hmc_sampler(
            theta,
            step_size_,
            log_post,
            grad_fn,
            num_leapfrogs_,
            inv_mass,
            rng
        );

        // Update model state with new parameters
        model.set_vectorized_parameters(result.state);

        return result;
    }

    double step_size_;
    double target_acceptance_;
    int num_leapfrogs_;
    int no_warmup_;
    int warmup_iteration_;
    arma::vec dual_avg_state_;
};
