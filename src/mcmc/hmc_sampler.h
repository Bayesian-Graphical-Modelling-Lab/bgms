#pragma once

#include "mcmc/adaptive_gradient_sampler.h"
#include "mcmc/mcmc_hmc.h"

/**
 * HMCSampler - Hamiltonian Monte Carlo
 *
 * Fixed-length leapfrog integration. Inherits warmup adaptation
 * (step size + diagonal mass matrix) from AdaptiveGradientSampler.
 */
class HMCSampler : public AdaptiveGradientSampler {
public:
    explicit HMCSampler(const SamplerConfig& config)
        : AdaptiveGradientSampler(config.initial_step_size, config.target_acceptance, config.no_warmup),
          num_leapfrogs_(config.num_leapfrogs)
    {}

protected:
    SamplerResult do_gradient_step(BaseModel& model) override {
        arma::vec theta = model.get_vectorized_parameters();
        arma::vec inv_mass = model.get_active_inv_mass();
        SafeRNG& rng = model.get_rng();

        auto log_post = [&model](const arma::vec& params) -> double {
            return model.logp_and_gradient(params).first;
        };
        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };

        SamplerResult result = hmc_sampler(
            theta, step_size_, log_post, grad_fn,
            num_leapfrogs_, inv_mass, rng);

        model.set_vectorized_parameters(result.state);
        return result;
    }

private:
    int num_leapfrogs_;
};
