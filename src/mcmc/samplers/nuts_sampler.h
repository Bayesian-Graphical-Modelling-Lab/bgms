#pragma once

#include <utility>
#include "mcmc/samplers/sampler_base.h"
#include "mcmc/algorithms/nuts.h"

/**
 * NUTSSampler - No-U-Turn Sampler
 *
 * Adaptive tree-depth leapfrog integration. Inherits warmup adaptation
 * (step size + diagonal mass matrix) from GradientSamplerBase.
 */
class NUTSSampler : public GradientSamplerBase {
public:
    explicit NUTSSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : GradientSamplerBase(config.initial_step_size, config.target_acceptance, schedule),
          max_tree_depth_(config.max_tree_depth)
    {}

    bool has_nuts_diagnostics() const override { return true; }

protected:
    StepResult do_gradient_step(BaseModel& model) override {
        arma::vec theta = model.get_vectorized_parameters();
        SafeRNG& rng = model.get_rng();

        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(params);
        };

        arma::vec active_inv_mass = model.get_active_inv_mass();

        StepResult result = nuts_step(
            theta, step_size_, joint_fn,
            active_inv_mass, rng, max_tree_depth_
        );

        model.set_vectorized_parameters(result.state);
        return result;
    }

private:
    int max_tree_depth_;
};
