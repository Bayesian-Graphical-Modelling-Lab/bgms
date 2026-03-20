#pragma once

#include <utility>
#include "mcmc/samplers/sampler_base.h"
#include "mcmc/algorithms/nuts.h"
#include "mcmc/algorithms/leapfrog.h"

/**
 * NUTSSampler - No-U-Turn Sampler
 *
 * Adaptive tree-depth leapfrog integration. Inherits warmup adaptation
 * (step size + diagonal mass matrix) from GradientSamplerBase.
 *
 * When the model has constraints (edge selection enabled), uses the
 * RATTLE constrained integrator in full Cholesky space with position
 * and momentum projection.
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
        if (model.has_constraints()) {
            return do_constrained_step(model);
        }
        return do_unconstrained_step(model);
    }

private:
    StepResult do_unconstrained_step(BaseModel& model) {
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

    StepResult do_constrained_step(BaseModel& model) {
        arma::vec x = model.get_full_position();
        SafeRNG& rng = model.get_rng();

        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient_full(params);
        };

        ProjectFn project = [&model](arma::vec& pos, arma::vec& mom) {
            model.project_position(pos);
            model.project_momentum(mom, pos);
        };

        arma::vec inv_mass = model.get_inv_mass();

        StepResult result = nuts_step(
            x, step_size_, joint_fn,
            inv_mass, rng, max_tree_depth_,
            &project
        );

        model.set_full_position(result.state);
        return result;
    }

    int max_tree_depth_;
};
