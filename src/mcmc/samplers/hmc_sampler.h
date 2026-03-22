#pragma once

#include <utility>
#include "mcmc/samplers/sampler_base.h"
#include "mcmc/algorithms/hmc.h"
#include "mcmc/algorithms/leapfrog.h"

/**
 * HMCSampler - Hamiltonian Monte Carlo
 *
 * Fixed-length leapfrog integration. Inherits warmup adaptation
 * (step size + diagonal mass matrix) from GradientSamplerBase.
 *
 * For constrained models (edge selection), uses RATTLE integration:
 * full Cholesky space with position and momentum projection at each
 * leapfrog step.
 */
class HMCSampler : public GradientSamplerBase {
public:
    explicit HMCSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : GradientSamplerBase(config.initial_step_size, config.target_acceptance, schedule),
          num_leapfrogs_(config.num_leapfrogs)
    {}

protected:
    StepResult do_gradient_step(BaseModel& model) override {
        if (uses_constrained_integration(model)) {
            return do_constrained_step(model);
        }
        return do_unconstrained_step(model);
    }

private:
    StepResult do_unconstrained_step(BaseModel& model) {
        arma::vec theta = model.get_vectorized_parameters();
        arma::vec inv_mass = model.get_active_inv_mass();
        SafeRNG& rng = model.get_rng();

        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };
        auto joint_fn = [&model](const arma::vec& params) -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(params);
        };

        StepResult result = hmc_step(
            theta, step_size_, grad_fn, joint_fn,
            num_leapfrogs_, inv_mass, rng);

        model.set_vectorized_parameters(result.state);
        return result;
    }

    StepResult do_constrained_step(BaseModel& model) {
        model.reset_projection_cache();
        arma::vec x = model.get_full_position();
        SafeRNG& rng = model.get_rng();

        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient_full(params);
        };

        arma::vec inv_mass = model.get_inv_mass();

        ProjectPositionFn proj_pos = [&model, &inv_mass](arma::vec& pos) {
            model.project_position(pos, inv_mass);
        };

        ProjectMomentumFn proj_mom = [&model, &inv_mass](arma::vec& mom, const arma::vec& pos) {
            model.project_momentum(mom, pos, inv_mass);
        };

        StepResult result = hmc_step(
            x, step_size_, joint_fn,
            num_leapfrogs_, inv_mass,
            proj_pos, proj_mom, rng);

        model.set_full_position(result.state);
        return result;
    }

    int num_leapfrogs_;
};
