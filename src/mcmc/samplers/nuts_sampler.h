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
 * Two integration modes for constrained models (edge selection):
 * - RATTLE (default, sampler_type = "nuts"): Full Cholesky space with
 *   position and momentum projection at each leapfrog step.
 * - Null-space (sampler_type = "nuts-nullspace"): Reduced theta-space
 *   where constraints are satisfied by construction via the null-space
 *   basis. No projections needed.
 */
class NUTSSampler : public GradientSamplerBase {
public:
    explicit NUTSSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : GradientSamplerBase(config.initial_step_size, config.target_acceptance, schedule,
                              /* force_nullspace = */ config.sampler_type == "nuts-nullspace"),
          max_tree_depth_(config.max_tree_depth)
    {}

    bool has_nuts_diagnostics() const override { return true; }

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
        SafeRNG& rng = model.get_rng();

        if (has_dense_transform()) {
            // Dense mass transform: operate in z-space where z = L^{-1} theta
            recompute_dense_transform(model);
            arma::vec z = arma::solve(arma::trimatl(L_active_), theta);

            auto joint_fn_z = [this, &model](const arma::vec& z_param)
                -> std::pair<double, arma::vec> {
                arma::vec theta_param = L_active_ * z_param;
                auto [lp, g] = model.logp_and_gradient(theta_param);
                return {lp, L_active_.t() * g};
            };

            arma::vec unit_mass = arma::ones<arma::vec>(z.n_elem);

            StepResult result = nuts_step(
                z, step_size_, joint_fn_z,
                unit_mass, rng, max_tree_depth_
            );

            // Transform back to theta-space
            arma::vec theta_new = L_active_ * result.state;
            model.set_vectorized_parameters(theta_new);
            result.state = theta_new;
            return result;
        }

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

        arma::vec inv_mass = model.get_inv_mass();

        ProjectFn project = [&model, &inv_mass](arma::vec& pos, arma::vec& mom) {
            model.project_position(pos, inv_mass);
            model.project_momentum(mom, pos, inv_mass);
        };

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
