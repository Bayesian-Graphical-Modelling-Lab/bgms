#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <memory>
#include <utility>
#include "mcmc/execution/step_result.h"
#include "mcmc/algorithms/hmc.h"
#include "mcmc/samplers/hmc_adaptation.h"
#include "models/base_model.h"

// ---------------------------------------------------------------------------
// SamplerBase — abstract interface for all MCMC samplers
// ---------------------------------------------------------------------------

/**
 * SamplerBase - Abstract base class for MCMC samplers
 *
 * Provides a unified interface for all MCMC sampling algorithms:
 * - MetropolisSampler (component-wise random-walk Metropolis)
 * - HMCSampler (Hamiltonian Monte Carlo)
 * - NUTSSampler (No-U-Turn Sampler)
 *
 * The sampler internally decides whether to adapt based on the iteration
 * number and its warmup schedule reference.
 */
class SamplerBase {
public:
    virtual ~SamplerBase() = default;

    /**
     * Perform one MCMC step
     *
     * The sampler internally decides whether to adapt based on the
     * iteration number and its warmup schedule reference.
     *
     * @param model      The model to sample from
     * @param iteration  Current iteration (0-based, spans warmup + sampling)
     * @return StepResult with new state and diagnostics
     */
    virtual StepResult step(BaseModel& model, int iteration) = 0;

    /**
     * Initialize the sampler before the MCMC loop.
     * For gradient-based samplers, runs the step-size heuristic. Default no-op.
     */
    virtual void initialize(BaseModel& /*model*/) {}

    /**
     * Check if this sampler produces NUTS-style diagnostics
     * (tree depth, divergences, energy)
     */
    virtual bool has_nuts_diagnostics() const { return false; }

};

// ---------------------------------------------------------------------------
// GradientSamplerBase — shared base for gradient-based MCMC (HMC, NUTS)
// ---------------------------------------------------------------------------

/**
 * GradientSamplerBase - Base for gradient-based MCMC with warmup adaptation
 *
 * Uses HMCAdaptationController (from hmc_adaptation.h) with the shared
 * WarmupSchedule constructed by the runner.
 *
 * The adaptation controller handles:
 *  - Step-size dual averaging (Stages 1, 2, 3a, 3c)
 *  - Mass matrix estimation in doubling windows (Stage 2)
 *  - Step-size freezing at Stage 3b boundary
 */
class GradientSamplerBase : public SamplerBase {
public:
    GradientSamplerBase(double step_size, double target_acceptance,
                        WarmupSchedule& schedule)
        : step_size_(step_size),
          target_acceptance_(target_acceptance),
          schedule_(schedule),
          initialized_(false)
    {}

    StepResult step(BaseModel& model, int iteration) override {
        // Stage 3c boundary: edge selection just activated.
        // Restart dual averaging so adaptation can tune to the new
        // geometry (changed active parameters) quickly.
        if (schedule_.in_stage3c(iteration) && !stage3c_initialized_) {
            stage3c_initialized_ = true;
            if (uses_constrained_integration(model)) {
                // Re-run step-size heuristic with constrained integrator.
                // The step size from stages 1-2 was tuned for unconstrained
                // leapfrog; constraints are now active, so re-tune.
                SafeRNG& rng = model.get_rng();
                arma::vec x = model.get_full_position();
                arma::vec inv_mass = model.get_inv_mass();
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient_full(params);
                };
                ProjectPositionFn proj_pos = [&model, &inv_mass](arma::vec& pos) {
                    model.project_position(pos, inv_mass);
                };
                ProjectMomentumFn proj_mom = [&model, &inv_mass](arma::vec& mom, const arma::vec& pos) {
                    model.project_momentum(mom, pos, inv_mass);
                };
                double new_eps = heuristic_initial_step_size_constrained(
                    x, joint_fn, inv_mass, proj_pos, proj_mom, rng,
                    0.625, adapt_->current_step_size());
                adapt_->reinit_stepsize(new_eps);
            } else {
                adapt_->reinit_stepsize(adapt_->current_step_size());
            }
        }

        // Use adaptation controller's current step size for this iteration
        step_size_ = adapt_->current_step_size();

        StepResult result = do_gradient_step(model);

        // Let the adaptation controller handle step-size and mass-matrix logic.
        // For RATTLE (constrained) models, feed x-space samples so the mass
        // matrix is estimated in the same coordinate system NUTS operates in.
        arma::vec full_params = uses_constrained_integration(model)
            ? model.get_full_position()
            : model.get_full_vectorized_parameters();
        adapt_->update(full_params, result.accept_prob, iteration);

        // If mass matrix was just updated, apply it and re-run the step-size heuristic
        if (adapt_->mass_matrix_just_updated()) {
            arma::vec new_inv_mass = adapt_->inv_mass_diag();
            model.set_inv_mass(new_inv_mass);

            SafeRNG& rng = model.get_rng();

            if (uses_constrained_integration(model)) {
                arma::vec x = model.get_full_position();
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient_full(params);
                };
                ProjectPositionFn proj_pos = [&model, &new_inv_mass](arma::vec& pos) {
                    model.project_position(pos, new_inv_mass);
                };
                ProjectMomentumFn proj_mom = [&model, &new_inv_mass](arma::vec& mom, const arma::vec& pos) {
                    model.project_momentum(mom, pos, new_inv_mass);
                };
                double new_eps = heuristic_initial_step_size_constrained(
                    x, joint_fn, new_inv_mass, proj_pos, proj_mom, rng,
                    0.625, adapt_->current_step_size());
                adapt_->reinit_stepsize(new_eps);
            } else {
                arma::vec theta = model.get_vectorized_parameters();
                auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                    return model.logp_and_gradient(params).second;
                };
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient(params);
                };
                arma::vec active_inv_mass = model.get_active_inv_mass();
                double new_eps = heuristic_initial_step_size(
                    theta, grad_fn, joint_fn, active_inv_mass, rng,
                    0.625, adapt_->current_step_size());
                adapt_->reinit_stepsize(new_eps);
            }
        }

        // Update step_size_ from controller (may have changed due to mass update)
        step_size_ = adapt_->current_step_size();

        return result;
    }

    double get_step_size() const { return step_size_; }
    double get_averaged_step_size() const {
        return adapt_ ? adapt_->final_step_size() : step_size_;
    }
    const arma::vec& get_inv_mass() const { return adapt_->inv_mass_diag(); }

protected:
    virtual StepResult do_gradient_step(BaseModel& model) = 0;

    bool uses_constrained_integration(const BaseModel& model) const {
        return model.has_constraints();
    }

    double step_size_;
    double target_acceptance_;

public:
    /**
     * Initialize the adaptation controller and run the step-size heuristic.
     * Called by the runner before the MCMC loop.
     */
    void initialize(BaseModel& model) override {
        if (initialized_) return;
        do_initialize(model);
        initialized_ = true;
    }

private:
    void do_initialize(BaseModel& model) {
        int dim = static_cast<int>(model.full_parameter_dimension());
        SafeRNG& rng = model.get_rng();

        // Initialize inverse mass to ones
        arma::vec init_inv_mass = arma::ones<arma::vec>(dim);
        model.set_inv_mass(init_inv_mass);

        double init_eps;

        if (uses_constrained_integration(model)) {
            // Project initial position onto constraint manifold before
            // computing step size. The MLE initialization may violate
            // K_ij = 0 constraints for excluded edges.
            arma::vec x = model.get_full_position();
            arma::vec r_dummy = arma::zeros<arma::vec>(x.n_elem);
            model.project_position(x);
            model.project_momentum(r_dummy, x);
            model.set_full_position(x);

            x = model.get_full_position();
            auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                return model.logp_and_gradient_full(params).second;
            };
            auto joint_fn = [&model](const arma::vec& params)
                -> std::pair<double, arma::vec> {
                return model.logp_and_gradient_full(params);
            };
            init_eps = heuristic_initial_step_size(
                x, grad_fn, joint_fn, rng, target_acceptance_);
        } else {
            arma::vec theta = model.get_vectorized_parameters();
            auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                return model.logp_and_gradient(params).second;
            };
            auto joint_fn = [&model](const arma::vec& params)
                -> std::pair<double, arma::vec> {
                return model.logp_and_gradient(params);
            };
            init_eps = heuristic_initial_step_size(
                theta, grad_fn, joint_fn, rng, target_acceptance_);
        }

        step_size_ = init_eps;

        // Construct the adaptation controller with the shared schedule
        adapt_ = std::make_unique<HMCAdaptationController>(
            dim, init_eps, target_acceptance_, schedule_,
            /*learn_mass_matrix=*/true);
    }

    WarmupSchedule& schedule_;
    bool initialized_;
    bool stage3c_initialized_ = false;
    std::unique_ptr<HMCAdaptationController> adapt_;
};
