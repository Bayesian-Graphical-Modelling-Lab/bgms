#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <memory>
#include "mcmc/base_sampler.h"
#include "mcmc/hamiltonian.h"
#include "mcmc/adaptation.h"
#include "mcmc/sampler_config.h"
#include "models/base_model.h"

/**
 * AdaptiveGradientSampler - Base for gradient-based MCMC with warmup adaptation
 *
 * Uses HMCAdaptationController (from adaptation.h) with the shared
 * WarmupSchedule constructed by the runner.
 *
 * The adaptation controller handles:
 *  - Step-size dual averaging (Stages 1, 2, 3a, 3c)
 *  - Mass matrix estimation in doubling windows (Stage 2)
 *  - Step-size freezing at Stage 3b boundary
 */
class AdaptiveGradientSampler : public BaseSampler {
public:
    AdaptiveGradientSampler(double step_size, double target_acceptance,
                            WarmupSchedule& schedule)
        : step_size_(step_size),
          target_acceptance_(target_acceptance),
          schedule_(schedule),
          initialized_(false)
    {}

    SamplerResult step(BaseModel& model, int iteration) override {
        // Use adaptation controller's current step size for this iteration
        step_size_ = adapt_->current_step_size();

        SamplerResult result = do_gradient_step(model);

        // Let the adaptation controller handle step-size and mass-matrix logic
        arma::vec full_params = model.get_full_vectorized_parameters();
        adapt_->update(full_params, result.accept_prob, iteration);

        // If mass matrix was just updated, apply it and re-run the step-size heuristic
        if (adapt_->mass_matrix_just_updated()) {
            arma::vec new_inv_mass = adapt_->inv_mass_diag();
            model.set_inv_mass(new_inv_mass);

            arma::vec theta = model.get_vectorized_parameters();
            SafeRNG& rng = model.get_rng();
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
    virtual SamplerResult do_gradient_step(BaseModel& model) = 0;

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
        arma::vec theta = model.get_vectorized_parameters();
        SafeRNG& rng = model.get_rng();

        // Initialize inverse mass to ones
        arma::vec init_inv_mass = arma::ones<arma::vec>(dim);
        model.set_inv_mass(init_inv_mass);

        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };
        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(params);
        };

        double init_eps = heuristic_initial_step_size(
            theta, grad_fn, joint_fn, rng, target_acceptance_);
        step_size_ = init_eps;

        // Construct the adaptation controller with the shared schedule
        adapt_ = std::make_unique<HMCAdaptationController>(
            dim, init_eps, target_acceptance_, schedule_);
    }

    WarmupSchedule& schedule_;
    bool initialized_;
    std::unique_ptr<HMCAdaptationController> adapt_;
};
