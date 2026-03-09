#pragma once

#include <utility>
#include "mcmc/samplers/sampler_base.h"
#include "mcmc/algorithms/nuts.h"
#include "models/mixed/mixed_mrf_model.h"

/**
 * HybridNUTSSampler - NUTS for unconstrained block + MH for Kyy
 *
 * Designed for MixedMRFModel where the NUTS block covers (mux, Kxx, muy, Kxy)
 * and the SPD-constrained Kyy is updated via component-wise Metropolis.
 * Inherits warmup adaptation (step size, diagonal mass matrix) from
 * GradientSamplerBase for the NUTS block.  Kyy proposal SDs are adapted
 * via the embedded Robbins-Monro schedule inside the model.
 */
class HybridNUTSSampler : public GradientSamplerBase {
public:
    explicit HybridNUTSSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : GradientSamplerBase(config.initial_step_size, config.target_acceptance, schedule),
          max_tree_depth_(config.max_tree_depth),
          schedule_(schedule)
    {}

    bool has_nuts_diagnostics() const override { return true; }

    void initialize(BaseModel& model) override {
        // Initialize NUTS adaptation (step-size heuristic + mass matrix)
        GradientSamplerBase::initialize(model);

        // Initialize Kyy Metropolis adaptation (stores total_warmup_)
        model.init_metropolis_adaptation(schedule_);
    }

protected:
    StepResult do_gradient_step(BaseModel& model) override {
        // --- Phase 1: NUTS step for the unconstrained block ---
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

        // --- Phase 2: Kyy Metropolis step ---
        auto& mixed = static_cast<MixedMRFModel&>(model);
        mixed.do_pairwise_continuous_metropolis_step(current_iteration_);

        return result;
    }

private:
    int max_tree_depth_;
    WarmupSchedule& schedule_;

public:
    // The iteration counter is set by the overridden step() method via
    // the base class.  We store it so do_gradient_step() can pass it
    // to the Kyy Metropolis update for Robbins-Monro adaptation.
    StepResult step(BaseModel& model, int iteration) override {
        current_iteration_ = iteration;
        return GradientSamplerBase::step(model, iteration);
    }

private:
    int current_iteration_ = -1;
};
