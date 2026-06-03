#pragma once

#include "mcmc/samplers/sampler_base.h"
#include "mcmc/execution/sampler_config.h"
#include "mcmc/execution/warmup_schedule.h"
#include "models/base_model.h"

/**
 * MetropolisSampler - Metropolis-Hastings sampler
 *
 * Delegates to the model's component-wise Metropolis updates. Proposal-SD
 * adaptation is handled by the model via MetropolisAdaptationController
 * instances, initialized lazily on the first step.
 *
 * This is a thin wrapper that provides a uniform interface consistent
 * with other samplers (NUTS), but the actual sampling logic
 * (component-wise updates, Gibbs sweeps, etc.) is model-specific.
 */
class MetropolisSampler : public SamplerBase {
public:
    /**
     * Construct Metropolis sampler with configuration and warmup schedule
     * @param config    Sampler configuration
     * @param schedule  Shared warmup schedule
     */
    MetropolisSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : schedule_(schedule), initialized_(false) {
        (void)config;
    }

    /**
     * Eager initialization: sets up Metropolis adaptation controllers.
     * Called from the runner before the MCMC loop.
     */
    void initialize(BaseModel& model) override {
        if (initialized_) return;
        model.init_metropolis_adaptation(schedule_);
        initialized_ = true;
    }

    /**
     * Perform one Metropolis step with adaptation.
     *
     * Returns the mean acceptance probability over all components updated
     * in the sweep (read from the model). This value is stored as the
     * per-iteration trace surfaced as the am_accept_prob__ output column
     * (build_output.R). It does NOT feed fit$am_diag$accept_prob, which is
     * a separate per-parameter empirical move-rate matrix computed in R by
     * summarize_am_diagnostics() from the stored samples.
     */
    StepResult step(BaseModel& model, int iteration) override {
        if (!initialized_) {
            initialize(model);
        }

        model.do_one_metropolis_step(iteration);

        StepResult result;
        result.accept_prob = model.last_metropolis_mean_accept_prob();
        return result;
    }

private:
    WarmupSchedule& schedule_;
    bool initialized_;
};
