#pragma once

#include <RcppArmadillo.h>
#include "mcmc/base_sampler.h"
#include "mcmc/mcmc_utils.h"
#include "mcmc/sampler_config.h"
#include "mcmc/mcmc_adaptation.h"
#include "models/base_model.h"

/**
 * MHSampler - Metropolis-Hastings sampler
 *
 * Delegates to the model's component-wise MH updates. RWM proposal-SD
 * adaptation is handled by the model via RWMAdaptationController instances,
 * initialized lazily on the first step.
 *
 * This is a thin wrapper that provides a uniform interface consistent
 * with other samplers (NUTS, HMC), but the actual sampling logic
 * (component-wise updates, Gibbs sweeps, etc.) is model-specific.
 */
class MHSampler : public BaseSampler {
public:
    /**
     * Construct MH sampler with configuration and warmup schedule
     * @param config    Sampler configuration
     * @param schedule  Shared warmup schedule
     */
    MHSampler(const SamplerConfig& config, WarmupSchedule& schedule)
        : schedule_(schedule), initialized_(false) {
        (void)config;
    }

    /**
     * Eager initialization: sets up RWM adaptation controllers.
     * Called from the runner before the MCMC loop.
     */
    void initialize(BaseModel& model) override {
        if (initialized_) return;
        model.init_mh_adaptation(schedule_);
        initialized_ = true;
    }

    /**
     * Perform one MH step with RWM adaptation
     */
    SamplerResult step(BaseModel& model, int iteration) override {
        if (!initialized_) {
            initialize(model);
        }

        model.do_one_mh_step(iteration);

        SamplerResult result;
        result.accept_prob = 1.0;
        return result;
    }

private:
    WarmupSchedule& schedule_;
    bool initialized_;
};
