#pragma once

#include <RcppArmadillo.h>
#include "base_sampler.h"
#include "mcmc_utils.h"
#include "sampler_config.h"
#include "../base_model.h"

/**
 * MHSampler - Metropolis-Hastings sampler
 *
 * Delegates to the model's component-wise MH updates. The model
 * handles proposal adaptation internally during warmup.
 *
 * This is a thin wrapper that provides a uniform interface consistent
 * with other samplers (NUTS, HMC), but the actual sampling logic
 * (component-wise updates, Gibbs sweeps, etc.) is model-specific.
 */
class MHSampler : public BaseSampler {
public:
    /**
     * Construct MH sampler with configuration
     * @param config  Sampler configuration
     */
    explicit MHSampler(const SamplerConfig& config)
        : no_warmup_(config.no_warmup),
          warmup_iteration_(0)
    {}

    /**
     * Perform one MH step during warmup
     *
     * The model handles proposal adaptation internally.
     */
    SamplerResult warmup_step(BaseModel& model) override {
        warmup_iteration_++;
        model.do_one_mh_step();

        SamplerResult result;
        result.state = model.get_full_vectorized_parameters();
        result.accept_prob = 1.0;  // Not tracked for component-wise MH
        return result;
    }

    /**
     * Finalize warmup phase
     *
     * Nothing to do - model handles adaptation internally.
     */
    void finalize_warmup() override {
        // Model handles proposal finalization internally
    }

    /**
     * Perform one MH step during sampling
     */
    SamplerResult sample_step(BaseModel& model) override {
        model.do_one_mh_step();

        SamplerResult result;
        result.state = model.get_full_vectorized_parameters();
        result.accept_prob = 1.0;
        return result;
    }

private:
    int no_warmup_;
    int warmup_iteration_;
};
