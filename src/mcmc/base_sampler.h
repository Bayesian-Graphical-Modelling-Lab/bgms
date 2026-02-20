#pragma once

#include <RcppArmadillo.h>
#include "mcmc/mcmc_utils.h"
#include "mcmc/sampler_config.h"
#include "models/base_model.h"

/**
 * BaseSampler - Abstract base class for MCMC samplers
 *
 * Provides a unified interface for all MCMC sampling algorithms:
 * - Random-walk Metropolis-Hastings
 * - Hamiltonian Monte Carlo (HMC)
 * - No-U-Turn Sampler (NUTS)
 *
 * The sampler internally decides whether to adapt based on the iteration
 * number and its warmup schedule reference.
 */
class BaseSampler {
public:
    virtual ~BaseSampler() = default;

    /**
     * Perform one MCMC step
     *
     * The sampler internally decides whether to adapt based on the
     * iteration number and its warmup schedule reference.
     *
     * @param model      The model to sample from
     * @param iteration  Current iteration (0-based, spans warmup + sampling)
     * @return SamplerResult with new state and diagnostics
     */
    virtual SamplerResult step(BaseModel& model, int iteration) = 0;

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
