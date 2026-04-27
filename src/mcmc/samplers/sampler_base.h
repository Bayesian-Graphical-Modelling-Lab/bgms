#pragma once

#include "mcmc/execution/step_result.h"
#include "models/base_model.h"

// ---------------------------------------------------------------------------
// SamplerBase — abstract interface for all MCMC samplers
// ---------------------------------------------------------------------------

/**
 * SamplerBase - Abstract base class for MCMC samplers
 *
 * Provides a unified interface for all MCMC sampling algorithms:
 * - MetropolisSampler (component-wise random-walk Metropolis)
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
