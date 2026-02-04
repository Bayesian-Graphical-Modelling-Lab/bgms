#pragma once

#include <RcppArmadillo.h>
#include "mcmc_utils.h"
#include "sampler_config.h"
#include "../base_model.h"

/**
 * BaseSampler - Abstract base class for MCMC samplers
 *
 * Provides a unified interface for all MCMC sampling algorithms:
 * - Random-walk Metropolis-Hastings
 * - Hamiltonian Monte Carlo (HMC)
 * - No-U-Turn Sampler (NUTS)
 *
 * All samplers follow the same workflow:
 *   1. warmup_step() during warmup (may adapt parameters)
 *   2. finalize_warmup() after warmup completes
 *   3. sample_step() during sampling (fixed parameters)
 *
 * The sampler asks the model for logp/gradient evaluations but owns
 * the sampling algorithm logic.
 */
class BaseSampler {
public:
    virtual ~BaseSampler() = default;

    /**
     * Perform one step during warmup phase
     *
     * During warmup, samplers may adapt their parameters (step size,
     * proposal covariance, mass matrix, etc.)
     *
     * @param model  The model to sample from
     * @return SamplerResult with new state and diagnostics
     */
    virtual SamplerResult warmup_step(BaseModel& model) = 0;

    /**
     * Finalize warmup phase
     *
     * Called after all warmup iterations complete. Samplers should
     * fix their adapted parameters for the sampling phase.
     */
    virtual void finalize_warmup() {}

    /**
     * Perform one step during sampling phase
     *
     * Sampling steps use fixed parameters (no adaptation).
     *
     * @param model  The model to sample from
     * @return SamplerResult with new state and diagnostics
     */
    virtual SamplerResult sample_step(BaseModel& model) = 0;

    /**
     * Check if this sampler produces NUTS-style diagnostics
     * (tree depth, divergences, energy)
     */
    virtual bool has_nuts_diagnostics() const { return false; }
};
