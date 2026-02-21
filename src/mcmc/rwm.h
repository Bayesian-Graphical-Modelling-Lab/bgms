#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/sampler_result.h"
struct SafeRNG;


/**
 * Performs one step of Random Walk Metropolis sampling for a scalar parameter
 *
 * Proposes from a symmetric normal distribution centered at the current state,
 * then accepts or rejects via the Metropolis-Hastings criterion.
 *
 * @param current_state  Current scalar parameter value
 * @param step_size      Standard deviation of the Gaussian proposal
 * @param log_post       Log-posterior function (scalar to scalar)
 * @param rng            Thread-safe random number generator
 * @return SamplerResult with accepted state (1-element vector) and acceptance probability
 */
SamplerResult rwm_sampler(
    double current_state,
    double step_size,
    const std::function<double(double)>& log_post,
    SafeRNG& rng
);