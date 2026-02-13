#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_utils.h"
struct SafeRNG;


/**
 * HMC sampler using joint log_post+gradient function.
 *
 * Eliminates redundant probability computations by:
 *  - Using joint at θ₀ to get (log_post_0, grad_0)
 *  - Using grad-only for intermediate positions
 *  - Using joint at θ_L to get (log_post_L, grad_L)
 *
 * This avoids computing probabilities twice at endpoints.
 */
SamplerResult hmc_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng
);