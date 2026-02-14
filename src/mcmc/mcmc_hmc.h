#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_utils.h"
struct SafeRNG;


/**
 * Function: hmc_sampler
 *
 * Performs one iteration of Hamiltonian Monte Carlo sampling.
 *
 * Uses joint function at endpoints for log_post+gradient, grad-only for
 * intermediate steps, avoiding redundant probability computations.
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