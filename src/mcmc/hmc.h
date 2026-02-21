#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/sampler_result.h"
struct SafeRNG;


/**
 * Performs one iteration of Hamiltonian Monte Carlo sampling
 *
 * Proposes a new state by simulating Hamiltonian dynamics through leapfrog
 * integration, then accepts or rejects via the Metropolis criterion.
 * Uses joint function at endpoints for log_post+gradient, grad-only for
 * intermediate steps, avoiding redundant probability computations.
 *
 * @param init_theta     Initial parameter vector (position)
 * @param step_size      Leapfrog integration step size (epsilon)
 * @param grad           Gradient function
 * @param joint          Joint log-posterior + gradient function
 * @param num_leapfrogs  Number of leapfrog steps per proposal
 * @param inv_mass_diag  Diagonal of the inverse mass matrix
 * @param rng            Thread-safe random number generator
 * @return SamplerResult with accepted state and acceptance probability
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