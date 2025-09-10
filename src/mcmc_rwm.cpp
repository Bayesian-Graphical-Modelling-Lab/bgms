#include <RcppArmadillo.h>
#include <cmath>
#include <functional>
#include "mcmc_utils.h"
#include "mcmc_rwm.h"
#include "rng_utils.h"
using namespace Rcpp;



/**
 * Function: rwm_sampler
 *
 * Performs one step of the Random Walk Metropolis (RWM) sampler for a scalar parameter.
 *
 * Proposes a new value from a symmetric normal distribution centered at the current state,
 * evaluates the Metropolis-Hastings acceptance probability, and accepts or rejects the proposal.
 *
 * Inputs:
 *  - current_state: The current scalar parameter value.
 *  - step_size: Standard deviation of the Gaussian proposal distribution.
 *  - log_post: A function returning the log posterior at a given scalar value.
 *
 * Returns:
 *  - SamplerResult:
 *      * state: A 1-element vector containing the accepted or rejected value.
 *      * accept_prob: The acceptance probability (between 0 and 1).
 */
SamplerResult rwm_sampler(
    double current_state,
    double step_size,
    const std::function<double(double)>& log_post,
    SafeRNG& rng
) {
  double proposed_state = rnorm(rng, current_state, step_size);
  double log_accept = log_post(proposed_state) - log_post(current_state);
  double accept_prob = std::min(1.0, std::exp(log_accept));

  double state = (runif(rng) < accept_prob) ? proposed_state : current_state;

  arma::vec State(1);
  State[0] = state;

  return {State, accept_prob};
}

