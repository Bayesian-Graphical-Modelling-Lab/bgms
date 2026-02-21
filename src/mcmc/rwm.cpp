#include <RcppArmadillo.h>
#include <cmath>
#include <functional>
#include "mcmc/sampler_result.h"
#include "mcmc/rwm.h"
#include "rng/rng_utils.h"


SamplerResult rwm_sampler(
    double current_state,
    double step_size,
    const std::function<double(double)>& log_post,
    SafeRNG& rng
) {
  double proposed_state = rnorm(rng, current_state, step_size);
  double log_accept = log_post(proposed_state) - log_post(current_state);
  double accept_prob = std::min(1.0, MY_EXP(log_accept));

  double state = (runif(rng) < accept_prob) ? proposed_state : current_state;

  arma::vec State(1);
  State[0] = state;

  return {State, accept_prob};
}

