#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_hmc.h"
#include "mcmc/mcmc_leapfrog.h"
#include "mcmc/mcmc_utils.h"
#include "rng/rng_utils.h"



SamplerResult hmc_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng
) {
  // Sample initial momentum
  arma::vec init_r = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, init_theta.n_elem);

  // Get log_post and gradient at initial position using joint
  auto [log_post_0, grad_0] = joint(init_theta);

  // Leapfrog with joint - passes initial grad, uses joint at final position
  LeapfrogJointResult result = leapfrog(
    init_theta, init_r, step_size, grad, joint, num_leapfrogs,
    inv_mass_diag, &grad_0
  );

  // Hamiltonians using the values from joint calls
  double current_H = -log_post_0 + kinetic_energy(init_r, inv_mass_diag);
  double proposed_H = -result.log_post + kinetic_energy(result.r, inv_mass_diag);
  double log_accept_prob = current_H - proposed_H;

  arma::vec state = (MY_LOG(runif(rng)) < log_accept_prob) ? result.theta : init_theta;

  double accept_prob = std::min(1.0, MY_EXP(log_accept_prob));

  return {state, accept_prob};
}
