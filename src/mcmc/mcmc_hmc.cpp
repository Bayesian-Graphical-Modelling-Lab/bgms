#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_hmc.h"
#include "mcmc/mcmc_leapfrog.h"
#include "mcmc/mcmc_utils.h"
#include "rng/rng_utils.h"



SamplerResult hmc_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<double(const arma::vec&)>& log_post,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng
) {
  arma::vec theta = init_theta;
  arma::vec init_r = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, theta.n_elem);
  arma::vec r = init_r;

  std::tie(theta, r) = leapfrog(
    theta, r, step_size, grad, num_leapfrogs, inv_mass_diag
  );

  // If leapfrog produced NaN/Inf, reject immediately
  if (theta.has_nan() || theta.has_inf() || r.has_nan() || r.has_inf()) {
    return {init_theta, 0.0};
  }

  // Hamiltonians
  double current_H = -log_post(init_theta) + kinetic_energy(init_r, inv_mass_diag);
  double proposed_H = -log_post(theta) + kinetic_energy(r, inv_mass_diag);
  double log_accept_prob = current_H - proposed_H;

  // NaN guard: treat non-finite Hamiltonian as rejection
  if (!std::isfinite(log_accept_prob)) {
    return {init_theta, 0.0};
  }

  arma::vec state = (MY_LOG(runif(rng)) < log_accept_prob) ? theta : init_theta;
  double accept_prob = std::min(1.0, MY_EXP(log_accept_prob));

  return {state, accept_prob};
}
