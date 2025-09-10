#include <RcppArmadillo.h>
#include <functional>
#include "mcmc_hmc.h"
#include "mcmc_leapfrog.h"
#include "mcmc_utils.h"
#include "rng_utils.h"
using namespace Rcpp;



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

  // Hamiltonians
  double current_H = -log_post(init_theta) + kinetic_energy(init_r, inv_mass_diag);
  double proposed_H = -log_post(theta) + kinetic_energy(r, inv_mass_diag);
  double log_accept_prob = current_H - proposed_H;

  arma::vec state = (std::log(runif(rng)) < log_accept_prob) ? theta : init_theta;

  double accept_prob = std::min(1.0, std::exp(log_accept_prob));

  return {state, accept_prob};
}
