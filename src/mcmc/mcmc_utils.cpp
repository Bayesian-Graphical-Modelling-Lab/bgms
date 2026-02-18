#include <RcppArmadillo.h>
#include <cmath>
#include <functional>
#include "mcmc/mcmc_leapfrog.h"
#include "mcmc/mcmc_utils.h"
#include "rng/rng_utils.h"



/**
 * Function: kinetic_energy
 *
 * Computes the kinetic energy of a momentum vector r, assuming a standard multivariate normal distribution.
 * This is used as part of the Hamiltonian energy in Hamiltonian Monte Carlo.
 *
 * Inputs:
 *  - r: The momentum vector.
 *  - inv_mass_diag: Diagonal of Inverse Mass Matrix
 * Returns:
 *  - The scalar kinetic energy value (0.5 * r^T * r).
 */
double kinetic_energy(const arma::vec& r, const arma::vec& inv_mass_diag) {
  return 0.5 * arma::dot(r % inv_mass_diag, r);
}



/**
 * Function: heuristic_initial_step_size
 *
 * Finds a reasonable initial step size for HMC by iteratively doubling or
 * halving until the acceptance probability is near the target. Uses identity
 * mass matrix. This overload delegates to the full version with inv_mass_diag
 * set to ones.
 *
 * Inputs:
 *  - theta: Initial parameter vector.
 *  - grad: Function returning gradient of log-posterior.
 *  - joint: Function returning both log-posterior and gradient.
 *  - rng: Thread-safe random number generator.
 *  - target_acceptance: Target acceptance probability (default 0.65).
 *  - init_step: Starting step size to search from.
 *  - max_attempts: Maximum doubling/halving iterations.
 *
 * Returns:
 *  - A step size yielding acceptance probability near target_acceptance.
 */
double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    SafeRNG& rng,
    double target_acceptance,
    double init_step,
    int max_attempts
) {
  arma::vec inv_mass_diag = arma::ones<arma::vec>(theta.n_elem);
  return heuristic_initial_step_size(
    theta, grad, joint, inv_mass_diag, rng, target_acceptance, init_step, max_attempts
  );
}


/**
 * Function: heuristic_initial_step_size
 *
 * Finds a reasonable initial step size for HMC by iteratively doubling or
 * halving until the acceptance probability is near the target. Performs single
 * leapfrog steps from the initial position, adjusting step size based on
 * whether the Hamiltonian change exceeds log(0.5) or log(2).
 *
 * Inputs:
 *  - theta: Initial parameter vector.
 *  - grad: Function returning gradient of log-posterior.
 *  - joint: Function returning both log-posterior and gradient.
 *  - inv_mass_diag: Diagonal inverse mass matrix.
 *  - rng: Thread-safe random number generator.
 *  - target_acceptance: Target acceptance probability (default 0.65).
 *  - init_step: Starting step size to search from.
 *  - max_attempts: Maximum doubling/halving iterations.
 *
 * Returns:
 *  - A step size yielding acceptance probability near target_acceptance.
 */
double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    double target_acceptance,
    double init_step,
    int max_attempts
) {
  double eps = init_step;

  // Use joint at initial position (compute once - position doesn't change)
  auto [logp0, grad0] = joint(theta);

  // Sample initial momentum and evaluate
  arma::vec r = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, theta.n_elem);
  double kin0 = kinetic_energy(r, inv_mass_diag);
  double H0 = logp0 - kin0;

  // One leapfrog step using leapfrog
  LeapfrogJointResult result = leapfrog(
    theta, r, eps, grad, joint, 1, inv_mass_diag, &grad0
  );

  double kin1 = kinetic_energy(result.r, inv_mass_diag);
  double H1 = result.log_post - kin1;

  int direction = 2 * (H1 - H0 > MY_LOG(0.5)) - 1;  // +1 or -1

  int attempts = 0;
  while (direction * (H1 - H0) > -direction * MY_LOG(2.0) && attempts < max_attempts) {
    eps = (direction == 1) ? 2.0 * eps : 0.5 * eps;

    // Resample momentum on each iteration for step size search
    r = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, theta.n_elem);
    kin0 = kinetic_energy(r, inv_mass_diag);
    H0 = logp0 - kin0;

    // One leapfrog step from original position with new momentum
    result = leapfrog(theta, r, eps, grad, joint, 1, inv_mass_diag, &grad0);

    // Evaluate Hamiltonian
    kin1 = kinetic_energy(result.r, inv_mass_diag);
    H1 = result.log_post - kin1;

    attempts++;
  }

  return eps;
}
