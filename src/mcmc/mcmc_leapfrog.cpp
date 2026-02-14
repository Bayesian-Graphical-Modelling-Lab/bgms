#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_leapfrog.h"
#include "mcmc/mcmc_memoization.h"



/**
 * Function: leapfrog_memo
 *
 * Performs a single leapfrog step using memoization for the gradient evaluations.
 * This improves computational efficiency when the same positions are revisited.
 *
 * Inputs:
 *  - theta: Current position (parameter vector).
 *  - r: Current momentum vector.
 *  - eps: Step size for integration.
 *  - memo: Memoizer object caching gradient evaluations.
 *
 * Returns:
 *  - A pair containing:
 *      - Updated position vector.
 *      - Updated momentum vector.
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  auto grad1 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad1;
  theta_new += eps * (inv_mass_diag % r_half);
  auto grad2 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}


/**
 * Function: leapfrog
 *
 * Performs leapfrog integration for Hamiltonian Monte Carlo. Simulates
 * Hamiltonian dynamics by alternating half-step momentum updates with
 * full-step position updates over the specified number of steps.
 *
 * Inputs:
 *  - theta_init: Initial parameter vector (position).
 *  - r_init: Initial momentum vector.
 *  - eps: Step size for integration.
 *  - grad: Function returning gradient of log-posterior at a position.
 *  - joint: Function returning both log-posterior and gradient at a position.
 *  - num_leapfrogs: Number of leapfrog steps to perform.
 *  - inv_mass_diag: Diagonal inverse mass matrix.
 *  - init_grad: Optional pre-computed gradient at theta_init (avoids recomputation).
 *
 * Returns:
 *  - LeapfrogJointResult containing:
 *      - theta: Final position vector.
 *      - r: Final momentum vector.
 *      - log_post: Log-posterior at final position.
 *      - grad: Gradient at final position.
 *
 * Notes:
 *  - Uses grad-only function for intermediate positions (efficiency).
 *  - Uses joint function at final position to obtain both log_post and gradient.
 */
LeapfrogJointResult leapfrog(
    const arma::vec& theta_init,
    const arma::vec& r_init,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    const arma::vec* init_grad
) {
  arma::vec r = r_init;
  arma::vec theta = theta_init;

  // Use provided initial gradient or compute it
  arma::vec grad_theta = init_grad ? *init_grad : grad(theta_init);

  // All steps except the last one
  for (int step = 0; step < num_leapfrogs - 1; step++) {
    // Half-step momentum
    r += 0.5 * eps * grad_theta;

    // Full step position
    theta += eps * (inv_mass_diag % r);

    // Update gradient (intermediate position - only need grad)
    grad_theta = grad(theta);

    // Final half-step momentum
    r += 0.5 * eps * grad_theta;
  }

  // Final step: use joint to get both log_post and gradient
  if (num_leapfrogs >= 1) {
    // Half-step momentum
    r += 0.5 * eps * grad_theta;

    // Full step position
    theta += eps * (inv_mass_diag % r);

    // Use joint at final position
    auto [log_post_final, grad_final] = joint(theta);

    // Final half-step momentum
    r += 0.5 * eps * grad_final;

    return {theta, r, log_post_final, grad_final};
  }

  // Edge case: num_leapfrogs == 0 (shouldn't happen in practice)
  auto [log_post, grad_vec] = joint(theta);
  return {theta, r, log_post, grad_vec};
}