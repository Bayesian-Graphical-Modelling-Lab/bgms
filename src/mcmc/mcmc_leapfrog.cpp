#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_leapfrog.h"
#include "mcmc/mcmc_memoization.h"



 /**
 * Function: leapfrog
 *
 * Performs a num_leapfrog leapfrog integration steps using the standard gradient function.
 * Used to simulate Hamiltonian dynamics in HMC/NUTS.
 *
 * Inputs:
 *  - theta: Current position (parameter vector).
 *  - r: Current momentum vector.
 *  - eps: Step size for integration.
 *  - grad: Gradient function of the log posterior.
 *  - num_leapfrogs: Number of leapfrog steps
 *
 * Returns:
 *  - A pair containing:
 *      - Updated position vector.
 *      - Updated momentum vector.
 */
std::pair<arma::vec, arma::vec> leapfrog(
    const arma::vec& theta_init,
    const arma::vec& r_init,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag
) {
  arma::vec r = r_init;
  arma::vec theta = theta_init;
  arma::vec grad_theta = grad(theta_init);

  for(int step = 0; step < num_leapfrogs; step++) {
    // Half-step momentum
    r += 0.5 * eps * grad_theta;

    // Full step position
    theta += eps * (inv_mass_diag % r);

    // Update gradient
    grad_theta = grad(theta);

    // Final half-step momentum
    r += 0.5 * eps * grad_theta;
  }

  return {theta, r};
}



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
  if (grad1.n_elem != theta_new.n_elem) {
      Rcpp::Rcout << "LEAPFROG: grad1.n_elem=" << grad1.n_elem
                  << " theta.n_elem=" << theta_new.n_elem << std::endl;
  }
  r_half += 0.5 * eps * grad1;

  if (inv_mass_diag.n_elem != r_half.n_elem) {
      Rcpp::Rcout << "LEAPFROG: inv_mass_diag.n_elem=" << inv_mass_diag.n_elem
                  << " r_half.n_elem=" << r_half.n_elem << std::endl;
  }
  theta_new += eps * (inv_mass_diag % r_half);
  auto grad2 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}