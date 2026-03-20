#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc/algorithms/leapfrog.h"


std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  const arma::vec& grad1 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad1;

  theta_new += eps * (inv_mass_diag % r_half);
  const arma::vec& grad2 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}


std::pair<arma::vec, arma::vec> leapfrog_constrained(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    const ProjectFn& project
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  // Half-step momentum
  const arma::vec& grad1 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad1;

  // Full-step position
  theta_new += eps * (inv_mass_diag % r_half);

  // --- RATTLE position-constraint step ---
  // Save pre-projection position to compute the correction
  arma::vec theta_pre = theta_new;

  // Project position only (discard the momentum projection from this call)
  arma::vec r_temp = r_half;
  project(theta_new, r_temp);

  // RATTLE momentum correction: couples position and momentum updates
  // In standard RATTLE, the Lagrange multiplier from the position constraint
  // also appears in the momentum half-step. The correction is:
  //   Δr = M · Δx / ε = Δx / (ε · M^{-1})
  arma::vec delta_x = theta_new - theta_pre;
  r_half += delta_x / (eps * inv_mass_diag);

  // Note: mid-step momentum is NOT projected to the tangent space.
  // Only the final momentum gets the velocity-constraint projection.

  memo.invalidate();

  // Second half-step momentum (re-evaluates gradient at projected position)
  const arma::vec& grad2 = memo.cached_grad(theta_new);
  r_half += 0.5 * eps * grad2;

  // --- RATTLE velocity-constraint step ---
  // Project final momentum onto cotangent space
  // (position projection is a no-op since theta_new is already on manifold)
  project(theta_new, r_half);

  return {theta_new, r_half};
}


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
