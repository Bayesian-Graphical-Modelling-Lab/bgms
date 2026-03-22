#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc/algorithms/leapfrog.h"
#include "mcmc/profiler.h"


std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  BGMS_PROF_START(_t_g1);
  const arma::vec& grad1 = memo.cached_grad(theta_new);
  BGMS_PROF_RECORD("lf.grad1", _t_g1);

  r_half += 0.5 * eps * grad1;
  theta_new += eps * (inv_mass_diag % r_half);

  BGMS_PROF_START(_t_g2);
  const arma::vec& grad2 = memo.cached_grad(theta_new);
  BGMS_PROF_RECORD("lf.grad2", _t_g2);

  r_half += 0.5 * eps * grad2;

  return {theta_new, r_half};
}


std::pair<arma::vec, arma::vec> leapfrog_constrained(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    const ProjectPositionFn& project_position,
    const ProjectMomentumFn& project_momentum
) {
  arma::vec r_half = r;
  arma::vec theta_new = theta;

  // --- Step 1: Half-step momentum ---
  BGMS_PROF_START(_t1);
  const arma::vec& grad1 = memo.cached_grad(theta_new);
  BGMS_PROF_RECORD("rlf.grad1", _t1);

  BGMS_PROF_START(_t2);
  r_half += 0.5 * eps * grad1;
  BGMS_PROF_RECORD("rlf.half1", _t2);

  // --- Step 2: Project momentum onto cotangent space ---
  BGMS_PROF_START(_t2b);
  project_momentum(r_half, theta_new);
  BGMS_PROF_RECORD("rlf.proj_mom_pre", _t2b);

  // --- Step 3: Full-step position ---
  BGMS_PROF_START(_t3);
  theta_new += eps * (inv_mass_diag % r_half);
  BGMS_PROF_RECORD("rlf.pos", _t3);

  // --- Step 4: SHAKE — position-only projection ---
  BGMS_PROF_START(_t4);
  arma::vec theta_pre = theta_new;
  project_position(theta_new);
  BGMS_PROF_RECORD("rlf.proj_pos", _t4);

  // --- Step 5: Momentum correction for constraint forces ---
  BGMS_PROF_START(_t5);
  arma::vec delta_x = theta_new - theta_pre;
  r_half += delta_x / (eps * inv_mass_diag);
  BGMS_PROF_RECORD("rlf.mom_correct", _t5);

  memo.invalidate();

  // --- Step 6: Second half-step momentum ---
  BGMS_PROF_START(_t6);
  const arma::vec& grad2 = memo.cached_grad(theta_new);
  BGMS_PROF_RECORD("rlf.grad2", _t6);

  BGMS_PROF_START(_t7);
  r_half += 0.5 * eps * grad2;
  BGMS_PROF_RECORD("rlf.half2", _t7);

  // --- Step 7: Project momentum onto cotangent space ---
  BGMS_PROF_START(_t8);
  project_momentum(r_half, theta_new);
  BGMS_PROF_RECORD("rlf.proj_mom", _t8);

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
