#include <RcppArmadillo.h>
#include <algorithm>
#include <functional>
#include <memory>
#include <utility>
#include "mcmc/algorithms/leapfrog.h"
#include "mcmc/algorithms/nuts.h"
#include "mcmc/algorithms/hmc.h"
#include "rng/rng_utils.h"


// The generalized U-turn criterion used here is described in Betancourt (2017).
// The implementation follows the approach in STAN's base_nuts.hpp (BSD-3-Clause license).
//
// References:
//   Betancourt, M. (2017). A Conceptual Introduction to Hamiltonian Monte Carlo.
//     arXiv preprint arXiv:1701.02434.
//   Stan Development Team. base_nuts.hpp.
//     https://github.com/stan-dev/stan/blob/develop/src/stan/mcmc/hmc/nuts/base_nuts.hpp


// Computes the generalized U-turn criterion for the NUTS algorithm
//
// @param p_sharp_minus  Sharp momentum (M^{-1} p) at backward end
// @param p_sharp_plus   Sharp momentum (M^{-1} p) at forward end
// @param rho            Sum of momenta along the trajectory
// @return true if criterion satisfied (continue), false if U-turn detected (stop)
bool compute_criterion(const arma::vec& p_sharp_minus,
                       const arma::vec& p_sharp_plus,
                       const arma::vec& rho) {
  return arma::dot(p_sharp_plus, rho) > 0 && arma::dot(p_sharp_minus, rho) > 0;
}



// Recursively builds a binary tree of leapfrog steps in the NUTS algorithm
//
// Explores forward or backward in time, evaluating trajectory termination
// criteria. Based on Algorithm 6 in Hoffman & Gelman (2014).
//
// @param theta      Current position at the base of the tree
// @param r          Current momentum at the base of the tree
// @param log_u      Log slice variable for accept/reject decision
// @param v          Direction of expansion (-1 backward, +1 forward)
// @param j          Current tree depth
// @param step_size  Step size used in leapfrog integration
// @param theta_0    Initial position at the start of sampling
// @param r0         Initial momentum at the start of sampling
// @param logp0      Log posterior at initial position
// @param kin0       Kinetic energy at initial momentum
// @param memo       Memoizer object for caching evaluations
// @return BuildTreeResult with updated endpoints, candidate sample, and diagnostics
BuildTreeResult build_tree(
    const arma::vec& theta,
    const arma::vec& r,
    double log_u,
    int v,
    int j,
    double step_size,
    const arma::vec& theta_0,
    const arma::vec& r0,
    const double logp0,
    const double kin0,
    Memoizer& memo,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    const ProjectPositionFn* project_position,
    const ProjectMomentumFn* project_momentum,
    bool reverse_check,
    double reverse_check_tol
) {
  constexpr double Delta_max = 1000.0;

  if (j == 0) {
    // Base case: take a single leapfrog step
    arma::vec theta_new, r_new;
    bool non_reversible = false;
    if (project_position && project_momentum) {
      // Always run the checked variant so we can observe reversibility.
      // The reverse_check flag controls whether we ACT on the result
      // (i.e. terminate the tree). Observation is always on.
      auto checked = leapfrog_constrained_checked(
        theta, r, v * step_size, memo, inv_mass_diag,
        *project_position, *project_momentum,
        reverse_check_tol
      );
      theta_new = std::move(checked.theta);
      r_new = std::move(checked.r);
      non_reversible = !checked.reversible;
    } else {
      std::tie(theta_new, r_new) = leapfrog_memo(
        theta, r, v * step_size, memo, inv_mass_diag
      );
    }

    // Non-reversible step terminates the tree only when reverse_check is on.
    // During warmup, we record but don't act.
    if (reverse_check && non_reversible) {
      arma::vec p_sharp = inv_mass_diag % r_new;
      BuildTreeResult result;
      result.theta_min = theta_new;
      result.theta_plus = theta_new;
      result.r_min = r_new;
      result.r_plus = r_new;
      result.rho = r_new;
      result.p_beg = r_new;
      result.p_end = r_new;
      result.r_prime = std::move(r_new);
      result.theta_prime = std::move(theta_new);
      result.p_sharp_beg = p_sharp;
      result.p_sharp_end = std::move(p_sharp);
      result.n_prime = 0;
      result.s_prime = 0;
      result.alpha = 0.0;
      result.n_alpha = 1;
      result.divergent = false;
      result.non_reversible = true;
      return result;
    }

    auto logp = memo.cached_log_post(theta_new);
    double kin = kinetic_energy(r_new, inv_mass_diag);
    int n_new = 1 * (log_u <= logp - kin);
    int s_new = 1 * (log_u <= Delta_max + logp - kin);
    bool divergent = (s_new == 0);
    double alpha = std::min(1.0, MY_EXP(logp - kin - logp0 + kin0));

    // Sharp momentum (velocity): M^{-1} * p
    arma::vec p_sharp = inv_mass_diag % r_new;

    BuildTreeResult result;
    result.theta_min = theta_new;
    result.theta_plus = theta_new;
    result.r_min = r_new;
    result.r_plus = r_new;
    result.rho = r_new;
    result.p_beg = r_new;
    result.p_end = r_new;
    result.r_prime = std::move(r_new);
    result.theta_prime = std::move(theta_new);
    result.p_sharp_beg = p_sharp;
    result.p_sharp_end = std::move(p_sharp);
    result.n_prime = n_new;
    result.s_prime = s_new;
    result.alpha = alpha;
    result.n_alpha = 1;
    result.divergent = divergent;
    result.non_reversible = non_reversible;  // record even when not acting
    return result;

  } else {
    // Recursion: build the first subtree
    BuildTreeResult init_result = build_tree(
      theta, r, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0, memo,
      inv_mass_diag, rng, project_position, project_momentum, reverse_check,
      reverse_check_tol
    );

    if (init_result.s_prime == 0) {
      // First subtree is invalid, return early
      return init_result;
    }

    bool divergent = init_result.divergent;
    bool non_reversible = init_result.non_reversible;

    // Extract values from init subtree (move — init_result not used again)
    arma::vec theta_min = std::move(init_result.theta_min);
    arma::vec r_min = std::move(init_result.r_min);
    arma::vec theta_plus = std::move(init_result.theta_plus);
    arma::vec r_plus = std::move(init_result.r_plus);
    arma::vec theta_prime = std::move(init_result.theta_prime);
    arma::vec r_prime = std::move(init_result.r_prime);
    arma::vec rho_init = std::move(init_result.rho);
    arma::vec p_sharp_init_beg = std::move(init_result.p_sharp_beg);
    arma::vec p_sharp_init_end = std::move(init_result.p_sharp_end);
    arma::vec p_init_beg = std::move(init_result.p_beg);
    arma::vec p_init_end = std::move(init_result.p_end);
    int n_prime = init_result.n_prime;
    double alpha_prime = init_result.alpha;
    int n_alpha_prime = init_result.n_alpha;

    // Build the second subtree in the same direction
    BuildTreeResult final_result;
    if (v == -1) {
      final_result = build_tree(
        theta_min, r_min, log_u, v, j - 1, step_size, theta_0, r0, logp0,
        kin0, memo, inv_mass_diag, rng, project_position, project_momentum,
        reverse_check, reverse_check_tol
      );
      // Update backward boundary
      theta_min = std::move(final_result.theta_min);
      r_min = std::move(final_result.r_min);
    } else {
      final_result = build_tree(
        theta_plus, r_plus, log_u, v, j - 1, step_size, theta_0, r0, logp0,
        kin0, memo, inv_mass_diag, rng, project_position, project_momentum,
        reverse_check, reverse_check_tol
      );
      // Update forward boundary
      theta_plus = std::move(final_result.theta_plus);
      r_plus = std::move(final_result.r_plus);
    }

    if (final_result.s_prime == 0) {
      // Second subtree is invalid - return early with s_prime=0
      BuildTreeResult result;
      result.theta_min = std::move(theta_min);
      result.r_min = std::move(r_min);
      result.theta_plus = std::move(theta_plus);
      result.r_plus = std::move(r_plus);
      result.theta_prime = std::move(theta_prime);
      result.r_prime = std::move(r_prime);
      rho_init += final_result.rho;
      result.rho = std::move(rho_init);
      result.p_sharp_beg = std::move(p_sharp_init_beg);
      result.p_sharp_end = std::move(final_result.p_sharp_end);
      result.p_beg = std::move(p_init_beg);
      result.p_end = std::move(final_result.p_end);
      result.n_prime = n_prime + final_result.n_prime;
      result.s_prime = 0;
      result.alpha = alpha_prime + final_result.alpha;
      result.n_alpha = n_alpha_prime + final_result.n_alpha;
      result.divergent = divergent || final_result.divergent;
      result.non_reversible = non_reversible || final_result.non_reversible;
      return result;
    }

    // Extract values from final subtree (move — final_result not used again)
    arma::vec rho_final = std::move(final_result.rho);
    arma::vec p_sharp_final_beg = std::move(final_result.p_sharp_beg);
    arma::vec p_sharp_final_end = std::move(final_result.p_sharp_end);
    arma::vec p_final_beg = std::move(final_result.p_beg);
    arma::vec p_final_end = std::move(final_result.p_end);
    int n_double_prime = final_result.n_prime;
    double alpha_double_prime = final_result.alpha;
    int n_alpha_double_prime = final_result.n_alpha;
    divergent = divergent || final_result.divergent;
    non_reversible = non_reversible || final_result.non_reversible;

    // Multinomial sampling from the combined subtree
    double denom = static_cast<double>(n_prime + n_double_prime);
    double prob = static_cast<double>(n_double_prime) / denom;

    if (runif(rng) < prob) {
      theta_prime = std::move(final_result.theta_prime);
      r_prime = std::move(final_result.r_prime);
    }

    alpha_prime += alpha_double_prime;
    n_alpha_prime += n_alpha_double_prime;
    n_prime += n_double_prime;

    // Combine rho from both subtrees
    arma::vec rho_subtree = rho_init + rho_final;

    // Generalized U-turn criterion (three checks like STAN)

    // 1. Check criterion around merged subtrees
    bool persist_criterion = compute_criterion(p_sharp_init_beg, p_sharp_final_end, rho_subtree);

    // 2. Check between subtrees: from start of tree to start of final subtree
    arma::vec rho_extended = rho_init + p_final_beg;
    persist_criterion = persist_criterion &&
      compute_criterion(p_sharp_init_beg, p_sharp_final_beg, rho_extended);

    // 3. Check between subtrees: from end of init subtree to end of tree
    rho_extended = rho_final + p_init_end;
    persist_criterion = persist_criterion &&
      compute_criterion(p_sharp_init_end, p_sharp_final_end, rho_extended);

    int s_prime = persist_criterion ? 1 : 0;

    BuildTreeResult result;
    result.theta_min = std::move(theta_min);
    result.r_min = std::move(r_min);
    result.theta_plus = std::move(theta_plus);
    result.r_plus = std::move(r_plus);
    result.theta_prime = std::move(theta_prime);
    result.r_prime = std::move(r_prime);
    result.rho = std::move(rho_subtree);
    result.p_sharp_beg = std::move(p_sharp_init_beg);
    result.p_sharp_end = std::move(p_sharp_final_end);
    result.p_beg = std::move(p_init_beg);
    result.p_end = std::move(p_final_end);
    result.n_prime = n_prime;
    result.s_prime = s_prime;
    result.alpha = alpha_prime;
    result.n_alpha = n_alpha_prime;
    result.divergent = divergent;
    result.non_reversible = non_reversible;
    return result;
  }
}


StepResult nuts_step(
    const arma::vec& init_theta,
    double step_size,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    int max_depth,
    const ProjectPositionFn* project_position,
    const ProjectMomentumFn* project_momentum,
    bool reverse_check,
    double reverse_check_tol
) {
  // Create Memoizer with joint function
  Memoizer memo(joint);
  bool any_divergence = false;
  bool any_non_reversible = false;

  arma::vec r0 = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, init_theta.n_elem);

  // Project initial momentum onto cotangent space (momentum-only)
  if (project_momentum) {
    (*project_momentum)(r0, init_theta);
  }

  auto logp0 = memo.cached_log_post(init_theta);
  double kin0 = kinetic_energy(r0, inv_mass_diag);
  double joint0 = logp0 - kin0;
  double log_u = log(runif(rng)) + joint0;

  arma::vec theta_min = init_theta, r_min = r0;
  arma::vec theta_plus = init_theta, r_plus = r0;
  arma::vec theta = init_theta;
  arma::vec r = r0;

  arma::vec p_sharp_bck_bck = inv_mass_diag % r0;
  arma::vec p_sharp_fwd_fwd = p_sharp_bck_bck;
  arma::vec p_fwd_bck = r0;
  arma::vec p_sharp_fwd_bck = p_sharp_bck_bck;
  arma::vec p_bck_fwd = r0;
  arma::vec p_sharp_bck_fwd = p_sharp_bck_bck;
  arma::vec rho = r0;

  int j = 0;
  int n = 1, s = 1;
  double alpha = 0.5;
  int n_alpha = 1;

  while (s == 1 && j < max_depth) {
    int v = runif(rng) < 0.5 ? -1 : 1;
    arma::vec rho_fwd, rho_bck;

    BuildTreeResult result;
    if (v == -1) {
      rho_fwd = rho;
      result = build_tree(
        theta_min, r_min, log_u, v, j, step_size, init_theta, r0, logp0, kin0, memo,
        inv_mass_diag, rng, project_position, project_momentum, reverse_check,
        reverse_check_tol
      );
      theta_min = std::move(result.theta_min);
      r_min = std::move(result.r_min);
      rho_bck = std::move(result.rho);
      p_sharp_bck_bck = std::move(result.p_sharp_beg);
      p_bck_fwd = std::move(result.p_end);
      p_sharp_bck_fwd = std::move(result.p_sharp_end);
    } else {
      rho_bck = rho;
      result = build_tree(
        theta_plus, r_plus, log_u, v, j, step_size, init_theta, r0, logp0, kin0, memo,
        inv_mass_diag, rng, project_position, project_momentum, reverse_check,
        reverse_check_tol
      );
      theta_plus = std::move(result.theta_plus);
      r_plus = std::move(result.r_plus);
      rho_fwd = std::move(result.rho);
      p_sharp_fwd_fwd = std::move(result.p_sharp_end);
      p_fwd_bck = std::move(result.p_beg);
      p_sharp_fwd_bck = std::move(result.p_sharp_beg);
    }

    any_divergence = any_divergence || result.divergent;
    any_non_reversible = any_non_reversible || result.non_reversible;
    alpha = result.alpha;
    n_alpha = result.n_alpha;

    if (result.s_prime == 1) {
      double prob = static_cast<double>(result.n_prime) / static_cast<double>(n);
      if (runif(rng) < prob) {
        theta = std::move(result.theta_prime);
        r = std::move(result.r_prime);
      }
    }

    rho = rho_bck + rho_fwd;
    bool persist_criterion = true;

    if (result.s_prime == 1) {
      persist_criterion = compute_criterion(p_sharp_bck_bck, p_sharp_fwd_fwd, rho);
      arma::vec rho_extended = rho_bck + p_fwd_bck;
      persist_criterion = persist_criterion &&
        compute_criterion(p_sharp_bck_bck, p_sharp_fwd_bck, rho_extended);
      rho_extended = rho_fwd + p_bck_fwd;
      persist_criterion = persist_criterion &&
        compute_criterion(p_sharp_bck_fwd, p_sharp_fwd_fwd, rho_extended);
    }

    s = result.s_prime * (persist_criterion ? 1 : 0);
    n += result.n_prime;
    j++;
  }

  double accept_prob = alpha / static_cast<double>(n_alpha);
  auto logp_final = memo.cached_log_post(theta);
  double kin_final = kinetic_energy(r, inv_mass_diag);
  double energy = -logp_final + kin_final;

  auto diag = std::make_shared<NUTSDiagnostics>();
  diag->tree_depth = j;
  diag->divergent = any_divergence;
  diag->non_reversible = any_non_reversible;
  diag->energy = energy;

  return {theta, accept_prob, diag};
}
