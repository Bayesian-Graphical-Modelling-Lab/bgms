#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/leapfrog.h"
#include "mcmc/memoization.h"
#include "mcmc/nuts.h"
#include "mcmc/hamiltonian.h"
#include "rng/rng_utils.h"


/**
 * The generalized U-turn criterion used here is described in Betancourt (2017).
 * The implementation follows the approach in STAN's base_nuts.hpp (BSD-3-Clause license).
 *
 * References:
 *   Betancourt, M. (2017). A Conceptual Introduction to Hamiltonian Monte Carlo.
 *     arXiv preprint arXiv:1701.02434.
 *   Stan Development Team. base_nuts.hpp.
 *     https://github.com/stan-dev/stan/blob/develop/src/stan/mcmc/hmc/nuts/base_nuts.hpp
 */


/**
 * Computes the generalized U-turn criterion for the NUTS algorithm
 *
 * @param p_sharp_minus  Sharp momentum (M^{-1} p) at backward end
 * @param p_sharp_plus   Sharp momentum (M^{-1} p) at forward end
 * @param rho            Sum of momenta along the trajectory
 * @return true if criterion satisfied (continue), false if U-turn detected (stop)
 */
bool compute_criterion(const arma::vec& p_sharp_minus,
                       const arma::vec& p_sharp_plus,
                       const arma::vec& rho) {
  return arma::dot(p_sharp_plus, rho) > 0 && arma::dot(p_sharp_minus, rho) > 0;
}



/**
 * Recursively builds a binary tree of leapfrog steps in the NUTS algorithm
 *
 * Explores forward or backward in time, evaluating trajectory termination
 * criteria. Based on Algorithm 6 in Hoffman & Gelman (2014).
 *
 * @param theta      Current position at the base of the tree
 * @param r          Current momentum at the base of the tree
 * @param log_u      Log slice variable for accept/reject decision
 * @param v          Direction of expansion (-1 backward, +1 forward)
 * @param j          Current tree depth
 * @param step_size  Step size used in leapfrog integration
 * @param theta_0    Initial position at the start of sampling
 * @param r0         Initial momentum at the start of sampling
 * @param logp0      Log posterior at initial position
 * @param kin0       Kinetic energy at initial momentum
 * @param memo       Memoizer object for caching evaluations
 * @return BuildTreeResult with updated endpoints, candidate sample, and diagnostics
 */
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
    SafeRNG& rng
) {
  constexpr double Delta_max = 1000.0;

  if (j == 0) {
    // Base case: take a single leapfrog step
    arma::vec theta_new, r_new;
    std::tie(theta_new, r_new) = leapfrog_memo(
      theta, r, v * step_size, memo, inv_mass_diag
    );

    auto logp = memo.cached_log_post(theta_new);
    double kin = kinetic_energy(r_new, inv_mass_diag);
    int n_new = 1 * (log_u <= logp - kin);
    int s_new = 1 * (log_u <= Delta_max + logp - kin);
    bool divergent = (s_new == 0);
    double alpha = std::min(1.0, MY_EXP(logp - kin - logp0 + kin0));

    // Initialize rho with the momentum at this point
    arma::vec rho = r_new;
    // Sharp momentum (velocity): M^{-1} * p
    arma::vec p_sharp = inv_mass_diag % r_new;

    BuildTreeResult result;
    result.theta_min = theta_new;
    result.r_min = r_new;
    result.theta_plus = theta_new;
    result.r_plus = r_new;
    result.theta_prime = theta_new;
    result.r_prime = r_new;
    result.rho = rho;
    result.p_sharp_beg = p_sharp;
    result.p_sharp_end = p_sharp;
    result.p_beg = r_new;
    result.p_end = r_new;
    result.n_prime = n_new;
    result.s_prime = s_new;
    result.alpha = alpha;
    result.n_alpha = 1;
    result.divergent = divergent;
    return result;

  } else {
    // Recursion: build the first subtree
    BuildTreeResult init_result = build_tree(
      theta, r, log_u, v, j - 1, step_size, theta_0, r0, logp0, kin0, memo,
      inv_mass_diag, rng
    );

    if (init_result.s_prime == 0) {
      // First subtree is invalid, return early
      return init_result;
    }

    bool divergent = init_result.divergent;

    // Extract values from init subtree
    arma::vec theta_min = init_result.theta_min;
    arma::vec r_min = init_result.r_min;
    arma::vec theta_plus = init_result.theta_plus;
    arma::vec r_plus = init_result.r_plus;
    arma::vec theta_prime = init_result.theta_prime;
    arma::vec r_prime = init_result.r_prime;
    arma::vec rho_init = init_result.rho;
    arma::vec p_sharp_init_beg = init_result.p_sharp_beg;
    arma::vec p_sharp_init_end = init_result.p_sharp_end;
    arma::vec p_init_beg = init_result.p_beg;
    arma::vec p_init_end = init_result.p_end;
    int n_prime = init_result.n_prime;
    double alpha_prime = init_result.alpha;
    int n_alpha_prime = init_result.n_alpha;

    // Build the second subtree in the same direction
    BuildTreeResult final_result;
    if (v == -1) {
      final_result = build_tree(
        theta_min, r_min, log_u, v, j - 1, step_size, theta_0, r0, logp0,
        kin0, memo, inv_mass_diag, rng
      );
      // Update backward boundary
      theta_min = final_result.theta_min;
      r_min = final_result.r_min;
    } else {
      final_result = build_tree(
        theta_plus, r_plus, log_u, v, j - 1, step_size, theta_0, r0, logp0,
        kin0, memo, inv_mass_diag, rng
      );
      // Update forward boundary
      theta_plus = final_result.theta_plus;
      r_plus = final_result.r_plus;
    }

    if (final_result.s_prime == 0) {
      // Second subtree is invalid - return early with s_prime=0
      // p_sharp/p_beg/p_end values follow same convention as valid case
      // (init_beg and final_end), even though they won't be used for criterion checks
      BuildTreeResult result;
      result.theta_min = theta_min;
      result.r_min = r_min;
      result.theta_plus = theta_plus;
      result.r_plus = r_plus;
      result.theta_prime = theta_prime;
      result.r_prime = r_prime;
      result.rho = rho_init + final_result.rho;
      result.p_sharp_beg = p_sharp_init_beg;
      result.p_sharp_end = final_result.p_sharp_end;
      result.p_beg = p_init_beg;
      result.p_end = final_result.p_end;
      result.n_prime = n_prime + final_result.n_prime;
      result.s_prime = 0;
      result.alpha = alpha_prime + final_result.alpha;
      result.n_alpha = n_alpha_prime + final_result.n_alpha;
      result.divergent = divergent || final_result.divergent;
      return result;
    }

    // Extract values from final subtree
    arma::vec rho_final = final_result.rho;
    arma::vec p_sharp_final_beg = final_result.p_sharp_beg;
    arma::vec p_sharp_final_end = final_result.p_sharp_end;
    arma::vec p_final_beg = final_result.p_beg;
    arma::vec p_final_end = final_result.p_end;
    int n_double_prime = final_result.n_prime;
    double alpha_double_prime = final_result.alpha;
    int n_alpha_double_prime = final_result.n_alpha;
    divergent = divergent || final_result.divergent;

    // Multinomial sampling from the combined subtree
    double denom = static_cast<double>(n_prime + n_double_prime);
    double prob = static_cast<double>(n_double_prime) / denom;

    if (runif(rng) < prob) {
      theta_prime = final_result.theta_prime;
      r_prime = final_result.r_prime;
    }

    alpha_prime += alpha_double_prime;
    n_alpha_prime += n_alpha_double_prime;
    n_prime += n_double_prime;

    // Combine rho from both subtrees
    arma::vec rho_subtree = rho_init + rho_final;

    // Determine the sharp momenta at the boundaries of the combined subtree
    // Following STAN convention: "beg" = first visited (in build direction), "end" = last visited
    // This is the same regardless of direction - init subtree is always first, final is last
    arma::vec p_sharp_beg = p_sharp_init_beg;    // First visited in combined tree
    arma::vec p_sharp_end = final_result.p_sharp_end;  // Last visited in combined tree
    arma::vec p_beg = p_init_beg;
    arma::vec p_end = final_result.p_end;

    // However, theta_min/theta_plus track ABSOLUTE positions (backward/forward boundaries)
    // These DO depend on direction:
    // - For v=-1: we're extending backward, so theta_min (backward boundary) gets updated
    // - For v=+1: we're extending forward, so theta_plus (forward boundary) gets updated

    // Generalized U-turn criterion (three checks like STAN)
    // The tree structure is always: init_subtree -> final_subtree (in direction v)
    // So the junction is always between p_init_end and p_final_beg

    // 1. Check criterion around merged subtrees
    bool persist_criterion = compute_criterion(p_sharp_beg, p_sharp_end, rho_subtree);

    // 2. Check between subtrees: from start of tree to start of final subtree
    arma::vec rho_extended = rho_init + p_final_beg;
    persist_criterion = persist_criterion &&
      compute_criterion(p_sharp_beg, p_sharp_final_beg, rho_extended);

    // 3. Check between subtrees: from end of init subtree to end of tree
    rho_extended = rho_final + p_init_end;
    persist_criterion = persist_criterion &&
      compute_criterion(p_sharp_init_end, p_sharp_end, rho_extended);

    int s_prime = persist_criterion ? 1 : 0;

    BuildTreeResult result;
    result.theta_min = theta_min;
    result.r_min = r_min;
    result.theta_plus = theta_plus;
    result.r_plus = r_plus;
    result.theta_prime = theta_prime;
    result.r_prime = r_prime;
    result.rho = rho_subtree;
    result.p_sharp_beg = p_sharp_beg;
    result.p_sharp_end = p_sharp_end;
    result.p_beg = p_beg;
    result.p_end = p_end;
    result.n_prime = n_prime;
    result.s_prime = s_prime;
    result.alpha = alpha_prime;
    result.n_alpha = n_alpha_prime;
    result.divergent = divergent;
    return result;
  }
}


SamplerResult nuts_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    int max_depth
) {
  // Create Memoizer with joint function
  Memoizer memo(joint);
  bool any_divergence = false;

  arma::vec r0 = arma::sqrt(1.0 / inv_mass_diag) % arma_rnorm_vec(rng, init_theta.n_elem);
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
        inv_mass_diag, rng
      );
      theta_min = result.theta_min;
      r_min = result.r_min;
      rho_bck = result.rho;
      p_sharp_bck_bck = result.p_sharp_beg;
      p_bck_fwd = result.p_end;
      p_sharp_bck_fwd = result.p_sharp_end;
    } else {
      rho_bck = rho;
      result = build_tree(
        theta_plus, r_plus, log_u, v, j, step_size, init_theta, r0, logp0, kin0, memo,
        inv_mass_diag, rng
      );
      theta_plus = result.theta_plus;
      r_plus = result.r_plus;
      rho_fwd = result.rho;
      p_sharp_fwd_fwd = result.p_sharp_end;
      p_fwd_bck = result.p_beg;
      p_sharp_fwd_bck = result.p_sharp_beg;
    }

    any_divergence = any_divergence || result.divergent;
    alpha = result.alpha;
    n_alpha = result.n_alpha;

    if (result.s_prime == 1) {
      double prob = static_cast<double>(result.n_prime) / static_cast<double>(n);
      if (runif(rng) < prob) {
        theta = result.theta_prime;
        r = result.r_prime;
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
  diag->energy = energy;

  return {theta, accept_prob, diag};
}