#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc/mcmc_leapfrog.h"
#include "mcmc/mcmc_memoization.h"
#include "mcmc/mcmc_utils.h"
struct SafeRNG;


/**
 * Struct: BuildTreeResult
 *
 * Holds the return values of the recursive build_tree function used in the NUTS algorithm.
 * Each call to build_tree expands the sampling path and may return new candidate samples
 * or indicate when the trajectory should terminate.
 *
 * Fields:
 *  - theta_min: Leftmost position in the trajectory.
 *  - r_min: Corresponding momentum at theta_min.
 *  - theta_plus: Rightmost position in the trajectory.
 *  - r_plus: Corresponding momentum at theta_plus.
 *  - theta_prime: The current proposed sample (to possibly accept).
 *  - r_prime: Momentum corresponding to theta_prime (for energy diagnostics).
 *  - rho: Sum of momenta along the subtree trajectory (for generalized U-turn).
 *  - p_sharp_beg: Sharp momentum (M^{-1} * p) at beginning of subtree.
 *  - p_sharp_end: Sharp momentum (M^{-1} * p) at end of subtree.
 *  - p_beg: Momentum at beginning of subtree.
 *  - p_end: Momentum at end of subtree.
 *  - n_prime: Number of valid proposals from this subtree.
 *  - s_prime: Stop flag (1 = continue, 0 = stop expansion).
 *  - alpha: Sum of acceptance probabilities in the subtree.
 *  - n_alpha: Number of proposals contributing to alpha.
 *  - divergent: Is the transition diverging?
 */
struct BuildTreeResult {
  arma::vec theta_min;
  arma::vec r_min;
  arma::vec theta_plus;
  arma::vec r_plus;
  arma::vec theta_prime;
  arma::vec r_prime;
  arma::vec rho;
  arma::vec p_sharp_beg;
  arma::vec p_sharp_end;
  arma::vec p_beg;
  arma::vec p_end;
  int n_prime;
  int s_prime;
  double alpha;
  int n_alpha;
  bool divergent;
};



/**
 * Function: nuts_sampler
 *
 * Executes the No-U-Turn Sampler algorithm (NUTS).
 * Takes a joint log_post+gradient function for efficient memoization.
 *
 * The joint function computes both log-posterior and gradient together,
 * which is more efficient when they share common computations (e.g.,
 * normalization constants). The Memoizer caches both values together.
 *
 * Inputs:
 *  - init_theta: Initial position (parameter vector).
 *  - step_size: Step size for leapfrog integration.
 *  - joint: Function returning (log_post, gradient) pair.
 *  - inv_mass_diag: Diagonal inverse mass matrix.
 *  - rng: Random number generator.
 *  - max_depth: Maximum tree depth (default = 10).
 *
 * Returns:
 *  - SamplerResult with final position, acceptance probability, and diagnostics.
 */
SamplerResult nuts_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    int max_depth = 10);