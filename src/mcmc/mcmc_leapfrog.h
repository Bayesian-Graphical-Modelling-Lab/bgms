// mcmc_leapfrog.h
#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_memoization.h"



/**
 * Function: leapfrog_memo
 *
 * Performs a leapfrog step using a memoization wrapper to avoid redundant gradient evaluations.
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
);


/**
 * Struct: LeapfrogJointResult
 *
 * Return type for leapfrog, containing the final position, momentum,
 * log-posterior, and gradient at the final position.
 */
struct LeapfrogJointResult {
  arma::vec theta;      // Final position
  arma::vec r;          // Final momentum
  double log_post;      // Log-posterior at final position
  arma::vec grad;       // Gradient at final position
};


/**
 * Function: leapfrog
 *
 * Performs leapfrog integration using a joint log_post+gradient function.
 *
 * Key optimizations:
 *  - Accepts optional initial gradient to avoid recomputing grad(θ₀)
 *  - Uses joint function at the FINAL position to get both log_post and grad
 *  - Returns (theta_L, r_L, log_post_L, grad_L)
 *
 * This avoids redundant probability computations at endpoints when both
 * log_post and grad are needed (e.g., for Hamiltonian evaluation in HMC).
 *
 * Inputs:
 *  - theta: Initial position
 *  - r: Initial momentum
 *  - eps: Step size
 *  - grad: Gradient-only function (for intermediate steps)
 *  - joint: Joint function returning (log_post, grad) pair
 *  - num_leapfrogs: Number of leapfrog steps
 *  - inv_mass_diag: Diagonal inverse mass matrix
 *  - init_grad: Optional pre-computed gradient at theta (nullptr to compute)
 *
 * Returns:
 *  - LeapfrogJointResult with final position, momentum, log_post, and gradient
 */
LeapfrogJointResult leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    const arma::vec* init_grad = nullptr
);