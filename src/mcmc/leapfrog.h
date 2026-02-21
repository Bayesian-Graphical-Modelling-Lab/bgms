#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/memoization.h"



/**
 * Performs a single leapfrog step with memoized gradient evaluation
 *
 * @param theta          Current position (parameter vector)
 * @param r              Current momentum vector
 * @param eps            Step size for integration
 * @param memo           Memoizer caching gradient evaluations
 * @param inv_mass_diag  Diagonal of the inverse mass matrix
 * @return Pair of (updated position, updated momentum)
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
);


/**
 * LeapfrogJointResult - Return type for leapfrog integration
 *
 * Contains the final position, momentum, log-posterior, and gradient.
 */
struct LeapfrogJointResult {
  arma::vec theta;      ///< Final position
  arma::vec r;          ///< Final momentum
  double log_post;      ///< Log-posterior at final position
  arma::vec grad;       ///< Gradient at final position
};


/**
 * Function: leapfrog
 *
 * Performs leapfrog integration using a joint log_post+gradient function.
 *
 * Uses joint function at final position for both log_post and gradient,
 * grad-only at intermediate steps. Accepts optional pre-computed initial
 * gradient to avoid recomputation.
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