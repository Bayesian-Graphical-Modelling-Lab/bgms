#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <utility>
#include "mcmc/execution/step_result.h"
#include "mcmc/algorithms/leapfrog.h"
struct SafeRNG;


/**
 * Kinetic energy for Hamiltonian Monte Carlo
 *
 * Computes 0.5 * r^T * M^{-1} * r where M^{-1} is a diagonal mass matrix.
 *
 * @param r              Momentum vector
 * @param inv_mass_diag  Diagonal of the inverse mass matrix
 * @return Scalar kinetic energy
 */
double kinetic_energy(const arma::vec& r, const arma::vec& inv_mass_diag);


/**
 * Heuristic initial step size for HMC/NUTS (identity mass)
 *
 * Iteratively doubles or halves a candidate step size until a single leapfrog
 * step yields an acceptance probability near the target. Delegates to the
 * mass-matrix overload with inv_mass_diag = ones.
 *
 * @param theta             Initial parameter vector
 * @param grad              Gradient function
 * @param joint             Joint log-posterior + gradient function
 * @param rng               Random number generator
 * @param target_acceptance Target acceptance probability
 * @param init_step         Starting step size
 * @param max_attempts      Maximum doubling/halving iterations
 * @return Step size yielding acceptance probability near target
 */
double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    SafeRNG& rng,
    double target_acceptance = 0.625,
    double init_step = 1.0,
    int max_attempts = 20
);


/**
 * Heuristic initial step size for HMC/NUTS (with mass matrix)
 *
 * Same algorithm as the identity-mass overload, but samples momentum from
 * N(0, M) and evaluates kinetic energy with the supplied diagonal M^{-1}.
 *
 * @param theta             Initial parameter vector
 * @param grad              Gradient function
 * @param joint             Joint log-posterior + gradient function
 * @param inv_mass_diag     Diagonal of the inverse mass matrix
 * @param rng               Random number generator
 * @param target_acceptance Target acceptance probability
 * @param init_step         Starting step size
 * @param max_attempts      Maximum doubling/halving iterations
 * @return Step size yielding acceptance probability near target
 */
double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng,
    double target_acceptance = 0.625,
    double init_step = 1.0,
    int max_attempts = 20
);


/**
 * Heuristic initial step size for constrained HMC/NUTS (RATTLE)
 *
 * Same doubling/halving algorithm as the unconstrained overload, but uses
 * leapfrog_constrained for the trial step and projects initial momentum
 * onto the cotangent space before computing kinetic energy. This finds a
 * step size appropriate for the constrained manifold geometry.
 *
 * @param theta             Initial parameter vector (on constraint manifold)
 * @param joint             Joint log-posterior + gradient function
 * @param inv_mass_diag     Diagonal of the inverse mass matrix
 * @param project_position  SHAKE position projection callback
 * @param project_momentum  RATTLE momentum projection callback
 * @param rng               Random number generator
 * @param target_acceptance Target acceptance probability
 * @param init_step         Starting step size
 * @param max_attempts      Maximum doubling/halving iterations
 * @return Step size yielding acceptance probability near target
 */
double heuristic_initial_step_size_constrained(
    const arma::vec& theta,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const arma::vec& inv_mass_diag,
    const ProjectPositionFn& project_position,
    const ProjectMomentumFn& project_momentum,
    SafeRNG& rng,
    double target_acceptance = 0.625,
    double init_step = 1.0,
    int max_attempts = 20
);


/**
 * Performs one iteration of Hamiltonian Monte Carlo sampling
 *
 * Proposes a new state by simulating Hamiltonian dynamics through leapfrog
 * integration, then accepts or rejects via the Metropolis criterion.
 * Uses joint function at endpoints for log_post+gradient, grad-only for
 * intermediate steps, avoiding redundant probability computations.
 *
 * @param init_theta     Initial parameter vector (position)
 * @param step_size      Leapfrog integration step size (epsilon)
 * @param grad           Gradient function
 * @param joint          Joint log-posterior + gradient function
 * @param num_leapfrogs  Number of leapfrog steps per proposal
 * @param inv_mass_diag  Diagonal of the inverse mass matrix
 * @param rng            Thread-safe random number generator
 * @return StepResult with accepted state and acceptance probability
 */
StepResult hmc_step(
    const arma::vec& init_theta,
    double step_size,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng
);


/**
 * Constrained HMC step using RATTLE integration
 *
 * Same Metropolis accept/reject as the unconstrained overload, but uses
 * leapfrog_constrained (SHAKE + RATTLE projections) at each step to
 * keep the trajectory on the constraint manifold.
 *
 * @param init_theta        Initial parameter vector (on constraint manifold)
 * @param step_size         Leapfrog integration step size (epsilon)
 * @param joint             Joint log-posterior + gradient function
 * @param num_leapfrogs     Number of leapfrog steps per proposal
 * @param inv_mass_diag     Diagonal of the inverse mass matrix
 * @param project_position  SHAKE position projection callback
 * @param project_momentum  RATTLE momentum projection callback
 * @param rng               Thread-safe random number generator
 * @param reverse_check     Enable runtime reversibility check
 * @param reverse_check_tol Factor for eps²-scaled reversibility tolerance
 * @return StepResult with accepted state and acceptance probability
 */
StepResult hmc_step(
    const arma::vec& init_theta,
    double step_size,
    const std::function<std::pair<double, arma::vec>(const arma::vec&)>& joint,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    const ProjectPositionFn& project_position,
    const ProjectMomentumFn& project_momentum,
    SafeRNG& rng,
    bool reverse_check = true,
    double reverse_check_tol = 0.5
);
