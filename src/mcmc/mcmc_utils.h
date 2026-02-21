#pragma once

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include "math/explog_switch.h"
struct SafeRNG;


/**
 * DiagnosticsBase - Polymorphic base for per-iteration sampler diagnostics
 *
 * SamplerResult holds a shared_ptr<DiagnosticsBase> so the same return type
 * works for every sampler. Samplers that collect diagnostics (currently only
 * NUTS) define a derived struct; consumers downcast with dynamic_pointer_cast.
 * Samplers without diagnostics (MH, HMC) leave the pointer null.
 */
struct DiagnosticsBase {
  virtual ~DiagnosticsBase() = default;
};



/**
 * NUTSDiagnostics - Per-iteration NUTS diagnostics (derives from DiagnosticsBase)
 */
struct NUTSDiagnostics : public DiagnosticsBase {
  int tree_depth;    ///< Depth of the trajectory tree
  bool divergent;    ///< Whether a divergence occurred
  double energy;     ///< Final Hamiltonian (-log posterior + kinetic energy)
};



/**
 * SamplerResult - Outcome of one MCMC iteration
 *
 * Returned by BaseSampler::step(). The diagnostics pointer is non-null only
 * for NUTS (holds NUTSDiagnostics); MH and HMC leave it null.
 */
struct SamplerResult {
  arma::vec state;                              ///< Accepted parameter vector
  double accept_prob;                           ///< Acceptance probability
  std::shared_ptr<DiagnosticsBase> diagnostics; ///< NUTS diagnostics, or null
};



/**
 * Robbins-Monro update for MH proposal standard deviations
 *
 * Adjusts the proposal SD toward a target acceptance rate:
 *   sd += (observed_acceptance - target) * weight
 * The result is clamped to [0.001, 2.0]. NaN values are reset to 1.0.
 *
 * @param current_sd                          Current proposal standard deviation
 * @param observed_log_acceptance_probability Log acceptance probability from MH step
 * @param rm_weight                           Robbins-Monro weight (e.g. iteration^{-0.75})
 * @param target_acceptance                   Target acceptance rate
 * @return Updated proposal SD, clamped to [0.001, 2.0]
 */
inline double update_proposal_sd_with_robbins_monro (
    const double current_sd,
    const double observed_log_acceptance_probability,
    const double rm_weight,
    const double target_acceptance
) {
  constexpr double rm_lower_bound = 0.001;
  constexpr double rm_upper_bound = 2.0;

  double observed_acceptance_probability = 1.0;
  if (observed_log_acceptance_probability < 0.0) {
    observed_acceptance_probability = MY_EXP (observed_log_acceptance_probability);
  }

  double updated_sd = current_sd +
    (observed_acceptance_probability - target_acceptance) * rm_weight;

  if (std::isnan (updated_sd)) {
    updated_sd = 1.0;
  }

  return std::clamp (updated_sd, rm_lower_bound, rm_upper_bound);
}



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