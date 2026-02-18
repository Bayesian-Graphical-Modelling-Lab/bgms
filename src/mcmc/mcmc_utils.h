#pragma once

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include "math/explog_switch.h"
struct SafeRNG;

// (only if <algorithm> didn’t already provide it under C++17)
#if __cplusplus < 201703L
namespace std {
  template <class T>
  const T& clamp(const T& v, const T& lo, const T& hi) {
    return v < lo ? lo : hi < v ? hi : v;
  }
}
#endif



// -----------------------------------------------------------------------------
// MCMC Output Structures
// -----------------------------------------------------------------------------



/**
 * Struct: DiagnosticsBase
 *
 * Abstract base class for sampler-specific diagnostics.
 *
 * Allows storing runtime sampler diagnostics (e.g., tree depth, divergence) in
 * a polymorphic, type-safe way. Each sampler defines its own derived struct
 * inheriting from this base.
 *
 * Notes:
 *  - Must be used via std::shared_ptr.
 *  - Enables dynamic_pointer_cast for safe access to sampler-specific fields.
 *  - The virtual destructor ensures proper cleanup.
 */
struct DiagnosticsBase {
  virtual ~DiagnosticsBase() = default;
};



/**
 * Struct: NUTSDiagnostics
 *
 * Diagnostics collected during one iteration of the No-U-Turn Sampler (NUTS).
 *
 * Fields:
 *  - tree_depth: Depth of the final trajectory tree used during this iteration.
 *  - divergent: Whether a divergence occurred during trajectory simulation.
 *  - energy: Final Hamiltonian (−log posterior + kinetic energy) of accepted state.
 *
 * These diagnostics are used to assess performance and identify issues such as
 * poor geometry (e.g., divergences or saturated tree depth).
 */
struct NUTSDiagnostics : public DiagnosticsBase {
  int tree_depth;
  bool divergent;
  double energy;
};



/**
 * Struct: SamplerResult
 *
 * Represents the final outcome of one iteration of the NUTS sampler.
 *
 * Fields:
 *  - state: Final accepted position (parameter vector).
 *  - accept_prob: Acceptance probability.
 *  - diagnostics:
 */
struct SamplerResult {
  arma::vec state;
  double accept_prob;
  std::shared_ptr<DiagnosticsBase> diagnostics;
};



// -----------------------------------------------------------------------------
// MCMC Adaptation Utilities
// -----------------------------------------------------------------------------



/**
 * Function: update_proposal_sd_with_robbins_monro
 * Purpose: Performs Robbins-Monro updates for proposal standard deviations.
 *
 * Inputs:
 *  - current_sd: Current standard deviation of the proposal.
 *  - observed_log_acceptance_probability: Log acceptance probability from the Metropolis-Hastings step.
 *  - rm_weight: Robbins-Monro adaptation weight (e.g. iteration^{-0.75}).
 *
 * Returns:
 *  - Updated proposal standard deviation, clamped within bounds.
 */
inline double update_proposal_sd_with_robbins_monro (
    const double current_sd,
    const double observed_log_acceptance_probability,
    const double rm_weight,
    const double target_acceptance
) {
  constexpr double rm_lower_bound = 0.001;
  constexpr double rm_upper_bound = 2.0;

  // Normalize the acceptance probability
  double observed_acceptance_probability = 1.0;
  if (observed_log_acceptance_probability < 0.0) {
    observed_acceptance_probability = MY_EXP (observed_log_acceptance_probability);
  }

  // Robbins-Monro update step
  double updated_sd = current_sd +
    (observed_acceptance_probability - target_acceptance) * rm_weight;

  // Handle NaNs robustly
  if (std::isnan (updated_sd)) {
    updated_sd = 1.0;
  }

  return std::clamp (updated_sd, rm_lower_bound, rm_upper_bound);
}



double kinetic_energy(const arma::vec& r, const arma::vec& inv_mass_diag);



/**
 * Step size heuristic using joint log_post+gradient function.
 *
 * Eliminates redundant probability computations by using joint at both
 * the initial and proposed positions (1 leapfrog step each iteration).
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
 * Step size heuristic with mass matrix, using joint function.
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