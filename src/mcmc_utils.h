#pragma once

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include "explog_switch.h"
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
 * Adapts the log step size using dual averaging during MCMC burn-in.
 *
 * Implements the dual averaging algorithm to adaptively update the MALA step size
 * toward a target acceptance rate. Only used during burn-in.
 *
 * Inputs:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration number (starting from 1).
 *  - state: Vector of length 3, modified in-place:
 *      * state[0] = current log step size
 *      * state[1] = running average of log step size
 *      * state[2] = running average of acceptance error
 *
 * Modifies:
 *  - `state` vector in-place to update the log step size parameters.
 */
inline void update_step_size_with_dual_averaging (
    const double initial_step_size,
    const double acceptance_probability,
    const int iteration,
    arma::vec& state,
    const double target_acceptance
) {
  const double target_log_step_size = MY_LOG (10.0 * initial_step_size);
  constexpr int stabilization_offset = 10;

  double& log_step_size = state[0];
  double& log_step_size_avg = state[1];
  double& acceptance_error_avg = state[2];

  const double adjusted_iter = iteration + stabilization_offset;
  const double error = target_acceptance - acceptance_probability;

  acceptance_error_avg = (1.0 - 1.0 / adjusted_iter) * acceptance_error_avg +
    (1.0 / adjusted_iter) * error;

  log_step_size = target_log_step_size - std::sqrt (static_cast<double> (iteration)) / 0.05 * acceptance_error_avg;

  const double weight = std::pow (static_cast<double> (iteration), -0.75);
  log_step_size_avg = weight * log_step_size + (1.0 - weight) * log_step_size_avg;
}



/**
 * Function: update_step_size_with_robbins_monro
 *
 * Performs Robbins-Monro adaptation of the step size on the log scale.
 * Applies exponential decay to gradually reduce adaptation influence over time.
 *
 * Inputs:
 *  - acceptance_probability: Observed acceptance rate at the current iteration.
 *  - iteration: Current iteration number (starting from 1).
 *  - step_size_mala: Step size to update (in-place).
 *
 * Notes:
 *  - Decay rate is 0.75 for stable adaptation.
 */
inline void update_step_size_with_robbins_monro (
    const double acceptance_probability,
    const int iteration,
    double& step_size,
    const double target_acceptance
) {
  constexpr double decay_rate = 0.75;

  const double error = acceptance_probability - target_acceptance;
  const double decay = std::pow (static_cast<double> (iteration), -decay_rate);

  double log_step_size = MY_LOG (step_size);
  log_step_size += error * decay;
  step_size = MY_EXP (log_step_size);
}



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



double heuristic_initial_step_size(
    const arma::vec& theta,
    const std::function<double(const arma::vec&)>& log_post,
    const std::function<arma::vec(const arma::vec&)>& grad,
    SafeRNG& rng,
    double target_acceptance = 0.8,
    double init_step = 1.0,
    int max_attempts = 100
);