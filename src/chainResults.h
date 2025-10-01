#pragma once

#include <RcppArmadillo.h>
#include <string>

/**
 * Container for the result of a single MCMC chain.
 *
 * Fields:
 *  - error: True if the chain terminated with an error, false otherwise.
 *  - error_msg: Error message in case of failure (empty if no error).
 *  - chain_id: Integer identifier for the chain (1-based).
 *  - result: Rcpp::List containing the chain’s outputs (samples, diagnostics, etc.).
 *
 * Usage:
 *  - Used in parallel samplers to collect per-chain results.
 *  - Checked after execution to propagate errors or assemble outputs into R.
 */
struct ChainResult {
  bool error;
  std::string error_msg;
  int chain_id;
  bool userInterrupt;
  arma::mat main_effect_samples;
  arma::mat pairwise_effect_samples;
  arma::ivec treedepth_samples;
  arma::ivec divergent_samples;
  arma::vec  energy_samples;
  arma::imat  indicator_samples;
  arma::ivec allocation_samples;
};
