#ifndef SAMPLEROUTPUT_H
#define SAMPLEROUTPUT_H

#include <RcppArmadillo.h>



/**
 * Container for the output of a single MCMC chain.
 *
 * Stores posterior samples of main and pairwise effects, optional
 * inclusion indicators, and diagnostics for HMC/NUTS runs.
 *
 * Members:
 *  - main_samples:    [iter × (#main × groups)] matrix of main-effect samples.
 *  - pairwise_samples:[iter × (#pair × groups)] matrix of pairwise-effect samples.
 *  - indicator_samples:[iter × (#edges + #variables)] indicator samples (if used).
 *  - treedepth_samples:[iter] tree depth diagnostics (NUTS only).
 *  - divergent_samples:[iter] divergent transition flags (NUTS only).
 *  - energy_samples:   [iter] energy diagnostic (NUTS only).
 *  - chain_id:         Identifier of the chain.
 *  - has_indicator:    True if indicator samples are stored.
 */
struct SamplerOutput {
  arma::mat main_samples;
  arma::mat pairwise_samples;
  arma::imat indicator_samples;
  arma::ivec treedepth_samples;
  arma::ivec divergent_samples;
  arma::vec energy_samples;
  int chain_id;
  bool has_indicator;
  bool userInterrupt;
};

#endif
