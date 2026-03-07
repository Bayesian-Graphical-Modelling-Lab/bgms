#ifndef SAMPLEROUTPUT_H
#define SAMPLEROUTPUT_H

#include <RcppArmadillo.h>



/**
 * Container for the output of a single MCMC chain.
 *
 * Stores posterior samples of main and pairwise effects, optional
 * inclusion indicators, and diagnostics for HMC/NUTS runs.
 */
struct bgmCompareOutput {
  /// Main-effect samples [iter x (#main x groups)].
  arma::mat main_samples;
  /// Pairwise-effect samples [iter x (#pair x groups)].
  arma::mat pairwise_samples;
  /// Inclusion indicator samples [iter x (#edges + #variables)] (if used).
  arma::imat indicator_samples;

  /// Tree depth diagnostics [iter] (NUTS only).
  arma::ivec treedepth_samples;
  /// Divergent transition flags [iter] (NUTS only).
  arma::ivec divergent_samples;
  /// Energy diagnostic [iter] (NUTS only).
  arma::vec energy_samples;

  /// Identifier of the chain.
  int chain_id;
  /// True if indicator samples are stored.
  bool has_indicator;
  /// True if the chain was interrupted by the user.
  bool userInterrupt;
};

#endif
