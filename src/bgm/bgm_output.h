#pragma once
#include <RcppArmadillo.h>

struct bgmOutput {
  // required
  arma::mat main_samples;
  arma::mat pairwise_samples;

  // optional (only if edge_selection)
  arma::imat indicator_samples;
  arma::imat allocation_samples;

  // optional (only if NUTS)
  arma::ivec treedepth_samples;
  arma::ivec divergent_samples;
  arma::vec energy_samples;

  // metadata
  int  chain_id = -1;
  bool userInterrupt = false;
};
