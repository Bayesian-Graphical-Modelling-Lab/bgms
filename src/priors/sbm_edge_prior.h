#pragma once

#include <RcppArmadillo.h>
struct SafeRNG;



// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// ----------------------------------------------------------------------------|
arma::uvec block_allocations_mfm_sbm(arma::uvec cluster_assign,
                                                arma::uword no_variables,
                                                arma::vec log_Vn,
                                                arma::mat block_probs,
                                                arma::umat indicator,
                                                arma::uword dirichlet_alpha,
                                                double beta_bernoulli_alpha,
                                                double beta_bernoulli_beta,
                                                double beta_bernoulli_alpha_between,
                                                double beta_bernoulli_beta_between,
                                                bool sbm_singleton_boost,
                                                SafeRNG& rng);

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// ----------------------------------------------------------------------------|
arma::mat block_probs_mfm_sbm(arma::uvec cluster_assign,
                                        arma::umat indicator,
                                        arma::uword no_variables,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta,
                                        double beta_bernoulli_alpha_between,
                                        double beta_bernoulli_beta_between,
                                        SafeRNG& rng);