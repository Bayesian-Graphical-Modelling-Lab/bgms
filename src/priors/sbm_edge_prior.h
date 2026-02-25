#pragma once

/**
 * @file sbm_edge_prior.h
 * @brief Gibbs samplers for the Mixture of Finite Mixtures
 *        Stochastic Block Model (MFM-SBM) edge prior.
 *
 * Implements the block-allocation and block-probability updates
 * described in Geng, Bhattacharya & Pati (2019, JASA 114:526).
 * Called by StochasticBlockEdgePrior::update() in edge_prior.h.
 */

#include <RcppArmadillo.h>
struct SafeRNG;


/**
 * Sample block allocations for the MFM-SBM.
 *
 * Reassigns each variable to a cluster via a collapsed Gibbs step,
 * integrating over block-level inclusion probabilities.
 *
 * @param cluster_assign   Current cluster assignment (length p).
 * @param no_variables      Number of variables p.
 * @param log_Vn            Log partition coefficients from compute_Vn_mfm_sbm().
 * @param block_probs       Current block-level inclusion probability matrix.
 * @param indicator         Edge indicator matrix (p x p, upper-triangular).
 * @param dirichlet_alpha   Dirichlet concentration parameter for cluster sizes.
 * @param beta_bernoulli_alpha   Beta-Bernoulli alpha for within-block edges.
 * @param beta_bernoulli_beta    Beta-Bernoulli beta for within-block edges.
 * @param beta_bernoulli_alpha_between  Beta-Bernoulli alpha for between-block edges.
 * @param beta_bernoulli_beta_between   Beta-Bernoulli beta for between-block edges.
 * @param rng               Random number generator.
 * @return Updated cluster assignment vector (length p).
 */
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
                                                SafeRNG& rng);

/**
 * Sample block-level inclusion probabilities for the MFM-SBM.
 *
 * Draws within-block probabilities from Beta(alpha + included, beta + excluded)
 * and between-block probabilities from separate Beta hyperparameters.
 *
 * @param cluster_assign   Current cluster assignment (length p).
 * @param indicator         Edge indicator matrix (p x p, upper-triangular).
 * @param no_variables      Number of variables p.
 * @param beta_bernoulli_alpha   Beta-Bernoulli alpha for within-block edges.
 * @param beta_bernoulli_beta    Beta-Bernoulli beta for within-block edges.
 * @param beta_bernoulli_alpha_between  Beta-Bernoulli alpha for between-block edges.
 * @param beta_bernoulli_beta_between   Beta-Bernoulli beta for between-block edges.
 * @return Block-level inclusion probability matrix (K x K).
 */
arma::mat block_probs_mfm_sbm(arma::uvec cluster_assign,
                                        arma::umat indicator,
                                        arma::uword no_variables,
                                        double beta_bernoulli_alpha,
                                        double beta_bernoulli_beta,
                                        double beta_bernoulli_alpha_between,
                                        double beta_bernoulli_beta_between,
                                        SafeRNG& rng);