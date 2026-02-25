/**
 * @file sbm_edge_prior_interface.h
 * @brief R-callable function to compute MFM-SBM partition coefficients.
 */

#include <RcppArmadillo.h>
#include "math/explog_macros.h"


/**
 * Compute log partition coefficients V_n(t) for the MFM-SBM.
 *
 * Uses the log-sum-exp identity for numerical stability. Based on
 * Algorithm 1 from Geng, Bhattacharya & Pati (2019, JASA 114:526).
 *
 * @param num_variables   Number of variables p.
 * @param dirichlet_alpha Dirichlet concentration parameter.
 * @param t_max           Maximum number of clusters to evaluate.
 * @param lambda          Poisson rate for the prior on the number of clusters.
 * @return Vector of length t_max with log V_n(t) values.
 */
arma::vec compute_Vn_mfm_sbm(arma::uword num_variables,
                             double dirichlet_alpha,
                             arma::uword t_max,
                             double lambda);