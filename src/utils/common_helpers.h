#pragma once
#include <RcppArmadillo.h>



/**
 * Count the total number of main-effect parameters across all variables.
 *
 * Each variable contributes:
 *  - Ordinal variable: `num_categories[v]` parameters (one per category).
 *  - Blume–Capel variable: always 2 parameters (linear and quadratic terms).
 *
 * Inputs:
 *  - num_categories: Number of categories per variable [V].
 *  - is_ordinal_variable: Indicator of whether each variable is ordinal [V].
 *
 * Returns:
 *  - Total number of main-effect parameters (int).
 */
inline int count_num_main_effects(const arma::ivec& num_categories,
                                  const arma::uvec& is_ordinal_variable) {
  int n_params = 0;
  for (arma::uword i = 0; i < num_categories.n_elem; ++i) {
    n_params += is_ordinal_variable[i] ? num_categories[i] : 2;
  }
  return n_params;
}

/// MCMC update method.
enum UpdateMethod {
  adaptive_metropolis, ///< Robbins-Monro adaptive Metropolis-Hastings
  hamiltonian_mc,      ///< Fixed-trajectory Hamiltonian Monte Carlo
  nuts                 ///< No-U-Turn Sampler
};

/**
 * Convert a string identifier to an UpdateMethod enum value.
 *
 * @param update_method  One of "adaptive-metropolis", "hamiltonian-mc", "nuts".
 * @return Corresponding UpdateMethod value.
 */
inline UpdateMethod update_method_from_string(const std::string& update_method) {
  if (update_method == "adaptive-metropolis")
    return adaptive_metropolis;

  if (update_method == "hamiltonian-mc")
    return hamiltonian_mc;

  if (update_method == "nuts")
    return nuts;

  throw std::invalid_argument("Invalid update_method: " + update_method);
}

/// Edge inclusion prior type.
enum EdgePrior {
  Stochastic_Block, ///< MFM Stochastic Block Model
  Beta_Bernoulli,   ///< Shared Beta-Bernoulli prior
  Bernoulli,        ///< Fixed Bernoulli prior
  Not_Applicable    ///< No edge prior (all edges included)
};

/**
 * Convert a string identifier to an EdgePrior enum value.
 *
 * @param edge_prior  One of "Stochastic-Block", "Beta-Bernoulli",
 *                    "Bernoulli", "Not Applicable".
 * @return Corresponding EdgePrior value.
 */
inline EdgePrior edge_prior_from_string(const std::string& edge_prior) {
  if (edge_prior == "Stochastic-Block")
    return Stochastic_Block;

  if (edge_prior == "Beta-Bernoulli")
    return Beta_Bernoulli;

  if (edge_prior == "Bernoulli")
    return Bernoulli;

  if (edge_prior == "Not Applicable")
    return Not_Applicable;

  throw std::invalid_argument("Invalid edge_prior: " + edge_prior);
}