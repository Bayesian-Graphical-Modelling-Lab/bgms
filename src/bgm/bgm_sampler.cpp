#include <RcppArmadillo.h>
#include "bgm/bgm_helper.h"
#include "bgm/bgm_logp_and_grad.h"
#include "bgm/bgm_sampler.h"
#include "bgm/bgm_output.h"
#include "mcmc/mcmc_adaptation.h"
#include "mcmc/mcmc_hmc.h"
#include "mcmc/mcmc_leapfrog.h"
#include "mcmc/mcmc_nuts.h"
#include "mcmc/mcmc_rwm.h"
#include "mcmc/mcmc_utils.h"
#include "priors/sbm_edge_prior.h"
#include "sbm_edge_prior_interface.h"
#include "rng/rng_utils.h"
#include "utils/common_helpers.h"
#include "utils/progress_manager.h"



/**
 * Imputes missing values in the observation matrix for the bgm model.
 *
 * Missing entries are replaced by draws from their conditional
 * pseudo-likelihood, given the current main and pairwise effects.
 * After imputation, sufficient statistics and residuals are updated.
 *
 * Procedure:
 *  - For each missing (person, variable):
 *    * Compute unnormalized probabilities across categories using
 *      the variable’s type (ordinal vs. Blume–Capel).
 *    * Sample a new value by inverse transform.
 *    * If the value changes, update:
 *        - observations
 *        - counts_per_category (ordinal)
 *        - blume_capel_stats (Blume–Capel)
 *        - residual_matrix (all variables for that person).
 *  - Finally, recompute pairwise_stats as XᵀX from updated observations.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths.
 *  - main_effects: Matrix of main-effect parameters (variables × categories).
 *  - observations: Matrix of categorical scores (persons × variables), updated in place.
 *  - counts_per_category: Category counts per variable (updated for ordinal vars).
 *  - blume_capel_stats: Sufficient statistics (updated for Blume–Capel vars).
 *  - num_categories: Number of categories per variable.
 *  - residual_matrix: Residual scores (persons × variables), updated in place.
 *  - missing_index: Matrix of missing entries (rows = [person, variable]).
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - pairwise_stats: Pairwise sufficient statistics, recomputed at the end.
 *  - rng: Random number generator.
 *
 * Notes:
 *  - Sampling is done via inverse-transform using cumulative probabilities.
 *  - For Blume–Capel variables, both linear and quadratic terms are used.
 *  - Residuals are updated incrementally for efficiency.
 *  - Pairwise statistics are recomputed in full at the end (XᵀX); this could
 *    be optimized further to update incrementally.
 */
void impute_missing_bgm (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& observations,
    arma::imat& counts_per_category,
    arma::imat& blume_capel_stats,
    const arma::ivec& num_categories,
    arma::mat& residual_matrix,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    arma::imat& pairwise_stats,
    SafeRNG& rng
) {
  const int num_variables = observations.n_cols;
  const int num_missings = missing_index.n_rows;
  const int max_num_categories = num_categories.max ();

  arma::vec category_probabilities (max_num_categories + 1);

  for (int miss = 0; miss < num_missings; miss++) {
    const int person = missing_index (miss, 0);
    const int variable = missing_index (miss, 1);

    const double residual_score = residual_matrix (person, variable);
    const int num_cats = num_categories (variable);
    const bool is_ordinal = is_ordinal_variable (variable);

    double cumsum = 0.0;

    if (is_ordinal) {
      // Compute cumulative unnormalized probabilities for ordinal variable
      cumsum = 1.0;
      category_probabilities[0] = cumsum;
      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const double exponent = main_effects (variable, cat) + score * residual_score;
        cumsum += MY_EXP (exponent);
        category_probabilities[score] = cumsum;
      }
    } else {
      // Compute probabilities for Blume-Capel variable
      const int ref = baseline_category (variable);

      cumsum = MY_EXP (
        main_effects (variable, 0) * ref + main_effects (variable, 1) * ref * ref
      );
      category_probabilities[0] = cumsum;

      for (int cat = 0; cat <= num_cats; cat++) {
        const int score = cat - ref;
        const double exponent =
          main_effects (variable, 0) * score +
          main_effects (variable, 1) * score * score +
          score * residual_score;
        cumsum += MY_EXP (exponent);
        category_probabilities[cat] = cumsum;
      }
    }

    // Sample from categorical distribution via inverse transform
    const double u = runif (rng) * cumsum;
    int sampled_score = 0;
    while (u > category_probabilities[sampled_score]) {
      sampled_score++;
    }

    int new_value = sampled_score;
    if(!is_ordinal)
      new_value -= baseline_category (variable);
    const int old_value = observations(person, variable);

    if (new_value != old_value) {
      // Update observation matrix
      observations(person, variable) = new_value;

      if (is_ordinal) {
        counts_per_category(old_value, variable)--;
        counts_per_category(new_value, variable)++;
      } else {
        const int delta = new_value - old_value;
        const int delta_sq = new_value * new_value - old_value * old_value;

        blume_capel_stats(0, variable) += delta;
        blume_capel_stats(1, variable) += delta_sq;
      }
      // Update residuals across all variables
      for (int var = 0; var < num_variables; var++) {
        const double delta_score = (new_value - old_value) * pairwise_effects(var, variable);
        residual_matrix(person, var) += delta_score;
      }
    }
  }

  // Update sufficient statistics for pairwise effects
  pairwise_stats = observations.t() * observations;
  // This could be done more elegantly....
}



/**
 * Heuristically finds an initial step size for HMC/NUTS (bgm model).
 *
 * The routine follows the scheme of Hoffman & Gelman (2014), Algorithm 4:
 * it searches for a step size ε such that the acceptance probability is close
 * to 0.5 under a single leapfrog step. This serves as a good starting point
 * for dual-averaging adaptation during warmup.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters.
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Hyperparameters for Beta priors.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - target_acceptance: Target acceptance probability (default ~0.65 in practice).
 *  - pairwise_stats: Pairwise sufficient statistics.
 *  - rng: Random number generator.
 *
 * Returns:
 *  - A scalar initial step size ε to be used in HMC/NUTS warmup.
 *
 * Notes:
 *  - Uses `vectorize_model_parameters_bgm()` and its inverse to map between
 *    parameter matrices and flat vectors.
 *  - Calls `log_pseudoposterior()` and `gradient_log_pseudoposterior()` internally.
 *  - This function is typically called once before adaptation starts.
 */
double find_initial_stepsize_bgm(
    arma::mat main_effects,
    arma::mat pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    const double target_acceptance,
    const arma::imat& pairwise_stats,
    SafeRNG& rng
) {
  arma::vec theta = vectorize_model_parameters_bgm(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );

  arma::vec grad_obs;
  arma::imat index_matrix;

  std::tie(grad_obs, index_matrix) = gradient_observed_active(
      inclusion_indicator, observations, num_categories, counts_per_category,
      blume_capel_stats, baseline_category, is_ordinal_variable, pairwise_stats
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgm(theta_vec, current_main, current_pair, inclusion_indicator,
                                     num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;

    return gradient_log_pseudoposterior (
        current_main, current_pair, inclusion_indicator, observations,
        num_categories, baseline_category, is_ordinal_variable, main_alpha,
        main_beta, pairwise_scale, pairwise_scaling_factors, rm, index_matrix, grad_obs
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgm(theta_vec, current_main, current_pair,
                                     inclusion_indicator,
                                     num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;
    return log_pseudoposterior(
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, counts_per_category, blume_capel_stats,
      baseline_category, is_ordinal_variable, main_alpha, main_beta,
      pairwise_scale, pairwise_scaling_factors, pairwise_stats, rm
    );
  };

  return heuristic_initial_step_size(theta, log_post, grad, rng, target_acceptance);
}



/**
 * Updates all main-effect parameters using random-walk Metropolis (bgm model).
 *
 * For each variable:
 *  - Ordinal variables: update each threshold parameter separately.
 *  - Blume–Capel variables: update the linear and quadratic parameters.
 *
 * Each parameter is updated by proposing a new value from a normal
 * distribution (centered at the current value, with proposal SD given by
 * `proposal_sd_main`), evaluating the log-pseudoposterior, and applying
 * the Metropolis acceptance step.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (updated in place).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - num_persons: Number of observations (not used directly here).
 *  - main_alpha, main_beta: Hyperparameters for Beta priors.
 *  - residual_matrix: Residual scores (persons × variables).
 *  - proposal_sd_main: Proposal standard deviations for each parameter (updated adaptively).
 *  - adapter: Random-walk adaptation controller (updated each iteration).
 *  - iteration: Current iteration number.
 *  - rng: Random number generator.
 *
 * Outputs:
 *  - main_effects: Updated parameter values after Metropolis steps.
 *  - proposal_sd_main: Updated proposal SDs if adaptation is active.
 *  - adapter: Updated with acceptance probabilities.
 *
 * Notes:
 *  - The log-pseudoposterior contribution for each parameter is evaluated via
 *    `log_pseudoposterior_main_effects_component()`.
 *  - Acceptance probabilities are tracked in `accept_prob_main` and passed to the adapter.
 *  - Currently `num_persons` is unused and could be removed from the signature.
 */
void update_main_effects_metropolis_bgm (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const int num_persons,
    const double main_alpha,
    const double main_beta,
    const arma::mat& residual_matrix,
    arma::mat& proposal_sd_main,
    RWMAdaptationController& adapter,
    const int iteration,
    SafeRNG& rng
) {
  const int num_vars = observations.n_cols;
  arma::umat index_mask_main = arma::ones<arma::umat>(proposal_sd_main.n_rows, proposal_sd_main.n_cols);
  arma::mat accept_prob_main = arma::ones<arma::mat>(proposal_sd_main.n_rows, proposal_sd_main.n_cols);

  for(int variable = 0; variable < num_vars; variable++) {
    const int num_cats = num_categories(variable);
    if(is_ordinal_variable[variable] == true) {
      for (int category = 0; category < num_cats; category++) {
        double& current = main_effects(variable, category);
        double proposal_sd = proposal_sd_main(variable, category);

        auto log_post = [&](double theta) {
          main_effects(variable, category) = theta;
          return log_pseudoposterior_main_effects_component(
            main_effects, residual_matrix, num_categories, counts_per_category,
            blume_capel_stats, baseline_category, is_ordinal_variable,
            main_alpha, main_beta, variable, category, -1
          );
        };

        SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);

        current = result.state[0];
        accept_prob_main(variable, category) = result.accept_prob;
      }
    } else {
      for (int parameter = 0; parameter < 2; parameter++) {
        double& current = main_effects(variable, parameter);
        double proposal_sd = proposal_sd_main(variable, parameter);

        auto log_post = [&](double theta) {
          main_effects(variable, parameter) = theta;
          return log_pseudoposterior_main_effects_component(
            main_effects, residual_matrix, num_categories, counts_per_category,
            blume_capel_stats, baseline_category, is_ordinal_variable,
            main_alpha, main_beta,
            variable, -1, parameter
          );
        };

        SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);
        current = result.state[0];
        accept_prob_main(variable, parameter) = result.accept_prob;
      }
    }
  }

  adapter.update(index_mask_main, accept_prob_main, iteration);
}



/**
 * Updates all active pairwise interaction parameters using random-walk Metropolis (bgm model).
 *
 * For each active pair (variable1, variable2):
 *  - Propose a new interaction value from a normal distribution centered at the current value,
 *    with proposal SD given by `proposal_sd_pairwise_effects`.
 *  - Evaluate the log-pseudoposterior with
 *    `log_pseudoposterior_interactions_component()`.
 *  - Accept or reject the proposal according to the Metropolis criterion.
 *  - If accepted, update the residual matrix incrementally for efficiency.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters (updated in place).
 *  - main_effects: Matrix of main-effect parameters.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd_pairwise_effects: Proposal SDs for each pair (updated adaptively).
 *  - adapter: Random-walk adaptation controller.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - num_persons: Number of observations (not used directly).
 *  - num_variables: Number of variables.
 *  - residual_matrix: Residual scores (persons × variables), updated in place.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - iteration: Current iteration number.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - rng: Random number generator.
 *
 * Outputs:
 *  - pairwise_effects: Updated with accepted proposals.
 *  - residual_matrix: Incrementally updated when pairwise_effects change.
 *  - proposal_sd_pairwise_effects: Updated if adaptation is active.
 *  - adapter: Updated with acceptance probabilities.
 *
 * Notes:
 *  - Only pairs with `inclusion_indicator(i,j) == 1` are updated.
 *  - Symmetry is enforced: (i,j) and (j,i) are always set to the same value.
 *  - `num_persons` is unused and could be removed from the signature.
 */
void update_pairwise_effects_metropolis_bgm (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    arma::mat& proposal_sd_pairwise_effects,
    RWMAdaptationController& adapter,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    const int num_variables,
    arma::mat& residual_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int iteration,
    const arma::imat& pairwise_stats,
    SafeRNG& rng
) {
  arma::mat accept_prob_pairwise = arma::zeros<arma::mat>(num_variables, num_variables);
  arma::umat index_mask_pairwise = arma::zeros<arma::umat>(num_variables, num_variables);

  for (int variable1 = 0; variable1 < num_variables - 1; variable1++) {
    for (int variable2 = variable1 + 1; variable2 < num_variables; variable2++) {
      if (inclusion_indicator(variable1, variable2) == 1) {
        double& value = pairwise_effects(variable1, variable2);
        double proposal_sd = proposal_sd_pairwise_effects(variable1, variable2);
        double current = value;

        auto log_post = [&](double theta) {
          pairwise_effects(variable1, variable2) = theta;
          pairwise_effects(variable2, variable1) = theta;

          return log_pseudoposterior_interactions_component(
            pairwise_effects, main_effects, observations, num_categories,
            inclusion_indicator, is_ordinal_variable, baseline_category,
            pairwise_scale, pairwise_scaling_factors, pairwise_stats, variable1, variable2
          );
        };

        SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);

        value = result.state[0];
        pairwise_effects(variable2, variable1) = value;

        if(current != value) {
          double delta = value - current;
          residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
          residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
        }

        accept_prob_pairwise(variable1, variable2) = result.accept_prob;
        index_mask_pairwise(variable1, variable2) = 1;
      }
    }
  }

  adapter.update(index_mask_pairwise, accept_prob_pairwise, iteration);
}



/**
 * Performs one Hamiltonian Monte Carlo (HMC) update of main and pairwise parameters (bgm model).
 *
 * Procedure:
 *  - Flatten parameters into a vector with `vectorize_model_parameters_bgm()`.
 *  - Define log-pseudoposterior and gradient functions using
 *    `log_pseudoposterior()` and `gradient_log_pseudoposterior_active()`.
 *  - Run the HMC leapfrog integrator via `hmc_sampler()`.
 *  - Unpack the accepted state back into `main_effects` and `pairwise_effects`.
 *  - Recompute the residual matrix and update the adaptation controller.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (updated in place).
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths (updated in place).
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - residual_matrix: Residual scores (persons × variables), updated in place.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - num_leapfrogs: Number of leapfrog steps for HMC integration.
 *  - iteration: Current iteration number.
 *  - adapt: HMC adaptation controller (step size, mass matrix).
 *  - learn_mass_matrix: If true, adapt the mass matrix during warmup.
 *  - selection: If true, restrict gradient/mass matrix to active parameters.
 *  - rng: Random number generator.
 *
 * Outputs:
 *  - main_effects, pairwise_effects: Updated if the HMC proposal is accepted.
 *  - residual_matrix: Recomputed as observations × pairwise_effects.
 *  - adapt: Updated with acceptance probability and new state.
 *
 * Notes:
 *  - Uses `inv_mass_active()` to extract the relevant diagonal mass matrix.
 *  - This update is called within the Gibbs sampling loop when HMC is selected.
 */
void update_hmc_bgm(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    arma::mat& residual_matrix,
    const arma::imat& pairwise_stats,
    const int num_leapfrogs,
    const int iteration,
    HMCAdaptationController& adapt,
    const bool learn_mass_matrix,
    const bool selection,
    SafeRNG& rng
) {
  arma::vec current_state = vectorize_model_parameters_bgm(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );

  arma::vec grad_obs;
  arma::imat index_matrix;

  std::tie(grad_obs, index_matrix) = gradient_observed_active(
      inclusion_indicator, observations, num_categories, counts_per_category,
      blume_capel_stats, baseline_category, is_ordinal_variable, pairwise_stats
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgm(theta_vec, current_main, current_pair, inclusion_indicator,
                                 num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;

    return gradient_log_pseudoposterior (
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, baseline_category, is_ordinal_variable, main_alpha,
      main_beta, pairwise_scale, pairwise_scaling_factors, rm, index_matrix, grad_obs
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgm(theta_vec, current_main, current_pair, inclusion_indicator,
                                 num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;
    return log_pseudoposterior (
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, counts_per_category, blume_capel_stats,
      baseline_category, is_ordinal_variable, main_alpha, main_beta,
      pairwise_scale, pairwise_scaling_factors, pairwise_stats, rm
    );
  };

  arma::vec active_inv_mass = inv_mass_active(
    adapt.inv_mass_diag(), inclusion_indicator, num_categories,
    is_ordinal_variable, selection
  );

  SamplerResult result = hmc_sampler(
    current_state, adapt.current_step_size(), log_post, grad, num_leapfrogs,
    active_inv_mass, rng
  );

  current_state = result.state;
  unvectorize_model_parameters_bgm(
    current_state, main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );
  residual_matrix = observations * pairwise_effects;

  adapt.update(current_state, result.accept_prob, iteration);

  // If mass matrix was just updated, re-run the heuristic to find a good
  // step size for the new mass matrix (following STAN's approach).
  // STAN uses the current step size as the starting point for the heuristic.
  if (adapt.mass_matrix_just_updated()) {
    arma::vec new_inv_mass = inv_mass_active(
      adapt.inv_mass_diag(), inclusion_indicator, num_categories,
      is_ordinal_variable, selection
    );
    double current_eps = adapt.current_step_size();
    double new_eps = heuristic_initial_step_size(
      current_state, log_post, grad, new_inv_mass, rng,
      0.625,        // target_acceptance
      current_eps   // init_step: use current step size as starting point
    );
    adapt.reinit_stepsize(new_eps);
  }
}



/**
 * Performs one No-U-Turn Sampler (NUTS) update of main and pairwise parameters (bgm model).
 *
 * Procedure:
 *  - Flatten parameters into a vector with `vectorize_model_parameters_bgm()`.
 *  - Define log-pseudoposterior and gradient functions using
 *    `log_pseudoposterior()` and `gradient_log_pseudoposterior_active()`.
 *  - Run the NUTS sampler via `nuts_sampler()`, building a trajectory
 *    up to the maximum tree depth.
 *  - Unpack the accepted state back into `main_effects` and `pairwise_effects`.
 *  - Recompute the residual matrix and update the adaptation controller.
 *
 * Inputs:
 *  - main_effects: Matrix of main-effect parameters (updated in place).
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths (updated in place).
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - counts_per_category: Category counts per variable (ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - residual_matrix: Residual scores (persons × variables), updated in place.
 *  - nuts_max_depth: Maximum tree depth for the NUTS trajectory.
 *  - iteration: Current iteration number.
 *  - adapt: HMC/NUTS adaptation controller (step size, mass matrix).
 *  - learn_mass_matrix: If true, adapt the mass matrix during warmup.
 *  - selection: If true, restrict gradient/mass matrix to active parameters.
 *  - rng: Random number generator.
 *
 * Returns:
 *  - SamplerResult containing:
 *      * state: Accepted parameter vector.
 *      * accept_prob: Acceptance probability for the proposal.
 *
 * Notes:
 *  - Uses `inv_mass_active()` to extract the relevant diagonal mass matrix.
 *  - The step size is managed by the adaptation controller (`adapt`).
 *  - This update is called within the Gibbs sampling loop when NUTS is selected.
 */
SamplerResult update_nuts_bgm(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const arma::ivec& baseline_category,
    const arma::uvec& is_ordinal_variable,
    const double main_alpha,
    const double main_beta,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    const arma::imat& pairwise_stats,
    arma::mat& residual_matrix,
    const int nuts_max_depth,
    const int iteration,
    HMCAdaptationController& adapt,
    const bool learn_mass_matrix,
    const bool selection,
    SafeRNG& rng
) {
  arma::vec current_state = vectorize_model_parameters_bgm(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );

  arma::vec grad_obs;
  arma::imat index_matrix;

  std::tie(grad_obs, index_matrix) = gradient_observed_active(
      inclusion_indicator, observations, num_categories, counts_per_category,
      blume_capel_stats, baseline_category, is_ordinal_variable, pairwise_stats
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgm(theta_vec, current_main, current_pair, inclusion_indicator,
                                     num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;

    return gradient_log_pseudoposterior (
        current_main, current_pair, inclusion_indicator, observations,
        num_categories, baseline_category, is_ordinal_variable, main_alpha,
        main_beta, pairwise_scale, pairwise_scaling_factors, rm, index_matrix, grad_obs
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters_bgm(theta_vec, current_main, current_pair,
                                 inclusion_indicator, num_categories,
                                 is_ordinal_variable);
    arma::mat rm = observations * current_pair;
    return log_pseudoposterior(
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, counts_per_category, blume_capel_stats,
      baseline_category, is_ordinal_variable, main_alpha, main_beta,
      pairwise_scale, pairwise_scaling_factors, pairwise_stats, rm
    );
  };

  arma::vec active_inv_mass = inv_mass_active(
    adapt.inv_mass_diag(), inclusion_indicator, num_categories,
    is_ordinal_variable, selection
  );

  SamplerResult result = nuts_sampler(
    current_state, adapt.current_step_size(), log_post, grad,
    active_inv_mass, rng, nuts_max_depth
  );

  current_state = result.state;
  unvectorize_model_parameters_bgm(
    current_state, main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );
  residual_matrix = observations * pairwise_effects;

  adapt.update(current_state, result.accept_prob, iteration);

  // If mass matrix was just updated, re-run the heuristic to find a good
  // step size for the new mass matrix (following STAN's approach).
  // STAN uses the current step size as the starting point for the heuristic.
  if (adapt.mass_matrix_just_updated()) {
    arma::vec new_inv_mass = inv_mass_active(
      adapt.inv_mass_diag(), inclusion_indicator, num_categories,
      is_ordinal_variable, selection
    );
    double current_eps = adapt.current_step_size();
    double new_eps = heuristic_initial_step_size(
      current_state, log_post, grad, new_inv_mass, rng,
      0.625,        // target_acceptance
      current_eps   // init_step: use current step size as starting point
    );
    adapt.reinit_stepsize(new_eps);
  }

  return result;
}



/**
 * Adapts proposal standard deviations for pairwise effects during warmup (bgm model).
 *
 * For each pairwise effect (variable1, variable2):
 *  - Propose a new value via random-walk Metropolis using the current proposal SD.
 *  - Accept/reject the proposal and update the parameter and residuals if accepted.
 *  - Update the proposal SD using a Robbins–Monro adaptation step to target
 *    the desired acceptance probability.
 *
 * Adaptation is only performed if the warmup schedule (`sched`) allows it
 * at the current iteration.
 *
 * Inputs:
 *  - proposal_sd_pairwise_effects: Current proposal SDs (symmetric matrix, updated in place).
 *  - pairwise_effects: Symmetric matrix of pairwise interaction parameters (updated in place).
 *  - main_effects: Matrix of main-effect parameters.
 *  - inclusion_indicator: Symmetric binary matrix of active pairwise effects.
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - residual_matrix: Residual scores (persons × variables), updated in place.
 *  - num_categories: Number of categories per variable.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - iteration: Current iteration number.
 *  - sched: Warmup schedule controller (determines if adaptation is active).
 *  - rng: Random number generator.
 *  - target_accept: Target acceptance probability (default 0.44 for RWM).
 *  - rm_decay: Robbins–Monro decay exponent (default 0.75).
 *
 * Outputs:
 *  - proposal_sd_pairwise_effects: Updated with adapted proposal SDs.
 *  - pairwise_effects: Updated with accepted proposals.
 *  - residual_matrix: Updated incrementally when proposals are accepted.
 *
 * Notes:
 *  - Adaptation stops once warmup is finished (`sched.adapt_proposal_sd()` returns false).
 *  - Symmetry is enforced: (i,j) and (j,i) entries of proposal SDs are always equal.
 *  - Uses Robbins–Monro stochastic approximation for stability.
 */
void tune_proposal_sd_bgm(
    arma::mat& proposal_sd_pairwise_effects,
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    const arma::imat& pairwise_stats,
    int iteration,
    const WarmupSchedule& sched,
    SafeRNG& rng,
    double target_accept = 0.44,
    double rm_decay = 0.75
)
{
  if (!sched.adapt_proposal_sd(iteration)) return;

  double t = iteration - sched.stage3b_start + 1;
  double rm_weight = std::pow(t, -rm_decay);

  const int num_variables = pairwise_effects.n_rows;

  for (int variable1 = 0; variable1 < num_variables - 1; variable1++) {
    for (int variable2 = variable1 + 1; variable2 < num_variables; variable2++) {
      double& value = pairwise_effects(variable1, variable2);
      double proposal_sd = proposal_sd_pairwise_effects(variable1, variable2);
      double current = value;

      auto log_post = [&](double theta) {
        pairwise_effects(variable1, variable2) = theta;
        pairwise_effects(variable2, variable1) = theta;

        return log_pseudoposterior_interactions_component(
          pairwise_effects, main_effects, observations, num_categories,
          inclusion_indicator, is_ordinal_variable, baseline_category,
          pairwise_scale, pairwise_scaling_factors, pairwise_stats, variable1, variable2
        );
      };

      SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);

      value = result.state[0];
      pairwise_effects(variable2, variable1) = value;

      if(current != value) {
        double delta = value - current;
        residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
        residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
      }

      proposal_sd = update_proposal_sd_with_robbins_monro(
        proposal_sd, MY_LOG(result.accept_prob), rm_weight, target_accept
      );

      proposal_sd_pairwise_effects(variable1, variable2) = proposal_sd;
      proposal_sd_pairwise_effects(variable2, variable1) = proposal_sd;
    }
  }
}



/**
 * Updates edge inclusion indicators and associated pairwise effects
 * using a Metropolis–Hastings step (bgm model).
 *
 * For each candidate interaction:
 *  - If the edge is currently absent, propose adding it with a random draw
 *    from a normal distribution centered at the current effect.
 *  - If the edge is currently present, propose removing it by setting the
 *    effect to zero.
 *  - Compute the log-acceptance ratio from:
 *      * log-pseudolikelihood ratio (`log_pseudolikelihood_ratio_interaction()`),
 *      * prior ratio (Cauchy prior on effects, Bernoulli prior on inclusion),
 *      * and proposal correction terms.
 *  - Accept or reject the proposal. If accepted:
 *      * Update the indicator matrix (symmetric),
 *      * Update the pairwise effect value (symmetric),
 *      * Incrementally update the residual matrix.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of pairwise interaction strengths (updated in place).
 *  - main_effects: Matrix of main-effect parameters.
 *  - indicator: Symmetric binary inclusion matrix (updated in place).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd: Proposal standard deviations for pairwise effects.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - index: Matrix mapping interaction index → (var1, var2).
 *  - num_interactions: Total number of candidate pairwise interactions.
 *  - num_persons: Number of observations.
 *  - residual_matrix: Residual scores (persons × variables), updated in place.
 *  - inclusion_probability: Prior inclusion probabilities for edges.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - rng: Random number generator.
 *
 * Outputs:
 *  - indicator: Updated with accepted edge inclusions/removals.
 *  - pairwise_effects: Updated pairwise effect values.
 *  - residual_matrix: Updated incrementally when changes occur.
 *
 * Notes:
 *  - Proposals are asymmetric: additions sample a new value, removals set to zero.
 *  - Inclusion prior is Bernoulli with edge-specific probability.
 *  - Residual updates are vectorized for efficiency.
 */
void update_indicator_edges_metropolis_bgm (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& proposal_sd,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    const arma::imat& index,
    const int num_interactions,
    const int num_persons,
    arma::mat& residual_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::imat& pairwise_stats,
    SafeRNG& rng
) {
  for (int cntr = 0; cntr < num_interactions; cntr++) {
    const int variable1 = index(cntr, 1);
    const int variable2 = index(cntr, 2);

    const double current_state = pairwise_effects(variable1, variable2);

    // Propose a new state: either add a new edge or remove an existing one
    const bool proposing_addition = (indicator(variable1, variable2) == 0);
    const double proposed_state = proposing_addition ? rnorm(rng, current_state, proposal_sd(variable1, variable2)) : 0.0;

    // Compute log pseudo-likelihood ratio
    double log_accept = log_pseudolikelihood_ratio_interaction (
      pairwise_effects, main_effects, observations, num_categories, num_persons,
      variable1, variable2, proposed_state, current_state, residual_matrix,
      is_ordinal_variable, baseline_category, pairwise_stats
    );

    // Add prior ratio and proposal correction
    const double inclusion_probability_ij = inclusion_probability(variable1, variable2);
    const double sd = proposal_sd(variable1, variable2);
    const double scaled_pairwise_scale = pairwise_scale * pairwise_scaling_factors(variable1, variable2);

    if (proposing_addition) {
      log_accept += R::dcauchy(proposed_state, 0.0, scaled_pairwise_scale, true);
      log_accept -= R::dnorm(proposed_state, current_state, sd, true);
      log_accept += MY_LOG (inclusion_probability_ij) - MY_LOG (1.0 - inclusion_probability_ij);
    } else {
      log_accept -= R::dcauchy(current_state, 0.0, scaled_pairwise_scale, true);
      log_accept += R::dnorm(current_state, proposed_state, sd, true);
      log_accept -= MY_LOG (inclusion_probability_ij) - MY_LOG (1.0 - inclusion_probability_ij);
    }

    // Metropolis-Hastings accept step
    if (MY_LOG (runif(rng)) < log_accept) {
      const int updated_indicator = 1 - indicator(variable1, variable2);
      indicator(variable1, variable2) = updated_indicator;
      indicator(variable2, variable1) = updated_indicator;

      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      const double delta = proposed_state - current_state;

      // Vectorized residual update
      residual_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
      residual_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
    }
  }
}



/**
 * Performs one Gibbs update step for the bgm model.
 *
 * The update sequence depends on the chosen sampling method and the current
 * stage of the warmup schedule:
 *
 *  - Step 0 (initialization):
 *    * If edge selection is enabled and the iteration marks the start of Stage 3c,
 *      initialize a random graph structure via `initialise_graph_bgm()`.
 *
 *  - Step 1 (edge selection):
 *    * If enabled, update edge inclusion indicators with
 *      `update_indicator_edges_metropolis_bgm()`.
 *
 *  - Step 2a (pairwise effects):
 *    * If update_method = "adaptive-metropolis", update interaction weights
 *      for active edges with `update_pairwise_effects_metropolis_bgm()`.
 *
 *  - Step 2b (main effects):
 *    * If update_method = "adaptive-metropolis", update main-effect parameters
 *      with `update_main_effects_metropolis_bgm()`.
 *
 *  - Step 2 (joint updates):
 *    * If update_method = "hamiltonian-mc", update all parameters jointly
 *      with `update_hmc_bgm()`.
 *    * If update_method = "nuts", update all parameters jointly with
 *      `update_nuts_bgm()`, and record diagnostics (tree depth, divergences,
 *      energy) after burn-in.
 *
 *  - Stage-3b proposal tuning:
 *    * Adapt proposal standard deviations for pairwise effects via
 *      `tune_proposal_sd_bgm()`, if the warmup schedule allows.
 *
 * Inputs:
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - proposal_sd_pairwise, proposal_sd_main: Proposal SDs for pairwise and main effects.
 *  - index: Interaction index matrix.
 *  - counts_per_category: Category counts per variable (ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors.
 *  - num_persons, num_variables, num_pairwise, num_main: Model dimensions.
 *  - inclusion_indicator: Symmetric binary inclusion matrix (updated in place).
 *  - pairwise_effects: Symmetric matrix of pairwise effects (updated in place).
 *  - main_effects: Matrix of main effects (updated in place).
 *  - residual_matrix: Residual scores (persons × variables), updated in place.
 *  - inclusion_probability: Prior inclusion probabilities for edges.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - iteration: Current iteration number.
 *  - update_method: Sampling method ("adaptive-metropolis", "hamiltonian-mc", "nuts").
 *  - pairwise_effect_indices: Indexing matrix for pairwise effects.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - hmc_num_leapfrogs: Number of leapfrog steps (HMC).
 *  - nuts_max_depth: Maximum tree depth (NUTS).
 *  - adapt, adapt_main, adapt_pairwise: Adaptation controllers for HMC/NUTS and RWM.
 *  - learn_mass_matrix: If true, adapt mass matrix during warmup.
 *  - schedule: Warmup schedule controller.
 *  - treedepth_samples, divergent_samples, energy_samples: Diagnostic storage (NUTS).
 *  - rng: Random number generator.
 *
 * Outputs:
 *  - main_effects, pairwise_effects, inclusion_indicator, residual_matrix:
 *    updated parameter states.
 *  - proposal_sd_pairwise, proposal_sd_main, adapt, adapt_main, adapt_pairwise:
 *    updated adaptation state.
 *  - treedepth_samples, divergent_samples, energy_samples:
 *    updated with diagnostics if NUTS is used and past burn-in.
 *
 * Notes:
 *  - The function serves as the central per-iteration update step in the Gibbs sampler.
 *  - Different update paths are taken depending on `update_method`.
 *  - Warmup schedule controls both edge selection and proposal-SD tuning.
 */
void gibbs_update_step_bgm (
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double pairwise_scale,
    const arma::mat& pairwise_scaling_factors,
    arma::mat& proposal_sd_pairwise,
    arma::mat& proposal_sd_main,
    const arma::imat& index,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    const double main_alpha,
    const double main_beta,
    const int num_persons,
    const int num_variables,
    const int num_pairwise,
    const int num_main,
    arma::imat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    arma::mat& residual_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const int iteration,
    const UpdateMethod update_method,
    const arma::imat& pairwise_effect_indices,
    arma::imat& pairwise_stats,
    const int hmc_num_leapfrogs,
    const int nuts_max_depth,
    HMCAdaptationController& adapt,
    RWMAdaptationController& adapt_main,
    RWMAdaptationController& adapt_pairwise,
    const bool learn_mass_matrix,
    WarmupSchedule const& schedule,
    arma::ivec& treedepth_samples,
    arma::ivec& divergent_samples,
    arma::vec& energy_samples,
    SafeRNG& rng
) {

  // Step 0: Initialise random graph structure when edge_selection = TRUE
  if (schedule.selection_enabled(iteration) && iteration == schedule.stage3c_start) {
    initialise_graph_bgm(
      inclusion_indicator, pairwise_effects, inclusion_probability, residual_matrix,
      observations, rng
    );
  }

  // Step 1: Edge selection via MH indicator updates (if enabled)
  if (schedule.selection_enabled(iteration)) {
    update_indicator_edges_metropolis_bgm (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise, pairwise_scale, pairwise_scaling_factors, index,
        num_pairwise, num_persons, residual_matrix, inclusion_probability,
        is_ordinal_variable, baseline_category, pairwise_stats,
        rng
    );
  }

  // Step 2a: Update interaction weights for active edges
  if (update_method == adaptive_metropolis) {
    update_pairwise_effects_metropolis_bgm (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise, adapt_pairwise, pairwise_scale,
        pairwise_scaling_factors, num_variables, residual_matrix, is_ordinal_variable, baseline_category,
        iteration, pairwise_stats, rng
    );
  }

  // Step 2b: Update main effect (main_effect) parameters
  if (update_method == adaptive_metropolis) {
    update_main_effects_metropolis_bgm (
        main_effects, observations, num_categories, counts_per_category,
        blume_capel_stats, baseline_category, is_ordinal_variable,
        num_persons, main_alpha, main_beta, residual_matrix,
        proposal_sd_main, adapt_main, iteration,
        rng
    );
  }

  // Step 2: Update joint parameters if applicable
  if (update_method == hamiltonian_mc) {
    update_hmc_bgm(
      main_effects, pairwise_effects, inclusion_indicator, observations,
      num_categories, counts_per_category, blume_capel_stats,
      baseline_category, is_ordinal_variable, main_alpha, main_beta,
      pairwise_scale, pairwise_scaling_factors, residual_matrix, pairwise_stats, hmc_num_leapfrogs,
      iteration, adapt, learn_mass_matrix, schedule.selection_enabled(iteration),
      rng
    );
  } else if (update_method == nuts) {
    SamplerResult result = update_nuts_bgm(
      main_effects, pairwise_effects, inclusion_indicator,
      observations, num_categories, counts_per_category, blume_capel_stats,
      baseline_category, is_ordinal_variable, main_alpha, main_beta,
      pairwise_scale, pairwise_scaling_factors, pairwise_stats, residual_matrix, nuts_max_depth,
      iteration, adapt, learn_mass_matrix, schedule.selection_enabled(iteration),
      rng
    );
    if (iteration >= schedule.total_warmup) {
      int sample_index = iteration - schedule.total_warmup;
      if (auto diag = std::dynamic_pointer_cast<NUTSDiagnostics>(result.diagnostics)) {
        treedepth_samples(sample_index) = diag->tree_depth;
        divergent_samples(sample_index) = diag->divergent ? 1 : 0;
        energy_samples(sample_index) = diag->energy;
      }
    }
  }

  /* --- 2b.  proposal-sd tuning during Stage-3b ------------------------------ */
  tune_proposal_sd_bgm(
    proposal_sd_pairwise, pairwise_effects, main_effects, inclusion_indicator,
    observations, residual_matrix, num_categories, is_ordinal_variable,
    baseline_category, pairwise_scale, pairwise_scaling_factors, pairwise_stats,
    iteration, schedule, rng
  );
}



/**
 * Runs a single Gibbs sampling chain for the bgm model.
 *
 * This function performs one full MCMC chain, including:
 *  - Initialization of parameters, proposal scales, and priors.
 *  - Warmup scheduling and adaptation of proposals or HMC/NUTS settings.
 *  - Gibbs update steps for main effects, pairwise effects, and edge indicators.
 *  - Optional missing data imputation.
 *  - Optional stochastic block model (SBM) prior on edge probabilities.
 *  - Storage of posterior samples and (for NUTS) diagnostic statistics.
 *
 * Inputs:
 *  - chain_id: Numeric identifier for this chain (1-based).
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  - pairwise_scale: Scale parameter of the Cauchy prior on interactions.
 *  - edge_prior: Prior type for edge inclusion ("Beta-Bernoulli" or "Stochastic-Block").
 *  - inclusion_probability: Matrix of prior inclusion probabilities (updated if SBM).
 *  - beta_bernoulli_alpha, beta_bernoulli_beta: Hyperparameters for Beta–Bernoulli edge prior.
 *  - dirichlet_alpha, lambda: Hyperparameters for SBM prior.
 *  - interaction_index_matrix: Index mapping for candidate interactions.
 *  - iter: Number of post-warmup iterations.
 *  - warmup: Number of warmup iterations.
 *  - counts_per_category: Category counts per variable (for ordinal variables).
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - na_impute: If true, perform missing data imputation each iteration.
 *  - missing_index: Locations of missing entries in observations.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - edge_selection: If true, update inclusion indicators during sampling.
 *  - update_method: Sampler type ("adaptive-metropolis", "hamiltonian-mc", "nuts").
 *  - pairwise_effect_indices: Indexing matrix for pairwise effects.
 *  - target_accept: Target acceptance rate for MH/HMC/NUTS updates.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - hmc_num_leapfrogs: Number of leapfrog steps (HMC).
 *  - nuts_max_depth: Maximum tree depth (NUTS).
 *  - learn_mass_matrix: If true, adapt the mass matrix during warmup.
 *  - rng: Random number generator.
 *
 * Returns:
 *  - Rcpp::List with MCMC outputs:
 *    * "main_samples": Matrix of main-effect samples (iter × parameters).
 *    * "pairwise_samples": Matrix of pairwise-effect samples (iter × pairs).
 *    * "indicator_samples": (if edge_selection) Binary inclusion indicators per iteration.
 *    * "allocations": (if SBM prior) Cluster allocations for variables.
 *    * "treedepth__", "divergent__", "energy__": (if NUTS) diagnostic samples.
 *    * "chain_id": ID of this chain.
 *
 * Notes:
 *  - The Gibbs update for each iteration is handled by `gibbs_update_step_bgm()`.
 *  - Proposal scales are adapted during warmup via Robbins–Monro updates.
 *  - Parallel execution across chains is handled by `run_bgm_parallel()`;
 *    this function is for one chain only.
 */
bgmOutput run_gibbs_sampler_bgm(
    int chain_id,
    arma::imat observations,
    const arma::ivec& num_categories,
    const double pairwise_scale,
    const EdgePrior edge_prior,
    arma::mat inclusion_probability,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double beta_bernoulli_alpha_between,
    const double beta_bernoulli_beta_between,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& interaction_index_matrix,
    const int iter,
    const int warmup,
    arma::imat counts_per_category,
    arma::imat blume_capel_stats,
    const double main_alpha,
    const double main_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool edge_selection,
    const UpdateMethod update_method,
    const arma::imat pairwise_effect_indices,
    const double target_accept,
    arma::imat pairwise_stats,
    const int hmc_num_leapfrogs,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    const arma::mat& pairwise_scaling_factors,
    SafeRNG& rng,
    ProgressManager& pm
) {
  // --- Setup: dimensions and storage structures
  const int num_variables = observations.n_cols;
  const int num_persons = observations.n_rows;
  const int max_num_categories = num_categories.max();
  const int num_pairwise = interaction_index_matrix.n_rows;

  // Initialize model parameter matrices
  arma::mat main_effects(num_variables, max_num_categories, arma::fill::zeros);
  arma::mat pairwise_effects(num_variables, num_variables, arma::fill::zeros);
  arma::imat inclusion_indicator(num_variables, num_variables, arma::fill::ones);

  // Residuals used in pseudo-likelihood computation
  arma::mat residual_matrix(num_persons, num_variables, arma::fill::zeros);

  // Allocate optional storage for MCMC samples
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::mat main_effect_samples(iter, num_main);
  arma::mat pairwise_effect_samples(iter, num_pairwise);
  arma::imat indicator_samples;
  arma::imat allocation_samples;

  if (edge_selection) {
    indicator_samples.set_size(iter, num_pairwise);
  }
  if (edge_selection && edge_prior == Stochastic_Block) {
    allocation_samples.set_size(iter, num_variables);
  }

  // For logging nuts performance
  arma::ivec treedepth_samples(iter, arma::fill::zeros);
  arma::ivec divergent_samples(iter, arma::fill::zeros);
  arma::vec energy_samples(iter, arma::fill::zeros);

  // Edge update shuffling setup
  arma::uvec v = arma::regspace<arma::uvec>(0, num_pairwise - 1);
  arma::uvec order(num_pairwise);
  arma::imat index(num_pairwise, 3);

  // SBM-specific structures
  arma::uvec K_values;
  arma::uvec cluster_allocations(num_variables);
  arma::mat cluster_prob(1, 1);
  arma::vec log_Vn(1);

  // --- Initialize SBM prior if applicable
  if (edge_prior == Stochastic_Block) {
    cluster_allocations[0] = 0;
    cluster_allocations[1] = 1;
    for (int i = 2; i < num_variables; i++) {
      cluster_allocations[i] = (runif(rng) > 0.5) ? 1 : 0;
    }

    cluster_prob = block_probs_mfm_sbm(
      cluster_allocations, arma::conv_to<arma::umat>::from(inclusion_indicator),
      num_variables, beta_bernoulli_alpha, beta_bernoulli_beta, beta_bernoulli_alpha_between,
      beta_bernoulli_beta_between, rng
    );

    for (int i = 0; i < num_variables - 1; i++) {
      for (int j = i + 1; j < num_variables; j++) {
        inclusion_probability(i, j) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
        inclusion_probability(j, i) = inclusion_probability(i, j);
      }
    }

    log_Vn = compute_Vn_mfm_sbm(num_variables, dirichlet_alpha, num_variables + 10, lambda);
  }

  // --- Initialize proposal SDs
  arma::mat proposal_sd_main(num_main, max_num_categories, arma::fill::ones);
  arma::mat proposal_sd_pairwise(num_variables, num_variables, arma::fill::ones);

  // --- Optional HMC/NUTS warmup stage
  double initial_step_size_joint = 1.0;
  if (update_method == hamiltonian_mc || update_method == nuts) {
    initial_step_size_joint = find_initial_stepsize_bgm(
      main_effects, pairwise_effects, inclusion_indicator, observations,
      num_categories, counts_per_category, blume_capel_stats,
      baseline_category, is_ordinal_variable, main_alpha, main_beta,
      pairwise_scale, pairwise_scaling_factors, target_accept, pairwise_stats, rng
    );
  }

  // --- Warmup scheduling + adaptation controller
  WarmupSchedule warmup_schedule(warmup, edge_selection, (update_method != adaptive_metropolis));

  HMCAdaptationController adapt_joint(
      num_main + num_pairwise, initial_step_size_joint, target_accept,
      warmup_schedule, learn_mass_matrix
  );
  RWMAdaptationController adapt_main(
      proposal_sd_main, warmup_schedule, target_accept
  );
  RWMAdaptationController adapt_pairwise(
      proposal_sd_pairwise, warmup_schedule, target_accept
  );

  const int total_iter = warmup_schedule.total_warmup + iter;

  bool userInterrupt = false;
  // --- Main Gibbs sampling loop
  for (int iteration = 0; iteration < total_iter; iteration++) {

    pm.update(chain_id - 1);
    if (pm.shouldExit()) {
      userInterrupt = true;
      break;
    }

    // Shuffle update order of edge indices
    order = arma_randperm(rng, num_pairwise);
    for (int i = 0; i < num_pairwise; i++) {
      index.row(i) = interaction_index_matrix.row(order(i));
    }

    // Optional imputation
    if (na_impute) {
      impute_missing_bgm (
          pairwise_effects, main_effects, observations, counts_per_category,
          blume_capel_stats, num_categories, residual_matrix,
          missing_index, is_ordinal_variable, baseline_category,
          pairwise_stats, rng
      );
    }

    // Main Gibbs update step for parameters
    gibbs_update_step_bgm (
        observations, num_categories, pairwise_scale, pairwise_scaling_factors, proposal_sd_pairwise,
        proposal_sd_main, index, counts_per_category, blume_capel_stats,
        main_alpha, main_beta, num_persons, num_variables, num_pairwise,
        num_main, inclusion_indicator, pairwise_effects, main_effects,
        residual_matrix, inclusion_probability, is_ordinal_variable,
        baseline_category, iteration, update_method, pairwise_effect_indices,
        pairwise_stats,
        hmc_num_leapfrogs, nuts_max_depth, adapt_joint, adapt_main, adapt_pairwise,
        learn_mass_matrix, warmup_schedule,
        treedepth_samples, divergent_samples, energy_samples, rng
    );

    // --- Update edge probabilities under the prior (if edge selection is active)
    if (warmup_schedule.selection_enabled(iteration)) {
      if (edge_prior == Beta_Bernoulli) {
        int num_edges_included = 0;
        for (int i = 0; i < num_variables - 1; i++)
          for (int j = i + 1; j < num_variables; j++)
            num_edges_included += inclusion_indicator(i, j);

        double prob = rbeta(rng,
          beta_bernoulli_alpha + num_edges_included,
          beta_bernoulli_beta + num_pairwise - num_edges_included
        );

        for (int i = 0; i < num_variables - 1; i++)
          for (int j = i + 1; j < num_variables; j++)
            inclusion_probability(i, j) = inclusion_probability(j, i) = prob;

      } else if (edge_prior == Stochastic_Block) {
        cluster_allocations = block_allocations_mfm_sbm(
          cluster_allocations, num_variables, log_Vn, cluster_prob,
          arma::conv_to<arma::umat>::from(inclusion_indicator), dirichlet_alpha,
          beta_bernoulli_alpha, beta_bernoulli_beta, beta_bernoulli_alpha_between,
          beta_bernoulli_beta_between, rng
        );

        cluster_prob = block_probs_mfm_sbm(
          cluster_allocations,
          arma::conv_to<arma::umat>::from(inclusion_indicator), num_variables,
          beta_bernoulli_alpha, beta_bernoulli_beta, beta_bernoulli_alpha_between,
          beta_bernoulli_beta_between, rng
        );

        for (int i = 0; i < num_variables - 1; i++) {
          for (int j = i + 1; j < num_variables; j++) {
            inclusion_probability(i, j) = inclusion_probability(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
          }
        }
      }
    }

    // --- Store states
    if (iteration >= warmup_schedule.total_warmup) {
      int sample_index = iteration - warmup_schedule.total_warmup;

      arma::vec vectorized_main = vectorize_main_effects_bgm(main_effects, num_categories, is_ordinal_variable);
      main_effect_samples.row(sample_index) = vectorized_main.t();

      for (int i = 0; i < num_pairwise; i++) {
        int v1 = interaction_index_matrix(i, 1);
        int v2 = interaction_index_matrix(i, 2);
        pairwise_effect_samples(sample_index, i) = pairwise_effects(v1, v2);
      }

      if (edge_selection) {
        for (int i = 0; i < num_pairwise; i++) {
          int v1 = interaction_index_matrix(i, 1);
          int v2 = interaction_index_matrix(i, 2);
          indicator_samples(sample_index, i) = inclusion_indicator(v1, v2);
        }
      }

      if (edge_selection && edge_prior == Stochastic_Block) {
        for (int v = 0; v < num_variables; v++) {
          allocation_samples(sample_index, v) = cluster_allocations[v] + 1;
        }
      }
    }
  }

  bgmOutput chain_result;
  chain_result.chain_id = chain_id;
  chain_result.userInterrupt = userInterrupt;

  chain_result.main_samples = main_effect_samples;
  chain_result.pairwise_samples = pairwise_effect_samples;
  if (update_method == nuts) {
    chain_result.treedepth_samples = treedepth_samples;
    chain_result.divergent_samples = divergent_samples;
    chain_result.energy_samples    = energy_samples;
  }
  if (edge_selection) {
    chain_result.indicator_samples = indicator_samples;
    if (edge_prior == Stochastic_Block)
      chain_result.allocation_samples = allocation_samples;
  }

  return chain_result;
}