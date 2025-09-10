#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bgm_helper.h"
#include "bgm_logp_and_grad.h"
#include "bgm_sampler.h"
#include "common_helpers.h"
#include "mcmc_adaptation.h"
#include "mcmc_hmc.h"
#include "mcmc_leapfrog.h"
#include "mcmc_nuts.h"
#include "mcmc_rwm.h"
#include "mcmc_utils.h"
#include "print_mutex.h"
#include "gibbs_functions_edge_prior.h"
#include "rng_utils.h"

using namespace Rcpp;


/**
 * Function: impute_missing_values_for_graphical_model
 *
 * Imputes missing values in the observation matrix using the current model parameters.
 * For each missing entry, a category is sampled from the posterior predictive distribution,
 * and if the value changes, all dependent structures are updated accordingly.
 *
 * Inputs:
 *  - pairwise_effects: matrix of interaction weights.
 *  - main_effects: matrix of threshold parameters (main effects).
 *  - missing_index: matrix of (person, variable) indices for missing entries.
 *  - num_categories: vector of category counts per variable.
 *  - is_ordinal_variable: logical vector (0/1) indicating ordinal vs. Blume-Capel.
 *  - reference_category: vector of reference categories per variable.
 *
 * Updates (in-place):
 *  - observations: matrix of imputed categorical scores.
 *  - num_obs_categories: counts of observed categories per variable.
 *  - sufficient_blume_capel: sufficient stats for Blume-Capel thresholds.
 *  - rest_matrix: updated linear predictors after imputation.
 */
void impute_missing_values_for_graphical_model (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& observations,
    arma::imat& num_obs_categories,
    arma::imat& sufficient_blume_capel,
    const arma::ivec& num_categories,
    arma::mat& rest_matrix,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    arma::imat& sufficient_pairwise,
    SafeRNG& rng
) {
  const int num_variables = observations.n_cols;
  const int num_missings = missing_index.n_rows;
  const int max_num_categories = num_categories.max ();

  arma::vec category_probabilities (max_num_categories + 1);

  for (int miss = 0; miss < num_missings; miss++) {
    const int person = missing_index (miss, 0);
    const int variable = missing_index (miss, 1);

    const double rest_score = rest_matrix (person, variable);
    const int num_cats = num_categories (variable);
    const bool is_ordinal = is_ordinal_variable (variable);

    double cumsum = 0.0;

    if (is_ordinal) {
      // Compute cumulative unnormalized probabilities for ordinal variable
      cumsum = 1.0;
      category_probabilities[0] = cumsum;
      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const double exponent = main_effects (variable, cat) + score * rest_score;
        cumsum += std::exp (exponent);
        category_probabilities[score] = cumsum;
      }
    } else {
      // Compute probabilities for Blume-Capel variable
      const int ref = reference_category (variable);

      cumsum = std::exp (main_effects (variable, 1) * ref * ref);
      category_probabilities[0] = cumsum;

      for (int cat = 0; cat < num_cats; cat++) {
        const int score = cat + 1;
        const int centered = score - ref;
        const double exponent =
          main_effects (variable, 0) * score +
          main_effects (variable, 1) * centered * centered +
          score * rest_score;
        cumsum += std::exp (exponent);
        category_probabilities[score] = cumsum;
      }
    }

    // Sample from categorical distribution via inverse transform
    const double u = rng () * cumsum;
    int sampled_score = 0;
    while (u > category_probabilities[sampled_score]) {
      sampled_score++;
    }

    const int new_value = sampled_score;
    const int old_value = observations(person, variable);

    if (new_value != old_value) {
      // Update observation matrix
      observations(person, variable) = new_value;

      if (is_ordinal) {
        num_obs_categories(old_value, variable)--;
        num_obs_categories(new_value, variable)++;
      } else {
        const int ref = reference_category(variable);
        const int delta = new_value - old_value;
        const int delta_sq =
          (new_value - ref) * (new_value - ref) -
          (old_value - ref) * (old_value - ref);

        sufficient_blume_capel(0, variable) += delta;
        sufficient_blume_capel(1, variable) += delta_sq;
      }
      // Update residuals across all variables
      for (int var = 0; var < num_variables; var++) {
        const double delta_score = (new_value - old_value) * pairwise_effects(var, variable);
        rest_matrix(person, var) += delta_score;
      }
    }
  }

  // Update sufficient statistics for pairwise effects
  sufficient_pairwise = observations.t() * observations;
  // This could be done more elegantly....
}



/**
 * Function: find_reasonable_initial_step_size
 *
 * Heuristically finds a reasonable initial step size for leapfrog-based MCMC algorithms
 * (such as HMC and NUTS), following the procedure described in:
 *
 *   Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn Sampler: Adaptively Setting
 *   Path Lengths in Hamiltonian Monte Carlo. Journal of Machine Learning Research, 15, 1593–1623.
 *   [Algorithm 4: Heuristic for Choosing an Initial Value of ε]
 *
 * The algorithm simulates a single leapfrog step and compares the resulting change
 * in the log joint density (Hamiltonian). If the proposal is too likely (acceptance too high),
 * the step size is increased. If it is too unlikely (acceptance too low), it is decreased.
 * The process repeats until the log acceptance probability is approximately −log(2),
 * corresponding to a target acceptance probability of ~0.5.
 *
 * Inputs:
 *  - main_effects: Matrix of main (threshold) parameters.
 *  - pairwise_effects: Matrix of pairwise interaction parameters.
 *  - inclusion_indicator: Binary matrix indicating active interactions.
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Vector of category counts per variable.
 *  - num_obs_categories: Observed category counts per variable (matrix).
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters for main effects.
 *  - interaction_scale: Scale of the Cauchy prior on interactions.
 *  - target_acceptance: Desired acceptance probability (typically ~0.65).
 *  - sufficient_pairwise: Sufficient statistics for pairwise interactions.
 *
 * Returns:
 *  - A scalar step size ε that yields roughly the target acceptance probability
 *    under a single leapfrog step.
 *
 * Note:
 *  - This function is suitable for both NUTS and standard HMC algorithms.
 *  - It is typically called once before warm-up/adaptation.
 */
double find_reasonable_initial_step_size(
    arma::mat main_effects,
    arma::mat pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const double interaction_scale,
    const double target_acceptance,
    const arma::imat& sufficient_pairwise,
    SafeRNG& rng
) {
  arma::vec theta = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(theta_vec, current_main, current_pair,
                                 inclusion_indicator,
                                 num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;
    return log_pseudoposterior(
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rm
    );
  };

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(theta_vec, current_main, current_pair, inclusion_indicator,
                                 num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;
    return gradient_log_pseudoposterior(
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rm
    );
  };

  return heuristic_initial_step_size(theta, log_post, grad, rng, target_acceptance);
}



/**
 * Function: update_main_effects_with_metropolis
 *
 * Performs an adaptive Metropolis update of all threshold and Blume-Capel parameters,
 * using Robbins-Monro tuning during the warmup phase.
 *
 * For ordinal variables, all threshold parameters are updated.
 * For categorical variables, two effect parameters are updated.
 * Acceptance probabilities are stored, and proposal SDs are tuned via Robbins-Monro.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold or effect parameters (updated in-place).
 *  - observations: Matrix of categorical responses.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category indicators.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable.
 *  - is_ordinal_variable: Indicator vector (1 = ordinal, 0 = categorical).
 *  - num_persons: Number of individuals.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - rest_matrix: Residual matrix used in log-posterior computation.
 *  - proposal_sd_main: Proposal SD matrix (updated in-place).
 *  - adapter: Robbins-Monro adapter used to update proposal SDs.
 *  - iteration: Current MCMC iteration.
 *
 * Modifies:
 *  - main_effects
 *  - proposal_sd_main
 *  - accept_prob_main
 */
void update_main_effects_with_metropolis (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const int num_persons,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& rest_matrix,
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
          return log_pseudoposterior_thresholds_component(
            main_effects, rest_matrix, num_categories, num_obs_categories,
            sufficient_blume_capel, reference_category, is_ordinal_variable,
            threshold_alpha, threshold_beta, variable, category, -1
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
          return log_pseudoposterior_thresholds_component(
            main_effects, rest_matrix, num_categories, num_obs_categories,
            sufficient_blume_capel, reference_category, is_ordinal_variable,
            threshold_alpha, threshold_beta,
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
 * Function: update_pairwise_effects_with_metropolis
 *
 * Performs adaptive Metropolis-Hastings updates for all active pairwise interactions,
 * using Robbins-Monro tuning during the warmup phase.
 *
 * For each pair (i, j) with inclusion_indicator(i, j) == 1, proposes a new interaction strength.
 * If accepted, updates the interaction matrix and residual matrix.
 * Acceptance probabilities are tracked per parameter and used to adapt proposal SDs.
 *
 * Inputs:
 *  - pairwise_effects: Matrix of interaction parameters (updated in-place).
 *  - main_effects: Matrix of threshold/main effect parameters.
 *  - inclusion_indicator: Binary matrix indicating which interactions are active.
 *  - observations: Matrix of observed categorical responses.
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd_pairwise_effects: Matrix of proposal standard deviations (updated in-place).
 *  - adapter: Robbins-Monro adapter for adaptive tuning (called internally).
 *  - interaction_scale: Scale parameter for the prior.
 *  - num_persons: Number of observations.
 *  - num_variables: Number of variables.
 *  - rest_matrix: Residual matrix used in log-posterior evaluation (updated in-place).
 *  - is_ordinal_variable: Indicator for ordinal variables.
 *  - reference_category: Reference category per variable.
 *  - iteration: Current MCMC iteration.
 *  - sufficient_pairwise: Sufficient statistics for pairwise terms.
 *
 * Modifies:
 *  - pairwise_effects
 *  - rest_matrix
 *  - proposal_sd_pairwise_effects
 *  - accept_prob_pairwise
 */
void update_pairwise_effects_with_metropolis (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    arma::mat& proposal_sd_pairwise_effects,
    RWMAdaptationController& adapter,
    const double interaction_scale,
    const int num_persons,
    const int num_variables,
    arma::mat& rest_matrix,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const int iteration,
    const arma::imat& sufficient_pairwise,
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
            inclusion_indicator, is_ordinal_variable, reference_category,
            interaction_scale, sufficient_pairwise, variable1, variable2
          );
        };

        SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);

        value = result.state[0];
        pairwise_effects(variable2, variable1) = value;

        if(current != value) {
          double delta = value - current;
          rest_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
          rest_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
        }

        accept_prob_pairwise(variable1, variable2) = result.accept_prob;
        index_mask_pairwise(variable1, variable2) = 1;
      }
    }
  }

  adapter.update(index_mask_pairwise, accept_prob_pairwise, iteration);
}



/**
 * Function: update_parameters_with_hmc
 *
 * Performs one HMC (Hamiltonian Monte Carlo) update of the main and pairwise effects
 * using vectorized leapfrog integration, with centralized step size and mass matrix adaptation.
 *
 * Step size adaptation is handled internally by the HMCAdaptationController,
 * which also manages diagonal mass matrix estimation across warmup stages.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - pairwise_effects: Matrix of interaction parameters (updated in-place).
 *  - inclusion_indicator: Binary matrix indicating which interactions are active.
 *  - observations: Matrix of categorical observations.
 *  - num_categories: Vector of category counts per variable.
 *  - num_obs_categories: Matrix of observed category counts.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel thresholds.
 *  - reference_category: Reference category per variable.
 *  - is_ordinal_variable: Indicator vector for ordinal vs. Blume-Capel variables.
 *  - threshold_alpha, threshold_beta: Hyperparameters for main effect priors.
 *  - interaction_scale: Scale parameter for the Cauchy prior on interactions.
 *  - rest_matrix: Matrix of residual scores (updated in-place).
 *  - sufficient_pairwise: Sufficient statistics for pairwise terms.
 *  - num_leapfrogs: Number of leapfrog steps per HMC proposal.
 *  - iteration: Current MCMC iteration.
 *  - adapt: HMC adaptation controller managing step size and mass matrix updates.
 *
 * Modifies:
 *  - main_effects, pairwise_effects: If the HMC proposal is accepted.
 *  - rest_matrix: Recomputed from updated pairwise_effects.
 *  - adapt: Updates step size and optionally the mass matrix during warmup.
 */
void update_parameters_with_hmc(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const double interaction_scale,
    arma::mat& rest_matrix,
    const arma::imat& sufficient_pairwise,
    const int num_leapfrogs,
    const int iteration,
    HMCAdaptationController& adapt,
    const bool learn_mass_matrix,
    const bool selection,
    SafeRNG& rng
) {
  arma::vec current_state = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(theta_vec, current_main, current_pair, inclusion_indicator,
                                 num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;

    return gradient_log_pseudoposterior_active (
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha,
      threshold_beta, interaction_scale, sufficient_pairwise, rm
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(theta_vec, current_main, current_pair, inclusion_indicator,
                                 num_categories, is_ordinal_variable);
    arma::mat rm = observations * current_pair;
    return log_pseudoposterior (
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rm
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
  unvectorize_model_parameters(
    current_state, main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );
  rest_matrix = observations * pairwise_effects;

  adapt.update(current_state, result.accept_prob, iteration);
}



/**
 * Function: update_parameters_with_nuts
 *
 * Performs one update of the main and pairwise effect parameters using
 * the No-U-Turn Sampler (NUTS), with centralized step size and mass matrix adaptation.
 *
 * Step size and mass matrix adaptation are handled via the HMCAdaptationController,
 * which manages warmup phases and dual averaging internally.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - pairwise_effects: Matrix of interaction parameters (updated in-place).
 *  - inclusion_indicator: Binary matrix indicating which interactions are active.
 *  - observations: Matrix of categorical observations.
 *  - num_categories: Vector of category counts per variable.
 *  - num_obs_categories: Matrix of observed category counts.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel thresholds.
 *  - reference_category: Reference category per variable.
 *  - is_ordinal_variable: Logical vector indicating ordinal vs. Blume-Capel (1 = ordinal).
 *  - threshold_alpha, threshold_beta: Hyperparameters for main effect priors.
 *  - interaction_scale: Scale parameter for the Cauchy prior on interactions.
 *  - sufficient_pairwise: Sufficient statistics for pairwise interactions.
 *  - rest_matrix: Matrix of residual scores (observations × variables), updated in-place.
 *  - nuts_max_depth: Maximum tree depth for the NUTS trajectory expansion.
 *  - iteration: Current iteration number.
 *  - adapt: Adaptation controller (step size + mass matrix).
 *
 * Modifies (in-place):
 *  - main_effects, pairwise_effects: Updated if NUTS proposal is accepted.
 *  - rest_matrix: Recomputed from the updated pairwise_effects.
 *  - adapt: Updated with step size and mass matrix changes if within warmup.
 */
SamplerResult update_parameters_with_nuts(
    arma::mat& main_effects,
    arma::mat& pairwise_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const double threshold_alpha,
    const double threshold_beta,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
    arma::mat& rest_matrix,
    const int nuts_max_depth,
    const int iteration,
    HMCAdaptationController& adapt,
    const bool learn_mass_matrix,
    const bool selection,
    SafeRNG& rng
) {
  arma::vec current_state = vectorize_model_parameters(
    main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );

  arma::mat current_main = main_effects;
  arma::mat current_pair = pairwise_effects;

  auto grad = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(theta_vec, current_main, current_pair,
                                 inclusion_indicator, num_categories,
                                 is_ordinal_variable);
    arma::mat rm = observations * current_pair;

    return gradient_log_pseudoposterior_active(
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha,
      threshold_beta, interaction_scale, sufficient_pairwise, rm
    );
  };

  auto log_post = [&](const arma::vec& theta_vec) {
    unvectorize_model_parameters(theta_vec, current_main, current_pair,
                                 inclusion_indicator, num_categories,
                                 is_ordinal_variable);
    arma::mat rm = observations * current_pair;
    return log_pseudoposterior(
      current_main, current_pair, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rm
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
  unvectorize_model_parameters(
    current_state, main_effects, pairwise_effects, inclusion_indicator,
    num_categories, is_ordinal_variable
  );
  rest_matrix = observations * pairwise_effects;

  adapt.update(current_state, result.accept_prob, iteration);

  return result;
}




/* ------------------------------------------------------------------ *
 Tune the Gaussian s.d. for *weight* proposals only.
 Runs   ONLY   when  schedule.adapt_proposal_sd(iter) == true,
 i.e.   Stage-3b,   indicator moves still *OFF*.
 * ------------------------------------------------------------------ */
void tune_pairwise_proposal_sd(
    arma::mat& proposal_sd_pairwise_effects,
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& inclusion_indicator,
    const arma::imat& observations,
    arma::mat& rest_matrix,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale,
    const arma::imat& sufficient_pairwise,
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
          inclusion_indicator, is_ordinal_variable, reference_category,
          interaction_scale, sufficient_pairwise, variable1, variable2
        );
      };

      SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng);

      value = result.state[0];
      pairwise_effects(variable2, variable1) = value;

      if(current != value) {
        double delta = value - current;
        rest_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
        rest_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
      }

      proposal_sd = update_proposal_sd_with_robbins_monro(
        proposal_sd, std::log(result.accept_prob), rm_weight, target_accept
      );

      proposal_sd_pairwise_effects(variable1, variable2) = proposal_sd;
      proposal_sd_pairwise_effects(variable2, variable1) = proposal_sd;
    }
  }
}



/**
 * Function: update_indicator_interaction_pair_with_metropolis
 *
 * Metropolis-Hastings update of pairwise inclusion indicators for a predefined set of edges.
 *
 * Inputs:
 *  - pairwise_effects: Matrix of interaction weights (updated in-place).
 *  - main_effects: Matrix of main effect (threshold) parameters.
 *  - indicator: Matrix of edge inclusion flags (updated in-place).
 *  - observations: Matrix of category scores.
 *  - num_categories: Number of categories per variable.
 *  - proposal_sd: Matrix of proposal standard deviations for pairwise effects.
 *  - interaction_scale: Scale parameter for the Cauchy prior.
 *  - index: List of interaction pairs to update.
 *  - num_interactions: Number of interaction pairs.
 *  - num_persons: Number of observations.
 *  - rest_matrix: Residual scores matrix (updated in-place).
 *  - inclusion_probability: Matrix of prior inclusion probabilities.
 *  - is_ordinal_variable: Logical vector indicating variable type.
 *  - reference_category: Reference category per variable (Blume-Capel).
 *
 * Modifies:
 *  - indicator
 *  - pairwise_effects
 *  - rest_matrix
 */
void update_indicator_interaction_pair_with_metropolis (
    arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    arma::imat& indicator,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& proposal_sd,
    const double interaction_scale,
    const arma::imat& index,
    const int num_interactions,
    const int num_persons,
    arma::mat& rest_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const arma::imat& sufficient_pairwise,
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
      variable1, variable2, proposed_state, current_state, rest_matrix,
      is_ordinal_variable, reference_category, sufficient_pairwise
    );

    // Add prior ratio and proposal correction
    const double inclusion_probability_ij = inclusion_probability(variable1, variable2);
    const double sd = proposal_sd(variable1, variable2);

    if (proposing_addition) {
      log_accept += R::dcauchy(proposed_state, 0.0, interaction_scale, true);
      log_accept -= R::dnorm(proposed_state, current_state, sd, true);
      log_accept += std::log (inclusion_probability_ij) - std::log (1.0 - inclusion_probability_ij);
    } else {
      log_accept -= R::dcauchy(current_state, 0.0, interaction_scale, true);
      log_accept += R::dnorm(current_state, proposed_state, sd, true);
      log_accept -= std::log (inclusion_probability_ij) - std::log (1.0 - inclusion_probability_ij);
    }

    // Metropolis-Hastings accept step
    if (std::log (runif(rng)) < log_accept) {
      const int updated_indicator = 1 - indicator(variable1, variable2);
      indicator(variable1, variable2) = updated_indicator;
      indicator(variable2, variable1) = updated_indicator;

      pairwise_effects(variable1, variable2) = proposed_state;
      pairwise_effects(variable2, variable1) = proposed_state;

      const double delta = proposed_state - current_state;

      // Vectorized residual update
      rest_matrix.col(variable1) += arma::conv_to<arma::vec>::from(observations.col(variable2)) * delta;
      rest_matrix.col(variable2) += arma::conv_to<arma::vec>::from(observations.col(variable1)) * delta;
    }
  }
}



/**
 * Performs a single iteration of the Gibbs sampler for graphical model parameters.
 *
 * This function performs a full Gibbs update sweep over:
 *   1. Inclusion indicators for pairwise interactions (if selection is enabled)
 *   2. Pairwise interaction coefficients (with MALA or adaptive Metropolis)
 *   3. Main effect (threshold) parameters (with MALA, Fisher-MALA, or Metropolis)
 *
 * The interaction updates support optional Fisher preconditioning, and the step sizes
 * are adapted using either dual averaging (during burn-in) or Robbins-Monro.
 *
 * Inputs:
 *  - observations: Matrix of observed categorical scores (persons × variables).
 *  - num_categories: Number of categories for each variable.
 *  - interaction_scale: Cauchy prior scale for pairwise interaction coefficients.
 *  - proposal_sd_pairwise: Proposal SDs for interaction updates (adaptive Metropolis).
 *  - proposal_sd_main: Proposal SDs for threshold updates (Blume-Capel variables).
 *  - index: List of candidate interaction pairs.
 *  - num_obs_categories: Number of observations per category per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - threshold_alpha, threshold_beta: Hyperparameters for main effect priors.
 *  - num_persons: Number of observations.
 *  - num_variables: Number of variables.
 *  - num_pairwise: Number of candidate interaction pairs.
 *  - num_main: Number of main effect parameters.
 *  - inclusion_indicator: Symmetric binary matrix of active interactions (updated).
 *  - pairwise_effects: Symmetric matrix of interaction strengths (updated).
 *  - main_effects: Matrix of threshold parameters (updated).
 *  - rest_matrix: Linear predictor matrix (updated if interaction/main effects change).
 *  - inclusion_probability: Matrix of prior inclusion probabilities.
 *  - rm_decay_rate: Robbins-Monro decay rate (e.g. 0.75).
 *  - is_ordinal_variable: Indicator vector for ordinal variables (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference categories for Blume-Capel variables.
 *  - edge_selection: Whether to update inclusion indicators this iteration.
 *  - step_size_main: Step size for MALA threshold updates (updated).
 *  - iteration: Current iteration number (starts at 0).
 *  - dual_averaging_main: Dual averaging state vector for main effect step size (updated).
 *  - total_burnin: Total number of burn-in iterations.
 *  - use_mala: Whether to use MALA (vs. Metropolis) for updates.
 *  - initial_step_size_main: Initial step size for MALA threshold updates.
 *  - sqrt_inv_fisher_main: Square root inverse Fisher matrix for threshold parameters (updated).
 *  - step_size_pairwise: Step size for interaction MALA updates (updated).
 *  - dual_averaging_pairwise: Dual averaging state for interaction step size (updated).
 *  - initial_step_size_pairwise: Initial step size for interaction MALA updates.
 *  - use_fisher_for_interactions: Whether to use Fisher-preconditioned MALA for interactions.
 *  - sqrt_inv_fisher_pairwise: Square root inverse Fisher matrix for interactions (updated).
 *
 * Updates (in-place):
 *  - inclusion_indicator
 *  - pairwise_effects
 *  - main_effects
 *  - rest_matrix
 *  - step_size_main
 *  - step_size_pairwise
 *  - dual_averaging_main
 *  - dual_averaging_pairwise
 *  - proposal_sd_main
 *  - proposal_sd_pairwise
 *  - sqrt_inv_fisher_main
 *  - sqrt_inv_fisher_pairwise
 */
void gibbs_update_step_for_graphical_model_parameters (
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    arma::mat& proposal_sd_pairwise,
    arma::mat& proposal_sd_main,
    const arma::imat& index,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const int num_persons,
    const int num_variables,
    const int num_pairwise,
    const int num_main,
    arma::imat& inclusion_indicator,
    arma::mat& pairwise_effects,
    arma::mat& main_effects,
    arma::mat& rest_matrix,
    const arma::mat& inclusion_probability,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const int iteration,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    arma::imat& sufficient_pairwise,
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
    initialise_graph(
      inclusion_indicator, pairwise_effects, inclusion_probability, rest_matrix,
      observations, rng
    );
  }

  // Step 1: Edge selection via MH indicator updates (if enabled)
  if (schedule.selection_enabled(iteration)) {
    update_indicator_interaction_pair_with_metropolis (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise, interaction_scale, index,
        num_pairwise, num_persons, rest_matrix, inclusion_probability,
        is_ordinal_variable, reference_category, sufficient_pairwise,
        rng
    );
  }

  // Step 2a: Update interaction weights for active edges
  if (update_method == "adaptive-metropolis") {
    update_pairwise_effects_with_metropolis (
        pairwise_effects, main_effects, inclusion_indicator, observations,
        num_categories, proposal_sd_pairwise, adapt_pairwise, interaction_scale,
        num_persons, num_variables, rest_matrix,
        is_ordinal_variable, reference_category, iteration, sufficient_pairwise,
        rng
    );
  }

  // Step 2b: Update main effect (threshold) parameters
  if (update_method == "adaptive-metropolis") {
    update_main_effects_with_metropolis (
        main_effects, observations, num_categories, num_obs_categories,
        sufficient_blume_capel, reference_category, is_ordinal_variable,
        num_persons, threshold_alpha, threshold_beta, rest_matrix,
        proposal_sd_main, adapt_main, iteration,
        rng
    );
  }

  // Step 2: Update joint parameters if applicable
  if (update_method == "hamiltonian-mc") {
    update_parameters_with_hmc(
      main_effects, pairwise_effects, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, rest_matrix, sufficient_pairwise, hmc_num_leapfrogs,
      iteration, adapt, learn_mass_matrix, schedule.selection_enabled(iteration),
      rng
    );
  } else if (update_method == "nuts") {
    SamplerResult result = update_parameters_with_nuts(
      main_effects, pairwise_effects, inclusion_indicator,
      observations, num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, sufficient_pairwise, rest_matrix, nuts_max_depth,
      iteration, adapt, learn_mass_matrix, schedule.selection_enabled(iteration),
      rng
    );
    if (iteration >= schedule.total_burnin) {
      int sample_index = iteration - schedule.total_burnin;
      if (auto diag = std::dynamic_pointer_cast<NUTSDiagnostics>(result.diagnostics)) {
        treedepth_samples(sample_index) = diag->tree_depth;
        divergent_samples(sample_index) = diag->divergent ? 1 : 0;
        energy_samples(sample_index) = diag->energy;
      }
    }
  }

  /* --- 2b.  proposal-sd tuning during Stage-3b ------------------------------ */
  tune_pairwise_proposal_sd(
    proposal_sd_pairwise, pairwise_effects, main_effects, inclusion_indicator,
    observations, rest_matrix, num_categories, is_ordinal_variable,
    reference_category, interaction_scale, sufficient_pairwise,
    iteration, schedule, rng
  );
}



/**
 * Function: run_gibbs_sampler_for_bgm
 *
 * Runs the full Gibbs sampler for a graphical model with ordinal and/or Blume-Capel variables.
 * Optionally performs edge selection using a Beta-Bernoulli or Stochastic Block prior.
 *
 * During each iteration, the algorithm:
 *  - Updates missing values (if imputation enabled)
 *  - Updates main effects and interactions using adaptive Metropolis, MALA, or Fisher-MALA
 *  - Updates inclusion indicators if edge selection is active
 *  - Adapts proposal variances using Robbins-Monro
 *  - Saves MCMC samples or updates running averages
 *
 * Inputs:
 *  - observations: Imputed categorical data.
 *  - num_categories: Number of categories per variable.
 *  - interaction_scale: Scale for Cauchy prior on pairwise interactions.
 *  - edge_prior: Type of edge prior ("Beta-Bernoulli" or "Stochastic-Block").
 *  - inclusion_probability: Matrix of edge inclusion probabilities (updated if SBM/Beta-Bernoulli).
 *  - beta_bernoulli_alpha, beta_bernoulli_beta: Beta prior parameters for edge inclusion.
 *  - dirichlet_alpha, lambda: SBM prior parameters.
 *  - interaction_index_matrix: Matrix of pairwise edge indices (E × 3).
 *  - iter: Number of post-burn-in iterations.
 *  - burnin: Number of burn-in iterations.
 *  - num_obs_categories: Category counts for each variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - threshold_alpha, threshold_beta: Prior parameters for logistic-Beta.
 *  - na_impute: Whether to impute missing values.
 *  - missing_index: List of missing value locations.
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - reference_category: Reference category for centering (for Blume-Capel).
 *  - save_main, save_pairwise, save_indicator: Whether to save MCMC samples.
 *  - display_progress: Show progress bar during sampling.
 *  - edge_selection: Whether to enable edge selection during burn-in.
 *  - update_method: One of "adaptive-metropolis", "adaptive-mala", or "fisher-mala".
 *    Controls which sampler is used for main effects and interactions.
 *
 * Returns:
 *  - List containing:
 *    - "main": Posterior mean of main effects
 *    - "pairwise": Posterior mean of interaction effects
 *    - "inclusion_indicator": Posterior mean of inclusion matrix
 *    - (optional) "main_samples", "pairwise_samples", "inclusion_indicator_samples"
 *    - (optional) "allocations": Cluster allocations (if SBM used)
 */
Rcpp::List run_gibbs_sampler_for_bgm(
    int chain_id,
    arma::imat observations,
    const arma::ivec& num_categories,
    const double interaction_scale,
    const std::string& edge_prior,
    arma::mat inclusion_probability,
    const double beta_bernoulli_alpha,
    const double beta_bernoulli_beta,
    const double dirichlet_alpha,
    const double lambda,
    const arma::imat& interaction_index_matrix,
    const int iter,
    const int burnin,
    arma::imat num_obs_categories,
    arma::imat sufficient_blume_capel,
    const double threshold_alpha,
    const double threshold_beta,
    const bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat pairwise_effect_indices,
    const double target_accept,
    arma::imat sufficient_pairwise,
    const int hmc_num_leapfrogs,
    const int nuts_max_depth,
    const bool learn_mass_matrix,
    SafeRNG& rng
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
  arma::mat rest_matrix(num_persons, num_variables, arma::fill::zeros);

  // Allocate optional storage for MCMC samples
  const int num_main = count_num_main_effects(num_categories, is_ordinal_variable);
  arma::mat main_effect_samples(iter, num_main);
  arma::mat pairwise_effect_samples(iter, num_pairwise);
  arma::imat indicator_samples;
  arma::imat allocation_samples;

  if (edge_selection) {
    indicator_samples.set_size(iter, num_pairwise);
  }
  if (edge_selection && edge_prior == "Stochastic-Block") {
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
  if (edge_prior == "Stochastic-Block") {
    cluster_allocations[0] = 0;
    cluster_allocations[1] = 1;
    for (int i = 2; i < num_variables; i++) {
      cluster_allocations[i] = (runif(rng) > 0.5) ? 1 : 0;
    }

    cluster_prob = block_probs_mfm_sbm(
      cluster_allocations, arma::conv_to<arma::umat>::from(inclusion_indicator),
      num_variables, beta_bernoulli_alpha, beta_bernoulli_beta, rng
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
  if (update_method == "hamiltonian-mc" || update_method == "nuts") {
    initial_step_size_joint = find_reasonable_initial_step_size(
      main_effects, pairwise_effects, inclusion_indicator, observations,
      num_categories, num_obs_categories, sufficient_blume_capel,
      reference_category, is_ordinal_variable, threshold_alpha, threshold_beta,
      interaction_scale, target_accept, sufficient_pairwise, rng
    );
  }

  // --- Warmup scheduling + adaptation controller
  WarmupSchedule warmup_schedule(burnin, edge_selection, (update_method != "adaptive-metropolis"));
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

  const int total_iter = warmup_schedule.total_burnin + iter;
  const int print_every = std::max(1, total_iter / 10);

  // --- Main Gibbs sampling loop
  for (int iteration = 0; iteration < total_iter; iteration++) {
    if (iteration % print_every == 0) {
      tbb::mutex::scoped_lock lock(get_print_mutex());
      std::cout
      // Rcpp::Rcout
      << "[bgm] chain " << chain_id
      << " iteration " << iteration
      << " / " << total_iter
      << std::endl;
    }

    // Shuffle update order of edge indices
    order = arma_randperm(rng, num_pairwise);
    for (int i = 0; i < num_pairwise; i++) {
      index.row(i) = interaction_index_matrix.row(order(i));
    }

    // Optional imputation
    if (na_impute) {
      impute_missing_values_for_graphical_model (
          pairwise_effects, main_effects, observations, num_obs_categories,
          sufficient_blume_capel, num_categories, rest_matrix,
          missing_index, is_ordinal_variable, reference_category,
          sufficient_pairwise, rng
      );
    }

    // Main Gibbs update step for parameters
    gibbs_update_step_for_graphical_model_parameters (
        observations, num_categories, interaction_scale, proposal_sd_pairwise,
        proposal_sd_main, index, num_obs_categories, sufficient_blume_capel,
        threshold_alpha, threshold_beta, num_persons, num_variables, num_pairwise,
        num_main, inclusion_indicator, pairwise_effects, main_effects,
        rest_matrix, inclusion_probability, is_ordinal_variable,
        reference_category, iteration, update_method, pairwise_effect_indices,
        sufficient_pairwise,
        hmc_num_leapfrogs, nuts_max_depth, adapt_joint, adapt_main, adapt_pairwise,
        learn_mass_matrix, warmup_schedule,
        treedepth_samples, divergent_samples, energy_samples, rng
    );

    // --- Update edge probabilities under the prior (if edge selection is active)
    if (warmup_schedule.selection_enabled(iteration)) {
      if (edge_prior == "Beta-Bernoulli") {
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

      } else if (edge_prior == "Stochastic-Block") {
        cluster_allocations = block_allocations_mfm_sbm(
          cluster_allocations, num_variables, log_Vn, cluster_prob,
          arma::conv_to<arma::umat>::from(inclusion_indicator), dirichlet_alpha,
          beta_bernoulli_alpha, beta_bernoulli_beta, rng
        );

        cluster_prob = block_probs_mfm_sbm(
          cluster_allocations,
          arma::conv_to<arma::umat>::from(inclusion_indicator), num_variables,
          beta_bernoulli_alpha, beta_bernoulli_beta, rng
        );

        for (int i = 0; i < num_variables - 1; i++) {
          for (int j = i + 1; j < num_variables; j++) {
            inclusion_probability(i, j) = inclusion_probability(j, i) = cluster_prob(cluster_allocations[i], cluster_allocations[j]);
          }
        }
      }
    }

    // --- Store states
    if (iteration >= warmup_schedule.total_burnin) {
      int sample_index = iteration - warmup_schedule.total_burnin;

      arma::vec vectorized_main = vectorize_thresholds(main_effects, num_categories, is_ordinal_variable);
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

      if (edge_selection && edge_prior == "Stochastic-Block") {
        for (int v = 0; v < num_variables; v++) {
          allocation_samples(sample_index, v) = cluster_allocations[v] + 1;
        }
      }
    }
  }

  Rcpp::List out;
  out["main_samples"] = main_effect_samples;
  out["pairwise_samples"] = pairwise_effect_samples;

  if (update_method == "nuts") {
    out["treedepth__"] = treedepth_samples;
    out["divergent__"] = divergent_samples;
    out["energy__"] = energy_samples;
  }

  if (edge_selection) {
    out["indicator_samples"] = indicator_samples;
  }

  if (edge_selection && edge_prior == "Stochastic-Block") {
    out["allocations"] = allocation_samples;
  }

  out["chain_id"] = chain_id;
  return out;
}