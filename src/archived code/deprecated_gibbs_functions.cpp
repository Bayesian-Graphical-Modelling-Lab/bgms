// Deprecated: Former function to update main effect parameters of ordinal variables.
// Retained for reference and archival purposes only. Not included in build.
/**
 * Function: update_regular_thresholds_with_metropolis
 *
 * Performs a Metropolis-Hastings update for each threshold of an ordinal variable.
 * Uses a generalized beta-prime proposal and logistic-Beta prior.
 *
 * Inputs:
 *  - main_effects: Matrix of thresholds (updated in-place).
 *  - observations: Matrix of category scores.
 *  - num_categories: Vector of number of categories per variable.
 *  - num_obs_categories: Count matrix of observed scores.
 *  - num_persons: Number of individuals.
 *  - variable: Index of the variable being updated.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - residual_matrix: Residual scores for each observation and variable.
 *
 * Modifies:
 *  - main_effects (only for the specified variable)
 */
void update_regular_thresholds_with_metropolis (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const int num_persons,
    const int variable,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix
) {
  arma::vec g(num_persons);
  arma::vec q(num_persons);
  const int num_cats = num_categories(variable);

  for (int category = 0; category < num_cats; category++) {
    double current = main_effects(variable, category);
    double exp_current = std::exp (current);
    double c = (threshold_alpha + threshold_beta) / (1.0 + exp_current);

    for (int person = 0; person < num_persons; person++) {
      double rest_score = residual_matrix(person, variable);
      double denom = 1.0;
      double numer = std::exp ((category + 1) * rest_score);

      for (int cat = 0; cat < num_cats; cat++) {
        if (cat != category) {
          denom += std::exp (main_effects(variable, cat) + (cat + 1) * rest_score);
        }
      }

      g(person) = denom;
      q(person) = numer;
      c += numer / (denom + numer * exp_current);
    }

    c /= (num_persons + threshold_alpha + threshold_beta - exp_current * c);

    // Sample from generalized beta-prime proposal
    double a = num_obs_categories(category + 1, variable) + threshold_alpha;
    double b = num_persons + threshold_beta - num_obs_categories(category + 1, variable);
    double tmp = R::rbeta(a, b);
    double proposed = std::log (tmp / (1.0 - tmp) / c);
    double exp_proposed = std::exp (proposed);

    // Compute MH acceptance probability
    double log_acceptance_probability = 0.0;
    for (int person = 0; person < num_persons; person++) {
      log_acceptance_probability += std::log (g(person) + q(person) * exp_current);
      log_acceptance_probability -= std::log (g(person) + q(person) * exp_proposed);
    }

    log_acceptance_probability -= (threshold_alpha + threshold_beta) * std::log1p(exp_proposed);
    log_acceptance_probability += (threshold_alpha + threshold_beta) * std::log1p(exp_current);
    log_acceptance_probability -= (a + b) * std::log1p(c * exp_current);
    log_acceptance_probability += (a + b) * std::log1p(c * exp_proposed);

    if (std::log (R::unif_rand()) < log_acceptance_probability) {
      main_effects(variable, category) = proposed;
    }
  }
}



// Deprecated: Former function to update main effect parameters of Blume-Capel variables.
// Retained for reference and archival purposes only. Not included in build.

/**
 * Function: update_blumecapel_thresholds_with_adaptive_metropolis
 *
 * Performs an adaptive Metropolis update of the Blume-Capel threshold parameters
 * (linear and quadratic) for a single variable, with Robbins-Monro tuning.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - observations: Matrix of categorical scores.
 *  - num_categories: Number of categories per variable.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - num_persons: Number of observations.
 *  - variable: Index of the variable being updated.
 *  - reference_category: Reference category per variable.
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *  - residual_matrix: Residual scores.
 *  - proposal_sd_blumecapel: Matrix of proposal SDs for each variable (updated in-place).
 *  - exp_neg_log_t_rm_adaptation_rate: Robbins-Monro adaptation weight.
 *
 * Modifies:
 *  - main_effects (for the given variable)
 *  - proposal_sd_blumecapel
 */
void update_blumecapel_thresholds_with_adaptive_metropolis (
    arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& sufficient_blume_capel,
    const int num_persons,
    const int variable,
    const arma::ivec& reference_category,
    const double threshold_alpha,
    const double threshold_beta,
    const arma::mat& residual_matrix,
    arma::mat& proposal_sd_blumecapel,
    const double exp_neg_log_t_rm_adaptation_rate
) {
  const int num_cats = num_categories(variable);
  const int ref = reference_category(variable);

  // --- Define helper for prior contribution
  auto log_beta_prior_diff = [&](double curr, double prop) {
    return (threshold_alpha + threshold_beta) *
      (std::log1p(std::exp (curr)) - std::log1p(std::exp (prop)));
  };

  // --- Update each threshold parameter: 0 = linear, 1 = quadratic
  for (int param = 0; param < 2; param++) {
    double& proposal_sd = proposal_sd_blumecapel(variable, param);
    double current = main_effects(variable, param);
    double proposed = R::rnorm(current, proposal_sd);
    double diff = proposed - current;

    arma::vec numer_current(num_cats + 1);
    arma::vec numer_proposed(num_cats + 1);

    // --- Step 1: Construct numerators for softmax (for all categories)
    for (int cat = 0; cat <= num_cats; cat++) {
      int centered = cat - ref;
      if (param == 0) {
        // Linear update
        double quad = main_effects(variable, 1) * centered * centered;
        numer_current(cat) = current * cat + quad;
        numer_proposed(cat) = proposed * cat + quad;
      } else {
        // Quadratic update
        double lin = main_effects(variable, 0) * cat;
        numer_current(cat) = current * centered * centered + lin;
        numer_proposed(cat) = proposed * centered * centered + lin;
      }
    }

    // --- Step 2: Compute lbound for numerical stability
    double max_curr = numer_current.max();
    double max_prop = numer_proposed.max();
    double lbound = (max_curr > 0.0 || max_prop > 0.0) ? std::max(max_curr, max_prop) : 0.0;

    // --- Step 3: Likelihood ratio
    // Accumulate log acceptance probability based on change in pseudo-likelihood

    // Contribution from sufficient statistics and prior
    double log_accept = diff * (threshold_alpha + sufficient_blume_capel(param, variable));

    // Vectorized likelihood ratio for all persons:
    //
    // For each person, compute:
    //   log p(y_i | proposed) - log p(y_i | current)
    //   using softmax-style normalization for categorical probabilities.
    //
    // The bound stabilizes exponentials across categories and persons.
    arma::vec rest_score = residual_matrix.col(variable);                       // Person-wise residuals
    arma::vec bound = arma::max(rest_score, arma::zeros<arma::vec>(num_persons)) * num_cats + lbound;

    arma::vec denom_curr = arma::exp (numer_current(0) - bound);                 // Score = 0 contribution
    arma::vec denom_prop = arma::exp (numer_proposed(0) - bound);

    for (int cat = 0; cat < num_cats; cat++) {
      arma::vec score_term = (cat + 1) * rest_score - bound;

      // Compute exponentials for each category and add to denominator
      denom_curr += arma::exp (numer_current(cat + 1) + score_term);
      denom_prop += arma::exp (numer_proposed(cat + 1) + score_term);
    }

    // Accumulate the person-wise log ratio contributions
    log_accept += arma::accu (arma::log (denom_curr) - arma::log (denom_prop));

    // --- Step 4: Add prior ratio
    log_accept += log_beta_prior_diff(current, proposed);

    // --- Step 5: Metropolis accept/reject
    if (std::log (R::unif_rand()) < log_accept) {
      main_effects(variable, param) = proposed;
    }

    // --- Step 6: Robbins-Monro proposal adaptation
    proposal_sd = update_proposal_sd_with_robbins_monro (
      proposal_sd, log_accept, exp_neg_log_t_rm_adaptation_rate
    );
  }
}



// Deprecated: Former function to compute gradient vector for the interactions.
// Retained for reference and archival purposes only. Not included in build.

/**
 * Function: gradient_log_pseudoposterior_interactions
 *
 * Computes the gradient of the log pseudoposterior with respect to the
 * active interaction parameters. This is used in MALA updates of the
 * pairwise interaction matrix.
 *
 * Inputs:
 *  - pairwise_effects: Symmetric matrix of interaction parameters [V × V].
 *  - main_effects: Matrix of main effect (threshold) parameters [V × max_categories].
 *  - observations: Matrix of ordinal and Blume-Capel scores (row = individual, col = variable).
 *  - num_categories: Vector of number of categories per variable.
 *  - inclusion_indicator: Symmetric binary matrix indicating which interactions are active.
 *  - is_ordinal_variable: Logical vector indicating ordinal (1) or Blume-Capel (0) variables.
 *  - reference_category: Vector of reference categories for BC variables.
 *  - interaction_scale: Cauchy prior scale on interaction weights.
 *
 * Returns:
 *  - Gradient vector for all pairwise interactions (in upper-triangle order).
 */
arma::vec gradient_log_pseudoposterior_interactions (
    const arma::mat& pairwise_effects,
    const arma::mat& main_effects,
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::imat& inclusion_indicator,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    const double interaction_scale
) {
  const int num_variables = observations.n_cols;
  const int num_observations = observations.n_rows;
  const int num_interactions = (num_variables * (num_variables - 1)) / 2;

  arma::vec gradient (num_interactions, arma::fill::zeros);
  int interaction_index = -1;

  for (int var1 = 0; var1 < num_variables - 1; var1++) {
    for (int var2 = var1 + 1; var2 < num_variables; var2++) {
      interaction_index++;

      if (inclusion_indicator (var1, var2) == 0)
        continue;

      // Convert observed scores to integer vectors
      const arma::ivec responses_var1 = arma::conv_to<arma::ivec>::from (observations.col (var1));
      const arma::ivec responses_var2 = arma::conv_to<arma::ivec>::from (observations.col (var2));


      // First-order gradient term from complete data
      gradient (interaction_index) = 2.0 * arma::dot (responses_var1, responses_var2);

      // --- Contribution from variable var1
      int num_categories_var1 = num_categories (var1);
      arma::vec rest_scores = observations * pairwise_effects.col (var1);
      arma::vec numerator = arma::zeros (num_observations);
      arma::vec denominator = arma::zeros (num_observations);
      arma::vec bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var1;

      if (is_ordinal_variable (var1)) {
        denominator += arma::exp ( -bounds );
        for (int category = 0; category < num_categories_var1; category++) {
          arma::vec exponent = main_effects (var1, category) + (category + 1) * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += (category + 1) * responses_var2 % weight;
        }
      } else {
        const int ref_cat = reference_category (var1);
        for (int category = 0; category <= num_categories_var1; category++) {
          int centered_cat = category - ref_cat;
          double lin_term = main_effects (var1, 0) * category;
          double quad_term = main_effects (var1, 1) * centered_cat * centered_cat;
          arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += category * responses_var2 % weight;
        }
      }

      gradient (interaction_index) -= arma::accu (numerator / denominator);

      // --- Contribution from variable var2
      int num_categories_var2 = num_categories (var2);
      rest_scores = observations * pairwise_effects.col (var2);
      numerator.zeros ();
      denominator.zeros ();
      bounds = arma::max (rest_scores, arma::zeros<arma::vec> (num_observations)) * num_categories_var2;

      if (is_ordinal_variable (var2)) {
        denominator += arma::exp ( -bounds );
        for (int category = 0; category < num_categories_var2; category++) {
          arma::vec exponent = main_effects (var2, category) + (category + 1) * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += (category + 1) * responses_var1 % weight;
        }
      } else {
        const int ref_cat = reference_category (var2);
        for (int category = 0; category <= num_categories_var2; category++) {
          int centered_cat = category - ref_cat;
          double lin_term = main_effects (var2, 0) * category;
          double quad_term = main_effects (var2, 1) * centered_cat * centered_cat;
          arma::vec exponent = lin_term + quad_term + category * rest_scores - bounds;
          arma::vec weight = arma::exp (exponent);
          denominator += weight;
          numerator += category * responses_var1 % weight;
        }
      }

      gradient (interaction_index) -= arma::accu (numerator / denominator);

      // ---- Gradient contribution from Cauchy prior
      const double effect = pairwise_effects (var1, var2);
      gradient (interaction_index) -= 2.0 * effect / (effect * effect + interaction_scale * interaction_scale);
    }
  }

  return gradient;
}



// Deprecated: Former threshold update function from Gibbs sampler
// Retained for reference and archival purposes only. Not included in build.

/**
 * Function: update_thresholds_with_adaptive_mala
 *
 * Performs a MALA update of threshold parameters with adaptive step size tuning.
 * Applies dual averaging during burn-in and Robbins-Monro afterward.
 *
 * Inputs:
 *  - main_effects: Matrix of threshold parameters (updated in-place).
 *  - step_size_mala: MALA step size (updated in-place).
 *  - residual_matrix: Residual scores for each observation and variable.
 *  - num_categories: Number of categories per variable.
 *  - num_obs_categories: Observed category count matrix.
 *  - sufficient_blume_capel: Sufficient statistics for Blume-Capel variables.
 *  - reference_category: Reference category per variable (for Blume-Capel).
 *  - is_ordinal_variable: Logical vector (1 = ordinal, 0 = Blume-Capel).
 *  - iteration: Current iteration number.
 *  - burnin: Total number of burn-in iterations.
 *  - dual_averaging_state: Dual averaging state vector [log_eps, log_eps_avg, H_bar] (updated in-place).
 *  - threshold_alpha, threshold_beta: Prior hyperparameters.
 *
 * Modifies:
 *  - main_effects
 *  - step_size_mala
 *  - dual_averaging_state
 */
void update_thresholds_with_adaptive_mala (
    arma::mat& main_effects,
    double& step_size_mala,
    const arma::mat& residual_matrix,
    const arma::ivec& num_categories,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    const arma::ivec& reference_category,
    const arma::uvec& is_ordinal_variable,
    const int iteration,
    const int burnin,
    arma::vec& dual_averaging_state,
    const double threshold_alpha,
    const double threshold_beta,
    const double initial_step_size_mala
) {
  // --- Step 1: Flatten current parameters and compute gradient & posterior
  arma::vec flat_theta = vectorize_thresholds (
    main_effects, num_categories, is_ordinal_variable
  );
  arma::vec grad = gradient_log_pseudoposterior_thresholds (
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );
  const double log_post_current = log_pseudoposterior_thresholds (
    main_effects, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // --- Step 2: Propose new parameters using MALA
  const double sqrt_step = std::sqrt(step_size_mala);
  arma::vec proposal = flat_theta + 0.5 * step_size_mala * grad + sqrt_step * arma::randn(flat_theta.n_elem);
  arma::mat proposed_thresholds = reconstruct_threshold_matrix (
    proposal, num_categories, is_ordinal_variable
  );

  // --- Step 3: Evaluate proposed state
  const double log_post_proposal = log_pseudoposterior_thresholds (
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );
  arma::vec grad_proposal = gradient_log_pseudoposterior_thresholds (
    proposed_thresholds, residual_matrix, num_categories, num_obs_categories,
    sufficient_blume_capel, reference_category, is_ordinal_variable,
    threshold_alpha, threshold_beta
  );

  // --- Step 4: Compute forward and backward proposal densities
  const arma::vec forward_mean = flat_theta + 0.5 * step_size_mala * grad;
  const arma::vec backward_mean = proposal + 0.5 * step_size_mala * grad_proposal;

  const double log_forward = -0.5 / step_size_mala * arma::accu(arma::square(proposal - forward_mean));
  const double log_backward = -0.5 / step_size_mala * arma::accu(arma::square(flat_theta - backward_mean));

  // --- Step 5: Accept/reject
  const double log_acceptance = log_post_proposal + log_backward - log_post_current - log_forward;
  if (std::log(R::unif_rand()) < log_acceptance) {
    main_effects = proposed_thresholds;
  }

  const double accept_prob = std::min(1.0, std::exp(log_acceptance));

  // --- Step 6: Adapt step size
  if (iteration <= burnin) {
    update_step_size_with_dual_averaging (
        initial_step_size_mala, accept_prob, iteration + 1, dual_averaging_state);
    step_size_mala = std::exp(dual_averaging_state[1]);
  } else {
    update_step_size_with_robbins_monro(accept_prob, iteration - burnin, step_size_mala);
  }
}

update_thresholds_with_adaptive_mala(
  main_effects, step_size_mala, residual_matrix, num_categories,
  num_obs_categories, sufficient_blume_capel, reference_category,
  is_ordinal_variable, iteration, total_burnin, dual_averaging_state,
  threshold_alpha, threshold_beta, initial_step_size_mala
);