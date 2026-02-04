#include <RcppArmadillo.h>
#include "omrf_model.h"
#include "adaptiveMetropolis.h"
#include "rng/rng_utils.h"
#include "mcmc/mcmc_hmc.h"
#include "mcmc/mcmc_nuts.h"
#include "mcmc/mcmc_rwm.h"
#include "mcmc/mcmc_utils.h"
#include "mcmc/mcmc_adaptation.h"
#include "mcmc/mcmc_runner.h"
#include "math/explog_switch.h"
#include "utils/common_helpers.h"
#include "utils/variable_helpers.h"


// =============================================================================
// Constructor
// =============================================================================

OMRFModel::OMRFModel(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    double main_alpha,
    double main_beta,
    double pairwise_scale,
    bool edge_selection
) :
    n_(observations.n_rows),
    p_(observations.n_cols),
    observations_(observations),
    num_categories_(num_categories),
    is_ordinal_variable_(is_ordinal_variable),
    baseline_category_(baseline_category),
    inclusion_probability_(inclusion_probability),
    main_alpha_(main_alpha),
    main_beta_(main_beta),
    pairwise_scale_(pairwise_scale),
    edge_selection_(edge_selection),
    edge_selection_active_(false),
    proposal_(AdaptiveProposal(1, 500)),  // Will be resized later
    step_size_(0.1),
    has_missing_(false),
    gradient_cache_valid_(false)
{
    // Initialize parameter dimensions
    num_main_ = count_num_main_effects_internal();
    num_pairwise_ = (p_ * (p_ - 1)) / 2;

    // Initialize parameters
    int max_cats = num_categories_.max();
    main_effects_ = arma::zeros<arma::mat>(p_, max_cats);
    pairwise_effects_ = arma::zeros<arma::mat>(p_, p_);
    edge_indicators_ = initial_edge_indicators;

    // Initialize proposal SDs
    proposal_sd_main_ = arma::ones<arma::mat>(p_, max_cats) * 0.5;
    proposal_sd_pairwise_ = arma::ones<arma::mat>(p_, p_) * 0.5;

    // Initialize adaptive proposal
    proposal_ = AdaptiveProposal(num_main_ + num_pairwise_, 500);

    // Initialize mass matrix
    inv_mass_ = arma::ones<arma::vec>(num_main_ + num_pairwise_);

    // Pre-compute observations as double (for efficient matrix operations)
    observations_double_ = arma::conv_to<arma::mat>::from(observations_);

    // Compute sufficient statistics
    compute_sufficient_statistics();

    // Initialize residual matrix
    update_residual_matrix();

    // Build interaction index
    build_interaction_index();
}


// =============================================================================
// Copy constructor
// =============================================================================

OMRFModel::OMRFModel(const OMRFModel& other)
    : BaseModel(other),
      n_(other.n_),
      p_(other.p_),
      observations_(other.observations_),
      observations_double_(other.observations_double_),
      num_categories_(other.num_categories_),
      is_ordinal_variable_(other.is_ordinal_variable_),
      baseline_category_(other.baseline_category_),
      counts_per_category_(other.counts_per_category_),
      blume_capel_stats_(other.blume_capel_stats_),
      pairwise_stats_(other.pairwise_stats_),
      residual_matrix_(other.residual_matrix_),
      main_effects_(other.main_effects_),
      pairwise_effects_(other.pairwise_effects_),
      edge_indicators_(other.edge_indicators_),
      inclusion_probability_(other.inclusion_probability_),
      main_alpha_(other.main_alpha_),
      main_beta_(other.main_beta_),
      pairwise_scale_(other.pairwise_scale_),
      edge_selection_(other.edge_selection_),
      edge_selection_active_(other.edge_selection_active_),
      num_main_(other.num_main_),
      num_pairwise_(other.num_pairwise_),
      proposal_(other.proposal_),
      proposal_sd_main_(other.proposal_sd_main_),
      proposal_sd_pairwise_(other.proposal_sd_pairwise_),
      rng_(other.rng_),
      step_size_(other.step_size_),
      inv_mass_(other.inv_mass_),
      has_missing_(other.has_missing_),
      missing_index_(other.missing_index_),
      grad_obs_cache_(other.grad_obs_cache_),
      index_matrix_cache_(other.index_matrix_cache_),
      gradient_cache_valid_(other.gradient_cache_valid_),
      interaction_index_(other.interaction_index_)
{
}


// =============================================================================
// Sufficient statistics computation
// =============================================================================

void OMRFModel::compute_sufficient_statistics() {
    int max_cats = num_categories_.max();

    // Category counts for ordinal variables
    counts_per_category_ = arma::zeros<arma::imat>(max_cats + 1, p_);
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            for (size_t i = 0; i < n_; ++i) {
                int cat = observations_(i, v);
                if (cat >= 0 && cat <= num_categories_(v)) {
                    counts_per_category_(cat, v)++;
                }
            }
        }
    }

    // Blume-Capel statistics (linear and quadratic sums)
    blume_capel_stats_ = arma::zeros<arma::imat>(2, p_);
    for (size_t v = 0; v < p_; ++v) {
        if (!is_ordinal_variable_(v)) {
            int baseline = baseline_category_(v);
            for (size_t i = 0; i < n_; ++i) {
                int s = observations_(i, v) - baseline;
                blume_capel_stats_(0, v) += s;      // linear
                blume_capel_stats_(1, v) += s * s;  // quadratic
            }
        }
    }

    // Pairwise statistics (X^T X) - use pre-computed transformed observations
    arma::mat ps = observations_double_.t() * observations_double_;
    pairwise_stats_ = arma::conv_to<arma::imat>::from(ps);
}


size_t OMRFModel::count_num_main_effects_internal() const {
    size_t count = 0;
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            count += num_categories_(v);
        } else {
            count += 2;  // linear and quadratic for Blume-Capel
        }
    }
    return count;
}


void OMRFModel::build_interaction_index() {
    interaction_index_ = arma::zeros<arma::imat>(num_pairwise_, 3);
    int idx = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            interaction_index_(idx, 0) = idx;
            interaction_index_(idx, 1) = v1;
            interaction_index_(idx, 2) = v2;
            idx++;
        }
    }
}


void OMRFModel::update_residual_matrix() {
    // Use pre-computed transformed observations (computed once in constructor)
    residual_matrix_ = observations_double_ * pairwise_effects_;
}


void OMRFModel::set_pairwise_effects(const arma::mat& pairwise_effects) {
    pairwise_effects_ = pairwise_effects;
    update_residual_matrix();
    invalidate_gradient_cache();
}


// =============================================================================
// BaseModel interface implementation
// =============================================================================

size_t OMRFModel::parameter_dimension() const {
    // Count active parameters: main effects + included pairwise effects
    size_t active = num_main_;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                active++;
            }
        }
    }
    return active;
}


void OMRFModel::set_seed(int seed) {
    rng_ = SafeRNG(seed);
}


std::unique_ptr<BaseModel> OMRFModel::clone() const {
    return std::make_unique<OMRFModel>(*this);
}


void OMRFModel::set_adaptive_proposal(AdaptiveProposal proposal) {
    proposal_ = proposal;
}


// =============================================================================
// Parameter vectorization
// =============================================================================

arma::vec OMRFModel::vectorize_parameters() const {
    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    arma::vec param_vec(num_main_ + num_active);
    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                param_vec(offset++) = main_effects_(v, c);
            }
        } else {
            param_vec(offset++) = main_effects_(v, 0);  // linear
            param_vec(offset++) = main_effects_(v, 1);  // quadratic
        }
    }

    // Active pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                param_vec(offset++) = pairwise_effects_(v1, v2);
            }
        }
    }

    return param_vec;
}


void OMRFModel::unvectorize_parameters(const arma::vec& param_vec) {
    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                main_effects_(v, c) = param_vec(offset++);
            }
        } else {
            main_effects_(v, 0) = param_vec(offset++);  // linear
            main_effects_(v, 1) = param_vec(offset++);  // quadratic
        }
    }

    // Active pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                double val = param_vec(offset++);
                pairwise_effects_(v1, v2) = val;
                pairwise_effects_(v2, v1) = val;
            }
        }
    }

    update_residual_matrix();
    invalidate_gradient_cache();
}


arma::vec OMRFModel::get_vectorized_parameters() const {
    return vectorize_parameters();
}


void OMRFModel::set_vectorized_parameters(const arma::vec& parameters) {
    unvectorize_parameters(parameters);
}


arma::vec OMRFModel::get_full_vectorized_parameters() const {
    // Fixed-size vector: all main effects + ALL pairwise effects
    arma::vec param_vec(num_main_ + num_pairwise_);
    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                param_vec(offset++) = main_effects_(v, c);
            }
        } else {
            param_vec(offset++) = main_effects_(v, 0);  // linear
            param_vec(offset++) = main_effects_(v, 1);  // quadratic
        }
    }

    // ALL pairwise effects (zeros for inactive edges)
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            param_vec(offset++) = pairwise_effects_(v1, v2);
        }
    }

    return param_vec;
}


arma::ivec OMRFModel::get_vectorized_indicator_parameters() {
    arma::ivec indicators(num_pairwise_);
    int idx = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            indicators(idx++) = edge_indicators_(v1, v2);
        }
    }
    return indicators;
}


arma::vec OMRFModel::get_active_inv_mass() const {
    if (!edge_selection_active_) {
        return inv_mass_;
    }

    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    arma::vec active_inv_mass(num_main_ + num_active);
    active_inv_mass.head(num_main_) = inv_mass_.head(num_main_);

    int offset_full = num_main_;
    int offset_active = num_main_;

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                active_inv_mass(offset_active) = inv_mass_(offset_full);
                offset_active++;
            }
            offset_full++;
        }
    }

    return active_inv_mass;
}


void OMRFModel::vectorize_parameters_into(arma::vec& param_vec) const {
    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    // Resize if needed (should rarely happen after first call)
    size_t needed_size = num_main_ + num_active;
    if (param_vec.n_elem != needed_size) {
        param_vec.set_size(needed_size);
    }

    int offset = 0;

    // Main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                param_vec(offset++) = main_effects_(v, c);
            }
        } else {
            param_vec(offset++) = main_effects_(v, 0);  // linear
            param_vec(offset++) = main_effects_(v, 1);  // quadratic
        }
    }

    // Active pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                param_vec(offset++) = pairwise_effects_(v1, v2);
            }
        }
    }
}


void OMRFModel::get_active_inv_mass_into(arma::vec& active_inv_mass) const {
    if (!edge_selection_active_) {
        // No edge selection - just use full inv_mass
        if (active_inv_mass.n_elem != inv_mass_.n_elem) {
            active_inv_mass.set_size(inv_mass_.n_elem);
        }
        active_inv_mass = inv_mass_;
        return;
    }

    // Count active parameters
    int num_active = 0;
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                num_active++;
            }
        }
    }

    size_t needed_size = num_main_ + num_active;
    if (active_inv_mass.n_elem != needed_size) {
        active_inv_mass.set_size(needed_size);
    }

    active_inv_mass.head(num_main_) = inv_mass_.head(num_main_);

    int offset_full = num_main_;
    int offset_active = num_main_;

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                active_inv_mass(offset_active) = inv_mass_(offset_full);
                offset_active++;
            }
            offset_full++;
        }
    }
}


// =============================================================================
// Log-pseudoposterior computation
// =============================================================================

double OMRFModel::logp(const arma::vec& parameters) {
    // Unvectorize into temporary matrices (safe approach)
    arma::mat temp_main = main_effects_;
    arma::mat temp_pairwise = pairwise_effects_;

    // Unvectorize parameters into temporaries
    int offset = 0;
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                temp_main(v, c) = parameters(offset++);
            }
        } else {
            temp_main(v, 0) = parameters(offset++);
            temp_main(v, 1) = parameters(offset++);
        }
    }

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                temp_pairwise(v1, v2) = parameters(offset++);
                temp_pairwise(v2, v1) = temp_pairwise(v1, v2);
            }
        }
    }

    // Compute residual matrix from temp_pairwise
    arma::mat temp_residual = arma::conv_to<arma::mat>::from(observations_) * temp_pairwise;

    // Compute log-posterior with temporaries
    return log_pseudoposterior_with_state(temp_main, temp_pairwise, temp_residual);
}


double OMRFModel::log_pseudoposterior_with_state(
    const arma::mat& main_eff,
    const arma::mat& pairwise_eff,
    const arma::mat& residual_mat
) const {
    double log_post = 0.0;

    auto log_beta_prior = [this](double x) {
        return x * main_alpha_ - std::log1p(std::exp(x)) * (main_alpha_ + main_beta_);
    };

    // Main effect contributions (priors and sufficient statistics)
    for (size_t v = 0; v < p_; ++v) {
        int num_cats = num_categories_(v);

        if (is_ordinal_variable_(v)) {
            for (int c = 0; c < num_cats; ++c) {
                log_post += log_beta_prior(main_eff(v, c));
                log_post += main_eff(v, c) * counts_per_category_(c + 1, v);
            }
        } else {
            log_post += log_beta_prior(main_eff(v, 0));
            log_post += log_beta_prior(main_eff(v, 1));
            log_post += main_eff(v, 0) * blume_capel_stats_(0, v);
            log_post += main_eff(v, 1) * blume_capel_stats_(1, v);
        }
    }

    // Log-denominator contributions using vectorized helpers
    for (size_t v = 0; v < p_; ++v) {
        int num_cats = num_categories_(v);
        arma::vec residual_score = residual_mat.col(v);
        arma::vec bound = num_cats * residual_score;

        arma::vec denom(n_, arma::fill::zeros);
        if (is_ordinal_variable_(v)) {
            // Extract main effect parameters for this variable
            arma::vec main_effect_param = main_eff.row(v).cols(0, num_cats - 1).t();
            denom = compute_denom_ordinal(residual_score, main_effect_param, bound);
        } else {
            int ref = baseline_category_(v);
            double lin_effect = main_eff(v, 0);
            double quad_effect = main_eff(v, 1);
            // This updates bound in-place
            denom = compute_denom_blume_capel(residual_score, lin_effect, quad_effect, ref, num_cats, bound);
        }
        log_post -= arma::accu(bound + ARMA_MY_LOG(denom));
    }

    // Pairwise effect contributions: sufficient statistics + Cauchy prior
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                double effect = pairwise_eff(v1, v2);
                // Sufficient statistics term (data likelihood contribution)
                log_post += 2.0 * pairwise_stats_(v1, v2) * effect;
                // Cauchy prior using R's dcauchy for consistency
                log_post += R::dcauchy(effect, 0.0, pairwise_scale_, true);
            }
        }
    }

    return log_post;
}


double OMRFModel::log_pseudoposterior_internal() const {
    return log_pseudoposterior_with_state(main_effects_, pairwise_effects_, residual_matrix_);
}


double OMRFModel::log_pseudoposterior_main_component(int variable, int category, int parameter) const {
    double log_post = 0.0;

    // Lambda for Beta-prime prior on main effects (matches original implementation)
    // log p(theta) = alpha * theta - (alpha + beta) * log(1 + exp(theta))
    auto log_beta_prior = [this](double x) {
        return x * main_alpha_ - std::log1p(std::exp(x)) * (main_alpha_ + main_beta_);
    };

    int num_cats = num_categories_(variable);
    arma::vec bound = num_cats * residual_matrix_.col(variable);

    if (is_ordinal_variable_(variable)) {
        // Ordinal variable: use category
        log_post += log_beta_prior(main_effects_(variable, category));
        log_post += main_effects_(variable, category) * counts_per_category_(category + 1, variable);

        // Log-denominator contribution
        for (size_t i = 0; i < n_; ++i) {
            double max_val = 0.0;
            for (int c = 0; c < num_cats; ++c) {
                double val = main_effects_(variable, c) + (c + 1) * residual_matrix_(i, variable);
                if (val > max_val) max_val = val;
            }

            double denom = std::exp(-max_val);
            for (int c = 0; c < num_cats; ++c) {
                double val = main_effects_(variable, c) + (c + 1) * residual_matrix_(i, variable);
                denom += std::exp(val - max_val);
            }
            log_post -= (max_val + std::log(denom));
        }
    } else {
        // Blume-Capel: use parameter (0 = linear, 1 = quadratic)
        log_post += log_beta_prior(main_effects_(variable, parameter));
        log_post += main_effects_(variable, parameter) * blume_capel_stats_(parameter, variable);

        int baseline = baseline_category_(variable);
        for (size_t i = 0; i < n_; ++i) {
            double max_val = -std::numeric_limits<double>::infinity();
            for (int c = 0; c <= num_cats; ++c) {
                int s = c - baseline;
                double val = main_effects_(variable, 0) * s + main_effects_(variable, 1) * s * s + s * residual_matrix_(i, variable);
                if (val > max_val) max_val = val;
            }

            double denom = 0.0;
            for (int c = 0; c <= num_cats; ++c) {
                int s = c - baseline;
                double val = main_effects_(variable, 0) * s + main_effects_(variable, 1) * s * s + s * residual_matrix_(i, variable);
                denom += std::exp(val - max_val);
            }
            log_post -= (max_val + std::log(denom));
        }
    }

    return log_post;
}


double OMRFModel::log_pseudoposterior_pairwise_component(int var1, int var2) const {
    double log_post = 2.0 * pairwise_effects_(var1, var2) * pairwise_stats_(var1, var2);

    // Contribution from both variables' pseudo-likelihoods
    for (int var : {var1, var2}) {
        int num_cats = num_categories_(var);
        int other_var = (var == var1) ? var2 : var1;

        for (size_t i = 0; i < n_; ++i) {
            double max_val = -std::numeric_limits<double>::infinity();

            if (is_ordinal_variable_(var)) {
                max_val = 0.0;
                for (int c = 0; c < num_cats; ++c) {
                    double val = main_effects_(var, c) + (c + 1) * residual_matrix_(i, var);
                    if (val > max_val) max_val = val;
                }

                double denom = std::exp(-max_val);
                for (int c = 0; c < num_cats; ++c) {
                    double val = main_effects_(var, c) + (c + 1) * residual_matrix_(i, var);
                    denom += std::exp(val - max_val);
                }
                log_post -= (max_val + std::log(denom));
            } else {
                int baseline = baseline_category_(var);
                for (int c = 0; c <= num_cats; ++c) {
                    int s = c - baseline;
                    double val = main_effects_(var, 0) * s + main_effects_(var, 1) * s * s + s * residual_matrix_(i, var);
                    if (val > max_val) max_val = val;
                }

                double denom = 0.0;
                for (int c = 0; c <= num_cats; ++c) {
                    int s = c - baseline;
                    double val = main_effects_(var, 0) * s + main_effects_(var, 1) * s * s + s * residual_matrix_(i, var);
                    denom += std::exp(val - max_val);
                }
                log_post -= (max_val + std::log(denom));
            }
        }
    }

    // Cauchy prior if edge is included
    if (edge_indicators_(var1, var2) == 1) {
        log_post += R::dcauchy(pairwise_effects_(var1, var2), 0.0, pairwise_scale_, true);
    }

    return log_post;
}


double OMRFModel::compute_log_likelihood_ratio_for_variable(
    int variable,
    const arma::vec& interacting_score,
    double proposed_state,
    double current_state
) const {
    double log_ratio = 0.0;
    int num_cats = num_categories_(variable);

    for (size_t i = 0; i < n_; ++i) {
        double rest_minus = residual_matrix_(i, variable) - current_state * interacting_score(i);
        double rest_prop = rest_minus + proposed_state * interacting_score(i);
        double rest_curr = rest_minus + current_state * interacting_score(i);

        double max_prop = -std::numeric_limits<double>::infinity();
        double max_curr = -std::numeric_limits<double>::infinity();

        if (is_ordinal_variable_(variable)) {
            max_prop = 0.0;
            max_curr = 0.0;
            for (int c = 0; c < num_cats; ++c) {
                double val_prop = main_effects_(variable, c) + (c + 1) * rest_prop;
                double val_curr = main_effects_(variable, c) + (c + 1) * rest_curr;
                if (val_prop > max_prop) max_prop = val_prop;
                if (val_curr > max_curr) max_curr = val_curr;
            }

            double denom_prop = std::exp(-max_prop);
            double denom_curr = std::exp(-max_curr);
            for (int c = 0; c < num_cats; ++c) {
                double val_prop = main_effects_(variable, c) + (c + 1) * rest_prop;
                double val_curr = main_effects_(variable, c) + (c + 1) * rest_curr;
                denom_prop += std::exp(val_prop - max_prop);
                denom_curr += std::exp(val_curr - max_curr);
            }

            log_ratio += (max_curr + std::log(denom_curr)) - (max_prop + std::log(denom_prop));
        } else {
            int baseline = baseline_category_(variable);
            for (int c = 0; c <= num_cats; ++c) {
                int s = c - baseline;
                double val_prop = main_effects_(variable, 0) * s + main_effects_(variable, 1) * s * s + s * rest_prop;
                double val_curr = main_effects_(variable, 0) * s + main_effects_(variable, 1) * s * s + s * rest_curr;
                if (val_prop > max_prop) max_prop = val_prop;
                if (val_curr > max_curr) max_curr = val_curr;
            }

            double denom_prop = 0.0;
            double denom_curr = 0.0;
            for (int c = 0; c <= num_cats; ++c) {
                int s = c - baseline;
                double val_prop = main_effects_(variable, 0) * s + main_effects_(variable, 1) * s * s + s * rest_prop;
                double val_curr = main_effects_(variable, 0) * s + main_effects_(variable, 1) * s * s + s * rest_curr;
                denom_prop += std::exp(val_prop - max_prop);
                denom_curr += std::exp(val_curr - max_curr);
            }

            log_ratio += (max_curr + std::log(denom_curr)) - (max_prop + std::log(denom_prop));
        }
    }

    return log_ratio;
}


double OMRFModel::log_pseudolikelihood_ratio_interaction(
    int variable1,
    int variable2,
    double proposed_state,
    double current_state
) const {
    double delta = proposed_state - current_state;
    double log_ratio = 2.0 * delta * pairwise_stats_(variable1, variable2);

    // For Blume-Capel variables, transform observations by subtracting baseline
    auto get_transformed_score = [this](int var) -> arma::vec {
        arma::vec score = arma::conv_to<arma::vec>::from(observations_.col(var));
        if (!is_ordinal_variable_(var)) {
            score -= static_cast<double>(baseline_category_(var));
        }
        return score;
    };

    // Contribution from variable1
    arma::vec interacting_score = get_transformed_score(variable2);
    log_ratio += compute_log_likelihood_ratio_for_variable(variable1, interacting_score, proposed_state, current_state);

    // Contribution from variable2
    interacting_score = get_transformed_score(variable1);
    log_ratio += compute_log_likelihood_ratio_for_variable(variable2, interacting_score, proposed_state, current_state);

    return log_ratio;
}


// =============================================================================
// Gradient computation
// =============================================================================

void OMRFModel::ensure_gradient_cache() {
    if (gradient_cache_valid_) return;

    // Compute observed gradient and index matrix (constant during MCMC)
    int num_active = 0;
    index_matrix_cache_ = arma::zeros<arma::imat>(p_, p_);

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                index_matrix_cache_(v1, v2) = num_main_ + num_active;
                num_active++;
            }
        }
    }

    grad_obs_cache_ = arma::zeros<arma::vec>(num_main_ + num_active);

    // Observed statistics for main effects
    int offset = 0;
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                grad_obs_cache_(offset++) = counts_per_category_(c + 1, v);
            }
        } else {
            grad_obs_cache_(offset++) = blume_capel_stats_(0, v);
            grad_obs_cache_(offset++) = blume_capel_stats_(1, v);
        }
    }

    // Observed statistics for pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                grad_obs_cache_(offset++) = 2.0 * pairwise_stats_(v1, v2);
            }
        }
    }

    gradient_cache_valid_ = true;
}


arma::vec OMRFModel::gradient(const arma::vec& parameters) {
    // Unvectorize into temporary matrices (safe approach)
    arma::mat temp_main = main_effects_;
    arma::mat temp_pairwise = pairwise_effects_;

    // Unvectorize parameters into temporaries
    int offset = 0;
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                temp_main(v, c) = parameters(offset++);
            }
        } else {
            temp_main(v, 0) = parameters(offset++);
            temp_main(v, 1) = parameters(offset++);
        }
    }

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                temp_pairwise(v1, v2) = parameters(offset++);
                temp_pairwise(v2, v1) = temp_pairwise(v1, v2);
            }
        }
    }

    // Compute residual matrix from temp_pairwise
    arma::mat temp_residual = arma::conv_to<arma::mat>::from(observations_) * temp_pairwise;

    return gradient_with_state(temp_main, temp_pairwise, temp_residual);
}


arma::vec OMRFModel::gradient_with_state(
    const arma::mat& main_eff,
    const arma::mat& pairwise_eff,
    const arma::mat& residual_mat
) const {
    // Start with cached observed gradient
    arma::vec gradient = grad_obs_cache_;

    // Expected statistics for main and pairwise effects
    int offset = 0;
    for (size_t v = 0; v < p_; ++v) {
        int num_cats = num_categories_(v);
        arma::vec residual_score = residual_mat.col(v);
        arma::vec bound = num_cats * residual_score;

        if (is_ordinal_variable_(v)) {
            // Extract main effect parameters for this variable
            arma::vec main_param = main_eff.row(v).cols(0, num_cats - 1).t();

            // Use optimized helper function
            arma::mat probs = compute_probs_ordinal(main_param, residual_score, bound, num_cats);

            // Main effects gradient
            for (int c = 0; c < num_cats; ++c) {
                gradient(offset + c) -= arma::accu(probs.col(c + 1));
            }

            // Pairwise effects gradient
            for (size_t j = 0; j < p_; ++j) {
                if (edge_indicators_(v, j) == 0 || v == j) continue;

                arma::vec expected_value = arma::zeros<arma::vec>(n_);
                for (int c = 1; c <= num_cats; ++c) {
                    expected_value += c * probs.col(c) % observations_double_.col(j);
                }

                int location = (v < j) ? index_matrix_cache_(v, j) : index_matrix_cache_(j, v);
                gradient(location) -= arma::accu(expected_value);
            }

            offset += num_cats;
        } else {
            int ref = baseline_category_(v);
            double lin_eff = main_eff(v, 0);
            double quad_eff = main_eff(v, 1);

            // Use optimized helper function (updates bound in-place)
            arma::mat probs = compute_probs_blume_capel(residual_score, lin_eff, quad_eff, ref, num_cats, bound);

            arma::vec score = arma::regspace<arma::vec>(0, num_cats) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            // Main effects gradient
            gradient(offset) -= arma::accu(probs * score);
            gradient(offset + 1) -= arma::accu(probs * sq_score);

            // Pairwise effects gradient
            for (size_t j = 0; j < p_; ++j) {
                if (edge_indicators_(v, j) == 0 || v == j) continue;

                arma::vec expected_value = arma::zeros<arma::vec>(n_);
                for (int c = 0; c <= num_cats; ++c) {
                    int s = c - ref;
                    expected_value += s * probs.col(c) % observations_double_.col(j);
                }

                int location = (v < j) ? index_matrix_cache_(v, j) : index_matrix_cache_(j, v);
                gradient(location) -= arma::accu(expected_value);
            }

            offset += 2;
        }
    }

    // Prior gradients for main effects (Beta-prime prior)
    offset = 0;
    for (size_t v = 0; v < p_; ++v) {
        int num_pars = is_ordinal_variable_(v) ? num_categories_(v) : 2;
        for (int c = 0; c < num_pars; ++c) {
            double x = main_eff(v, c);
            double prob = 1.0 / (1.0 + std::exp(-x));
            gradient(offset + c) += main_alpha_ - (main_alpha_ + main_beta_) * prob;
        }
        offset += num_pars;
    }

    // Cauchy prior gradient for pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                int idx = index_matrix_cache_(v1, v2);
                double x = pairwise_eff(v1, v2);
                gradient(idx) -= 2.0 * x / (pairwise_scale_ * pairwise_scale_ + x * x);
            }
        }
    }

    return gradient;
}


arma::vec OMRFModel::gradient_internal() const {
    return gradient_with_state(main_effects_, pairwise_effects_, residual_matrix_);
}


std::pair<double, arma::vec> OMRFModel::logp_and_gradient(const arma::vec& parameters) {
    // Ensure gradient cache is initialized
    ensure_gradient_cache();

    // Use the external-state versions
    arma::mat temp_main = main_effects_;
    arma::mat temp_pairwise = pairwise_effects_;

    // Unvectorize parameters into temporaries
    int offset = 0;
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                temp_main(v, c) = parameters(offset++);
            }
        } else {
            temp_main(v, 0) = parameters(offset++);
            temp_main(v, 1) = parameters(offset++);
        }
    }

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                temp_pairwise(v1, v2) = parameters(offset++);
                temp_pairwise(v2, v1) = temp_pairwise(v1, v2);
            }
        }
    }

    arma::mat temp_residual = observations_double_ * temp_pairwise;

    double lp = log_pseudoposterior_with_state(temp_main, temp_pairwise, temp_residual);
    arma::vec grad = gradient_with_state(temp_main, temp_pairwise, temp_residual);

    return {lp, grad};
}


// =============================================================================
// Metropolis-Hastings updates
// =============================================================================

void OMRFModel::update_main_effect_parameter(int variable, int category, int parameter) {
    double proposal_sd;
    double& current = is_ordinal_variable_(variable)
        ? main_effects_(variable, category)
        : main_effects_(variable, parameter);

    proposal_sd = is_ordinal_variable_(variable)
        ? proposal_sd_main_(variable, category)
        : proposal_sd_main_(variable, parameter);

    int cat_for_log = is_ordinal_variable_(variable) ? category : -1;
    int par_for_log = is_ordinal_variable_(variable) ? -1 : parameter;

    auto log_post = [&](double theta) {
        double old_val = current;
        current = theta;
        double lp = log_pseudoposterior_main_component(variable, cat_for_log, par_for_log);
        current = old_val;
        return lp;
    };

    SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng_);
    current = result.state[0];
    invalidate_gradient_cache();
}


void OMRFModel::update_pairwise_effect(int var1, int var2) {
    if (edge_indicators_(var1, var2) == 0) return;

    double& value = pairwise_effects_(var1, var2);
    double proposal_sd = proposal_sd_pairwise_(var1, var2);
    double current = value;

    auto log_post = [&](double theta) {
        pairwise_effects_(var1, var2) = theta;
        pairwise_effects_(var2, var1) = theta;
        update_residual_matrix();
        return log_pseudoposterior_pairwise_component(var1, var2);
    };

    SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng_);

    value = result.state[0];
    pairwise_effects_(var2, var1) = value;

    if (current != value) {
        update_residual_matrix();
    }
    invalidate_gradient_cache();
}


void OMRFModel::update_edge_indicator(int var1, int var2) {
    double current_state = pairwise_effects_(var1, var2);
    double proposal_sd = proposal_sd_pairwise_(var1, var2);

    bool proposing_addition = (edge_indicators_(var1, var2) == 0);
    double proposed_state = proposing_addition ? rnorm(rng_, current_state, proposal_sd) : 0.0;

    double log_accept = log_pseudolikelihood_ratio_interaction(var1, var2, proposed_state, current_state);

    double incl_prob = inclusion_probability_(var1, var2);

    if (proposing_addition) {
        log_accept += R::dcauchy(proposed_state, 0.0, pairwise_scale_, true);
        log_accept -= R::dnorm(proposed_state, current_state, proposal_sd, true);
        log_accept += MY_LOG(incl_prob) - MY_LOG(1.0 - incl_prob);
    } else {
        log_accept -= R::dcauchy(current_state, 0.0, pairwise_scale_, true);
        log_accept += R::dnorm(current_state, proposed_state, proposal_sd, true);
        log_accept -= MY_LOG(incl_prob) - MY_LOG(1.0 - incl_prob);
    }

    if (MY_LOG(runif(rng_)) < log_accept) {
        int updated = 1 - edge_indicators_(var1, var2);
        edge_indicators_(var1, var2) = updated;
        edge_indicators_(var2, var1) = updated;

        pairwise_effects_(var1, var2) = proposed_state;
        pairwise_effects_(var2, var1) = proposed_state;

        update_residual_matrix();
        invalidate_gradient_cache();
    }
}


// =============================================================================
// Main update methods
// =============================================================================

void OMRFModel::do_one_mh_step() {
    // Update main effects
    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                update_main_effect_parameter(v, c, -1);
            }
        } else {
            for (int p = 0; p < 2; ++p) {
                update_main_effect_parameter(v, -1, p);
            }
        }
    }

    // Update pairwise effects
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            update_pairwise_effect(v1, v2);
        }
    }

    // Update edge indicators if in selection phase
    if (edge_selection_active_) {
        update_edge_indicators();
    }

    proposal_.increment_iteration();
}


void OMRFModel::update_edge_indicators() {
    for (size_t idx = 0; idx < num_pairwise_; ++idx) {
        int var1 = interaction_index_(idx, 1);
        int var2 = interaction_index_(idx, 2);
        update_edge_indicator(var1, var2);
    }
}


void OMRFModel::initialize_graph() {
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            double p = inclusion_probability_(v1, v2);
            int draw = (runif(rng_) < p) ? 1 : 0;
            edge_indicators_(v1, v2) = draw;
            edge_indicators_(v2, v1) = draw;
            if (!draw) {
                pairwise_effects_(v1, v2) = 0.0;
                pairwise_effects_(v2, v1) = 0.0;
            }
        }
    }
    update_residual_matrix();
    invalidate_gradient_cache();
}



void OMRFModel::impute_missing() {
    if (!has_missing_) return;

    // For each missing value, sample from conditional distribution
    for (size_t m = 0; m < missing_index_.n_rows; ++m) {
        int person = missing_index_(m, 0);
        int variable = missing_index_(m, 1);
        int num_cats = num_categories_(variable);

        arma::vec log_probs;
        if (is_ordinal_variable_(variable)) {
            log_probs.set_size(num_cats + 1);
            log_probs(0) = 0.0;
            for (int c = 0; c < num_cats; ++c) {
                log_probs(c + 1) = main_effects_(variable, c) + (c + 1) * residual_matrix_(person, variable);
            }
        } else {
            int baseline = baseline_category_(variable);
            log_probs.set_size(num_cats + 1);
            for (int c = 0; c <= num_cats; ++c) {
                int s = c - baseline;
                log_probs(c) = main_effects_(variable, 0) * s + main_effects_(variable, 1) * s * s + s * residual_matrix_(person, variable);
            }
        }

        // Sample from categorical
        double max_val = log_probs.max();
        arma::vec probs = arma::exp(log_probs - max_val);
        probs /= arma::sum(probs);

        double u = runif(rng_);
        double cumsum = 0.0;
        int new_value = 0;
        for (size_t c = 0; c < probs.n_elem; ++c) {
            cumsum += probs(c);
            if (u < cumsum) {
                new_value = c;
                break;
            }
        }

        int old_value = observations_(person, variable);
        if (new_value != old_value) {
            observations_(person, variable) = new_value;
            // Update sufficient statistics
            compute_sufficient_statistics();
            update_residual_matrix();
        }
    }
}


// =============================================================================
// Factory function
// =============================================================================

OMRFModel createOMRFModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection
) {
    arma::imat observations = Rcpp::as<arma::imat>(inputFromR["observations"]);
    arma::ivec num_categories = Rcpp::as<arma::ivec>(inputFromR["num_categories"]);
    arma::uvec is_ordinal_variable = Rcpp::as<arma::uvec>(inputFromR["is_ordinal_variable"]);
    arma::ivec baseline_category = Rcpp::as<arma::ivec>(inputFromR["baseline_category"]);

    double main_alpha = inputFromR.containsElementNamed("main_alpha")
        ? Rcpp::as<double>(inputFromR["main_alpha"]) : 1.0;
    double main_beta = inputFromR.containsElementNamed("main_beta")
        ? Rcpp::as<double>(inputFromR["main_beta"]) : 1.0;
    double pairwise_scale = inputFromR.containsElementNamed("pairwise_scale")
        ? Rcpp::as<double>(inputFromR["pairwise_scale"]) : 2.5;

    return OMRFModel(
        observations,
        num_categories,
        inclusion_probability,
        initial_edge_indicators,
        is_ordinal_variable,
        baseline_category,
        main_alpha,
        main_beta,
        pairwise_scale,
        edge_selection
    );
}


// =============================================================================
// R interface: sample_omrf_classed
// =============================================================================

// [[Rcpp::export]]
Rcpp::List sample_omrf_classed(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const int no_iter,
    const int no_warmup,
    const bool edge_selection,
    const std::string& sampler_type,
    const int seed
) {
    // Create model from R input
    OMRFModel model = createOMRFModelFromR(
        inputFromR,
        prior_inclusion_prob,
        initial_edge_indicators,
        edge_selection
    );

    // Set random seed
    model.set_seed(seed);

    // Storage for samples - use FIXED size (all parameters)
    int full_dim = model.full_parameter_dimension();
    arma::mat samples(no_iter, full_dim);
    arma::imat indicator_samples;
    if (edge_selection) {
        int num_edges = (model.get_p() * (model.get_p() - 1)) / 2;
        indicator_samples.set_size(no_iter, num_edges);
    }

    // NUTS/HMC diagnostics
    arma::ivec treedepth_samples;
    arma::ivec divergent_samples;
    arma::vec energy_samples;
    bool use_nuts = (sampler_type == "nuts");
    bool use_hmc = (sampler_type == "hmc");
    if (use_nuts || use_hmc) {
        treedepth_samples.set_size(no_iter);
        divergent_samples.set_size(no_iter);
        energy_samples.set_size(no_iter);
    }

    // Create sampler configuration
    SamplerConfig config;
    config.sampler_type = sampler_type;
    config.initial_step_size = 0.1;
    config.target_acceptance = 0.8;
    config.max_tree_depth = 10;
    config.num_leapfrogs = 10;
    config.no_warmup = no_warmup;

    // Create appropriate sampler
    std::unique_ptr<BaseSampler> sampler = create_sampler(config);

    // Warmup phase
    Rcpp::Rcout << "Running warmup (" << no_warmup << " iterations)..." << std::endl;
    for (int iter = 0; iter < no_warmup; ++iter) {
        SamplerResult result = sampler->warmup_step(model);

        // Edge selection only after initial warmup period
        if (edge_selection && iter > no_warmup / 2) {
            model.update_edge_indicators();
        }

        // Check for user interrupt
        if ((iter + 1) % 100 == 0) {
            Rcpp::checkUserInterrupt();
            Rcpp::Rcout << "  Warmup iteration " << (iter + 1) << "/" << no_warmup
                        << " (step_size=" << model.get_step_size() << ")" << std::endl;
        }
    }

    // Use averaged step size for sampling
    sampler->finalize_warmup();
    Rcpp::Rcout << "Warmup complete." << std::endl;

    // Sampling phase
    Rcpp::Rcout << "Running sampling (" << no_iter << " iterations)..." << std::endl;
    for (int iter = 0; iter < no_iter; ++iter) {
        SamplerResult result = sampler->sample_step(model);

        // Extract NUTS/HMC diagnostics if available
        if (sampler->has_nuts_diagnostics()) {
            if (auto nuts_diag = std::dynamic_pointer_cast<NUTSDiagnostics>(result.diagnostics)) {
                treedepth_samples(iter) = nuts_diag->tree_depth;
                divergent_samples(iter) = nuts_diag->divergent ? 1 : 0;
                energy_samples(iter) = nuts_diag->energy;
            }
        }

        if (edge_selection) {
            model.update_edge_indicators();
        }

        // Store samples - use FULL vectorization (fixed size)
        samples.row(iter) = model.get_full_vectorized_parameters().t();

        if (edge_selection) {
            arma::imat indicators = model.get_edge_indicators();
            int idx = 0;
            for (int i = 0; i < static_cast<int>(model.get_p()) - 1; ++i) {
                for (int j = i + 1; j < static_cast<int>(model.get_p()); ++j) {
                    indicator_samples(iter, idx++) = indicators(i, j);
                }
            }
        }

        // Check for user interrupt
        if ((iter + 1) % 100 == 0) {
            Rcpp::checkUserInterrupt();
            Rcpp::Rcout << "  Sampling iteration " << (iter + 1) << "/" << no_iter << std::endl;
        }
    }

    // Build output list
    Rcpp::List output;
    output["samples"] = samples;

    if (edge_selection) {
        output["indicator_samples"] = indicator_samples;
        // Compute posterior mean of edge indicators
        arma::vec posterior_mean_indicator = arma::mean(arma::conv_to<arma::mat>::from(indicator_samples), 0).t();
        output["posterior_mean_indicator"] = posterior_mean_indicator;
    }

    if (use_nuts || use_hmc) {
        output["treedepth"] = treedepth_samples;
        output["divergent"] = divergent_samples;
        output["energy"] = energy_samples;
        // Get final step size from sampler (NUTSSampler and HMCSampler have get_step_size())
        output["final_step_size"] = 0.0;  // Could add getter to sampler if needed
    }

    output["sampler_type"] = sampler_type;
    output["no_iter"] = no_iter;
    output["no_warmup"] = no_warmup;
    output["edge_selection"] = edge_selection;
    output["num_variables"] = model.get_p();
    output["num_observations"] = model.get_n();

    return output;
}