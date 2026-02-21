#include <RcppArmadillo.h>
#include "models/omrf/omrf_model.h"
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
    proposal_sd_main_ = arma::ones<arma::mat>(p_, max_cats);
    proposal_sd_pairwise_ = arma::ones<arma::mat>(p_, p_);

    // Initialize per-pair scaling factors (default: ones)
    pairwise_scaling_factors_ = arma::ones<arma::mat>(p_, p_);

    // Initialize mass matrix
    inv_mass_ = arma::ones<arma::vec>(num_main_ + num_pairwise_);

    // Pre-compute observations as double (for efficient matrix operations)
    observations_double_ = arma::conv_to<arma::mat>::from(observations_);
    observations_double_t_ = observations_double_.t();

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
      observations_double_t_(other.observations_double_t_),
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
      pairwise_scaling_factors_(other.pairwise_scaling_factors_),
      edge_selection_(other.edge_selection_),
      edge_selection_active_(other.edge_selection_active_),
      num_main_(other.num_main_),
      num_pairwise_(other.num_pairwise_),
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
      interaction_index_(other.interaction_index_),
      shuffled_edge_order_(other.shuffled_edge_order_)
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
    residual_matrix_ = observations_double_ * pairwise_effects_;
}


void OMRFModel::update_residual_columns(int var1, int var2, double delta) {
    residual_matrix_.col(var1) += delta * observations_double_.col(var2);
    residual_matrix_.col(var2) += delta * observations_double_.col(var1);
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


void OMRFModel::init_mh_adaptation(const WarmupSchedule& schedule) {
    rwm_main_adapter_ = std::make_unique<RWMAdaptationController>(
        proposal_sd_main_, schedule);
    rwm_pairwise_adapter_ = std::make_unique<RWMAdaptationController>(
        proposal_sd_pairwise_, schedule);
}


void OMRFModel::tune_proposal_sd(int iteration, const WarmupSchedule& schedule) {
    if (!schedule.adapt_proposal_sd(iteration)) return;

    const double target_accept = 0.44;
    const double rm_decay = 0.75;
    double t = iteration - schedule.stage3b_start + 1;
    double rm_weight = std::pow(t, -rm_decay);

    const int num_variables = static_cast<int>(p_);

    for (int variable1 = 0; variable1 < num_variables - 1; variable1++) {
        for (int variable2 = variable1 + 1; variable2 < num_variables; variable2++) {
            double current = pairwise_effects_(variable1, variable2);
            double proposal_sd = proposal_sd_pairwise_(variable1, variable2);

            auto log_post = [&](double theta) {
                double delta = theta - current;
                return log_pseudoposterior_pairwise_at_delta(variable1, variable2, delta);
            };

            SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng_);

            double value = result.state[0];
            pairwise_effects_(variable1, variable2) = value;
            pairwise_effects_(variable2, variable1) = value;

            if (current != value) {
                double delta = value - current;
                residual_matrix_.col(variable1) += arma::conv_to<arma::vec>::from(observations_.col(variable2)) * delta;
                residual_matrix_.col(variable2) += arma::conv_to<arma::vec>::from(observations_.col(variable1)) * delta;
            }

            proposal_sd = update_proposal_sd_with_robbins_monro(
                proposal_sd, MY_LOG(result.accept_prob), rm_weight, target_accept);
            proposal_sd_pairwise_(variable1, variable2) = proposal_sd;
            proposal_sd_pairwise_(variable2, variable1) = proposal_sd;
        }
    }

    invalidate_gradient_cache();
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


void OMRFModel::unvectorize_to_temps(
    const arma::vec& parameters,
    arma::mat& temp_main,
    arma::mat& temp_pairwise,
    arma::mat& temp_residual
) const {
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

    temp_residual = observations_double_ * temp_pairwise;
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
    arma::mat temp_main = main_effects_;
    arma::mat temp_pairwise = pairwise_effects_;
    arma::mat temp_residual;
    unvectorize_to_temps(parameters, temp_main, temp_pairwise, temp_residual);
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

    // Log-normalizer using joint logZ+probs helpers
    for (size_t v = 0; v < p_; ++v) {
        int num_cats = num_categories_(v);
        arma::vec residual_score = residual_mat.col(v);
        arma::vec bound = num_cats * residual_score;

        if (is_ordinal_variable_(v)) {
            arma::vec main_param = main_eff.row(v).cols(0, num_cats - 1).t();
            LogZAndProbs result = compute_logZ_and_probs_ordinal(
                main_param, residual_score, bound, num_cats);
            log_post -= arma::accu(result.log_Z);
        } else {
            int ref = baseline_category_(v);
            double lin = main_eff(v, 0), quad = main_eff(v, 1);
            LogZAndProbs result = compute_logZ_and_probs_blume_capel(
                residual_score, lin, quad, ref, num_cats, bound);
            log_post -= arma::accu(result.log_Z);
        }
    }

    // Pairwise effect contributions: sufficient statistics + Cauchy prior
    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            if (edge_indicators_(v1, v2) == 1) {
                double effect = pairwise_eff(v1, v2);
                log_post += 2.0 * pairwise_stats_(v1, v2) * effect;
                double scaled_scale = pairwise_scale_ * pairwise_scaling_factors_(v1, v2);
                log_post += R::dcauchy(effect, 0.0, scaled_scale, true);
            }
        }
    }

    return log_post;
}


double OMRFModel::log_pseudoposterior_internal() const {
    return log_pseudoposterior_with_state(main_effects_, pairwise_effects_, residual_matrix_);
}


double OMRFModel::log_pseudoposterior_main_component(int variable, int category, int parameter) const {
    double log_posterior = 0.0;

    auto log_beta_prior = [&](double main_effect_param) {
        return main_effect_param * main_alpha_ - std::log1p(MY_EXP(main_effect_param)) * (main_alpha_ + main_beta_);
    };

    const int num_cats = num_categories_(variable);
    arma::vec bound = num_cats * residual_matrix_.col(variable);

    if (is_ordinal_variable_(variable)) {
        const double value = main_effects_(variable, category);
        log_posterior += value * counts_per_category_(category + 1, variable);
        log_posterior += log_beta_prior(value);

        arma::vec residual_score = residual_matrix_.col(variable);
        arma::vec main_effect_param = main_effects_.row(variable).cols(0, num_cats - 1).t();

        arma::vec denom = compute_denom_ordinal(residual_score, main_effect_param, bound);
        log_posterior -= arma::accu(bound + ARMA_MY_LOG(denom));
    } else {
        const double value = main_effects_(variable, parameter);
        const double linear_main_effect = main_effects_(variable, 0);
        const double quadratic_main_effect = main_effects_(variable, 1);
        const int ref = baseline_category_(variable);

        log_posterior += value * blume_capel_stats_(parameter, variable);
        log_posterior += log_beta_prior(value);

        arma::vec residual_score = residual_matrix_.col(variable);
        arma::vec denom(n_, arma::fill::zeros);

        denom = compute_denom_blume_capel(
            residual_score, linear_main_effect, quadratic_main_effect, ref, num_cats, bound
        );

        log_posterior -= arma::accu(bound + ARMA_MY_LOG(denom));
    }

    return log_posterior;
}


double OMRFModel::log_pseudoposterior_pairwise_component(int var1, int var2) const {
    double log_post = 2.0 * pairwise_effects_(var1, var2) * pairwise_stats_(var1, var2);

    for (int var : {var1, var2}) {
        int num_cats = num_categories_(var);
        arma::vec residual_score = residual_matrix_.col(var);

        if (is_ordinal_variable_(var)) {
            arma::vec main_param = main_effects_.row(var).cols(0, num_cats - 1).t();
            arma::vec bound = num_cats * residual_score;
            arma::vec denom = compute_denom_ordinal(residual_score, main_param, bound);
            log_post -= arma::accu(bound + ARMA_MY_LOG(denom));
        } else {
            arma::vec bound(n_);
            arma::vec denom = compute_denom_blume_capel(
                residual_score, main_effects_(var, 0), main_effects_(var, 1),
                baseline_category_(var), num_cats, bound
            );
            log_post -= arma::accu(bound + ARMA_MY_LOG(denom));
        }
    }

    if (edge_indicators_(var1, var2) == 1) {
        double scaled_scale = pairwise_scale_ * pairwise_scaling_factors_(var1, var2);
        log_post += R::dcauchy(pairwise_effects_(var1, var2), 0.0, scaled_scale, true);
    }

    return log_post;
}


double OMRFModel::compute_log_likelihood_ratio_for_variable(
    int variable,
    const arma::ivec& interacting_score,
    double proposed_state,
    double current_state
) const {
    // Convert interaction score vector to double precision
    arma::vec interaction = arma::conv_to<arma::vec>::from(interacting_score);

    const int num_persons = static_cast<int>(n_);
    const int num_cats = num_categories_(variable);

    // Compute adjusted linear predictors without the current interaction
    arma::vec residual_score = residual_matrix_.col(variable) - interaction * current_state;
    arma::vec bounds = residual_score * num_cats;

    arma::vec denom_current = arma::zeros(num_persons);
    arma::vec denom_proposed = arma::zeros(num_persons);

    if (is_ordinal_variable_(variable)) {
        arma::vec main_param = main_effects_.row(variable).cols(0, num_cats - 1).t();

        denom_current += compute_denom_ordinal(
            residual_score + interaction * current_state, main_param, bounds
        );
        denom_proposed += compute_denom_ordinal(
            residual_score + interaction * proposed_state, main_param, bounds
        );
    } else {
        const int ref_cat = baseline_category_(variable);

        denom_current = compute_denom_blume_capel(
            residual_score + interaction * current_state, main_effects_(variable, 0),
            main_effects_(variable, 1), ref_cat, num_cats, bounds
        );
        double log_ratio = arma::accu(ARMA_MY_LOG(denom_current) + bounds);

        denom_proposed = compute_denom_blume_capel(
            residual_score + interaction * proposed_state, main_effects_(variable, 0),
            main_effects_(variable, 1), ref_cat, num_cats, bounds
        );
        log_ratio -= arma::accu(ARMA_MY_LOG(denom_proposed) + bounds);

        return log_ratio;
    }

    // Accumulated log-likelihood difference across persons
    return arma::accu(ARMA_MY_LOG(denom_current) - ARMA_MY_LOG(denom_proposed));
}


double OMRFModel::log_pseudolikelihood_ratio_interaction(
    int variable1,
    int variable2,
    double proposed_state,
    double current_state
) const {
    double log_ratio = 0.0;
    const double delta = proposed_state - current_state;

    arma::ivec score1 = observations_.col(variable1);
    arma::ivec score2 = observations_.col(variable2);

    log_ratio += 2.0 * pairwise_stats_(variable1, variable2) * delta;

    log_ratio += compute_log_likelihood_ratio_for_variable(
        variable1, score2, proposed_state, current_state
    );

    log_ratio += compute_log_likelihood_ratio_for_variable(
        variable2, score1, proposed_state, current_state
    );

    return log_ratio;
}


double OMRFModel::log_pseudoposterior_pairwise_at_delta(int var1, int var2, double delta) const {
    const int num_observations = static_cast<int>(n_);
    const double proposed_value = pairwise_effects_(var1, var2) + delta;

    double log_pseudo_posterior = 2.0 * proposed_value * pairwise_stats_(var1, var2);

    arma::vec obs_var1 = arma::conv_to<arma::vec>::from(observations_.col(var1));
    arma::vec obs_var2 = arma::conv_to<arma::vec>::from(observations_.col(var2));

    for (int var : {var1, var2}) {
        int num_cats = num_categories_(var);
        const arma::vec& obs_other = (var == var1) ? obs_var2 : obs_var1;

        arma::vec residual_score = residual_matrix_.col(var) + obs_other * delta;
        arma::vec denominator = arma::zeros(num_observations);
        arma::vec bound = num_cats * residual_score;

        if (is_ordinal_variable_(var)) {
            arma::vec main_effect_param = main_effects_.row(var).cols(0, num_cats - 1).t();
            denominator += compute_denom_ordinal(residual_score, main_effect_param, bound);
        } else {
            const int ref = baseline_category_(var);
            denominator = compute_denom_blume_capel(
                residual_score, main_effects_(var, 0), main_effects_(var, 1), ref, num_cats, bound
            );
        }

        log_pseudo_posterior -= arma::accu(ARMA_MY_LOG(denominator));
        log_pseudo_posterior -= arma::accu(bound);
    }

    if (edge_indicators_(var1, var2) == 1) {
        const double scaled_pairwise_scale = pairwise_scale_ * pairwise_scaling_factors_(var1, var2);
        log_pseudo_posterior += R::dcauchy(proposed_value, 0.0, scaled_pairwise_scale, true);
    }

    return log_pseudo_posterior;
}


// =============================================================================
// Gradient computation
// =============================================================================

void OMRFModel::ensure_gradient_cache() {
    if (gradient_cache_valid_) return;

    const int num_variables = static_cast<int>(p_);
    const int num_main = static_cast<int>(num_main_);
    index_matrix_cache_.set_size(num_variables, num_variables);
    index_matrix_cache_.zeros();

    // Count active pairwise effects + build index map
    int num_active = 0;
    for (int i = 0; i < num_variables - 1; i++) {
        for (int j = i + 1; j < num_variables; j++) {
            if (edge_indicators_(i, j) == 1) {
                index_matrix_cache_(i, j) = num_main + num_active++;
                index_matrix_cache_(j, i) = index_matrix_cache_(i, j);
            }
        }
    }

    // Allocate gradient vector (main + active pairwise only)
    grad_obs_cache_.set_size(num_main + num_active);
    grad_obs_cache_.zeros();

    // Observed statistics for main effects
    int offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        if (is_ordinal_variable_(variable)) {
            const int num_cats = num_categories_(variable);
            for (int cat = 0; cat < num_cats; cat++) {
                grad_obs_cache_(offset + cat) = counts_per_category_(cat + 1, variable);
            }
            offset += num_cats;
        } else {
            grad_obs_cache_(offset) = blume_capel_stats_(0, variable);
            grad_obs_cache_(offset + 1) = blume_capel_stats_(1, variable);
            offset += 2;
        }
    }

    // Observed statistics for pairwise effects
    for (int i = 0; i < num_variables - 1; i++) {
        for (int j = i + 1; j < num_variables; j++) {
            if (edge_indicators_(i, j) == 0) continue;
            int location = index_matrix_cache_(i, j);
            grad_obs_cache_(location) = 2.0 * pairwise_stats_(i, j);
        }
    }

    gradient_cache_valid_ = true;
}

arma::vec OMRFModel::gradient(const arma::vec& parameters) {
    ensure_gradient_cache();

    arma::mat temp_main = main_effects_;
    arma::mat temp_pairwise = pairwise_effects_;
    arma::mat temp_residual;
    unvectorize_to_temps(parameters, temp_main, temp_pairwise, temp_residual);

    const int num_variables = static_cast<int>(p_);

    // Start with cached observed gradient
    arma::vec gradient = grad_obs_cache_;

    // ---- Expected statistics ----
    int offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        const int num_cats = num_categories_(variable);
        arma::vec residual_score = temp_residual.col(variable);
        arma::vec bound = num_cats * residual_score;

        if (is_ordinal_variable_(variable)) {
            arma::vec main_param = temp_main.row(variable).cols(0, num_cats - 1).t();
            arma::mat probs = compute_probs_ordinal(
                main_param, residual_score, bound, num_cats
            );

            // Main effects gradient
            for (int cat = 0; cat < num_cats; cat++) {
                gradient(offset + cat) -= arma::accu(probs.col(cat + 1));
            }

            // Pairwise gradient (vectorized using BLAS)
            arma::vec weights = arma::regspace<arma::vec>(1, num_cats);
            arma::vec E = probs.cols(1, num_cats) * weights;
            arma::vec pw_grad = observations_double_t_ * E;
            for (int j = 0; j < num_variables; j++) {
                if (edge_indicators_(variable, j) == 0 || variable == j) continue;
                int location = (variable < j) ? index_matrix_cache_(variable, j) : index_matrix_cache_(j, variable);
                gradient(location) -= pw_grad(j);
            }
            offset += num_cats;
        } else {
            const int ref = baseline_category_(variable);
            const double lin_eff = temp_main(variable, 0);
            const double quad_eff = temp_main(variable, 1);

            arma::mat probs = compute_probs_blume_capel(
                residual_score, lin_eff, quad_eff, ref, num_cats, bound
            );

            arma::vec score = arma::regspace<arma::vec>(0, num_cats) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            // Main effects gradient
            gradient(offset) -= arma::accu(probs * score);
            gradient(offset + 1) -= arma::accu(probs * sq_score);

            // Pairwise gradient (vectorized using BLAS)
            arma::vec E = probs * score;
            arma::vec pw_grad = observations_double_t_ * E;
            for (int j = 0; j < num_variables; j++) {
                if (edge_indicators_(variable, j) == 0 || variable == j) continue;
                int location = (variable < j)
                    ? index_matrix_cache_(variable, j)
                    : index_matrix_cache_(j, variable);
                gradient(location) -= pw_grad(j);
            }
            offset += 2;
        }
    }

    // ---- Priors ----
    offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        if (is_ordinal_variable_(variable)) {
            const int num_cats = num_categories_(variable);
            for (int cat = 0; cat < num_cats; cat++) {
                const double p = 1.0 / (1.0 + MY_EXP(-temp_main(variable, cat)));
                gradient(offset + cat) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            offset += num_cats;
        } else {
            for (int k = 0; k < 2; k++) {
                const double param = temp_main(variable, k);
                const double p = 1.0 / (1.0 + MY_EXP(-param));
                gradient(offset + k) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            offset += 2;
        }
    }
    for (int i = 0; i < num_variables - 1; i++) {
        for (int j = i + 1; j < num_variables; j++) {
            if (edge_indicators_(i, j) == 0) continue;
            int location = index_matrix_cache_(i, j);
            const double effect = temp_pairwise(i, j);
            const double scaled_scale = pairwise_scale_ * pairwise_scaling_factors_(i, j);
            gradient(location) -= 2.0 * effect / (effect * effect + scaled_scale * scaled_scale);
        }
    }

    return gradient;
}


std::pair<double, arma::vec> OMRFModel::logp_and_gradient(const arma::vec& parameters) {
    ensure_gradient_cache();

    arma::mat temp_main(main_effects_.n_rows, main_effects_.n_cols, arma::fill::none);
    arma::mat temp_pairwise(p_, p_, arma::fill::zeros);
    arma::mat temp_residual;
    unvectorize_to_temps(parameters, temp_main, temp_pairwise, temp_residual);

    const int num_variables = static_cast<int>(p_);

    double log_pp = 0.0;
    arma::vec gradient = grad_obs_cache_;

    auto log_beta_prior = [&](double main_effect_param) {
        return main_effect_param * main_alpha_ -
               std::log1p(MY_EXP(main_effect_param)) * (main_alpha_ + main_beta_);
    };

    // ---- Main effects: priors + sufficient statistics ----
    for (int variable = 0; variable < num_variables; variable++) {
        if (is_ordinal_variable_(variable)) {
            const int num_cats = num_categories_(variable);
            for (int cat = 0; cat < num_cats; cat++) {
                double value = temp_main(variable, cat);
                log_pp += counts_per_category_(cat + 1, variable) * value;
                log_pp += log_beta_prior(value);
            }
        } else {
            double value = temp_main(variable, 0);
            log_pp += log_beta_prior(value);
            log_pp += blume_capel_stats_(0, variable) * value;

            value = temp_main(variable, 1);
            log_pp += log_beta_prior(value);
            log_pp += blume_capel_stats_(1, variable) * value;
        }
    }

    // ---- Pairwise effects: priors + sufficient statistics ----
    for (int var1 = 0; var1 < num_variables - 1; var1++) {
        for (int var2 = var1 + 1; var2 < num_variables; var2++) {
            if (edge_indicators_(var1, var2) == 0) continue;

            double value = temp_pairwise(var1, var2);
            log_pp += 2.0 * pairwise_stats_(var1, var2) * value;
            const double scaled_pairwise_scale = pairwise_scale_ * pairwise_scaling_factors_(var1, var2);
            log_pp += R::dcauchy(value, 0.0, scaled_pairwise_scale, true);
        }
    }

    // ---- Per-variable: joint computation of log-normalizer and gradient ----
    int offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        const int num_cats = num_categories_(variable);
        arma::vec residual_score = temp_residual.col(variable);
        arma::vec bound = num_cats * residual_score;

        if (is_ordinal_variable_(variable)) {
            arma::vec main_param = temp_main.row(variable).cols(0, num_cats - 1).t();

            // Joint computation: get both log_Z and probs in one pass
            LogZAndProbs result = compute_logZ_and_probs_ordinal(
                main_param, residual_score, bound, num_cats
            );

            // Use log_Z for log-pseudoposterior
            log_pp -= arma::accu(result.log_Z);

            // Use probs for gradient
            for (int cat = 0; cat < num_cats; cat++) {
                gradient(offset + cat) -= arma::accu(result.probs.col(cat + 1));
            }

            // Pairwise gradient contributions (vectorized using BLAS)
            arma::vec weights = arma::regspace<arma::vec>(1, num_cats);
            arma::vec E = result.probs.cols(1, num_cats) * weights;
            arma::vec pw_grad = observations_double_t_ * E;
            for (int j = 0; j < num_variables; j++) {
                if (edge_indicators_(variable, j) == 0 || variable == j) continue;
                int location = (variable < j) ? index_matrix_cache_(variable, j) : index_matrix_cache_(j, variable);
                gradient(location) -= pw_grad(j);
            }
            offset += num_cats;
        } else {
            const int ref = baseline_category_(variable);
            const double lin_eff = temp_main(variable, 0);
            const double quad_eff = temp_main(variable, 1);

            // Joint computation: get both log_Z and probs in one pass
            LogZAndProbs result = compute_logZ_and_probs_blume_capel(
                residual_score, lin_eff, quad_eff, ref, num_cats, bound
            );

            // Use log_Z for log-pseudoposterior
            log_pp -= arma::accu(result.log_Z);

            // Use probs for gradient
            arma::vec score = arma::regspace<arma::vec>(0, num_cats) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            gradient(offset)     -= arma::accu(result.probs * score);
            gradient(offset + 1) -= arma::accu(result.probs * sq_score);

            // Pairwise gradient contributions (vectorized using BLAS)
            arma::vec E = result.probs * score;
            arma::vec pw_grad = observations_double_t_ * E;
            for (int j = 0; j < num_variables; j++) {
                if (edge_indicators_(variable, j) == 0 || variable == j) continue;
                int location = (variable < j) ? index_matrix_cache_(variable, j) : index_matrix_cache_(j, variable);
                gradient(location) -= pw_grad(j);
            }
            offset += 2;
        }
    }

    // ---- Priors: gradient contributions ----
    offset = 0;
    for (int variable = 0; variable < num_variables; variable++) {
        if (is_ordinal_variable_(variable)) {
            const int num_cats = num_categories_(variable);
            for (int cat = 0; cat < num_cats; cat++) {
                const double p = 1.0 / (1.0 + MY_EXP(-temp_main(variable, cat)));
                gradient(offset + cat) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            offset += num_cats;
        } else {
            for (int k = 0; k < 2; k++) {
                const double param = temp_main(variable, k);
                const double p = 1.0 / (1.0 + MY_EXP(-param));
                gradient(offset + k) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            offset += 2;
        }
    }
    for (int i = 0; i < num_variables - 1; i++) {
        for (int j = i + 1; j < num_variables; j++) {
            if (edge_indicators_(i, j) == 0) continue;
            int location = index_matrix_cache_(i, j);
            const double effect = temp_pairwise(i, j);
            const double scaled_scale = pairwise_scale_ * pairwise_scaling_factors_(i, j);
            gradient(location) -= 2.0 * effect / (effect * effect + scaled_scale * scaled_scale);
        }
    }

    return {log_pp, gradient};
}


// =============================================================================
// Metropolis-Hastings updates
// =============================================================================

double OMRFModel::update_main_effect_parameter(int variable, int category, int parameter) {
    double& current = is_ordinal_variable_(variable)
        ? main_effects_(variable, category)
        : main_effects_(variable, parameter);

    double proposal_sd = is_ordinal_variable_(variable)
        ? proposal_sd_main_(variable, category)
        : proposal_sd_main_(variable, parameter);

    auto log_post = [&](double theta) {
        current = theta;
        return log_pseudoposterior_main_component(variable, category, parameter);
    };

    SamplerResult result = rwm_sampler(current, proposal_sd, log_post, rng_);
    current = result.state[0];
    return result.accept_prob;
}


double OMRFModel::update_pairwise_effect(int var1, int var2) {
    if (edge_indicators_(var1, var2) == 0) return 1.0;

    double current_value = pairwise_effects_(var1, var2);
    double proposal_sd = proposal_sd_pairwise_(var1, var2);

    auto log_post = [&](double theta) {
        double delta = theta - current_value;
        return log_pseudoposterior_pairwise_at_delta(var1, var2, delta);
    };

    SamplerResult result = rwm_sampler(current_value, proposal_sd, log_post, rng_);

    double value = result.state[0];
    pairwise_effects_(var1, var2) = value;
    pairwise_effects_(var2, var1) = value;

    if (current_value != value) {
        double delta = value - current_value;
        residual_matrix_.col(var1) += arma::conv_to<arma::vec>::from(observations_.col(var2)) * delta;
        residual_matrix_.col(var2) += arma::conv_to<arma::vec>::from(observations_.col(var1)) * delta;
    }

    return result.accept_prob;
}


void OMRFModel::update_edge_indicator(int var1, int var2) {
    const double current_state = pairwise_effects_(var1, var2);

    const bool proposing_addition = (edge_indicators_(var1, var2) == 0);
    const double proposed_state = proposing_addition
        ? rnorm(rng_, current_state, proposal_sd_pairwise_(var1, var2))
        : 0.0;

    double log_accept = log_pseudolikelihood_ratio_interaction(
        var1, var2, proposed_state, current_state
    );

    const double inclusion_probability_ij = inclusion_probability_(var1, var2);
    const double sd = proposal_sd_pairwise_(var1, var2);
    const double scaled_pairwise_scale = pairwise_scale_ * pairwise_scaling_factors_(var1, var2);

    if (proposing_addition) {
        log_accept += R::dcauchy(proposed_state, 0.0, scaled_pairwise_scale, true);
        log_accept -= R::dnorm(proposed_state, current_state, sd, true);
        log_accept += MY_LOG(inclusion_probability_ij) - MY_LOG(1.0 - inclusion_probability_ij);
    } else {
        log_accept -= R::dcauchy(current_state, 0.0, scaled_pairwise_scale, true);
        log_accept += R::dnorm(current_state, proposed_state, sd, true);
        log_accept -= MY_LOG(inclusion_probability_ij) - MY_LOG(1.0 - inclusion_probability_ij);
    }

    if (MY_LOG(runif(rng_)) < log_accept) {
        const int updated_indicator = 1 - edge_indicators_(var1, var2);
        edge_indicators_(var1, var2) = updated_indicator;
        edge_indicators_(var2, var1) = updated_indicator;

        pairwise_effects_(var1, var2) = proposed_state;
        pairwise_effects_(var2, var1) = proposed_state;

        const double delta = proposed_state - current_state;
        residual_matrix_.col(var1) += arma::conv_to<arma::vec>::from(observations_.col(var2)) * delta;
        residual_matrix_.col(var2) += arma::conv_to<arma::vec>::from(observations_.col(var1)) * delta;
    }
}


// =============================================================================
// Main update methods
// =============================================================================

void OMRFModel::do_one_mh_step(int iteration) {
    // --- Pairwise effects sweep ---
    arma::mat accept_prob_pairwise = arma::zeros<arma::mat>(p_, p_);
    arma::umat index_mask_pairwise = arma::zeros<arma::umat>(p_, p_);

    for (size_t v1 = 0; v1 < p_ - 1; ++v1) {
        for (size_t v2 = v1 + 1; v2 < p_; ++v2) {
            double ap = update_pairwise_effect(v1, v2);
            if (edge_indicators_(v1, v2) == 1) {
                accept_prob_pairwise(v1, v2) = ap;
                index_mask_pairwise(v1, v2) = 1;
            }
        }
    }

    if (rwm_pairwise_adapter_) {
        rwm_pairwise_adapter_->update(index_mask_pairwise, accept_prob_pairwise, iteration);
    }

    // --- Main effects sweep ---
    arma::umat index_mask_main = arma::ones<arma::umat>(
        proposal_sd_main_.n_rows, proposal_sd_main_.n_cols);
    arma::mat accept_prob_main = arma::ones<arma::mat>(
        proposal_sd_main_.n_rows, proposal_sd_main_.n_cols);

    for (size_t v = 0; v < p_; ++v) {
        if (is_ordinal_variable_(v)) {
            int num_cats = num_categories_(v);
            for (int c = 0; c < num_cats; ++c) {
                accept_prob_main(v, c) = update_main_effect_parameter(v, c, -1);
            }
        } else {
            for (int p = 0; p < 2; ++p) {
                accept_prob_main(v, p) = update_main_effect_parameter(v, -1, p);
            }
        }
    }

    if (rwm_main_adapter_) {
        rwm_main_adapter_->update(index_mask_main, accept_prob_main, iteration);
    }

    invalidate_gradient_cache();
}


void OMRFModel::prepare_iteration() {
    // Shuffle edge order unconditionally to advance the RNG state consistently.
    shuffled_edge_order_ = arma_randperm(rng_, num_pairwise_);
}


void OMRFModel::update_edge_indicators() {
    for (size_t i = 0; i < num_pairwise_; ++i) {
        int idx = shuffled_edge_order_(i);
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

    const int num_variables = p_;
    const int num_missings = missing_index_.n_rows;
    const int max_num_categories = num_categories_.max();

    arma::vec category_probabilities(max_num_categories + 1);

    for (int miss = 0; miss < num_missings; miss++) {
        const int person = missing_index_(miss, 0);
        const int variable = missing_index_(miss, 1);

        const double residual_score = residual_matrix_(person, variable);
        const int num_cats = num_categories_(variable);
        const bool is_ordinal = is_ordinal_variable_(variable);

        double cumsum = 0.0;

        if (is_ordinal) {
            cumsum = 1.0;
            category_probabilities[0] = cumsum;
            for (int cat = 0; cat < num_cats; cat++) {
                const int score = cat + 1;
                const double exponent = main_effects_(variable, cat) + score * residual_score;
                cumsum += MY_EXP(exponent);
                category_probabilities[score] = cumsum;
            }
        } else {
            const int ref = baseline_category_(variable);

            cumsum = MY_EXP(
                main_effects_(variable, 0) * ref + main_effects_(variable, 1) * ref * ref
            );
            category_probabilities[0] = cumsum;

            for (int cat = 0; cat <= num_cats; cat++) {
                const int score = cat - ref;
                const double exponent =
                    main_effects_(variable, 0) * score +
                    main_effects_(variable, 1) * score * score +
                    score * residual_score;
                cumsum += MY_EXP(exponent);
                category_probabilities[cat] = cumsum;
            }
        }

        // Sample from categorical distribution via inverse transform
        const double u = runif(rng_) * cumsum;
        int sampled_score = 0;
        while (u > category_probabilities[sampled_score]) {
            sampled_score++;
        }

        int new_value = sampled_score;
        if (!is_ordinal)
            new_value -= baseline_category_(variable);
        const int old_value = observations_(person, variable);

        if (new_value != old_value) {
            observations_(person, variable) = new_value;
            observations_double_(person, variable) = static_cast<double>(new_value);

            if (is_ordinal) {
                counts_per_category_(old_value, variable)--;
                counts_per_category_(new_value, variable)++;
            } else {
                const int delta = new_value - old_value;
                const int delta_sq = new_value * new_value - old_value * old_value;
                blume_capel_stats_(0, variable) += delta;
                blume_capel_stats_(1, variable) += delta_sq;
            }

            // Incrementally update residuals across all variables
            for (int var = 0; var < num_variables; var++) {
                const double delta_score = (new_value - old_value) * pairwise_effects_(var, variable);
                residual_matrix_(person, var) += delta_score;
            }
        }
    }

    // Recompute pairwise sufficient statistics
    arma::mat ps = observations_double_.t() * observations_double_;
    pairwise_stats_ = arma::conv_to<arma::imat>::from(ps);
}


void OMRFModel::set_missing_data(const arma::imat& missing_index) {
    missing_index_ = missing_index;
    has_missing_ = (missing_index.n_rows > 0 && missing_index.n_cols == 2);
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


