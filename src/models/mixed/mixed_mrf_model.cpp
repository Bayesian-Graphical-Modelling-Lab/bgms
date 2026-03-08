#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "rng/rng_utils.h"
#include "mcmc/execution/warmup_schedule.h"


// =============================================================================
// Constructor
// =============================================================================

MixedMRFModel::MixedMRFModel(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection,
    const std::string& pseudolikelihood,
    double main_alpha,
    double main_beta,
    double pairwise_scale,
    int seed
) :
    n_(discrete_observations.n_rows),
    p_(discrete_observations.n_cols),
    q_(continuous_observations.n_cols),
    discrete_observations_(discrete_observations),
    continuous_observations_(continuous_observations),
    num_categories_(num_categories),
    is_ordinal_variable_(is_ordinal_variable),
    baseline_category_(baseline_category),
    edge_indicators_(initial_edge_indicators),
    inclusion_probability_(inclusion_probability),
    edge_selection_(edge_selection),
    edge_selection_active_(false),
    main_alpha_(main_alpha),
    main_beta_(main_beta),
    pairwise_scale_(pairwise_scale),
    use_marginal_pl_(pseudolikelihood == "marginal"),
    rng_(seed)
{
    // Dimension counts
    num_main_ = count_num_main_effects();
    num_pairwise_xx_ = (p_ * (p_ - 1)) / 2;
    num_pairwise_yy_ = (q_ * (q_ - 1)) / 2;
    num_cross_ = p_ * q_;

    max_cats_ = num_categories_.max();

    // Center Blume-Capel observations at baseline category so that all
    // downstream code operates in a shifted coordinate system where the
    // reference corresponds to zero (same convention as OMRFModel).
    for(size_t s = 0; s < p_; ++s) {
        if(!is_ordinal_variable_(s)) {
            discrete_observations_.col(s) -= baseline_category_(s);
        }
    }
    discrete_observations_dbl_ = arma::conv_to<arma::mat>::from(discrete_observations_);

    // Compute sufficient statistics
    compute_sufficient_statistics();

    // Initialize parameters to zero
    mux_ = arma::zeros<arma::mat>(p_, max_cats_);
    muy_ = arma::zeros<arma::vec>(q_);
    Kxx_ = arma::zeros<arma::mat>(p_, p_);
    Kyy_ = arma::eye<arma::mat>(q_, q_);
    Kxy_ = arma::zeros<arma::mat>(p_, q_);

    // Initialize proposal SDs
    prop_sd_mux_ = arma::ones<arma::mat>(p_, max_cats_);
    prop_sd_muy_ = arma::ones<arma::vec>(q_);
    prop_sd_Kxx_ = arma::ones<arma::mat>(p_, p_);
    prop_sd_Kyy_ = arma::ones<arma::mat>(q_, q_);
    prop_sd_Kxy_ = arma::ones<arma::mat>(p_, q_);

    // Initialize Kyy caches (Kyy starts as identity)
    Kyy_chol_ = arma::eye<arma::mat>(q_, q_);
    inv_cholesky_yy_ = arma::eye<arma::mat>(q_, q_);
    covariance_yy_ = arma::eye<arma::mat>(q_, q_);
    Kyy_log_det_ = 0.0;

    // Rank-1 Cholesky update workspace
    precision_yy_proposal_ = arma::mat(q_, q_, arma::fill::none);
    kyy_vf1_ = arma::zeros<arma::vec>(q_);
    kyy_vf2_ = arma::zeros<arma::vec>(q_);
    kyy_u1_ = arma::zeros<arma::vec>(q_);
    kyy_u2_ = arma::zeros<arma::vec>(q_);

    // Initialize conditional mean: μ_y' + 2 x Kxy Kyy_inv
    //   With Kxy = 0 and Kyy = I, this is just 1 * μ_y' = 0.
    conditional_mean_ = arma::zeros<arma::mat>(n_, q_);

    // Initialize Theta (marginal PL only): Kxx + 2 Kxy Kyy_inv Kxy'
    //   With Kxy = 0, Theta = Kxx = 0.
    if(use_marginal_pl_) {
        Theta_ = arma::zeros<arma::mat>(p_, p_);
    }

    // Initialize edge-order permutation vectors
    edge_order_xx_ = arma::regspace<arma::uvec>(0, num_pairwise_xx_ - 1);
    edge_order_yy_ = arma::regspace<arma::uvec>(0, num_pairwise_yy_ - 1);
    edge_order_xy_ = arma::regspace<arma::uvec>(0, num_cross_ - 1);
}


// =============================================================================
// Copy constructor
// =============================================================================

MixedMRFModel::MixedMRFModel(const MixedMRFModel& other)
    : BaseModel(other),
      n_(other.n_),
      p_(other.p_),
      q_(other.q_),
      num_main_(other.num_main_),
      num_pairwise_xx_(other.num_pairwise_xx_),
      num_pairwise_yy_(other.num_pairwise_yy_),
      num_cross_(other.num_cross_),
      discrete_observations_(other.discrete_observations_),
      discrete_observations_dbl_(other.discrete_observations_dbl_),
      continuous_observations_(other.continuous_observations_),
      num_categories_(other.num_categories_),
      max_cats_(other.max_cats_),
      is_ordinal_variable_(other.is_ordinal_variable_),
      baseline_category_(other.baseline_category_),
      missing_index_discrete_(other.missing_index_discrete_),
      missing_index_continuous_(other.missing_index_continuous_),
      has_missing_(other.has_missing_),
      counts_per_category_(other.counts_per_category_),
      blume_capel_stats_(other.blume_capel_stats_),
      mux_(other.mux_),
      muy_(other.muy_),
      Kxx_(other.Kxx_),
      Kyy_(other.Kyy_),
      Kxy_(other.Kxy_),
      edge_indicators_(other.edge_indicators_),
      inclusion_probability_(other.inclusion_probability_),
      edge_selection_(other.edge_selection_),
      edge_selection_active_(other.edge_selection_active_),
      main_alpha_(other.main_alpha_),
      main_beta_(other.main_beta_),
      pairwise_scale_(other.pairwise_scale_),
      prop_sd_mux_(other.prop_sd_mux_),
      prop_sd_muy_(other.prop_sd_muy_),
      prop_sd_Kxx_(other.prop_sd_Kxx_),
      prop_sd_Kyy_(other.prop_sd_Kyy_),
      prop_sd_Kxy_(other.prop_sd_Kxy_),
      total_warmup_(other.total_warmup_),
      Kyy_chol_(other.Kyy_chol_),
      inv_cholesky_yy_(other.inv_cholesky_yy_),
      covariance_yy_(other.covariance_yy_),
      Kyy_log_det_(other.Kyy_log_det_),
      Theta_(other.Theta_),
      conditional_mean_(other.conditional_mean_),
      kyy_constants_(other.kyy_constants_),
      precision_yy_proposal_(other.precision_yy_proposal_),
      kyy_v1_(other.kyy_v1_),
      kyy_v2_(other.kyy_v2_),
      kyy_vf1_(other.kyy_vf1_),
      kyy_vf2_(other.kyy_vf2_),
      kyy_u1_(other.kyy_u1_),
      kyy_u2_(other.kyy_u2_),
      gradient_cache_valid_(false),
      use_marginal_pl_(other.use_marginal_pl_),
      rng_(other.rng_),
      edge_order_xx_(other.edge_order_xx_),
      edge_order_yy_(other.edge_order_yy_),
      edge_order_xy_(other.edge_order_xy_)
{
}


// =============================================================================
// Sufficient statistics
// =============================================================================

void MixedMRFModel::compute_sufficient_statistics() {
    // Category counts for ordinal variables
    counts_per_category_ = arma::zeros<arma::imat>(max_cats_ + 1, p_);
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(size_t i = 0; i < n_; ++i) {
                int cat = discrete_observations_(i, s);
                if(cat >= 0 && cat <= num_categories_(s)) {
                    counts_per_category_(cat, s)++;
                }
            }
        }
    }

    // Blume-Capel statistics (linear and quadratic sums of centered obs)
    blume_capel_stats_ = arma::zeros<arma::imat>(2, p_);
    for(size_t s = 0; s < p_; ++s) {
        if(!is_ordinal_variable_(s)) {
            for(size_t i = 0; i < n_; ++i) {
                int val = discrete_observations_(i, s);  // already centered
                blume_capel_stats_(0, s) += val;
                blume_capel_stats_(1, s) += val * val;
            }
        }
    }
}


size_t MixedMRFModel::count_num_main_effects() const {
    size_t count = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            count += num_categories_(s);
        } else {
            count += 2;  // linear α and quadratic β
        }
    }
    return count;
}


// =============================================================================
// Cache maintenance
// =============================================================================

void MixedMRFModel::recompute_conditional_mean() {
    // conditional_mean_ = 1*μ_y' + 2 * discrete_obs * Kxy * covariance_yy_
    conditional_mean_ = arma::repmat(muy_.t(), n_, 1) +
                        2.0 * discrete_observations_dbl_ * Kxy_ * covariance_yy_;
}

void MixedMRFModel::recompute_Kyy_decomposition() {
    Kyy_chol_ = arma::chol(Kyy_);                // upper Cholesky: Kyy = R'R
    arma::inv(inv_cholesky_yy_, arma::trimatu(Kyy_chol_));
    covariance_yy_ = inv_cholesky_yy_ * inv_cholesky_yy_.t();
    Kyy_log_det_ = cholesky_helpers::get_log_det(Kyy_chol_);
}

void MixedMRFModel::recompute_Theta() {
    // Θ = Kxx + 2 Kxy covariance_yy_ Kxy'
    Theta_ = Kxx_ + 2.0 * Kxy_ * covariance_yy_ * Kxy_.t();
}


// =============================================================================
// Parameter vectorization
// =============================================================================

// NUTS vectorization order (excludes Kyy — sampled by MH separately):
//   1. mux_: per-variable (ordinal: C_s thresholds; BC: 2 coefficients)
//   2. Kxx_: upper-triangular, row-major  — p(p-1)/2
//   3. muy_: all q means
//   4. Kxy_: all p*q entries, row-major
//
// Storage vectorization order (includes Kyy):
//   1–4. Same as NUTS order
//   5. Kyy_: upper-triangle including diagonal — q(q+1)/2

size_t MixedMRFModel::parameter_dimension() const {
    if(!edge_selection_active_) {
        return full_parameter_dimension();
    }
    // Count active NUTS parameters only (no Kyy)
    size_t dim = num_main_ + q_;  // mux + muy always active

    // Active Kxx edges
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(gxx(i, j)) dim++;
        }
    }

    // Active Kxy edges
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(gxy(i, j)) dim++;
        }
    }

    return dim;
}

size_t MixedMRFModel::full_parameter_dimension() const {
    // NUTS block: mux + Kxx upper-tri + muy + Kxy full (no Kyy)
    return num_main_ + num_pairwise_xx_ + q_ + num_cross_;
}

size_t MixedMRFModel::storage_dimension() const {
    // All parameters including Kyy
    return num_main_ + num_pairwise_xx_ + q_ +
           (q_ * (q_ + 1)) / 2 + num_cross_;
}

arma::vec MixedMRFModel::get_vectorized_parameters() const {
    // Active NUTS parameters only (excludes Kyy, excludes inactive edges)
    arma::vec out(parameter_dimension());
    size_t idx = 0;

    // 1. mux_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                out(idx++) = mux_(s, c);
            }
        } else {
            out(idx++) = mux_(s, 0);
            out(idx++) = mux_(s, 1);
        }
    }

    // 2. Kxx_ upper-triangular (active edges only when selection is active)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(!edge_selection_active_ || gxx(i, j) == 1) {
                out(idx++) = Kxx_(i, j);
            }
        }
    }

    // 3. muy_
    for(size_t j = 0; j < q_; ++j) {
        out(idx++) = muy_(j);
    }

    // 4. Kxy_ row-major (active edges only when selection is active)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(!edge_selection_active_ || gxy(i, j) == 1) {
                out(idx++) = Kxy_(i, j);
            }
        }
    }

    return out;
}

arma::vec MixedMRFModel::get_full_vectorized_parameters() const {
    // All NUTS parameters, fixed size (inactive edges are 0, no Kyy)
    arma::vec out(full_parameter_dimension(), arma::fill::zeros);
    size_t idx = 0;

    // 1. mux_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                out(idx++) = mux_(s, c);
            }
        } else {
            out(idx++) = mux_(s, 0);
            out(idx++) = mux_(s, 1);
        }
    }

    // 2. Kxx_ upper-triangular (all entries, zeros for inactive)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            out(idx++) = Kxx_(i, j);
        }
    }

    // 3. muy_
    for(size_t j = 0; j < q_; ++j) {
        out(idx++) = muy_(j);
    }

    // 4. Kxy_ row-major (all entries, zeros for inactive)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            out(idx++) = Kxy_(i, j);
        }
    }

    return out;
}

arma::vec MixedMRFModel::get_storage_vectorized_parameters() const {
    // All parameters including Kyy, fixed size
    arma::vec out(storage_dimension(), arma::fill::zeros);
    size_t idx = 0;

    // 1. mux_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                out(idx++) = mux_(s, c);
            }
        } else {
            out(idx++) = mux_(s, 0);
            out(idx++) = mux_(s, 1);
        }
    }

    // 2. Kxx_ upper-triangular
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            out(idx++) = Kxx_(i, j);
        }
    }

    // 3. muy_
    for(size_t j = 0; j < q_; ++j) {
        out(idx++) = muy_(j);
    }

    // 4. Kxy_ row-major
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            out(idx++) = Kxy_(i, j);
        }
    }

    // 5. Kyy_ upper-triangle including diagonal
    for(size_t i = 0; i < q_; ++i) {
        for(size_t j = i; j < q_; ++j) {
            out(idx++) = Kyy_(i, j);
        }
    }

    return out;
}

void MixedMRFModel::set_vectorized_parameters(const arma::vec& params) {
    // Unpack NUTS block only (no Kyy)
    size_t idx = 0;

    // 1. mux_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                mux_(s, c) = params(idx++);
            }
        } else {
            mux_(s, 0) = params(idx++);
            mux_(s, 1) = params(idx++);
        }
    }

    // 2. Kxx_ upper-triangular (active edges only when selection is active)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(!edge_selection_active_ || gxx(i, j) == 1) {
                Kxx_(i, j) = params(idx);
                Kxx_(j, i) = params(idx);
                idx++;
            }
        }
    }

    // 3. muy_
    for(size_t j = 0; j < q_; ++j) {
        muy_(j) = params(idx++);
    }

    // 4. Kxy_ row-major (active edges only when selection is active)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(!edge_selection_active_ || gxy(i, j) == 1) {
                Kxy_(i, j) = params(idx++);
            }
        }
    }

    // Refresh caches (Kyy unchanged, so no Kyy decomposition update needed)
    recompute_conditional_mean();
    if(use_marginal_pl_) {
        recompute_Theta();
    }
}

arma::vec MixedMRFModel::get_active_inv_mass() const {
    if(!edge_selection_active_) {
        return inv_mass_;
    }

    arma::vec active(parameter_dimension());
    // Main effects: always active
    active.head(num_main_) = inv_mass_.head(num_main_);

    size_t offset_full = num_main_;
    size_t offset_active = num_main_;

    // Kxx active edges
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(gxx(i, j) == 1) {
                active(offset_active++) = inv_mass_(offset_full);
            }
            offset_full++;
        }
    }

    // muy: always active
    for(size_t j = 0; j < q_; ++j) {
        active(offset_active++) = inv_mass_(offset_full++);
    }

    // Kxy active edges
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(gxy(i, j) == 1) {
                active(offset_active++) = inv_mass_(offset_full);
            }
            offset_full++;
        }
    }

    return active;
}

arma::ivec MixedMRFModel::get_vectorized_indicator_parameters() {
    size_t total = num_pairwise_xx_ + num_pairwise_yy_ + num_cross_;
    arma::ivec out(total);
    size_t idx = 0;

    // 1. Upper-triangle of Gxx
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            out(idx++) = gxx(i, j);
        }
    }

    // 2. Upper-triangle of Gyy
    for(size_t i = 0; i < q_ - 1; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            out(idx++) = gyy(i, j);
        }
    }

    // 3. Full Gxy block row-major
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            out(idx++) = gxy(i, j);
        }
    }

    return out;
}


// =============================================================================
// Infrastructure
// =============================================================================

void MixedMRFModel::set_seed(int seed) {
    rng_ = SafeRNG(seed);
}

std::unique_ptr<BaseModel> MixedMRFModel::clone() const {
    return std::make_unique<MixedMRFModel>(*this);
}


// =============================================================================
// Missing data imputation
// =============================================================================

void MixedMRFModel::set_missing_data(const arma::imat& missing_discrete,
                                      const arma::imat& missing_continuous) {
    missing_index_discrete_ = missing_discrete;
    missing_index_continuous_ = missing_continuous;
    has_missing_ = (missing_index_discrete_.n_rows > 0 ||
                    missing_index_continuous_.n_rows > 0);
}

void MixedMRFModel::impute_missing() {
    if(!has_missing_) return;

    // --- Phase 1: Impute discrete entries ---
    const int num_disc_missing = missing_index_discrete_.n_rows;
    if(num_disc_missing > 0) {
        arma::vec category_probabilities(max_cats_ + 1);

        for(int miss = 0; miss < num_disc_missing; miss++) {
            const int person = missing_index_discrete_(miss, 0);
            const int variable = missing_index_discrete_(miss, 1);
            const int num_cats = num_categories_(variable);

            // Rest score: sum_t x_vt * Kxx(t,s) + 2 * sum_j y_vj * Kxy(s,j)
            // Kxx diagonal is zero, so no self-interaction subtraction needed
            double rest_v = 0.0;
            for(size_t t = 0; t < p_; t++) {
                rest_v += discrete_observations_dbl_(person, t) * Kxx_(t, variable);
            }
            for(size_t j = 0; j < q_; j++) {
                rest_v += 2.0 * continuous_observations_(person, j) * Kxy_(variable, j);
            }

            double cumsum = 0.0;

            if(is_ordinal_variable_(variable)) {
                // P(x=0) = 1, P(x=c) ∝ exp(c * rest + mux(s, c-1))
                cumsum = 1.0;
                category_probabilities(0) = cumsum;
                for(int c = 1; c <= num_cats; c++) {
                    double exponent = static_cast<double>(c) * rest_v +
                                      mux_(variable, c - 1);
                    cumsum += std::exp(exponent);
                    category_probabilities(c) = cumsum;
                }
            } else {
                // Blume-Capel: categories centered at baseline
                const int ref = baseline_category_(variable);
                double alpha = mux_(variable, 0);
                double beta = mux_(variable, 1);
                cumsum = 0.0;
                for(int cat = 0; cat <= num_cats; cat++) {
                    const int score = cat - ref;
                    double exponent = alpha * score +
                                      beta * score * score +
                                      score * rest_v;
                    cumsum += std::exp(exponent);
                    category_probabilities(cat) = cumsum;
                }
            }

            // Sample via inverse-transform
            double u = runif(rng_) * cumsum;
            int sampled = 0;
            while(u > category_probabilities(sampled)) {
                sampled++;
            }

            int new_value = sampled;
            if(!is_ordinal_variable_(variable)) {
                new_value -= baseline_category_(variable);
            }
            const int old_value = discrete_observations_(person, variable);

            if(new_value != old_value) {
                discrete_observations_(person, variable) = new_value;
                discrete_observations_dbl_(person, variable) =
                    static_cast<double>(new_value);

                if(is_ordinal_variable_(variable)) {
                    counts_per_category_(old_value, variable)--;
                    counts_per_category_(new_value, variable)++;
                } else {
                    blume_capel_stats_(0, variable) += (new_value - old_value);
                    blume_capel_stats_(1, variable) +=
                        (new_value * new_value - old_value * old_value);
                }
            }
        }
    }

    // --- Phase 2: Refresh conditional_mean_ (depends on discrete data) ---
    if(num_disc_missing > 0 && missing_index_continuous_.n_rows > 0) {
        recompute_conditional_mean();
    }

    // --- Phase 3: Impute continuous entries ---
    const int num_cont_missing = missing_index_continuous_.n_rows;
    if(num_cont_missing > 0) {
        for(int miss = 0; miss < num_cont_missing; miss++) {
            const int person = missing_index_continuous_(miss, 0);
            const int variable = missing_index_continuous_(miss, 1);

            // Conditional: y_vj | y_{v,-j}, x ~ N(mu*, 1/Kyy_jj)
            // mu* = M_vj - (1/Kyy_jj) * sum_{k!=j} Kyy_jk * (y_vk - M_vk)
            double cond_mean = conditional_mean_(person, variable);
            for(size_t k = 0; k < q_; k++) {
                if(k != static_cast<size_t>(variable)) {
                    cond_mean -= (Kyy_(variable, k) / Kyy_(variable, variable)) *
                        (continuous_observations_(person, k) -
                         conditional_mean_(person, k));
                }
            }
            double cond_sd = std::sqrt(1.0 / Kyy_(variable, variable));

            continuous_observations_(person, variable) =
                rnorm(rng_, cond_mean, cond_sd);
        }
    }

    // Invalidate gradient cache (observations changed)
    invalidate_gradient_cache();
}


// =============================================================================
// Stubs (to be implemented in later phases)
// =============================================================================

void MixedMRFModel::do_one_metropolis_step(int iteration) {
    // Step 1: Update all main effects (ordinal thresholds or BC α/β)
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c)
                update_main_effect(s, c, iteration);
        } else {
            update_main_effect(s, 0, iteration);  // linear α
            update_main_effect(s, 1, iteration);  // quadratic β
        }
    }

    // Step 2: Update all continuous means
    for(size_t j = 0; j < q_; ++j)
        update_continuous_mean(j, iteration);

    // Step 3: Update Kxx (upper triangle, edge-gated)
    for(size_t i = 0; i < p_ - 1; ++i)
        for(size_t j = i + 1; j < p_; ++j)
            if(!edge_selection_active_ || gxx(i, j) == 1)
                update_Kxx(i, j, iteration);

    // Step 4: Update Kyy (off-diag + diagonal, edge-gated)
    if(q_ >= 2) {
        for(size_t i = 0; i < q_ - 1; ++i)
            for(size_t j = i + 1; j < q_; ++j)
                if(!edge_selection_active_ || gyy(i, j) == 1)
                    update_Kyy_offdiag(i, j, iteration);
    }
    for(size_t i = 0; i < q_; ++i)
        update_Kyy_diag(i, iteration);

    // Step 5: Update Kxy (edge-gated)
    for(size_t i = 0; i < p_; ++i)
        for(size_t j = 0; j < q_; ++j)
            if(!edge_selection_active_ || gxy(i, j) == 1)
                update_Kxy(i, j, iteration);

    // Edge-indicator updates are handled by ChainRunner, not here.
    // (Matches the OMRF pattern; avoids double-counting indicator proposals.)
}

void MixedMRFModel::do_kyy_metropolis_step(int iteration) {
    // Off-diagonal Kyy (edge-gated)
    if(q_ >= 2) {
        for(size_t i = 0; i < q_ - 1; ++i)
            for(size_t j = i + 1; j < q_; ++j)
                if(!edge_selection_active_ || gyy(i, j) == 1)
                    update_Kyy_offdiag(i, j, iteration);
    }
    // Diagonal Kyy (always updated)
    for(size_t i = 0; i < q_; ++i)
        update_Kyy_diag(i, iteration);
}

void MixedMRFModel::update_edge_indicators() {
    if(!edge_selection_active_) return;

    invalidate_gradient_cache();

    // Discrete-discrete edges (shuffled order)
    for(size_t e = 0; e < num_pairwise_xx_; ++e) {
        size_t idx = edge_order_xx_(e);
        // Decode upper-triangle index to (i, j)
        size_t i = 0, j = 1;
        size_t count = 0;
        for(i = 0; i < p_ - 1; ++i) {
            size_t row_len = p_ - 1 - i;
            if(count + row_len > idx) {
                j = i + 1 + (idx - count);
                break;
            }
            count += row_len;
        }
        update_edge_indicator_Kxx(i, j);
    }

    // Continuous-continuous edges (shuffled order)
    for(size_t e = 0; e < num_pairwise_yy_; ++e) {
        size_t idx = edge_order_yy_(e);
        size_t i = 0, j = 1;
        size_t count = 0;
        for(i = 0; i < q_ - 1; ++i) {
            size_t row_len = q_ - 1 - i;
            if(count + row_len > idx) {
                j = i + 1 + (idx - count);
                break;
            }
            count += row_len;
        }
        update_edge_indicator_Kyy(i, j);
    }

    // Cross edges (shuffled order)
    for(size_t e = 0; e < num_cross_; ++e) {
        size_t idx = edge_order_xy_(e);
        size_t i = idx / q_;
        size_t j = idx % q_;
        update_edge_indicator_Kxy(i, j);
    }
}

void MixedMRFModel::initialize_graph() {
    // Draw initial graph from prior inclusion probabilities.
    // Zero out parameters for excluded edges.
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(runif(rng_) >= inclusion_probability_(i, j)) {
                set_gxx(i, j, 0);
                Kxx_(i, j) = 0.0;
                Kxx_(j, i) = 0.0;
            }
        }
    }

    for(size_t i = 0; i < q_ - 1; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            if(runif(rng_) >= inclusion_probability_(p_ + i, p_ + j)) {
                set_gyy(i, j, 0);
                Kyy_(i, j) = 0.0;
                Kyy_(j, i) = 0.0;
            }
        }
    }
    // Recompute Kyy decomposition after potential zeroing
    recompute_Kyy_decomposition();
    recompute_conditional_mean();

    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(runif(rng_) >= inclusion_probability_(i, p_ + j)) {
                set_gxy(i, j, 0);
                Kxy_(i, j) = 0.0;
            }
        }
    }
    recompute_conditional_mean();
    if(use_marginal_pl_) recompute_Theta();
}

void MixedMRFModel::prepare_iteration() {
    // Shuffle edge-update order to avoid order bias.
    // Always called, even when edge selection is off, to keep RNG consistent.
    edge_order_xx_ = arma_randperm(rng_, num_pairwise_xx_);
    edge_order_yy_ = arma_randperm(rng_, num_pairwise_yy_);
    edge_order_xy_ = arma_randperm(rng_, num_cross_);
}

void MixedMRFModel::init_metropolis_adaptation(const WarmupSchedule& schedule) {
    total_warmup_ = schedule.total_warmup;
}

void MixedMRFModel::tune_proposal_sd(int /*iteration*/, const WarmupSchedule& /*schedule*/) {
    // Robbins-Monro adaptation is embedded in each MH update function,
    // gated by iteration < total_warmup_. No separate tuning pass needed.
}
