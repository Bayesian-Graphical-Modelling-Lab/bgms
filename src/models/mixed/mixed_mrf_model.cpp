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
    Kyy_inv_ = arma::eye<arma::mat>(q_, q_);
    Kyy_log_det_ = 0.0;

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
      Kyy_inv_(other.Kyy_inv_),
      Kyy_chol_(other.Kyy_chol_),
      Kyy_log_det_(other.Kyy_log_det_),
      Theta_(other.Theta_),
      conditional_mean_(other.conditional_mean_),
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
    // conditional_mean_ = 1*μ_y' + 2 * discrete_obs * Kxy * Kyy_inv
    conditional_mean_ = arma::repmat(muy_.t(), n_, 1) +
                        2.0 * discrete_observations_dbl_ * Kxy_ * Kyy_inv_;
}

void MixedMRFModel::recompute_Kyy_decomposition() {
    Kyy_chol_ = arma::chol(Kyy_);                // upper Cholesky
    Kyy_inv_ = arma::inv_sympd(Kyy_);
    Kyy_log_det_ = 2.0 * arma::sum(arma::log(Kyy_chol_.diag()));
}

void MixedMRFModel::recompute_Theta() {
    // Θ = Kxx + 2 Kxy Kyy_inv Kxy'
    Theta_ = Kxx_ + 2.0 * Kxy_ * Kyy_inv_ * Kxy_.t();
}


// =============================================================================
// Parameter vectorization
// =============================================================================

// Vectorization order (free parameters):
//   1. mux_: per-variable (ordinal: C_s thresholds; BC: 2 coefficients)
//   2. Kxx_: upper-triangular, row-major  — p(p-1)/2
//   3. muy_: all q means
//   4. Kyy_: upper-triangle including diagonal — q(q+1)/2
//   5. Kxy_: all p*q entries, row-major

size_t MixedMRFModel::parameter_dimension() const {
    if(!edge_selection_active_) {
        return full_parameter_dimension();
    }
    // Count active parameters only
    size_t dim = num_main_ + q_;  // mux + muy always active

    // Active Kxx edges
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(gxx(i, j)) dim++;
        }
    }

    // Kyy diagonal always active; off-diagonal gated by edge indicators
    dim += q_;  // diagonal
    for(size_t i = 0; i < q_ - 1; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            if(gyy(i, j)) dim++;
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
    // mux + Kxx upper-tri + muy + Kyy upper-tri-with-diag + Kxy full
    return num_main_ + num_pairwise_xx_ + q_ +
           (q_ * (q_ + 1)) / 2 + num_cross_;
}

arma::vec MixedMRFModel::get_vectorized_parameters() const {
    arma::vec out(full_parameter_dimension(), arma::fill::zeros);
    size_t idx = 0;

    // 1. mux_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                out(idx++) = mux_(s, c);
            }
        } else {
            out(idx++) = mux_(s, 0);  // linear α
            out(idx++) = mux_(s, 1);  // quadratic β
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

    // 4. Kyy_ upper-triangle including diagonal
    for(size_t i = 0; i < q_; ++i) {
        for(size_t j = i; j < q_; ++j) {
            out(idx++) = Kyy_(i, j);
        }
    }

    // 5. Kxy_ row-major
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            out(idx++) = Kxy_(i, j);
        }
    }

    return out;
}

arma::vec MixedMRFModel::get_full_vectorized_parameters() const {
    return get_vectorized_parameters();
}

void MixedMRFModel::set_vectorized_parameters(const arma::vec& params) {
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

    // 2. Kxx_ upper-triangular (mirror to lower)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            Kxx_(i, j) = params(idx);
            Kxx_(j, i) = params(idx);
            idx++;
        }
    }

    // 3. muy_
    for(size_t j = 0; j < q_; ++j) {
        muy_(j) = params(idx++);
    }

    // 4. Kyy_ upper-triangle including diagonal (mirror off-diag)
    for(size_t i = 0; i < q_; ++i) {
        for(size_t j = i; j < q_; ++j) {
            Kyy_(i, j) = params(idx);
            if(i != j) Kyy_(j, i) = params(idx);
            idx++;
        }
    }

    // 5. Kxy_ row-major
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            Kxy_(i, j) = params(idx++);
        }
    }

    // Refresh all caches
    recompute_Kyy_decomposition();
    recompute_conditional_mean();
    if(use_marginal_pl_) {
        recompute_Theta();
    }
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
// Stubs (to be implemented in later phases)
// =============================================================================

void MixedMRFModel::do_one_metropolis_step(int /*iteration*/) {
    // Phase B/C: MH updates for all parameter groups
}

void MixedMRFModel::update_edge_indicators() {
    // Phase D: reversible-jump edge birth/death moves
}

void MixedMRFModel::initialize_graph() {
    // Phase D: random graph initialization for edge selection
    // For now, start with all edges included (matching initial_edge_indicators)
}

void MixedMRFModel::prepare_iteration() {
    // Shuffle edge-update order to avoid order bias.
    // Always called, even when edge selection is off, to keep RNG consistent.
    edge_order_xx_ = arma_randperm(rng_, num_pairwise_xx_);
    edge_order_yy_ = arma_randperm(rng_, num_pairwise_yy_);
    edge_order_xy_ = arma_randperm(rng_, num_cross_);
}

void MixedMRFModel::init_metropolis_adaptation(const WarmupSchedule& /*schedule*/) {
    // Phase F: initialize Robbins-Monro controllers
}

void MixedMRFModel::tune_proposal_sd(int /*iteration*/, const WarmupSchedule& /*schedule*/) {
    // Phase F: Robbins-Monro proposal-SD tuning
}
