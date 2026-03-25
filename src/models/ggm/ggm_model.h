#pragma once

#include <array>
#include <memory>
#include "models/base_model.h"
#include "math/cholesky_helpers.h"
#include "rng/rng_utils.h"


/**
 * GGMModel - Gaussian Graphical Model
 *
 * Bayesian inference on the precision matrix (inverse covariance) of a
 * multivariate Gaussian via element-wise Metropolis-Hastings. Edge
 * selection uses a spike-and-slab prior with Cauchy slab.
 *
 * The Cholesky factor of the precision matrix is maintained incrementally
 * through rank-1 updates/downdates after each element change.
 */
class GGMModel : public BaseModel {
public:

    /**
     * Construct from raw observations.
     *
     * Computes the sufficient-statistic matrix S = X'X from the raw data.
     * When na_impute is true, the observation matrix is retained for
     * full-conditional imputation of missing entries.
     *
     * @param observations          Raw data matrix (n x p)
     * @param inclusion_probability Prior inclusion probabilities for each edge
     * @param initial_edge_indicators Initial edge inclusion indicators
     * @param edge_selection        Enable edge selection (spike-and-slab)
     * @param pairwise_scale        Scale parameter of Cauchy slab prior
     * @param na_impute             Retain observations for missing-data imputation
     */
    GGMModel(
            const arma::mat& observations,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection,
            const double pairwise_scale,
            const bool na_impute = false
    ) : n_(observations.n_rows),
        p_(observations.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(observations.t() * observations),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        pairwise_scale_(pairwise_scale),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        proposal_sds_(arma::vec(dim_, arma::fill::ones) * 0.25),
        observations_(na_impute ? observations : arma::mat()),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none))
    {}

    /**
     * Construct from sufficient statistics.
     *
     * Bypasses raw data storage; useful when only X'X and n are available.
     * Missing-data imputation is not supported with this constructor.
     *
     * @param n                     Number of observations
     * @param suf_stat              Sufficient-statistic matrix X'X (p x p)
     * @param inclusion_probability Prior inclusion probabilities for each edge
     * @param initial_edge_indicators Initial edge inclusion indicators
     * @param edge_selection        Enable edge selection (spike-and-slab)
     * @param pairwise_scale        Scale parameter of Cauchy slab prior
     */
    GGMModel(
            const int n,
            const arma::mat& suf_stat,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection,
            const double pairwise_scale
    ) : n_(n),
        p_(suf_stat.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(suf_stat),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        pairwise_scale_(pairwise_scale),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        proposal_sds_(arma::vec(dim_, arma::fill::ones) * 0.25),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none))
    {}

    /** Copy constructor for cloning (required for parallel chains). */
    GGMModel(const GGMModel& other)
        : BaseModel(other),
          n_(other.n_),
          p_(other.p_),
          dim_(other.dim_),
          suf_stat_(other.suf_stat_),
          inclusion_probability_(other.inclusion_probability_),
          edge_selection_(other.edge_selection_),
          pairwise_scale_(other.pairwise_scale_),
          precision_matrix_(other.precision_matrix_),
          cholesky_of_precision_(other.cholesky_of_precision_),
          inv_cholesky_of_precision_(other.inv_cholesky_of_precision_),
          covariance_matrix_(other.covariance_matrix_),
          edge_indicators_(other.edge_indicators_),
          vectorized_parameters_(other.vectorized_parameters_),
          vectorized_indicator_parameters_(other.vectorized_indicator_parameters_),
          proposal_sds_(other.proposal_sds_),
          total_warmup_(other.total_warmup_),
          rng_(other.rng_),
          observations_(other.observations_),
          has_missing_(other.has_missing_),
          missing_index_(other.missing_index_),
          precision_proposal_(other.precision_proposal_)
    {}

    /** @return false (GGM has no gradient implementation). */
    bool has_gradient()        const override { return false; }
    /** @return true (GGM supports adaptive Metropolis). */
    bool has_adaptive_metropolis()     const override { return true; }
    /** @return true when edge selection is enabled. */
    bool has_edge_selection()  const override { return edge_selection_; }
    /** @return true when missing-data imputation is active. */
    bool has_missing_data()    const override { return has_missing_; }

    /** Impute missing entries from full-conditional normal distributions. */
    void impute_missing() override;

    /**
     * Register missing-data locations.
     *
     * @param missing_index  M x 2 matrix of 0-based (row, col) indices
     * @throws std::logic_error if the model was constructed without na_impute
     */
    void set_missing_data(const arma::imat& missing_index) {
        if (observations_.n_elem == 0) {
            throw std::logic_error(
                "set_missing_data() called but observations_ is empty. "
                "The model must be constructed with na_impute=true to retain observations.");
        }
        missing_index_ = missing_index;
        has_missing_ = (missing_index.n_rows > 0 && missing_index.n_cols == 2);
    }

    /**
     * Enable or disable edge-selection proposals.
     * @param active  true to enable edge add-delete moves
     */
    void set_edge_selection_active(bool active) override {
        edge_selection_active_ = active;
    }

    /** Draw initial edge states from prior inclusion probabilities. */
    void initialize_graph() override;
    /** Store warmup length for Robbins-Monro proposal-SD adaptation. */
    void init_metropolis_adaptation(const WarmupSchedule& schedule) override;

    /** No-op: GGM handles edge indicator updates inside do_one_metropolis_step(). */
    void update_edge_indicators() override {}

    /**
     * Compute the Gaussian log-likelihood for a given precision matrix.
     * @param omega  Precision matrix
     */
    double log_likelihood(const arma::mat& omega) const { return log_density_impl(omega,  arma::chol(omega)); };
    /** Compute the Gaussian log-likelihood at the current precision matrix. */
    double log_likelihood()                       const { return log_density_impl(precision_matrix_, cholesky_of_precision_); }

    /**
     * Perform one full Metropolis sweep.
     *
     * Iterates over all off-diagonal entries (edge updates), all diagonal
     * entries, and (when active) all edge indicator add-delete moves.
     *
     * @param iteration  Current iteration index (for Robbins-Monro adaptation)
     */
    void do_one_metropolis_step(int iteration = -1) override;

    /** @return Number of upper-triangle elements p(p+1)/2. */
    size_t parameter_dimension() const override { return dim_; }
    /** @return Same as parameter_dimension() for GGM. */
    size_t full_parameter_dimension() const override { return dim_; }

    /**
     * Set random seed for reproducibility.
     * @param seed  Integer seed value
     */
    void set_seed(int seed) override {
        rng_ = SafeRNG(seed);
    }

    /** @return Upper triangle of the precision matrix as a vector. */
    arma::vec get_vectorized_parameters() const override {
        return extract_upper_triangle();
    }

    /** @return Same as get_vectorized_parameters() for GGM. */
    arma::vec get_full_vectorized_parameters() const override {
        return extract_upper_triangle();
    }

    /** @return Upper triangle of the edge-indicator matrix as an integer vector. */
    arma::ivec get_vectorized_indicator_parameters() override {
        size_t e = 0;
        for (size_t i = 0; i < p_; ++i) {
            for (size_t j = i; j < p_; ++j) {
                vectorized_indicator_parameters_(e) = edge_indicators_(i, j);
                ++e;
            }
        }
        return vectorized_indicator_parameters_;
    }

    /** @return Reference to the model's random number generator. */
    SafeRNG& get_rng() override { return rng_; }

    /** @return Current edge-indicator matrix. */
    const arma::imat& get_edge_indicators() const override {
        return edge_indicators_;
    }

    /** @return Mutable reference to the prior inclusion-probability matrix. */
    arma::mat& get_inclusion_probability() override {
        return inclusion_probability_;
    }

    /** @return Number of variables (p). */
    int get_num_variables() const override {
        return static_cast<int>(p_);
    }

    /** @return Number of unique off-diagonal pairs p(p-1)/2. */
    int get_num_pairwise() const override {
        return static_cast<int>(p_ * (p_ - 1) / 2);
    }

    /** @return Deep copy of this model. */
    std::unique_ptr<BaseModel> clone() const override {
        return std::make_unique<GGMModel>(*this);
    }

private:

    /** Extract upper triangle of the precision matrix into a vector. */
    arma::vec extract_upper_triangle() const {
        arma::vec result(dim_);
        size_t e = 0;
        for (size_t i = 0; i < p_; ++i) {
            for (size_t j = i; j < p_; ++j) {
                result(e) = precision_matrix_(i, j);
                ++e;
            }
        }
        return result;
    }

    /// Number of observations.
    size_t n_;
    /// Number of variables.
    size_t p_;
    /// Number of upper-triangle elements: p(p+1)/2.
    size_t dim_;
    /// Sufficient-statistic matrix X'X (p x p).
    arma::mat suf_stat_;
    /// Prior inclusion probabilities (p x p, symmetric).
    arma::mat inclusion_probability_;
    /// Whether the model was constructed with edge selection.
    bool edge_selection_;
    /// Whether edge add-delete proposals are currently active.
    bool edge_selection_active_ = false;
    /// Scale parameter of the Cauchy slab prior on off-diagonal elements.
    double pairwise_scale_;

    /// Precision matrix Omega, its Cholesky factor R (Omega = R'R),
    /// inverse Cholesky factor, and covariance matrix.
    arma::mat precision_matrix_, cholesky_of_precision_, inv_cholesky_of_precision_, covariance_matrix_;
    /// Current edge-indicator matrix (p x p, symmetric, 0/1).
    arma::imat edge_indicators_;
    /// Pre-allocated storage returned by get_vectorized_parameters().
    arma::vec vectorized_parameters_;
    /// Pre-allocated storage returned by get_vectorized_indicator_parameters().
    arma::ivec vectorized_indicator_parameters_;

    /// Proposal standard deviations for Metropolis updates (one per element).
    arma::vec proposal_sds_;
    /// Total number of warmup iterations (for Robbins-Monro adaptation).
    int total_warmup_ = 0;
    /// Random number generator.
    SafeRNG rng_;

    /// Raw observation matrix (n x p), only populated when na_impute=true.
    arma::mat observations_;
    /// Whether missing-data imputation is active.
    bool has_missing_ = false;
    /// M x 2 matrix of 0-based (row, col) indices of missing entries.
    arma::imat missing_index_;

    /**
     * Incrementally adjust S = X'X after replacing one observation value.
     *
     * @param variable  Column index of the changed variable
     * @param person    Row index of the changed observation
     * @param delta     Change in value (new - old)
     */
    void update_suf_stat_for_imputation(int variable, int person, double delta);

    /// Scratch matrix for proposed precision values.
    arma::mat precision_proposal_;

    /**
     * Workspace for conditional precision reparameterization.
     *
     * - [0] Phi_q1q
     * - [1] Phi_q1q1
     * - [2] omega_ij - Phi_q1q * Phi_q1q1
     * - [3] Phi_q1q1
     * - [4] omega_jj - Phi_q1q^2
     * - [5] constrained diagonal at x = 0
     */
    std::array<double, 6> constants_{};

    /**
     * Work vectors for rank-2 Cholesky update.
     *
     * A symmetric rank-2 update  A + vf1*vf2' + vf2*vf1'  is decomposed
     * into two rank-1 updates via  u1 = (vf1+vf2)/sqrt(2),
     * u2 = (vf1-vf2)/sqrt(2).
     */
    arma::vec v1_ = {0, -1};
    arma::vec v2_ = {0, 0};
    arma::vec vf1_ = arma::zeros<arma::vec>(p_);
    arma::vec vf2_ = arma::zeros<arma::vec>(p_);
    arma::vec u1_ = arma::zeros<arma::vec>(p_);
    arma::vec u2_ = arma::zeros<arma::vec>(p_);

    /**
     * Propose a new off-diagonal precision entry via a normal perturbation
     * on an unconstrained reparameterization. Accepts or rejects with a
     * Metropolis ratio using the Gaussian likelihood and Cauchy prior.
     *
     * @param i          Row index (i < j)
     * @param j          Column index
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void update_edge_parameter(size_t i, size_t j, int iteration);

    /**
     * Propose a new diagonal precision entry on the log scale.
     * Accepts or rejects with a Metropolis ratio using the Gaussian
     * likelihood, a Gamma(1,1) prior, and a Jacobian correction.
     *
     * @param i          Diagonal index
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void update_diagonal_parameter(size_t i, int iteration);

    /**
     * Metropolis-Hastings add-delete move for an edge indicator.
     *
     * If the edge is on, proposes deletion; if off, proposes a new value
     * from a scaled normal. Acceptance combines the likelihood ratio,
     * Bernoulli prior odds, Cauchy slab, and proposal density.
     *
     * @param i  Row index (i < j)
     * @param j  Column index
     */
    void update_edge_indicator_parameter_pair(size_t i, size_t j);

    /**
     * Precompute reparameterization constants for the (i, j) element.
     *
     * Derives six values from the cofactor structure of the inverse
     * precision matrix that allow off-diagonal proposals on an
     * unconstrained scale while deterministically satisfying the
     * positive-definiteness constraint on the diagonal.
     *
     * @param i  Row index
     * @param j  Column index
     */
    void get_constants(size_t i, size_t j);



    /**
     * Return the diagonal value omega_jj required to keep the precision
     * matrix positive definite after changing the off-diagonal element to x.
     *
     * @param x  Proposed off-diagonal value omega_ij
     * @return   Constrained diagonal value omega_jj
     */
    double constrained_diagonal(const double x) const;

    /**
     * Full Gaussian log-likelihood: n/2 * (log|Omega| - tr(Omega S) / n).
     *
     * @param omega  Precision matrix
     * @param phi    Upper-triangular Cholesky factor of omega
     */
    double log_density_impl(const arma::mat& omega, const arma::mat& phi) const;

    /**
     * Log-likelihood ratio for a proposed off-diagonal element change,
     * computed via the matrix-determinant lemma (rank-2 update).
     *
     * @param i  Row index of the changed element
     * @param j  Column index of the changed element
     */
    double log_density_impl_edge(size_t i, size_t j) const;

    /**
     * Log-likelihood ratio for a proposed diagonal element change,
     * computed via the matrix-determinant lemma (rank-1 update).
     *
     * @param j  Index of the changed diagonal element
     */
    double log_density_impl_diag(size_t j) const;



    /**
     * Update the Cholesky factor after changing an off-diagonal element.
     *
     * Decomposes the rank-2 change into two rank-1 updates and
     * recomputes the inverse Cholesky factor and covariance matrix.
     *
     * @param omega_ij_old  Previous value of omega(i,j)
     * @param omega_jj_old  Previous value of omega(j,j)
     * @param i             Row index
     * @param j             Column index
     */
    void cholesky_update_after_edge(double omega_ij_old, double omega_jj_old, size_t i, size_t j);

    /**
     * Update the Cholesky factor after changing a diagonal element.
     *
     * Applies a rank-1 update and recomputes the inverse Cholesky
     * factor and covariance matrix.
     *
     * @param omega_ii_old  Previous value of omega(i,i)
     * @param i             Diagonal index
     */
    void cholesky_update_after_diag(double omega_ii_old, size_t i);
};

/**
 * Construct a GGMModel from an R list.
 *
 * Dispatches to the sufficient-statistics constructor (when the list
 * contains `n` and `suf_stat`) or the raw-data constructor (when the
 * list contains `X`).
 *
 * @param inputFromR              R list with data (either `X` or `n` + `suf_stat`)
 * @param inclusion_probability   Prior inclusion probabilities for each edge
 * @param initial_edge_indicators Initial edge inclusion indicators
 * @param edge_selection          Enable edge selection (spike-and-slab)
 * @param pairwise_scale          Scale parameter of Cauchy slab prior
 * @param na_impute               Retain observations for missing-data imputation
 * @return Fully constructed GGMModel
 */
GGMModel createGGMModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection,
    const double pairwise_scale,
    const bool na_impute = false
);
