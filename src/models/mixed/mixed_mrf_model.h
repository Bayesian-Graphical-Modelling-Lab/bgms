#pragma once

#include <array>
#include <memory>
#include "models/base_model.h"
#include "math/cholesky_helpers.h"
#include "models/ggm/cholupdate.h"
#include "rng/rng_utils.h"

/**
 * MixedMRFModel - Mixed Markov Random Field Model
 *
 * Joint model for p discrete (ordinal or Blume-Capel) variables x and
 * q continuous variables y.  The joint density is:
 *
 *   log f(x, y) ∝ Σ_s μ_{x,s}(x_s) + x' Kxx x
 *                  - ½ (y - μ_y)' Kyy (y - μ_y) + 2 x' Kxy y
 *
 * Supports both conditional and marginal pseudo-likelihood, with and
 * without edge selection via spike-and-slab priors.
 *
 * Discrete variables are either ordinal (free category thresholds, category
 * 0 as reference) or Blume-Capel (linear α + quadratic β, user-specified
 * reference).  Blume-Capel observations are centered at their baseline
 * category in the constructor, matching the OMRFModel convention.
 *
 * Inherits from BaseModel for compatibility with the generic MCMC framework
 * (ChainRunner, MetropolisSampler, WarmupSchedule).
 */
class MixedMRFModel : public BaseModel {
public:

    // =========================================================================
    // Construction
    // =========================================================================

    /**
     * Construct from raw observations.
     *
     * @param discrete_observations   Integer matrix of discrete observations (n × p, 0-based)
     * @param continuous_observations  Continuous observations (n × q)
     * @param num_categories       Number of categories per discrete variable (p-vector)
     * @param is_ordinal_variable  1 = ordinal, 0 = Blume-Capel (p-vector)
     * @param baseline_category    Reference category per discrete variable (p-vector)
     * @param inclusion_probability Prior inclusion probabilities ((p+q) × (p+q))
     * @param initial_edge_indicators Initial edge inclusion matrix ((p+q) × (p+q))
     * @param edge_selection       Enable edge selection (spike-and-slab)
     * @param pseudolikelihood     "conditional" or "marginal"
     * @param main_alpha           Beta prior hyperparameter α for main effects
     * @param main_beta            Beta prior hyperparameter β for main effects
     * @param pairwise_scale       Scale parameter of Cauchy prior on interactions
     * @param seed                 RNG seed for reproducibility
     */
    MixedMRFModel(
        const arma::imat& discrete_observations,
        const arma::mat& continuous_observations,
        const arma::ivec& num_categories,
        const arma::uvec& is_ordinal_variable,
        const arma::ivec& baseline_category,
        const arma::mat& inclusion_probability,
        const arma::imat& initial_edge_indicators,
        bool edge_selection,
        const std::string& pseudolikelihood,
        double main_alpha = 1.0,
        double main_beta = 1.0,
        double pairwise_scale = 2.5,
        int seed = 1
    );

    /** Copy constructor for cloning (required for parallel chains). */
    MixedMRFModel(const MixedMRFModel& other);

    // =========================================================================
    // Capability queries
    // =========================================================================

    /** @return true (MixedMRFModel supports NUTS gradient for the unconstrained block). */
    bool has_gradient() const override { return true; }
    /** @return true (supports adaptive Metropolis via Robbins-Monro). */
    bool has_adaptive_metropolis() const override { return true; }
    /** @return true when edge selection is enabled. */
    bool has_edge_selection() const override { return edge_selection_; }
    /** @return true when missing-data imputation is active. */
    bool has_missing_data() const override { return has_missing_; }

    // =========================================================================
    // Core sampling methods
    // =========================================================================

    /**
     * Compute gradient of the log pseudo-posterior (NUTS block only).
     * @param parameters  NUTS-dimension parameter vector
     * @return Gradient vector (same size as parameters)
     */
    arma::vec gradient(const arma::vec& parameters) override;

    /**
     * Combined log pseudo-posterior and gradient evaluation.
     * @param parameters  NUTS-dimension parameter vector
     * @return Pair of (log-pseudo-posterior, gradient)
     */
    std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) override;

    /**
     * Perform one full Metropolis sweep over all parameter groups.
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void do_one_metropolis_step(int iteration = -1) override;

    /**
     * Update only Kyy parameters via Metropolis (used by hybrid NUTS+MH).
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void do_kyy_metropolis_step(int iteration = -1);

    /**
     * Initialize Metropolis adaptation controllers for proposal-SD tuning.
     * Called before warmup begins.
     */
    void init_metropolis_adaptation(const WarmupSchedule& schedule) override;

    /**
     * Tune proposal SDs via Robbins-Monro (Stage 3b).
     * Called every iteration; checks schedule internally.
     */
    void tune_proposal_sd(int iteration, const WarmupSchedule& schedule) override;

    /**
     * Shuffle edge update order at the start of each iteration.
     * Advances the RNG state consistently even when edge selection is off.
     */
    void prepare_iteration() override;

    // =========================================================================
    // Edge selection
    // =========================================================================

    /** Perform one sweep of reversible-jump edge birth/death moves. */
    void update_edge_indicators() override;

    /** Initialize a random graph structure for starting edge selection. */
    void initialize_graph() override;

    /**
     * Enable or disable edge-selection proposals.
     * @param active  true to enable edge birth/death moves
     */
    void set_edge_selection_active(bool active) override {
        edge_selection_active_ = active;
    }

    // =========================================================================
    // Parameter vectorization
    // =========================================================================

    /**
     * Dimensionality of the active NUTS parameter space (excludes Kyy).
     * When edge selection is active, excludes parameters for inactive edges.
     */
    size_t parameter_dimension() const override;

    /**
     * Full NUTS-block dimension (all NUTS params, regardless of edge state).
     * Excludes Kyy. Used for mass-matrix sizing and adaptation.
     */
    size_t full_parameter_dimension() const override;

    /**
     * Storage dimension (all parameters including Kyy, regardless of edge state).
     * Used for fixed-size sample storage.
     */
    size_t storage_dimension() const override;

    /** Get active NUTS parameters as a flat vector (excludes Kyy). */
    arma::vec get_vectorized_parameters() const override;

    /** Get all NUTS parameters (inactive edges are 0, excludes Kyy). */
    arma::vec get_full_vectorized_parameters() const override;

    /** Get all parameters including Kyy for sample storage. */
    arma::vec get_storage_vectorized_parameters() const override;

    /** Set NUTS parameters from a flat vector (does not touch Kyy). */
    void set_vectorized_parameters(const arma::vec& params) override;

    /** Get vectorized edge indicators (Gxx upper-tri, Gyy upper-tri, Gxy full). */
    arma::ivec get_vectorized_indicator_parameters() override;

    /** Get active subset of inverse mass diagonal (NUTS params only, excludes Kyy). */
    arma::vec get_active_inv_mass() const override;

    // =========================================================================
    // Infrastructure
    // =========================================================================

    /** Set random seed for reproducibility. */
    void set_seed(int seed) override;

    /** Clone the model for parallel execution. */
    std::unique_ptr<BaseModel> clone() const override;

    /** @return Reference to the model's random number generator. */
    SafeRNG& get_rng() override { return rng_; }

    /** @return Current edge-indicator matrix ((p+q) × (p+q)). */
    const arma::imat& get_edge_indicators() const override {
        return edge_indicators_;
    }

    /** @return Mutable reference to the prior inclusion-probability matrix. */
    arma::mat& get_inclusion_probability() override {
        return inclusion_probability_;
    }

    /** @return Total number of variables (p + q). */
    int get_num_variables() const override {
        return static_cast<int>(p_ + q_);
    }

    /**
     * Number of unique off-diagonal pairs in the (p+q) × (p+q) indicator
     * matrix: p(p-1)/2 + q(q-1)/2 + p*q.
     */
    int get_num_pairwise() const override {
        return static_cast<int>(num_pairwise_xx_ + num_pairwise_yy_ + num_cross_);
    }

    // =========================================================================
    // Missing data
    // =========================================================================

    /** Impute missing entries from full-conditional distributions. */
    void impute_missing() override;

    /**
     * Register missing-data locations for discrete and continuous sub-matrices.
     *
     * @param missing_discrete  M_d x 2 matrix of 0-based (row, col) indices into discrete_observations_
     * @param missing_continuous M_c x 2 matrix of 0-based (row, col) indices into continuous_observations_
     */
    void set_missing_data(const arma::imat& missing_discrete,
                          const arma::imat& missing_continuous);

private:

    // =========================================================================
    // Counts and dimensions
    // =========================================================================

    size_t n_;                          ///< Number of observations
    size_t p_;                          ///< Number of discrete variables
    size_t q_;                          ///< Number of continuous variables
    size_t num_main_;                   ///< Total main-effect params (sum C_s for ord + 2 per BC)
    size_t num_pairwise_xx_;            ///< p(p-1)/2
    size_t num_pairwise_yy_;            ///< q(q-1)/2
    size_t num_cross_;                  ///< p * q

    // =========================================================================
    // Data
    // =========================================================================

    arma::imat discrete_observations_;   ///< Discrete observations (n x p), BC columns centered
    arma::mat discrete_observations_dbl_; ///< Double version (post-centering)
    arma::mat continuous_observations_;  ///< Continuous observations (n x q)
    arma::ivec num_categories_;         ///< Categories per discrete variable (p-vector)
    int max_cats_;                      ///< max(num_categories)
    arma::uvec is_ordinal_variable_;    ///< 1 = ordinal, 0 = Blume-Capel (p-vector)
    arma::ivec baseline_category_;      ///< Reference category per discrete variable (p-vector)

    // =========================================================================
    // Missing data
    // =========================================================================

    arma::imat missing_index_discrete_;   ///< M_d x 2 (row, col) for missing discrete entries
    arma::imat missing_index_continuous_; ///< M_c x 2 (row, col) for missing continuous entries
    bool has_missing_ = false;            ///< Whether imputation is active

    // =========================================================================
    // Sufficient statistics
    // =========================================================================

    arma::imat counts_per_category_;    ///< (max_cats+1) x p category counts (ordinal only)
    arma::imat blume_capel_stats_;      ///< 2 x p linear/quadratic sums (BC only)

    // =========================================================================
    // Parameters
    // =========================================================================

    arma::mat mux_;                     ///< p x max_cats main effects (thresholds or alpha/beta)
    arma::vec muy_;                     ///< q-vector continuous means
    arma::mat Kxx_;                     ///< p x p discrete interactions (symmetric, zero diag)
    arma::mat Kyy_;                     ///< q x q SPD precision matrix
    arma::mat Kxy_;                     ///< p x q cross-type interactions

    // =========================================================================
    // Edge indicators
    // =========================================================================

    /// Combined (p+q) x (p+q) indicator matrix.
    /// Gxx block: rows [0,p), cols [0,p) -- symmetric, zero diag.
    /// Gyy block: rows [p,p+q), cols [p,p+q) -- symmetric, zero diag.
    /// Gxy block: rows [0,p), cols [p,p+q) -- full p x q rectangle.
    arma::imat edge_indicators_;
    arma::mat inclusion_probability_;   ///< Prior inclusion probabilities
    bool edge_selection_;               ///< Enable edge selection
    bool edge_selection_active_;        ///< Currently in edge selection phase

    // =========================================================================
    // Priors
    // =========================================================================

    double main_alpha_;                 ///< Beta prior alpha for main effects
    double main_beta_;                  ///< Beta prior beta for main effects
    double pairwise_scale_;             ///< Cauchy scale for interaction priors

    // =========================================================================
    // Proposal SDs (Robbins-Monro adapted)
    // =========================================================================

    arma::mat prop_sd_mux_;             ///< p x max_cats
    arma::vec prop_sd_muy_;             ///< q-vector
    arma::mat prop_sd_Kxx_;             ///< p x p
    arma::mat prop_sd_Kyy_;             ///< q x q
    arma::mat prop_sd_Kxy_;             ///< p x q
    int total_warmup_ = 0;              ///< Stored by init_metropolis_adaptation

    // =========================================================================
    // Cached quantities
    // =========================================================================

    arma::mat Kyy_chol_;                ///< q x q upper Cholesky of Kyy (Kyy = R'R)
    arma::mat inv_cholesky_yy_;         ///< q x q R^{-1} (upper triangular)
    arma::mat covariance_yy_;           ///< q x q Kyy^{-1} = R^{-1} R^{-T}
    double Kyy_log_det_;                ///< log|Kyy|
    arma::mat Theta_;                   ///< p x p Kxx + 2 Kxy covariance_yy_ Kxy' (marginal PL)
    arma::mat conditional_mean_;        ///< n x q mu_y' + 2 discrete_obs Kxy covariance_yy_

    // Rank-1 Cholesky update workspace
    std::array<double, 6> kyy_constants_{};  ///< Reparameterization constants
    arma::mat precision_yy_proposal_;        ///< q x q scratch for proposed Kyy
    arma::vec kyy_v1_ = {0, -1};             ///< Rank-2 decomposition helper 1
    arma::vec kyy_v2_ = {0, 0};              ///< Rank-2 decomposition helper 2
    arma::vec kyy_vf1_;                      ///< q-vector, zeroed between uses
    arma::vec kyy_vf2_;                      ///< q-vector, zeroed between uses
    arma::vec kyy_u1_;                       ///< q-vector workspace
    arma::vec kyy_u2_;                       ///< q-vector workspace

    // =========================================================================
    // Gradient cache (populated by ensure_gradient_cache)
    // =========================================================================

    arma::mat discrete_observations_dbl_t_; ///< p x n transpose (BLAS gradient)
    arma::vec grad_obs_cache_;          ///< Cached observed-data gradient component
    arma::imat kxx_index_cache_;        ///< p x p map from (i,j) to gradient index
    arma::imat kxy_index_cache_;        ///< p x q map from (i,j) to gradient index
    int muy_grad_offset_ = 0;           ///< Offset of muy block in gradient vector
    bool gradient_cache_valid_ = false; ///< Whether gradient cache is current

    // =========================================================================
    // Configuration
    // =========================================================================

    bool use_marginal_pl_;              ///< true = marginal, false = conditional

    // =========================================================================
    // RNG and edge-update order
    // =========================================================================

    SafeRNG rng_;                       ///< Per-chain random number generator
    arma::uvec edge_order_xx_;          ///< Shuffled xx-edge pair indices
    arma::uvec edge_order_yy_;          ///< Shuffled yy-edge pair indices
    arma::uvec edge_order_xy_;          ///< Shuffled xy-edge pair indices

    // =========================================================================
    // Private helpers
    // =========================================================================

    /** Count total main-effect parameters across all discrete variables. */
    size_t count_num_main_effects() const;

    /** Compute category counts and BC sufficient statistics from discrete_observations_. */
    void compute_sufficient_statistics();

    /** Recompute conditional_mean_ from muy_, Kxy_, covariance_yy_. */
    void recompute_conditional_mean();

    /** Recompute Kyy_chol_, inv_cholesky_yy_, covariance_yy_, Kyy_log_det_ from Kyy_. */
    void recompute_Kyy_decomposition();

    /** Recompute Theta_ from Kxx_, Kxy_, covariance_yy_ (marginal PL only). */
    void recompute_Theta();

    // =========================================================================
    // Gradient helpers (implemented in mixed_mrf_gradient.cpp)
    // =========================================================================

    /** Rebuild gradient index maps after edge-indicator changes. */
    void ensure_gradient_cache();

    /** Mark gradient cache as stale (call after edge-indicator changes). */
    void invalidate_gradient_cache();

    /** Unpack NUTS-vector into temporary parameter matrices (no model mutation). */
    void unvectorize_nuts_to_temps(
        const arma::vec& params,
        arma::mat& temp_mux,
        arma::mat& temp_Kxx,
        arma::vec& temp_muy,
        arma::mat& temp_Kxy
    ) const;

    // =========================================================================
    // Likelihood functions (implemented in mixed_mrf_likelihoods.cpp)
    // =========================================================================

    /** Conditional OMRF pseudolikelihood for discrete variable s, summed over all n. */
    double log_conditional_omrf(int s) const;

    /** Marginal OMRF pseudolikelihood for discrete variable s, using Theta_. */
    double log_marginal_omrf(int s) const;

    /** Conditional GGM log-likelihood: log f(y | x), using cached decomposition. */
    double log_conditional_ggm() const;

    // =========================================================================
    // MH update functions (implemented in mixed_mrf_metropolis.cpp)
    // =========================================================================

    // --- Rank-1 Kyy proposal helpers (permutation-free) ---

    // Extract reparameterization constants for the (i,j) off-diagonal Kyy update.
    // Populates kyy_constants_[0..5] from Kyy_chol_ and covariance_yy_.
    void get_kyy_constants(int i, int j);

    // Constrained diagonal value for a proposed off-diagonal Kyy element.
    double kyy_constrained_diagonal(double x) const;

    // Log-likelihood ratio for a proposed off-diagonal Kyy change (rank-2).
    // Assumes precision_yy_proposal_ is already filled by the caller.
    double log_ggm_ratio_edge(int i, int j) const;

    // Log-likelihood ratio for a proposed diagonal Kyy change (rank-1).
    // Assumes precision_yy_proposal_ is already filled by the caller.
    double log_ggm_ratio_diag(int i) const;

    // Rank-1 Cholesky update after accepting an off-diagonal Kyy change.
    void cholesky_update_after_kyy_edge(double old_ij, double old_jj, int i, int j);

    // Rank-1 Cholesky update after accepting a diagonal Kyy change.
    void cholesky_update_after_kyy_diag(double old_ii, int i);

    // --- Parameter update sweeps ---

    /** Update one main-effect: mux_(s, c). Ordinal threshold or BC α/β. */
    void update_main_effect(int s, int c, int iteration);

    /** Update one continuous mean: muy_(j). */
    void update_continuous_mean(int j, int iteration);

    /** Update one discrete interaction: Kxx_(i, j). Symmetric. */
    void update_Kxx(int i, int j, int iteration);

    /** Update one off-diagonal precision element: Kyy_(i, j). Cholesky-based. */
    void update_Kyy_offdiag(int i, int j, int iteration);

    /** Update one diagonal precision element: Kyy_(i, i). Log-scale Cholesky. */
    void update_Kyy_diag(int i, int iteration);

    /** Update one cross interaction: Kxy_(i, j). */
    void update_Kxy(int i, int j, int iteration);

    // --- Edge-indicator update sweeps (Phase D) ---

    /** Reversible-jump birth/death for one Kxx edge (discrete-discrete). */
    void update_edge_indicator_Kxx(int i, int j);

    /** Reversible-jump birth/death for one Kyy edge (continuous-continuous). */
    void update_edge_indicator_Kyy(int i, int j);

    /** Reversible-jump birth/death for one Kxy edge (cross-type). */
    void update_edge_indicator_Kxy(int i, int j);

    // =========================================================================
    // Edge-indicator accessor helpers
    // =========================================================================

    int gxx(int i, int j) const { return edge_indicators_(i, j); }
    int gyy(int i, int j) const { return edge_indicators_(p_ + i, p_ + j); }
    int gxy(int i, int j) const { return edge_indicators_(i, p_ + j); }

    void set_gxx(int i, int j, int val) {
        edge_indicators_(i, j) = val;
        edge_indicators_(j, i) = val;
    }
    void set_gyy(int i, int j, int val) {
        edge_indicators_(p_ + i, p_ + j) = val;
        edge_indicators_(p_ + j, p_ + i) = val;
    }
    void set_gxy(int i, int j, int val) {
        edge_indicators_(i, p_ + j) = val;
    }
};
