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

    // Test helpers need access to private likelihood functions
    friend Rcpp::List test_mixed_mrf_likelihoods(
        const arma::imat&, const arma::mat&, const arma::ivec&,
        const arma::uvec&, const arma::ivec&, const arma::mat&,
        const arma::imat&, bool, const std::string&, const arma::vec&, int
    );

    // Test helper for rank-1 Cholesky correctness (T28/T29)
    friend Rcpp::List test_mixed_mrf_cholesky(
        const arma::imat&, const arma::mat&, const arma::ivec&,
        const arma::uvec&, const arma::ivec&, const arma::mat&,
        const arma::imat&, const arma::vec&, int, int, int
    );

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

    /** @return false (MixedMRFModel is Metropolis-only, no gradient). */
    bool has_gradient() const override { return false; }
    /** @return true (supports adaptive Metropolis via Robbins-Monro). */
    bool has_adaptive_metropolis() const override { return true; }
    /** @return true when edge selection is enabled. */
    bool has_edge_selection() const override { return edge_selection_; }
    /** @return false (missing data not yet supported). */
    bool has_missing_data() const override { return false; }

    // =========================================================================
    // Core sampling methods
    // =========================================================================

    /**
     * Perform one full Metropolis sweep over all parameter groups.
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void do_one_metropolis_step(int iteration = -1) override;

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
     * Dimensionality of the active parameter space.
     * When edge selection is active, excludes parameters for inactive edges.
     */
    size_t parameter_dimension() const override;

    /**
     * Full parameter dimension (all parameters, regardless of edge state).
     * Used for fixed-size sample storage.
     */
    size_t full_parameter_dimension() const override;

    /** Get active parameters as a flat vector. */
    arma::vec get_vectorized_parameters() const override;

    /** Get all parameters as a flat vector (inactive edges are 0). */
    arma::vec get_full_vectorized_parameters() const override;

    /** Set parameters from a flat vector. */
    void set_vectorized_parameters(const arma::vec& params) override;

    /** Get vectorized edge indicators (Gxx upper-tri, Gyy upper-tri, Gxy full). */
    arma::ivec get_vectorized_indicator_parameters() override;

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
    // Missing data (not yet supported)
    // =========================================================================

    /** No-op: missing data not supported in this PR. */
    void impute_missing() override {}

private:

    // =========================================================================
    // Counts and dimensions
    // =========================================================================

    size_t n_;                          // Number of observations
    size_t p_;                          // Number of discrete variables
    size_t q_;                          // Number of continuous variables
    size_t num_main_;                   // Total main-effect params (Σ C_s for ord + 2 per BC)
    size_t num_pairwise_xx_;            // p(p-1)/2
    size_t num_pairwise_yy_;            // q(q-1)/2
    size_t num_cross_;                  // p * q

    // =========================================================================
    // Data
    // =========================================================================

    arma::imat discrete_observations_;   // Discrete observations (n × p)
                                        //   BC columns centered at baseline_category_ in ctor.
    arma::mat discrete_observations_dbl_; // Double version (post-centering)
    arma::mat continuous_observations_;  // Continuous observations (n × q)
    arma::ivec num_categories_;         // Categories per discrete variable (p-vector)
    int max_cats_;                      // max(num_categories)
    arma::uvec is_ordinal_variable_;    // 1 = ordinal, 0 = Blume-Capel (p-vector)
    arma::ivec baseline_category_;      // Reference category per discrete variable (p-vector)

    // =========================================================================
    // Sufficient statistics
    // =========================================================================

    arma::imat counts_per_category_;    // (max_cats+1) × p  category counts (ordinal only)
    arma::imat blume_capel_stats_;      // 2 × p  linear/quadratic sums (BC only)

    // =========================================================================
    // Parameters
    // =========================================================================

    arma::mat mux_;                     // p × max_cats  main effects
                                        //   Ordinal: mux_(s,c) = threshold for category c+1;
                                        //     category 0 is reference (fixed at 0).
                                        //   Blume-Capel: mux_(s,0) = linear α_s,
                                        //     mux_(s,1) = quadratic β_s.
    arma::vec muy_;                     // q-vector  continuous means
    arma::mat Kxx_;                     // p × p  discrete interactions (symmetric, zero diag)
    arma::mat Kyy_;                     // q × q  SPD precision matrix
    arma::mat Kxy_;                     // p × q  cross-type interactions

    // =========================================================================
    // Edge indicators
    // =========================================================================

    // Combined (p+q) × (p+q) indicator matrix:
    //   Gxx block : rows [0,p),    cols [0,p)    — symmetric, zero diag
    //   Gyy block : rows [p,p+q),  cols [p,p+q)  — symmetric, zero diag
    //   Gxy block : rows [0,p),    cols [p,p+q)   — full p×q rectangle
    arma::imat edge_indicators_;
    arma::mat inclusion_probability_;
    bool edge_selection_;
    bool edge_selection_active_;

    // =========================================================================
    // Priors
    // =========================================================================

    double main_alpha_;                 // Beta prior α for main effects
    double main_beta_;                  // Beta prior β for main effects
    double pairwise_scale_;             // Cauchy scale for interaction priors

    // =========================================================================
    // Proposal SDs (Robbins-Monro adapted)
    // =========================================================================

    arma::mat prop_sd_mux_;             // p × max_cats
    arma::vec prop_sd_muy_;             // q-vector
    arma::mat prop_sd_Kxx_;             // p × p
    arma::mat prop_sd_Kyy_;             // q × q
    arma::mat prop_sd_Kxy_;             // p × q

    // =========================================================================
    // Cached quantities
    // =========================================================================

    arma::mat Kyy_chol_;                // q × q  upper Cholesky of Kyy (Kyy = R'R)
    arma::mat inv_cholesky_yy_;         // q × q  R^{-1} (upper triangular)
    arma::mat covariance_yy_;           // q × q  Kyy^{-1} = R^{-1} R^{-T}
    double Kyy_log_det_;                // log|Kyy|
    arma::mat Theta_;                   // p × p  Kxx + 2 Kxy covariance_yy_ Kxy' (marginal PL only)
    arma::mat conditional_mean_;        // n × q  μ_y' + 2 discrete_obs Kxy covariance_yy_

    // Rank-1 Cholesky update workspace
    std::array<double, 6> kyy_constants_{};  // reparameterization constants
    arma::mat precision_yy_proposal_;        // q × q scratch for proposed Kyy
    arma::vec kyy_v1_ = {0, -1};             // rank-2 decomposition helpers
    arma::vec kyy_v2_ = {0, 0};
    arma::vec kyy_vf1_;                      // q-vectors, zeroed between uses
    arma::vec kyy_vf2_;
    arma::vec kyy_u1_;
    arma::vec kyy_u2_;

    // =========================================================================
    // Configuration
    // =========================================================================

    bool use_marginal_pl_;              // true = marginal, false = conditional

    // =========================================================================
    // RNG and edge-update order
    // =========================================================================

    SafeRNG rng_;
    arma::uvec edge_order_xx_;          // Shuffled xx-edge pair indices
    arma::uvec edge_order_yy_;          // Shuffled yy-edge pair indices
    arma::uvec edge_order_xy_;          // Shuffled xy-edge pair indices

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
    // Likelihood functions (implemented in mixed_mrf_likelihoods.cpp)
    // =========================================================================

    /** Conditional OMRF pseudolikelihood for discrete variable s, summed over all n. */
    double log_conditional_omrf(int s) const;

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
    void update_main_effect(int s, int c);

    /** Update one continuous mean: muy_(j). */
    void update_continuous_mean(int j);

    /** Update one discrete interaction: Kxx_(i, j). Symmetric. */
    void update_Kxx(int i, int j);

    /** Update one off-diagonal precision element: Kyy_(i, j). Cholesky-based. */
    void update_Kyy_offdiag(int i, int j);

    /** Update one diagonal precision element: Kyy_(i, i). Log-scale Cholesky. */
    void update_Kyy_diag(int i);

    /** Update one cross interaction: Kxy_(i, j). */
    void update_Kxy(int i, int j);

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
