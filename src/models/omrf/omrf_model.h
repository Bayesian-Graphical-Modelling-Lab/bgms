#pragma once

#include <memory>
#include <functional>
#include "models/base_model.h"
#include "mcmc/adaptation.h"
#include "rng/rng_utils.h"
#include "mcmc/sampler_result.h"
#include "utils/common_helpers.h"

/**
 * OMRFModel - Ordinal Markov Random Field Model
 *
 * A class-based implementation of the OMRF model for Bayesian inference on
 * ordinal and Blume-Capel variables. This class encapsulates:
 *   - Parameter storage (main effects, pairwise effects, edge indicators)
 *   - Sufficient statistics computation
 *   - Log-pseudoposterior and gradient evaluations
 *   - Adaptive Metropolis-Hastings updates for individual parameters
 *   - NUTS/HMC updates for joint parameter sampling
 *   - Edge selection (spike-and-slab) with asymmetric proposals
 *
 * Inherits from BaseModel for compatibility with the generic MCMC framework.
 */
class OMRFModel : public BaseModel {
public:

    /**
     * Constructor from raw observations
     *
     * @param observations        Integer matrix of categorical observations (persons × variables)
     * @param num_categories      Number of categories per variable
     * @param inclusion_probability Prior inclusion probabilities for edges
     * @param initial_edge_indicators Initial edge inclusion matrix
     * @param is_ordinal_variable Indicator (1 = ordinal, 0 = Blume-Capel)
     * @param baseline_category   Reference categories for Blume-Capel variables
     * @param main_alpha          Beta prior hyperparameter α for main effects
     * @param main_beta           Beta prior hyperparameter β for main effects
     * @param pairwise_scale      Scale parameter of Cauchy prior on interactions
     * @param edge_selection      Enable edge selection (spike-and-slab)
     */
    OMRFModel(
        const arma::imat& observations,
        const arma::ivec& num_categories,
        const arma::mat& inclusion_probability,
        const arma::imat& initial_edge_indicators,
        const arma::uvec& is_ordinal_variable,
        const arma::ivec& baseline_category,
        double main_alpha = 1.0,
        double main_beta = 1.0,
        double pairwise_scale = 2.5,
        bool edge_selection = true
    );

    /**
     * Copy constructor for cloning (required for parallel chains)
     */
    OMRFModel(const OMRFModel& other);

    // =========================================================================
    // BaseModel interface implementation
    // =========================================================================

    bool has_gradient()    const override { return true; }
    bool has_adaptive_mh() const override { return true; }
    bool has_edge_selection() const override { return edge_selection_; }
    bool has_missing_data() const override { return has_missing_; }

    /**
     * Compute log-pseudoposterior for given parameter vector
     */
    double logp(const arma::vec& parameters) override;

    /**
     * Compute gradient of log-pseudoposterior
     */
    arma::vec gradient(const arma::vec& parameters) override;

    /**
     * Combined log-posterior and gradient evaluation (more efficient)
     */
    std::pair<double, arma::vec> logp_and_gradient(const arma::vec& parameters) override;

    /**
     * Perform one adaptive MH step (updates all parameters)
     * @param iteration  Current iteration (for Robbins-Monro adaptation)
     */
    void do_one_mh_step(int iteration = -1) override;

    /**
     * Initialize RWM adaptation controllers for proposal-SD tuning
     * Must be called before warmup begins (e.g., by MHSampler on first step)
     */
    void init_mh_adaptation(const WarmupSchedule& schedule) override;

    /**
     * Stage 3b: tune pairwise proposal SDs via Robbins-Monro
     * Called every iteration from the runner; checks schedule internally
     */
    void tune_proposal_sd(int iteration, const WarmupSchedule& schedule) override;

    /**
     * Return dimensionality of active parameter space
     */
    size_t parameter_dimension() const override;

    /**
     * Set random seed for reproducibility
     */
    void set_seed(int seed) override;

    /**
     * Get vectorized parameters (main effects + active pairwise effects)
     */
    arma::vec get_vectorized_parameters() const override;

    /**
     * Set parameters from vectorized form
     */
    void set_vectorized_parameters(const arma::vec& parameters) override;

    /**
     * Get vectorized edge indicators
     */
    arma::ivec get_vectorized_indicator_parameters() override;

    /**
     * Clone the model for parallel execution
     */
    std::unique_ptr<BaseModel> clone() const override;

    /**
     * Get RNG for samplers
     */
    SafeRNG& get_rng() override { return rng_; }

    // =========================================================================
    // OMRF-specific methods
    // =========================================================================



    /**
     * Shuffle edge update order at the start of each iteration.
     * The shuffled order is stored in shuffled_edge_order_ for use by
     * update_edge_indicators(). Called unconditionally to advance the
     * RNG state consistently.
     */
    void prepare_iteration() override;

    /**
     * Update edge indicators via Metropolis-Hastings
     */
    void update_edge_indicators() override;

    /**
     * Initialize random graph structure (for starting edge selection)
     */
    void initialize_graph() override;

    /**
     * Impute missing values (if any)
     */
    void impute_missing() override;

    /**
     * Set missing data information
     */
    void set_missing_data(const arma::imat& missing_index);

    // =========================================================================
    // Accessors
    // =========================================================================

    const arma::mat& get_main_effects() const { return main_effects_; }
    const arma::mat& get_pairwise_effects() const { return pairwise_effects_; }
    const arma::imat& get_edge_indicators() const override { return edge_indicators_; }
    arma::mat& get_inclusion_probability() override { return inclusion_probability_; }
    const arma::mat& get_residual_matrix() const { return residual_matrix_; }

    void set_main_effects(const arma::mat& main_effects) { main_effects_ = main_effects; }
    void set_pairwise_effects(const arma::mat& pairwise_effects);
    void set_edge_indicators(const arma::imat& edge_indicators) { edge_indicators_ = edge_indicators; }

    int get_num_variables() const override { return static_cast<int>(p_); }
    int get_num_pairwise() const override { return static_cast<int>(num_pairwise_); }
    size_t num_variables() const { return p_; }
    size_t num_observations() const { return n_; }
    size_t num_main_effects() const { return num_main_; }
    size_t num_pairwise_effects() const { return num_pairwise_; }

    // Shorthand accessors (for interface compatibility)
    size_t get_p() const { return p_; }
    size_t get_n() const { return n_; }

    // Adaptation control
    void set_step_size(double step_size) override { step_size_ = step_size; }
    double get_step_size() const override { return step_size_; }
    void set_inv_mass(const arma::vec& inv_mass) override { inv_mass_ = inv_mass; }
    const arma::vec& get_inv_mass() const override { return inv_mass_; }

    /**
     * Get full dimension (main + ALL pairwise, regardless of edge indicators)
     * Used for fixed-size sample storage
     */
    size_t full_parameter_dimension() const override { return num_main_ + num_pairwise_; }

    /**
     * Get all parameters in a fixed-size vector (inactive edges are 0)
     * Used for sample storage to avoid dimension changes
     */
    arma::vec get_full_vectorized_parameters() const override;

    // Proposal SD access (for external adaptation)
    arma::mat& get_proposal_sd_main() { return proposal_sd_main_; }
    arma::mat& get_proposal_sd_pairwise() { return proposal_sd_pairwise_; }

    // Pairwise scaling factors
    void set_pairwise_scaling_factors(const arma::mat& sf) { pairwise_scaling_factors_ = sf; }

    // Control edge selection phase
    void set_edge_selection_active(bool active) override { edge_selection_active_ = active; }
    bool is_edge_selection_active() const { return edge_selection_active_; }

private:
    // =========================================================================
    // Data members
    // =========================================================================

    // Data
    size_t n_;                          // Number of observations
    size_t p_;                          // Number of variables
    arma::imat observations_;           // Categorical observations (n × p)
    arma::mat observations_double_;     // Observations as double (for efficient matrix ops)
    arma::mat observations_double_t_;   // Transposed observations (for BLAS pairwise gradient)
    arma::ivec num_categories_;         // Categories per variable
    arma::uvec is_ordinal_variable_;    // 1 = ordinal, 0 = Blume-Capel
    arma::ivec baseline_category_;      // Reference category for Blume-Capel

    // Sufficient statistics
    arma::imat counts_per_category_;    // Category counts (max_cats+1 × p)
    arma::imat blume_capel_stats_;      // [linear_sum, quadratic_sum] for BC vars (2 × p)
    arma::imat pairwise_stats_;         // X^T X
    arma::mat residual_matrix_;         // X * pairwise_effects (n × p)

    // Parameters
    arma::mat main_effects_;            // Main effect parameters (p × max_cats)
    arma::mat pairwise_effects_;        // Pairwise interaction strengths (p × p, symmetric)
    arma::imat edge_indicators_;        // Edge inclusion indicators (p × p, symmetric binary)

    // Priors
    arma::mat inclusion_probability_;   // Prior inclusion probabilities
    double main_alpha_;                 // Beta prior α
    double main_beta_;                  // Beta prior β
    double pairwise_scale_;             // Cauchy scale for pairwise effects
    arma::mat pairwise_scaling_factors_; // Per-pair scaling factors for Cauchy prior

    // Model configuration
    bool edge_selection_;               // Enable edge selection
    bool edge_selection_active_;        // Currently in edge selection phase

    // Dimension tracking
    size_t num_main_;                   // Total number of main effect parameters
    size_t num_pairwise_;               // Number of possible pairwise effects

    // Proposal SDs (adapted by RWMAdaptationController during warmup)
    arma::mat proposal_sd_main_;
    arma::mat proposal_sd_pairwise_;

    // RWM adaptation controllers (created by init_mh_adaptation)
    std::unique_ptr<RWMAdaptationController> rwm_main_adapter_;
    std::unique_ptr<RWMAdaptationController> rwm_pairwise_adapter_;

    // RNG
    SafeRNG rng_;

    // NUTS/HMC settings
    double step_size_;
    arma::vec inv_mass_;

    // Missing data handling
    bool has_missing_;
    arma::imat missing_index_;

    // Cached gradient components
    arma::vec grad_obs_cache_;
    arma::imat index_matrix_cache_;
    bool gradient_cache_valid_;

    // Interaction indexing (for edge updates)
    arma::imat interaction_index_;
    arma::uvec shuffled_edge_order_;  // Pre-shuffled order (set in prepare_iteration)

    // =========================================================================
    // Private helper methods
    // =========================================================================

    /**
     * Compute sufficient statistics from observations
     */
    void compute_sufficient_statistics();

    /**
     * Count total number of main effect parameters
     */
    size_t count_num_main_effects_internal() const;

    /**
     * Build interaction index matrix
     */
    void build_interaction_index();

    /**
     * Update residual matrix after pairwise effects change
     */
    void update_residual_matrix();

    /**
     * Incrementally update two residual columns after a single pairwise effect change
     */
    void update_residual_columns(int var1, int var2, double delta);

    /**
     * Invalidate gradient cache (call after parameter changes)
     */
    void invalidate_gradient_cache() { gradient_cache_valid_ = false; }

    /**
     * Ensure gradient cache is valid
     */
    void ensure_gradient_cache();

    // -------------------------------------------------------------------------
    // Log-posterior components
    // -------------------------------------------------------------------------

    /**
     * Full log-pseudoposterior (internal, uses current state)
     */
    double log_pseudoposterior_internal() const;

    /**
     * Full log-pseudoposterior with external state (avoids modifying model)
     */
    double log_pseudoposterior_with_state(
        const arma::mat& main_eff,
        const arma::mat& pairwise_eff,
        const arma::mat& residual_mat
    ) const;

    /**
     * Log-posterior for single main effect component
     */
    double log_pseudoposterior_main_component(int variable, int category, int parameter) const;

    /**
     * Log-posterior for single pairwise interaction
     */
    double log_pseudoposterior_pairwise_component(int var1, int var2) const;

    /**
     * Log-likelihood ratio for variable update
     */
    double compute_log_likelihood_ratio_for_variable(
        int variable,
        const arma::ivec& interacting_score,
        double proposed_state,
        double current_state
    ) const;

    /**
     * Log-pseudolikelihood ratio for interaction update
     */
    double log_pseudolikelihood_ratio_interaction(
        int variable1,
        int variable2,
        double proposed_state,
        double current_state
    ) const;

    /**
     * Log-pseudoposterior of a pairwise interaction at proposed = current + delta.
     * Used by tune_proposal_sd() for Robbins-Monro adaptation.
     */
    double log_pseudoposterior_pairwise_at_delta(int var1, int var2, double delta) const;

    // -------------------------------------------------------------------------
    // Parameter vectorization
    // -------------------------------------------------------------------------

    /**
     * Flatten parameters to vector
     */
    arma::vec vectorize_parameters() const;

    /**
     * Flatten parameters into pre-allocated vector (avoids allocation)
     */
    void vectorize_parameters_into(arma::vec& param_vec) const;

    /**
     * Unflatten vector to parameter matrices
     */
    void unvectorize_parameters(const arma::vec& param_vec);

    /**
     * Unvectorize a parameter vector into temporary main/pairwise matrices,
     * then compute the corresponding residual matrix.
     */
    void unvectorize_to_temps(
        const arma::vec& parameters,
        arma::mat& temp_main,
        arma::mat& temp_pairwise,
        arma::mat& temp_residual
    ) const;

    /**
     * Extract active inverse mass (only for included edges)
     */
    arma::vec get_active_inv_mass() const override;

    /**
     * Extract active inverse mass into pre-allocated vector (avoids allocation)
     */
    void get_active_inv_mass_into(arma::vec& active_inv_mass) const;

    // -------------------------------------------------------------------------
    // Metropolis updates
    // -------------------------------------------------------------------------

    /**
     * Update single main effect parameter via RWM
     * @return acceptance probability (for RWM adaptation)
     */
    double update_main_effect_parameter(int variable, int category, int parameter);

    /**
     * Update single pairwise effect via RWM
     * @return acceptance probability (for RWM adaptation)
     */
    double update_pairwise_effect(int var1, int var2);

    /**
     * Update single edge indicator (spike-and-slab)
     */
    void update_edge_indicator(int var1, int var2);
};


/**
 * Factory function to create OMRFModel from R inputs
 */
OMRFModel createOMRFModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection = true
);
