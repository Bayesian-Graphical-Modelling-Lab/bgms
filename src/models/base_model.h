#pragma once

#include <RcppArmadillo.h>
#include <stdexcept>
#include <memory>

// Forward declarations
struct StepResult;
struct SafeRNG;
struct WarmupSchedule;

/**
 * BaseModel — Abstract interface for all graphical models.
 *
 * Defines the virtual methods that the MCMC framework (MetropolisSampler,
 * NUTSSampler, ChainRunner) calls during sampling. Most methods are pure
 * virtual (`= 0`) so the compiler enforces implementation in every
 * subclass. The exceptions are gradient-only methods (gradient,
 * logp_and_gradient, set_vectorized_parameters) which throw at runtime
 * and are guarded by the has_gradient() capability query.
 *
 * Subclass hierarchy:
 *   - GGMModel  — Gaussian Graphical Model (precision matrix, Metropolis only)
 *   - OMRFModel — Ordinal Markov Random Field (Metropolis + NUTS/HMC)
 *
 * Methods fall into several groups:
 *   - **Capability queries** (has_gradient, has_adaptive_metropolis, etc.)
 *     — samplers inspect these to choose the right algorithm.
 *   - **Sampling steps** (do_one_metropolis_step, gradient, logp_and_gradient)
 *     — called by the sampler each iteration.
 *   - **Edge selection** (update_edge_indicators, initialize_graph)
 *     — spike-and-slab structure learning.
 *   - **Parameter access** (get/set_vectorized_parameters, get_full_vectorized_parameters)
 *     — used by NUTS for momentum-based proposals and by the runner for output.
 *   - **Adaptation** (set_step_size, set_inv_mass, tune_proposal_sd)
 *     — tuned during warmup.
 *   - **Missing data** (has_missing_data, impute_missing)
 *     — full-conditional imputation between iterations.
 */
class BaseModel {
public:
    virtual ~BaseModel() = default;

    // =========================================================================
    // Capability queries
    // =========================================================================

    /** @return true if the model provides a gradient (enables NUTS/HMC). */
    virtual bool has_gradient() const { return false; }
    /** @return true if the model supports adaptive Metropolis. */
    virtual bool has_adaptive_metropolis() const { return false; }
    /** @return true if the model supports NUTS (default: same as has_gradient). */
    virtual bool has_nuts() const { return has_gradient(); }
    /** @return true if the model supports edge selection (spike-and-slab). */
    virtual bool has_edge_selection() const { return false; }

    // =========================================================================
    // Core sampling methods
    // =========================================================================

    /**
     * Compute the gradient of the log-(pseudo)posterior.
     * @param parameters  Vectorized parameters at which to evaluate
     * @return Gradient vector
     */
    virtual arma::vec gradient(const arma::vec& parameters) {
        throw std::runtime_error("Gradient not implemented for this model");
    }

    /**
     * Combined log-(pseudo)posterior and gradient evaluation.
     *
     * More efficient than calling log-posterior and gradient() separately
     * because intermediate computations can be shared.
     *
     * @param parameters  Vectorized parameters at which to evaluate
     * @return Pair of (log-posterior value, gradient vector)
     */
    virtual std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) {
        throw std::runtime_error("logp_and_gradient not implemented for this model");
    }

    /**
     * Perform one full Metropolis sweep over all parameters.
     *
     * The model handles its own parameter grouping (e.g. off-diagonal,
     * diagonal, edge indicators in GGM; main, pairwise in OMRF).
     *
     * @param iteration  Current iteration index (for Robbins-Monro adaptation)
     */
    virtual void do_one_metropolis_step(int iteration = -1) = 0;

    /**
     * Initialize Metropolis adaptation controllers.
     *
     * Called once before the MCMC loop begins. Subclasses store the
     * warmup schedule for later use by Robbins-Monro adaptation.
     *
     * @param schedule  Warmup schedule defining adaptation phases
     */
    virtual void init_metropolis_adaptation(const WarmupSchedule& /*schedule*/) {}

    /**
     * Tune pairwise proposal SDs via Robbins-Monro (warmup Stage 3b).
     *
     * Called every iteration from the runner; the implementation checks
     * the schedule internally to decide whether to adapt.
     *
     * @param iteration  Current iteration index
     * @param schedule   Warmup schedule
     */
    virtual void tune_proposal_sd(int /*iteration*/, const WarmupSchedule& /*schedule*/) {}

    /**
     * Called at the start of every iteration before edge selection and sampling.
     *
     * Subclasses may use this to shuffle edge update order or advance the
     * RNG state consistently.
     */
    virtual void prepare_iteration() {}

    // =========================================================================
    // Edge selection
    // =========================================================================

    /**
     * Update edge indicators via Metropolis-Hastings add-delete moves.
     *
     * Only meaningful when has_edge_selection() returns true. GGMModel
     * handles this inside do_one_metropolis_step() instead.
     */
    virtual void update_edge_indicators() = 0;

    // =========================================================================
    // Parameter vectorization
    // =========================================================================

    /** @return Active parameters as a flat vector (dimension may change with edge selection). */
    virtual arma::vec get_vectorized_parameters() const = 0;

    /**
     * Set parameters from a flat vector (inverse of get_vectorized_parameters).
     * @param parameters  Vectorized parameter values
     */
    virtual void set_vectorized_parameters(const arma::vec& parameters) {
        throw std::runtime_error("set_vectorized_parameters method must be implemented in derived class");
    }

    /** @return Edge indicators as a flat integer vector. */
    virtual arma::ivec get_vectorized_indicator_parameters() = 0;

    /**
     * @return Full parameter dimension (fixed size, includes inactive parameters).
     *
     * Used by GradientSamplerBase for mass-matrix sizing and adaptation.
     * For most models this equals the storage dimension. For models where
     * some parameters are not sampled by NUTS (e.g., MixedMRFModel's
     * continuous precision),
     * this returns the NUTS-block dimension.
     * Defaults to parameter_dimension().
     */
    virtual size_t full_parameter_dimension() const {
        return parameter_dimension();
    }

    /**
     * @return All parameters in a fixed-size vector (inactive edges are 0).
     *
     * Used by GradientSamplerBase for adaptation (online covariance).
     * Dimension must match full_parameter_dimension().
     */
    virtual arma::vec get_full_vectorized_parameters() const = 0;

    /** @return Dimensionality of the active parameter space. Pure virtual. */
    virtual size_t parameter_dimension() const = 0;

    /**
     * @return Dimension for sample storage (includes all parameters).
     *
     * For most models this equals full_parameter_dimension(). Override
     * when storage needs more entries than the NUTS block (e.g., continuous
     * precision parameters in MixedMRFModel).
     */
    virtual size_t storage_dimension() const {
        return full_parameter_dimension();
    }

    /**
     * @return All parameters in a fixed-size vector for sample storage.
     *
     * Dimension must match storage_dimension(). Default delegates to
     * get_full_vectorized_parameters().
     */
    virtual arma::vec get_storage_vectorized_parameters() const {
        return get_full_vectorized_parameters();
    }

    // =========================================================================
    // Infrastructure
    // =========================================================================

    /**
     * Set the random seed for reproducibility.
     * @param seed  Integer seed value
     */
    virtual void set_seed(int seed) = 0;

    /** @return Deep copy of this model (for parallel chains). */
    virtual std::unique_ptr<BaseModel> clone() const = 0;

    /** @return Reference to the model's random number generator. */
    virtual SafeRNG& get_rng() = 0;

    // =========================================================================
    // NUTS/HMC adaptation
    // =========================================================================

    /**
     * Set the leapfrog step size for gradient-based samplers.
     * @param step_size  New step size
     */
    virtual void set_step_size(double step_size) { step_size_ = step_size; }
    /** @return Current leapfrog step size. */
    virtual double get_step_size() const { return step_size_; }

    /**
     * Set the inverse mass matrix diagonal for HMC/NUTS.
     * @param inv_mass  Diagonal elements of the inverse mass matrix
     */
    virtual void set_inv_mass(const arma::vec& inv_mass) { inv_mass_ = inv_mass; }
    /** @return Current inverse mass matrix diagonal. */
    virtual const arma::vec& get_inv_mass() const { return inv_mass_; }

    /**
     * @return Active subset of the inverse mass diagonal.
     *
     * For models with edge selection, this may return only the entries
     * corresponding to included edges. Default: returns the full diagonal.
     */
    virtual arma::vec get_active_inv_mass() const { return inv_mass_; }

    // =========================================================================
    // Edge selection control
    // =========================================================================

    /**
     * Enable or disable edge-selection proposals.
     * @param active  true to enable edge add-delete moves
     */
    virtual void set_edge_selection_active(bool active) {
        (void)active;
    }

    /** Draw initial edge states from prior inclusion probabilities. */
    virtual void initialize_graph() {}

    // =========================================================================
    // Missing data
    // =========================================================================

    /** @return true when missing-data imputation is active. */
    virtual bool has_missing_data() const { return false; }

    /** Impute missing entries from full-conditional distributions. */
    virtual void impute_missing() {}

    // =========================================================================
    // Edge prior support
    // =========================================================================

    /** @return Current edge-indicator matrix. */
    virtual const arma::imat& get_edge_indicators() const = 0;

    /** @return Mutable reference to the prior inclusion-probability matrix. */
    virtual arma::mat& get_inclusion_probability() = 0;

    /** @return Number of variables (p). */
    virtual int get_num_variables() const = 0;

    /** @return Number of unique off-diagonal pairs p(p-1)/2. */
    virtual int get_num_pairwise() const = 0;

protected:
    BaseModel() = default;
    /// Leapfrog step size for NUTS/HMC.
    double step_size_ = 0.1;
    /// Inverse mass matrix diagonal for NUTS/HMC.
    arma::vec inv_mass_;
};
