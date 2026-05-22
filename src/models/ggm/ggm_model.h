#pragma once

#include <array>
#include <memory>
#include "models/base_model.h"
#include "math/cholesky_helpers.h"
#include "rng/rng_utils.h"
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "priors/parameter_prior.h"
#include "mcmc/samplers/metropolis_adaptation.h"


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
            std::unique_ptr<BaseParameterPrior> interaction_prior,
            std::unique_ptr<BaseParameterPrior> diagonal_prior,
            const bool na_impute = false
    ) : n_(observations.n_rows - 1),  // centered data has n-1 effective df
        p_(observations.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(observations.t() * observations),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        interaction_prior_(std::move(interaction_prior)),
        diagonal_prior_(std::move(diagonal_prior)),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        proposal_sds_(arma::mat(dim_, 1, arma::fill::ones) * 0.25),
        num_pairwise_(p_ * (p_ - 1) / 2),
        observations_(na_impute ? observations : arma::mat()),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none))
    {
        int num_edges = arma::accu(edge_indicators_) / 2;
        int max_edges = static_cast<int>(p_ * (p_ - 1) / 2);
        has_sparse_graph_ = !edge_selection_ && (num_edges < max_edges);
        initialize_precision_from_mle();
    }

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
     * @param interaction_prior_type Type of interaction prior (Cauchy or Normal)
     */
    GGMModel(
            const int n,
            const arma::mat& suf_stat,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection,
            std::unique_ptr<BaseParameterPrior> interaction_prior,
            std::unique_ptr<BaseParameterPrior> diagonal_prior
    ) : n_(n),
        p_(suf_stat.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(suf_stat),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        interaction_prior_(std::move(interaction_prior)),
        diagonal_prior_(std::move(diagonal_prior)),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        proposal_sds_(arma::mat(dim_, 1, arma::fill::ones) * 0.25),
        num_pairwise_(p_ * (p_ - 1) / 2),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none))
    {
        int num_edges = arma::accu(edge_indicators_) / 2;
        int max_edges = static_cast<int>(p_ * (p_ - 1) / 2);
        has_sparse_graph_ = !edge_selection_ && (num_edges < max_edges);
        initialize_precision_from_mle();
    }

    /** Copy constructor for cloning (required for parallel chains). */
    GGMModel(const GGMModel& other)
        : BaseModel(other),
          target_accept_(other.target_accept_),
          determinant_tilt_(other.determinant_tilt_),
          use_gg_prior_(other.use_gg_prior_),
          V_ij_(other.V_ij_),
          t_(other.t_),
          gg_hyperprior_(other.gg_hyperprior_),
          gg_tcch_a_(other.gg_tcch_a_),
          gg_tcch_b_(other.gg_tcch_b_),
          gg_tcch_r_(other.gg_tcch_r_),
          gg_tcch_s_(other.gg_tcch_s_),
          gg_tcch_u_(other.gg_tcch_u_),
          n_at_gg_setup_(other.n_at_gg_setup_),
          gg_log_g_proposal_sd_(other.gg_log_g_proposal_sd_),
          gg_log_g_n_accept_(other.gg_log_g_n_accept_),
          gg_log_g_n_total_(other.gg_log_g_n_total_),
          tcch_warned_(other.tcch_warned_),
          n_(other.n_),
          p_(other.p_),
          dim_(other.dim_),
          suf_stat_(other.suf_stat_),
          inclusion_probability_(other.inclusion_probability_),
          edge_selection_(other.edge_selection_),
          has_sparse_graph_(other.has_sparse_graph_),
          interaction_prior_(other.interaction_prior_->clone()),
          diagonal_prior_(other.diagonal_prior_->clone()),
          precision_matrix_(other.precision_matrix_),
          cholesky_of_precision_(other.cholesky_of_precision_),
          inv_cholesky_of_precision_(other.inv_cholesky_of_precision_),
          covariance_matrix_(other.covariance_matrix_),
          edge_indicators_(other.edge_indicators_),
          vectorized_parameters_(other.vectorized_parameters_),
          vectorized_indicator_parameters_(other.vectorized_indicator_parameters_),
          proposal_sds_(other.proposal_sds_),
          shuffled_edge_order_(other.shuffled_edge_order_),
          num_pairwise_(other.num_pairwise_),
          rng_(other.rng_),
          observations_(other.observations_),
          has_missing_(other.has_missing_),
          missing_index_(other.missing_index_),
          precision_proposal_(other.precision_proposal_),
          constraint_structure_(other.constraint_structure_),
          gradient_engine_(other.gradient_engine_),
          constraint_dirty_(other.constraint_dirty_),
          theta_valid_(other.theta_valid_),
          theta_(other.theta_)
    {
        // The cloned interaction_prior_ / diagonal_prior_ may be
        // GraphicalG-typed and still carry pointers into the source
        // model's V_ij_ / t_. Rebind them to this clone's state.
        gg_rebind_priors_();
    }

    /** @return true (GGM supports NUTS via free-element Cholesky gradient). */
    bool has_gradient()        const override { return true; }
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

    /**
     * Set the Robbins-Monro target acceptance rate used by the
     * adaptive-Metropolis updates of this GGM. Honoured by all
     * Metropolis sweeps (off-diagonal and diagonal).
     */
    void set_metropolis_target_accept(double target) override {
        target_accept_ = target;
    }

    /**
     * Construct Robbins-Monro adaptation controller for the per-iteration
     * MH proposal SDs. Called once by MetropolisSampler before warmup;
     * under NUTS this is never called and the controller stays null.
     */
    void init_metropolis_adaptation(const WarmupSchedule& schedule) override;

    /**
     * Set the determinant-tilt exponent delta. Adds delta * log|K| to the
     * log-prior, pushing the chain away from the PD-cone boundary. delta = 0
     * (default) recovers the untilted target. Currently consumed only by
     * the NUTS gradient engine; the MH path is unchanged. Triggers an engine
     * rebuild on the next gradient call.
     */
    void set_determinant_tilt(double delta) {
        determinant_tilt_ = delta;
        constraint_dirty_ = true;
    }

    /**
     * G-prior hyperprior on g = t² for the Graphical G-prior. Selects the
     * update kernel used at end of sweep (see `gg_update_t_()`).
     *
     *   Fixed         : g held at the user-supplied value (no update).
     *   ConjugateGamma: Gamma(a₀, b₀) prior on t = √g; closed-form Gibbs
     *                   update at sweep end (see implementation plan §5a.3).
     *   ZellnerSiow   : inverse-Gamma(½, ½) on g, induces a Cauchy marginal
     *                   per edge after marginalising g. MH update.
     *   HyperG        : Beta hyper-g family (a − 2)·(1 + g/n)^{−a/2}.
     *   HyperGOverN   : hyper-g / n re-scaling.
     *   TCCH          : umbrella family parameterised by (a, b, r, s, u).
     *
     * Default `Fixed` with `g = 1` reproduces a constant-scale Normal slab.
     */
    enum class GGHyperprior {
        Fixed, ConjugateGamma, ZellnerSiow, HyperG, HyperGOverN, TCCH
    };

    /**
     * Enable the Graphical G-prior. Replaces the slab interaction_prior and
     * the diagonal_prior with GraphicalGPrior / GraphicalGDiag instances
     * bound to the model's internal V_ij table and t state.
     *
     * - `V_ij` is built once at this call time from the model's sufficient
     *   statistic: V_ij = n / (4 · suf_stat_ii · suf_stat_jj). This is the
     *   data-adaptive Fisher-info-based per-edge scale of IDEA.md §5.1.
     * - `t` is initialised to `sqrt(g_init)` and is updated each sweep
     *   under the chosen hyperprior.
     * - The user-supplied `interaction_prior` and `diagonal_prior` are
     *   discarded — the GG-prior fixes both (paired by the scale-matching
     *   identity that the joint Z(Γ, g; δ) = g^{q δ/2} · Z̃(Γ) collapses).
     *
     * Called from `sample_ggm` after the model has been constructed.
     */
    void enable_gg_prior(
        GGHyperprior hyperprior,
        double g_init       = 1.0,
        double tcch_a       = 1.0,
        double tcch_b       = 1.0,
        double tcch_r       = 0.0,
        double tcch_s       = 0.0,
        double tcch_u       = 1.0);

    bool gg_prior_enabled() const { return use_gg_prior_; }
    double gg_current_t()  const { return t_; }
    double gg_current_g()  const { return t_ * t_; }

    /**
     * Mute the data likelihood: set n_ = 0 and suf_stat_ = 0, then rebuild
     * the gradient engine. Used to run a prior-only chain (samples from
     * the joint prior p(K, Gamma) under the original V_ij computed from
     * the data). Call AFTER enable_gg_prior() so the GG-prior V_ij has
     * already been cached from the original n_ and suf_stat_.
     */
    void set_prior_only();
    bool prior_only() const { return prior_only_; }

    /** Per-iteration GG-prior diagnostic getters used by the chain
     *  runner. Active iff the GG-prior is enabled on this model. */
    bool   has_gg_diagnostics() const override { return use_gg_prior_; }
    double current_gg_t()       const override { return t_; }
    const arma::mat& gg_V_ij() const { return V_ij_; }
    /** Read-only view of the current precision matrix. Used by test
     *  harnesses to snapshot the chain state directly. */
    const arma::mat& get_precision_matrix() const { return precision_matrix_; }

    /** Shuffle edge visit order (random scan). */
    void prepare_iteration() override;

    /** Sweep over edges in shuffled order, proposing add/remove moves. */
    void update_edge_indicators() override;

    /**
     * Element-wise MH updates for proposal-SD tuning during stage 3b.
     *
     * Runs off-diagonal and diagonal Metropolis updates with
     * Robbins-Monro adaptation, following the OMRF pattern.
     */
    void tune_proposal_sd(int iteration, const WarmupSchedule& schedule) override;

    /**
     * Combined log-posterior and gradient for NUTS.
     *
     * Uses the free-element Cholesky parameterization:
     * theta = (psi_1, f_2, psi_2, ..., f_p, psi_p) where psi_q = log(phi_qq)
     * and x_q = N_q f_q gives the off-diagonal Cholesky entries.
     *
     * @param parameters  Active theta vector (dimension = p + |E|)
     * @return (log-posterior, gradient) pair
     */
    std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) override;

    /**
     * Set model state from a theta vector (inverse of get_vectorized_parameters).
     *
     * Runs the forward map theta -> Phi -> K and updates all internal
     * matrices (precision, Cholesky, inverse Cholesky, covariance).
     *
     * @param parameters  Active theta vector (dimension = p + |E|)
     */
    void set_vectorized_parameters(const arma::vec& parameters) override;

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

    /**
     * @return Active theta dimension: p + |E| (diagonals + included edges).
     *
     * Changes when edge indicators toggle. Used by NUTS for leapfrog
     * integration.
     */
    size_t parameter_dimension() const override;

    /**
     * @return Full theta dimension: p + p(p-1)/2 (all possible off-diag slots).
     *
     * Fixed across all graphs. Used by the adaptation controller for
     * mass-matrix sizing.
     */
    size_t full_parameter_dimension() const override;

    /**
     * @return Storage dimension for sample output: p(p+1)/2 (upper triangle of K).
     *
     * Preserves the existing output contract: downstream R code expects
     * the upper triangle of the precision matrix.
     */
    size_t storage_dimension() const override { return dim_; }

    /**
     * Set random seed for reproducibility.
     * @param seed  Integer seed value
     */
    void set_seed(int seed) override {
        rng_ = SafeRNG(seed);
    }

    /**
     * @return Active theta vector: (psi_1, f_2, psi_2, ..., f_p, psi_p).
     *
     * Dimension = parameter_dimension() = p + |E|. Used by NUTS as the
     * current state. Recomputed lazily from Phi when stale.
     */
    arma::vec get_vectorized_parameters() const override;

    /**
     * @return Full (zero-padded) theta vector for mass-matrix adaptation.
     *
     * Dimension = full_parameter_dimension() = p + p(p-1)/2. Inactive
     * edges have their f_q slots set to zero.
     */
    arma::vec get_full_vectorized_parameters() const override;

    /**
     * @return Upper triangle of the precision matrix for sample storage.
     *
     * Preserves the existing output contract.
     */
    arma::vec get_storage_vectorized_parameters() const override {
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

    /**
     * @return Active subset of the inverse mass diagonal.
     *
     * Filters the full inv_mass_ (dimension p + p(p-1)/2) to active
     * parameters only (dimension p + |E|). For columns where N_q != I,
     * rotates the per-Cholesky-entry variances into f_q coordinates.
     */
    arma::vec get_active_inv_mass() const override;

    // GGMModel uses the theta-space (free-element Cholesky) NUTS path
    // exclusively. RATTLE projection (project_position/project_momentum,
    // full-position get/set, full-space gradient) is not implemented;
    // the BaseModel defaults are sufficient — they are never reached
    // because has_constraints() defaults to false. See nuts_sampler.h
    // for the do_unconstrained_step path actually taken.

    /** @return Deep copy of this model. */
    std::unique_ptr<BaseModel> clone() const override {
        return std::make_unique<GGMModel>(*this);
    }

private:

    // Robbins-Monro target acceptance rate for adaptive-Metropolis
    // proposal-SD tuning. Set via set_metropolis_target_accept(); defaults
    // to 0.44 (componentwise random-walk Metropolis optimum).
    double target_accept_ = 0.44;

    /// Per-iteration adaptation controller (MH mode only — under NUTS this
    /// stays null and the stage-3b path in tune_proposal_sd is used instead).
    std::unique_ptr<MetropolisAdaptationController> metropolis_adapter_;

    // Determinant-tilt exponent (see set_determinant_tilt). Forwarded to
    // GGMGradientEngine on every rebuild.
    double determinant_tilt_ = 0.0;

    // ---- Graphical G-prior state (see enable_gg_prior) -----------------
    // V_ij_(i, j) = n / (4 · suf_stat_(i,i) · suf_stat_(j,j))  for i ≠ j
    //             (diagonal entries unused; left zero)
    // t_         = √g, the shared scale hyperparameter
    // The GraphicalGPrior / GraphicalGDiag instances stored in
    // interaction_prior_ / diagonal_prior_ hold non-owning pointers into
    // these two members; they MUST be rebound when the model is cloned.
    bool         use_gg_prior_ = false;
    /// When true, n_ and suf_stat_ have been zeroed so the chain targets the
    /// prior alone. The GG-prior V_ij was cached from the original data at
    /// enable_gg_prior() time, so the prior over (K, Gamma) is unchanged.
    bool         prior_only_   = false;
    arma::mat    V_ij_;
    double       t_            = 1.0;
    GGHyperprior gg_hyperprior_ = GGHyperprior::Fixed;
    // tCCH-family hyperparameters (see IDEA.md §2 and Li & Clyde 2018 for
    // the parameterisation). For non-tCCH hyperpriors only a subset is
    // used; the defaults reproduce Zellner-Siow at a=b=½ etc.
    double       gg_tcch_a_ = 1.0;
    double       gg_tcch_b_ = 1.0;
    double       gg_tcch_r_ = 0.0;
    double       gg_tcch_s_ = 0.0;
    double       gg_tcch_u_ = 1.0;
    /// Original sample size captured at enable_gg_prior() time. Used by the
    /// hyper-g/n hyperprior, which scales by n; under prior-only mode n_ is
    /// zeroed but this remains the data-defined n for the prior on g.
    int          n_at_gg_setup_ = 0;
    /// MH-on-log(g) proposal SD (Gaussian random walk on log scale). Used by
    /// the ZS / hyper-g / hyper-g/n / tCCH hyperprior branches that lack a
    /// closed-form Gibbs update. Default 1.0 in log space gives reasonable
    /// mixing for the prior families implemented here.
    double       gg_log_g_proposal_sd_ = 1.0;
    /// Accepted / total MH-on-log(g) proposals (diagnostic only).
    long long    gg_log_g_n_accept_ = 0;
    long long    gg_log_g_n_total_  = 0;
    /// One-shot flag so the tCCH not-yet-implemented warning fires only
    /// once per chain instead of once per sweep.
    bool         tcch_warned_       = false;

    /// One-shot V_ij table construction from suf_stat_. Called by
    /// enable_gg_prior(); idempotent.
    void compute_gg_V_ij_();
    /// End-of-sweep update of t under the current hyperprior. Includes
    /// the Cholesky / covariance rescaling that follows from rescaling
    /// Θ ← (t_new/t_old) · Θ in the (η, t)-parameterisation. Called from
    /// prepare_iteration() / do_one_metropolis_step() (TBD).
    void gg_update_t_();
    /// MH on log(g) under a scale-matched joint proposal (Θ, g) →
    /// (α · Θ, g_new). Drives the ZS / hyper-g / hyper-g/n branches that
    /// lack a closed-form Gibbs update. The scale-matching collapses the
    /// slab contribution to the MH ratio, leaving only the hyperprior,
    /// diagonal Gamma, det-tilt, likelihood, and Jacobian terms.
    void gg_mh_update_log_g_();
    /// Evaluates log π(g) for the current hyperprior family. Returns 0.0
    /// for families without an MH update (Fixed, ConjugateGamma).
    double gg_log_hyperprior_(double g) const;
    /// Rebind GraphicalG{Prior,Diag} pointers after the model is cloned.
    /// Called from the copy constructor; safe to call when the priors
    /// aren't GraphicalG (no-ops).
    void gg_rebind_priors_();

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
    /// Whether the initial graph excludes any edges (triggers RATTLE).
    bool has_sparse_graph_ = false;
    /// Prior on off-diagonal precision elements (interactions).
    std::unique_ptr<BaseParameterPrior> interaction_prior_;
    /// Prior on diagonal precision elements (scale).
    std::unique_ptr<BaseParameterPrior> diagonal_prior_;

    /// Precision matrix Omega, its Cholesky factor R (Omega = R'R),
    /// inverse Cholesky factor, and covariance matrix.
    arma::mat precision_matrix_, cholesky_of_precision_, inv_cholesky_of_precision_, covariance_matrix_;
    /// Current edge-indicator matrix (p x p, symmetric, 0/1).
    arma::imat edge_indicators_;
    /// Pre-allocated storage returned by get_vectorized_parameters().
    arma::vec vectorized_parameters_;
    /// Pre-allocated storage returned by get_vectorized_indicator_parameters().
    arma::ivec vectorized_indicator_parameters_;

    /// Proposal standard deviations for Metropolis updates (one per element,
    /// stored as a (dim_, 1) matrix so it can be wrapped by
    /// MetropolisAdaptationController).
    arma::mat proposal_sds_;

    /// Shuffled edge visit order for random-scan edge selection.
    arma::uvec shuffled_edge_order_;
    /// Number of unique off-diagonal pairs: p(p-1)/2.
    size_t num_pairwise_ = 0;
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
     * @param i  Row index (i < j)
     * @param j  Column index
     * @return   Metropolis acceptance probability min(1, exp(ln_alpha)),
     *           or 0.0 if the edge is inactive (caller masks it out).
     */
    double update_edge_parameter(size_t i, size_t j);

    /**
     * Propose a new diagonal precision entry on the log scale.
     * Accepts or rejects with a Metropolis ratio using the Gaussian
     * likelihood, a Gamma(1,1) prior, and a Jacobian correction.
     *
     * @param i  Diagonal index
     * @return   Metropolis acceptance probability min(1, exp(ln_alpha)).
     */
    double update_diagonal_parameter(size_t i);

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
     * log|K_prop| - log|K_curr| for a rank-2 off-diagonal proposal at (i, j),
     * computed via the matrix-determinant lemma in O(p). Reads
     * precision_matrix_, precision_proposal_, and covariance_matrix_; assumes
     * precision_proposal_ has already been set up at (i, j), (j, i), (j, j).
     * Used to add the determinant-tilt term delta * (log|K_prop| - log|K_curr|)
     * to MH ratios.
     */
    double log_det_ratio_edge(size_t i, size_t j) const;

    /**
     * log|K_prop| - log|K_curr| for a rank-1 diagonal proposal at j.
     * Computed via the matrix-determinant lemma in O(1). Reads the same
     * cached state as log_det_ratio_edge.
     */
    double log_det_ratio_diag(size_t j) const;



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

    /**
     * Recompute Cholesky and its inverse from the precision matrix.
     *
     * Used as a fallback when accumulated rank-1 updates/downdates
     * cause numerical drift that makes the triangular inverse fail.
     * Resets both cholesky_of_precision_ and inv_cholesky_of_precision_
     * from precision_matrix_, then recomputes covariance_matrix_.
     */
    void refresh_cholesky();

    /**
     * Initialize precision matrix at the regularized MLE.
     *
     * Computes K = n * inv(S + delta * I) where delta provides
     * Ledoit-Wolf-style shrinkage toward identity. Gives NUTS a
     * starting point near the posterior mode, avoiding the step-size
     * instability that arises when starting from K = I far from the
     * mode.
     */
    void initialize_precision_from_mle();

    // =================================================================
    // NUTS gradient support
    // =================================================================

    /// Graph constraint structure (rebuilt when edge indicators change).
    GraphConstraintStructure constraint_structure_;
    /// Gradient engine for the free-element Cholesky parameterization.
    GGMGradientEngine gradient_engine_;
    /// Whether the constraint structure needs rebuilding.
    bool constraint_dirty_ = true;
    /// Whether theta_ is in sync with cholesky_of_precision_.
    mutable bool theta_valid_ = false;
    /// Cached theta vector (active parameterization).
    mutable arma::vec theta_;

public:
    /**
     * Rebuild the constraint structure and gradient engine from current
     * edge indicators. Called lazily before gradient evaluation.
     */
    void ensure_constraint_structure();

    /**
     * Convert the current Cholesky factor to the theta parameterization.
     *
     * For each column q, computes psi_q = log(phi_qq) and
     * f_q = N_q^T x_q where x_q = Phi[0:q-1, q].
     */
    void recompute_theta() const;
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
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> diagonal_prior,
    const bool na_impute = false
);
