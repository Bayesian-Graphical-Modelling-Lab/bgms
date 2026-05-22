#pragma once

#include <array>
#include <memory>
#include <vector>
#include "models/base_model.h"
#include "math/cholesky_helpers.h"
#include "rng/rng_utils.h"
#include "models/ggm/graph_constraint_structure.h"
#include "models/ggm/ggm_gradient.h"
#include "models/ggm/degord_sampler.h"
#include "priors/parameter_prior.h"
#include "mcmc/samplers/metropolis_adaptation.h"


/**
 * Graph-prior specification for the GGM with edge selection.
 *
 *   Joint         (default): π_joint(K, Γ) ∝ slab·diag·|K|^δ·1{K∈M+(Γ)}·π(Γ).
 *                 Γ marginal is π(Γ)·Z(Γ).
 *   Hierarchical: π_hier(K, Γ)  ∝ slab·diag·|K|^δ·1{K∈M+(Γ)}/Z(Γ)·π(Γ).
 *                 Γ marginal is π(Γ) directly. Requires the Z(Γ) ratio to be
 *                 estimated unbiasedly per between-edge proposal; implemented
 *                 via the DEGORD-permuted V/RR estimator (Phase 2 + 3).
 *
 * Hierarchical mode requires the slab to be NormalPrior and the diagonal to
 * be GammaScalePrior (the closed-form log_Z_NLO_gamma machinery only
 * supports this prior family). Construction will throw if hierarchical is
 * requested under any other family.
 */
enum class GraphPriorSpec { Joint, Hierarchical };


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

    /** Copy constructor for cloning (required for parallel chains).
     *
     * Note on hierarchical-spec state: only the *configuration* (graph_prior_spec_,
     * v_M_inner_, v_kappa_, v_rho_) is copied. The lazy state (pools, chain aux,
     * log_Z_NLO_curr_) is reset so each cloned chain rebuilds it on first
     * ensure_hierarchical_state_() call - this is intentional because pools are
     * RNG-derived and would otherwise share random draws across chains.
     */
    GGMModel(const GGMModel& other)
        : BaseModel(other),
          target_accept_(other.target_accept_),
          determinant_tilt_(other.determinant_tilt_),
          graph_prior_spec_(other.graph_prior_spec_),
          hierarchical_state_built_(false),
          prior_params_extracted_(false),
          use_manuscript_nlo_(other.use_manuscript_nlo_),
          mh_U_(other.mh_U_),
          mh_U_local_K_(other.mh_U_local_K_),
          mh_U_local_K_global_freq_(other.mh_U_local_K_global_freq_),
          plug_in_nlo_(other.plug_in_nlo_),
          use_sd_between_step_(other.use_sd_between_step_),
          use_sd_lspace_(other.use_sd_lspace_),
          prior_only_(other.prior_only_),
          n_pd_reverts_(other.n_pd_reverts_),
          v_M_inner_(other.v_M_inner_),
          v_kappa_(other.v_kappa_),
          v_rho_(other.v_rho_),
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
    {}

    /** @return true (GGM supports NUTS via free-element Cholesky gradient). */
    bool has_gradient()        const override { return true; }
    /** @return true (GGM supports adaptive Metropolis). */
    bool has_adaptive_metropolis()     const override { return true; }
    /** @return true when edge selection is enabled. */
    bool has_edge_selection()  const override { return edge_selection_; }
    /** @return true when missing-data imputation is active. */
    bool has_missing_data()    const override { return has_missing_; }

    /** @return true under hierarchical graph_prior_spec — the only path
     *  where V(Γ, U) is computed and sign / log|V| are meaningful. */
    bool has_v_ratio_diagnostics() const override {
        return graph_prior_spec_ == GraphPriorSpec::Hierarchical;
    }
    int    current_sign_V()    const override { return current_sign_V_; }
    double current_log_abs_V() const override { return current_log_abs_V_; }
    int    current_K_depth()   const override { return v_K_depth_; }

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
     * Switch the chain to hierarchical-spec inference (default is Joint).
     * Validates the slab/diag prior family is (NormalPrior, GammaScalePrior)
     * — throws std::runtime_error if not. Lazy: state is built on first
     * use (next prepare_iteration or between-edge proposal).
     */
    void set_graph_prior_spec(GraphPriorSpec spec);

    /**
     * Configure the V/RR estimator tuning. Defaults: M_inner=100, kappa=1.0,
     * rho=0.5. Only consumed when graph_prior_spec_ == Hierarchical.
     */
    void set_z_ratio_tuning(int M_inner, double kappa, double rho);

    /**
     * Select the analytic NLO formula used as the RR centring in the
     * between-Γ MH ratio. `false` (default) keeps the pre-2026-05-21 bgms
     * formula (log_Z_NLO_gamma). `true` switches to the manuscript App C
     * NLO (eq:NLO-decomp; ~/SV/Z/notes/2026-05-21_message-to-bgms-companion-NLO-fix.md).
     * The manuscript form drops the `(ν_i + α) / M_v[i]` factor in B_e
     * compared to bgms's existing formula; at α=1, δ=0 the two coincide
     * bit-exactly, with divergence growing with δ. Only consumed under
     * Hierarchical graph_prior_spec and α=1; at α ≠ 1 the flag is a no-op
     * (manuscript Na/Nb/Nc terms are not ported yet).
     */
    void set_use_manuscript_nlo(bool on) { use_manuscript_nlo_ = on; }
    bool use_manuscript_nlo() const { return use_manuscript_nlo_; }

    /**
     * Switch the between-Γ move to the Savage-Dickey variant (handoff from
     * Z, 2026-05-22). Replaces the joint-spec / hierarchical-spec MH ratio
     * with the marginalised-γ form using a 1D Laplace+NLO conditional
     * density at K_ij = 0; on ADD accept, K_ij is drawn from the 1D
     * Laplace proposal; on DEL accept, K_ij = 0 deterministically.
     * Within-K Roverato and diagonal moves remain unchanged. Default
     * `false` keeps the existing between-Γ path.
     */
    void set_use_sd_between_step(bool on) { use_sd_between_step_ = on; }
    bool use_sd_between_step() const { return use_sd_between_step_; }

    /**
     * Switch the SD between-Γ move to the L-space (Cholesky) variant. The
     * chain dynamics shift from "K_ij ↔ 0 with K_jj fixed" (Z's K-space SD)
     * to "l_ji ↔ m_ij(K_{-ij}) with K_jj implicitly slaved through K = LL^T".
     * PD becomes automatic — l_ji ∈ ℝ unconstrained — and the per-edge
     * conditional is exactly Gaussian at α = 1 (closed-form Gibbs). The
     * MH ratio includes a per-step Jacobian log l_ii. Only active when
     * use_sd_between_step_ is also true. α = 1 only (α > 1 not implemented).
     */
    void set_use_sd_lspace(bool on) { use_sd_lspace_ = on; }
    bool use_sd_lspace() const { return use_sd_lspace_; }

    /**
     * Diagnostic: number of times the L-space SD between-step's PD-revert
     * defense fired during this chain. Should be zero when the Bunch
     * permutation path is numerically stable; non-zero is a canary for either
     * a bug in the extraction or genuine pathological K. Reset via
     * reset_pd_revert_count().
     */
    long n_pd_reverts() const { return n_pd_reverts_; }
    void reset_pd_revert_count() { n_pd_reverts_ = 0; }

    /**
     * Disable the likelihood contribution in all MH ratios. When true the
     * chain targets the prior π(Γ, K) instead of the posterior π(Γ, K | Y).
     * Used for stationarity / calibration tests: under prior-only mode the
     * empirical edge-inclusion frequencies should converge to inclusion_
     * probability_. Affects update_edge_parameter, update_diagonal_parameter,
     * and update_edge_indicator_parameter_pair_sd. Default `false`.
     */
    void set_prior_only(bool on) { prior_only_ = on; }
    bool prior_only() const { return prior_only_; }

    /**
     * Enable a Metropolis–Hastings step on the auxiliary U-pool at the start
     * of each iteration in place of the legacy "draw U fresh from μ".
     *
     * Background: the chain is a block PMMH on (Γ, K, U, N) targeting
     * π(Γ, K) · V(U, N) · μ(U) · P(N). Under that target, p(U | Γ, K) is
     * proportional to V(U; Γ, K) · μ(U) — NOT μ(U) alone. Refreshing U
     * by a fresh draw from μ skips the V-tilt and breaks invariance;
     * empirically this yields a small but systematic Γ-marginal bias
     * (~−0.001 nats at p=20, p_inc=0.05; 5/5 seeds same sign in our
     * tests).
     *
     * With this flag, the U-refresh is replaced by a proper MH step:
     * propose U_new ~ μ, N_new ~ P(N) and accept with log α =
     * log|V(Γ, U_new)| − log|V(Γ, U_old)| (μ, P(N) cancel by proposal
     * symmetry). Cross-implementation experiments confirm this collapses
     * the bias to MC noise (companion-AI delivery 2026-05-21).
     *
     * Default `false` preserves the pre-2026-05-21 chain and SBC-clean
     * baselines.
     */
    void set_mh_U(bool on) { mh_U_ = on; }
    bool mh_U() const { return mh_U_; }

    /**
     * Plug-in mNLO mode: replace the RR/V/U machinery with a deterministic
     * closed-form Z(Γ) ratio in the between-Γ MH. Trades exactness for
     * predictable cost (no K-tail, no pool work, flat per-toggle wall in p).
     *
     * Bias: smooth in Γ, of order |log Z(Γ) − log_Z_NLO(Γ)| per toggle —
     * bit-exact at δ=0, controlled at small δ, grows with δ. The chain
     * targets π(Γ) · exp(log_Z_NLO(Γ) − log Z(Γ)), i.e. the hierarchical
     * target distorted by the centring miss. Under good mNLO centring this
     * distortion is small and operationally negligible vs the K-scaling
     * cost of the exact RR variant in dense / large-p regimes.
     *
     * Mutually exclusive with the RR machinery: when set, refresh_auxiliary_u
     * early-returns and the between-Γ MH skips V_log_pair entirely. The
     * mh_U / mh_U_local_K flags are ignored.
     *
     * Default `false`: exact RR + (optional) mh_U + (optional) local-K path.
     */
    void set_plug_in_nlo(bool on) { plug_in_nlo_ = on; }
    bool plug_in_nlo() const { return plug_in_nlo_; }

    /**
     * Enable a local random-walk move on the RR truncation depth K instead
     * of (or alongside) the fresh-from-prior K proposal. Targeted fix for
     * the PMMH-on-RR K-dwell observed at p=20+, p_inc=0.05 (mean K drifts
     * 2.4×−2.7× above the Geom prior; long streaks at K=5).
     *
     * Proposal: with probability mh_U_local_K_global_freq_ do the global
     * fresh-from-prior step (keeps escape route alive); otherwise propose
     * K_new ∈ {K_old−1, K_old+1} with reflection at 0. MH ratio includes
     * the geometric-prior ratio ρ^(K_new−K_old) and the reflection boundary
     * correction (±log 2 at K∈{0,1}).
     *
     * Default `false` preserves the fresh-from-prior dynamic.
     */
    void set_mh_U_local_K(bool on) { mh_U_local_K_ = on; }
    bool mh_U_local_K() const { return mh_U_local_K_; }

    /// Mixture fraction of fresh-from-prior K refreshes inside the local-K
    /// regime. 0.02 (= 1-in-50) by default. Used only when mh_U_local_K_.
    void set_mh_U_local_K_global_freq(double f) {
        mh_U_local_K_global_freq_ = f;
    }
    double mh_U_local_K_global_freq() const {
        return mh_U_local_K_global_freq_;
    }

    /** Shuffle edge visit order (random scan). */
    void prepare_iteration() override;

    /** V/RR U-pool refresh — gated by the chain runner on
     *  WarmupSchedule::u_refresh_enabled(iter). See cpp for details. */
    void refresh_auxiliary_u() override;

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

    // ---- Hierarchical-spec inference (Phase 4) ----
    // Default is Joint (the existing chain semantics). Hierarchical mode
    // multiplies the between-edge MH ratio by V(Γ_curr)/V(Γ_star), an
    // unbiased estimator of Z(Γ_star)/Z(Γ_curr), to convert the joint-
    // marginal Γ target π(Γ)·Z(Γ) into the user-specified π(Γ).
    GraphPriorSpec graph_prior_spec_ = GraphPriorSpec::Joint;
    bool   hierarchical_state_built_ = false;
    bool   prior_params_extracted_   = false;
    // Toggle for the manuscript App C NLO closed form (see
    // set_use_manuscript_nlo). Default `false` preserves the pre-2026-05-21
    // bgms formula and the SBC-clean baselines built on it.
    bool   use_manuscript_nlo_ = false;
    // Toggle for the MH-on-U fix at sweep boundary (see set_mh_U). Default
    // `false` preserves the legacy fresh-from-μ refresh.
    bool   mh_U_ = false;
    // Toggle for local random-walk on K_depth (see set_mh_U_local_K).
    // Active only when mh_U_ is also true. Default `false` keeps the
    // fresh-from-prior K proposal.
    bool   mh_U_local_K_ = false;
    double mh_U_local_K_global_freq_ = 0.02;
    // Plug-in mNLO mode (see set_plug_in_nlo). When true, the RR/U/K machinery
    // is fully bypassed and the between-Γ MH uses the closed-form NLO ratio
    // directly. Default `false` keeps the exact PMMH path.
    bool   plug_in_nlo_ = false;
    // Savage-Dickey between-step variant (handoff 2026-05-22). When true,
    // update_edge_indicators() dispatches to the SD between-step instead of
    // the existing joint/hierarchical-spec path. See set_use_sd_between_step.
    bool   use_sd_between_step_ = false;
    // L-space SD variant. Only consumed when use_sd_between_step_ is true.
    // See set_use_sd_lspace.
    bool   use_sd_lspace_ = false;
    // Disables the likelihood contribution in all MH ratios when true; used
    // to verify Γ-marginal stationarity against the prior. See set_prior_only.
    bool   prior_only_ = false;
    // Count of PD-revert events inside update_edge_indicator_parameter_pair_sd
    // _lspace. Diagnostic only.
    long   n_pd_reverts_ = 0;
    int    v_M_inner_ = 100;
    double v_kappa_   = 1.0;
    double v_rho_     = 0.5;
    // Per-Γ_curr state (only relevant when graph_prior_spec_ == Hierarchical):
    //   chain_aux_degord_       : (alpha, beta, sigma, delta) constants for DEGORD.
    //   log_Z_NLO_curr_         : analytic centring at Γ_curr (incremented on accept).
    //   v_pools_t_, v_K_depth_  : current U = (K_depth, K pools), refreshed per iteration.
    degord::ChainAux       chain_aux_degord_;
    double                 log_Z_NLO_curr_ = 0.0;
    int                    v_K_depth_     = 0;
    std::vector<arma::mat> v_pools_t_;
    // Extracted prior-family params (cached at ensure_hierarchical_state_).
    double prior_sigma_ = 1.0;  // NormalPrior scale
    double prior_alpha_ = 1.0;  // GammaScalePrior shape
    double prior_beta_  = 1.0;  // GammaScalePrior rate

    // Running V-ratio state for Lyne-style sign-corrected ergodic averaging
    // (F2). Updated inside update_edge_indicator_parameter_pair on accept;
    // chain runner snapshots into ChainResult at end of each sampling
    // iteration. v_diag_initialized_ stays false until the first V_log_pair
    // call produces a finite value, so chain-runner reads NaN for any
    // iteration that hits no toggle proposals (degenerate).
    int    current_sign_V_     = 1;
    double current_log_abs_V_  = std::numeric_limits<double>::quiet_NaN();
    bool   v_diag_initialized_ = false;
    // Endpoints of the toggle whose DEGORD permutation π_{i,j} was last used
    // to evaluate the cached V state (current_log_abs_V_, current_sign_V_).
    // mh_on_U_step_ reuses this permutation so V_old can be read from the
    // cache rather than recomputed. Initialised to (0, 1) — the canonical
    // permutation — so mh_U is well-defined even if no toggle has fired yet.
    int    last_v_pi_i_ = 0;
    int    last_v_pi_j_ = 1;

    // ---- Hierarchical-spec auto-reject counters ----
    // Counters for the between-Γ MH step's auto-reject sentinel. Incremented
    // inside update_edge_indicator_parameter_pair whenever the V/RR machinery
    // forces ln_alpha = -inf. Each proposal is classified into one of three
    // failure modes:
    //   *_nonfinite : non-finite log|V| at Γ_curr or Γ_star (V = 0 or RR
    //                 underflow / overflow).
    //   *_signzero  : sign(V) == 0 sentinel at Γ_curr or Γ_star.
    //   *_signflip  : both signs are finite and non-zero but differ — the
    //                 deferred Lyne sign-corrected weighting would handle
    //                 this; the current chain auto-rejects.
    // Counters survive across iterations (cumulative); read once at end of
    // chain via get_diagnostics_summary().
    mutable long long n_hier_add_attempts_   = 0;
    mutable long long n_hier_add_nonfinite_  = 0;
    mutable long long n_hier_add_signzero_   = 0;
    mutable long long n_hier_add_signflip_   = 0;
    mutable long long n_hier_del_attempts_   = 0;
    mutable long long n_hier_del_nonfinite_  = 0;
    mutable long long n_hier_del_signzero_   = 0;
    mutable long long n_hier_del_signflip_   = 0;

    // MH-on-U counters (set_mh_U fix). Each iteration's start adds one to
    // n_mh_U_attempts_ (after lazy init), plus exactly one of accepts or
    // the failure-mode counters.
    mutable long long n_mh_U_attempts_   = 0;
    mutable long long n_mh_U_accepts_    = 0;
    mutable long long n_mh_U_nonfinite_  = 0;
    mutable long long n_mh_U_signzero_   = 0;
    mutable long long n_mh_U_signflip_   = 0;

    // Local-K diagnostic counters (only incremented when mh_U_local_K_).
    // Tracking up/down separately surfaces the signature of a healthy
    // local-K kernel: acc_K_down > acc_K_up because the geometric prior
    // pulls the chain back toward small K.
    mutable long long n_mh_U_local_up_attempts_   = 0;
    mutable long long n_mh_U_local_up_accepts_    = 0;
    mutable long long n_mh_U_local_down_attempts_ = 0;
    mutable long long n_mh_U_local_down_accepts_  = 0;
    mutable long long n_mh_U_local_global_steps_  = 0;  // fresh-from-prior fraction

  public:
    /// @inheritdoc
    Rcpp::List get_diagnostics_summary() const override {
        return Rcpp::List::create(
            Rcpp::Named("hier_add_attempts")  = static_cast<double>(n_hier_add_attempts_),
            Rcpp::Named("hier_add_nonfinite") = static_cast<double>(n_hier_add_nonfinite_),
            Rcpp::Named("hier_add_signzero")  = static_cast<double>(n_hier_add_signzero_),
            Rcpp::Named("hier_add_signflip")  = static_cast<double>(n_hier_add_signflip_),
            Rcpp::Named("hier_del_attempts")  = static_cast<double>(n_hier_del_attempts_),
            Rcpp::Named("hier_del_nonfinite") = static_cast<double>(n_hier_del_nonfinite_),
            Rcpp::Named("hier_del_signzero")  = static_cast<double>(n_hier_del_signzero_),
            Rcpp::Named("hier_del_signflip")  = static_cast<double>(n_hier_del_signflip_),
            Rcpp::Named("mh_U_attempts")      = static_cast<double>(n_mh_U_attempts_),
            Rcpp::Named("mh_U_accepts")       = static_cast<double>(n_mh_U_accepts_),
            Rcpp::Named("mh_U_nonfinite")     = static_cast<double>(n_mh_U_nonfinite_),
            Rcpp::Named("mh_U_signzero")      = static_cast<double>(n_mh_U_signzero_),
            Rcpp::Named("mh_U_signflip")      = static_cast<double>(n_mh_U_signflip_),
            Rcpp::Named("mh_U_local_up_att")  = static_cast<double>(n_mh_U_local_up_attempts_),
            Rcpp::Named("mh_U_local_up_acc")  = static_cast<double>(n_mh_U_local_up_accepts_),
            Rcpp::Named("mh_U_local_down_att")= static_cast<double>(n_mh_U_local_down_attempts_),
            Rcpp::Named("mh_U_local_down_acc")= static_cast<double>(n_mh_U_local_down_accepts_),
            Rcpp::Named("mh_U_local_global")  = static_cast<double>(n_mh_U_local_global_steps_)
        );
    }

  private:

    /// Lazy initialiser for the V/RR machinery. Validates prior family,
    /// builds chain_aux_degord_, computes log_Z_NLO_curr_ via full-recompute,
    /// draws the first U-pool. Idempotent (no-op when state is fresh).
    void ensure_hierarchical_state_();
    /// Lightweight prior-family validator and parameter extractor. Sets
    /// prior_sigma_, prior_alpha_, prior_beta_ from the configured slab and
    /// diagonal priors. Used by paths that need the closed-form constants
    /// (e.g. the SD between-step) but not the V/RR state. Idempotent.
    void ensure_prior_params_extracted_();
    /// Draw a fresh (K_depth, pools_t) U for the V estimator.
    void refresh_z_ratio_pool_();
    /// MH step on (U, K_depth) at sweep boundary. Dispatches between the
    /// fresh-from-prior kernel (global proposal) and the local random-walk
    /// on K kernel (when mh_U_local_K_ is set), mixing them at the
    /// configured global frequency.
    void mh_on_U_step_();
    /// Fresh-from-prior MH proposal on (U, K_depth). Accepts on
    /// log|V_new| − log|V_old|; μ and P(N) cancel by proposal symmetry.
    void mh_on_U_step_global_();
    /// Local random-walk MH proposal on K_depth with reflection at 0; U is
    /// extended/shrunk in product-form fashion. Includes the geometric-prior
    /// ratio ρ^(K_new−K_old) and ±log 2 boundary corrections at K∈{0,1}.
    void mh_on_U_step_local_K_();

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
     * Savage-Dickey variant of the between-Γ MH step (handoff 2026-05-22).
     * Uses a marginalised-γ MH ratio with K_ij collapsed and the per-step
     * SD posterior density at K_ij = 0 evaluated via 1D Laplace+NLO at the
     * chain's current K_{-ij}. On ADD accept K_ij is drawn from the 1D
     * Laplace proposal; on DEL accept K_ij is set to zero. K_jj is
     * unchanged in this step.
     *
     * Activated by use_sd_between_step_ = true. Honours prior_only_:
     * when true, the BF reduces to 1 and the MH ratio is the prior odds
     * alone (used for Γ-marginal stationarity tests).
     *
     * @param i  Row index (i < j)
     * @param j  Column index
     */
    void update_edge_indicator_parameter_pair_sd(size_t i, size_t j);

    /**
     * L-space variant of the SD between-Γ step. Closed-form Gibbs at α = 1:
     * derives the trailing-2×2 Cholesky factor of the DEGORD-permuted K from
     * the 2×2 corner of Σ, computes the slave value m_ij and conditional
     * posterior moments (τ_post, μ_post) of l_ji in closed form, evaluates the
     * SD log-BF as a Gaussian density at m_ij minus log l_ii minus
     * log_slab_at_0, and on accept updates K_ij and K_jj jointly through
     * Δl_ji. PD is automatic; no truncation; no status=1 mode.
     *
     * Restricted to α = 1. Asserts diagonal prior shape parameter = 1.
     */
    void update_edge_indicator_parameter_pair_sd_lspace(size_t i, size_t j);

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
     * End-of-sweep drift check on covariance_matrix_.
     *
     * The SMW rank-1/rank-2 updates that replace the per-accept
     * O(p^3) refresh in cholesky_update_after_{edge,diag} are
     * backward stable but still accumulate FP error in
     * covariance_matrix_ over a long chain. This helper measures
     * max_i |sum_k cov(i,k) * K(k,i) - 1| -- the worst-case
     * diagonal entry of cov*K minus the identity -- and triggers
     * refresh_cholesky() when it exceeds kCovDriftTol_.
     * Computed in O(p^2) so the check fits comfortably inside the
     * MH sweep budget.
     */
    void check_and_refresh_if_drift_();

    /// Absolute tolerance on max|diag(cov*K) - 1| before refresh.
    /// Set conservatively; on a clean refresh this quantity is
    /// O(p * eps * cond(K)).
    static constexpr double kCovDriftTol_ = 1e-8;

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
