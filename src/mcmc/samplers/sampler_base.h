#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <memory>
#include <utility>
#include "mcmc/execution/step_result.h"
#include "mcmc/algorithms/hmc.h"
#include "mcmc/samplers/hmc_adaptation.h"
#include "models/base_model.h"

// ---------------------------------------------------------------------------
// SamplerBase — abstract interface for all MCMC samplers
// ---------------------------------------------------------------------------

/**
 * SamplerBase - Abstract base class for MCMC samplers
 *
 * Provides a unified interface for all MCMC sampling algorithms:
 * - MetropolisSampler (component-wise random-walk Metropolis)
 * - HMCSampler (Hamiltonian Monte Carlo)
 * - NUTSSampler (No-U-Turn Sampler)
 *
 * The sampler internally decides whether to adapt based on the iteration
 * number and its warmup schedule reference.
 */
class SamplerBase {
public:
    virtual ~SamplerBase() = default;

    /**
     * Perform one MCMC step
     *
     * The sampler internally decides whether to adapt based on the
     * iteration number and its warmup schedule reference.
     *
     * @param model      The model to sample from
     * @param iteration  Current iteration (0-based, spans warmup + sampling)
     * @return StepResult with new state and diagnostics
     */
    virtual StepResult step(BaseModel& model, int iteration) = 0;

    /**
     * Initialize the sampler before the MCMC loop.
     * For gradient-based samplers, runs the step-size heuristic. Default no-op.
     */
    virtual void initialize(BaseModel& /*model*/) {}

    /**
     * Check if this sampler produces NUTS-style diagnostics
     * (tree depth, divergences, energy)
     */
    virtual bool has_nuts_diagnostics() const { return false; }

    /**
     * Return the final adapted step size (0 for non-gradient samplers)
     */
    virtual double get_final_step_size() const { return 0.0; }
};

// ---------------------------------------------------------------------------
// GradientSamplerBase — shared base for gradient-based MCMC (HMC, NUTS)
// ---------------------------------------------------------------------------

/**
 * GradientSamplerBase - Base for gradient-based MCMC with warmup adaptation
 *
 * Uses HMCAdaptationController (from hmc_adaptation.h) with the shared
 * WarmupSchedule constructed by the runner.
 *
 * The adaptation controller handles:
 *  - Step-size dual averaging (Stages 1, 2, 3a, 3c)
 *  - Mass matrix estimation in doubling windows (Stage 2)
 *  - Step-size freezing at Stage 3b boundary
 */
class GradientSamplerBase : public SamplerBase {
public:
    GradientSamplerBase(double step_size, double target_acceptance,
                        WarmupSchedule& schedule,
                        bool force_nullspace = false)
        : step_size_(step_size),
          target_acceptance_(target_acceptance),
          use_dense_mass_(force_nullspace),
          schedule_(schedule),
          initialized_(false),
          force_nullspace_(force_nullspace)
    {}

    StepResult step(BaseModel& model, int iteration) override {
        // Stage 3c boundary: edge selection just activated.
        // Restart dual averaging so adaptation can tune to the new
        // geometry (changed active parameters) quickly.
        if (schedule_.in_stage3c(iteration) && !stage3c_initialized_) {
            stage3c_initialized_ = true;
            adapt_->reinit_stepsize(adapt_->current_step_size());
        }

        // Use adaptation controller's current step size for this iteration
        step_size_ = adapt_->current_step_size();

        StepResult result = do_gradient_step(model);

        // Let the adaptation controller handle step-size and mass-matrix logic.
        // For RATTLE (constrained) models, feed x-space samples so the mass
        // matrix is estimated in the same coordinate system NUTS operates in.
        arma::vec full_params = uses_constrained_integration(model)
            ? model.get_full_position()
            : model.get_full_vectorized_parameters();
        adapt_->update(full_params, result.accept_prob, iteration);

        // If mass matrix was just updated, apply it and re-run the step-size heuristic
        if (adapt_->mass_matrix_just_updated()) {
            arma::vec new_inv_mass = adapt_->inv_mass_diag();
            model.set_inv_mass(new_inv_mass);

            SafeRNG& rng = model.get_rng();

            if (uses_constrained_integration(model)) {
                arma::vec x = model.get_full_position();
                auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                    return model.logp_and_gradient_full(params).second;
                };
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient_full(params);
                };
                double new_eps = heuristic_initial_step_size(
                    x, grad_fn, joint_fn, new_inv_mass, rng,
                    0.625, adapt_->current_step_size());
                adapt_->reinit_stepsize(new_eps);
            } else if (use_dense_mass_ && adapt_->has_dense_covariance()) {
                // Dense mass: recompute L_full from new covariance, then L_active
                L_full_.reset();  // Force recomputation of L_full
                recompute_dense_transform(model);
                arma::vec theta = model.get_vectorized_parameters();
                arma::vec z = arma::solve(arma::trimatl(L_active_), theta);
                auto grad_fn_z = [this, &model](const arma::vec& z_param) -> arma::vec {
                    arma::vec theta_param = L_active_ * z_param;
                    return L_active_.t() * model.logp_and_gradient(theta_param).second;
                };
                auto joint_fn_z = [this, &model](const arma::vec& z_param)
                    -> std::pair<double, arma::vec> {
                    arma::vec theta_param = L_active_ * z_param;
                    auto [lp, g] = model.logp_and_gradient(theta_param);
                    return {lp, L_active_.t() * g};
                };
                arma::vec unit_mass = arma::ones<arma::vec>(z.n_elem);
                double new_eps = heuristic_initial_step_size(
                    z, grad_fn_z, joint_fn_z, unit_mass, rng,
                    0.625, adapt_->current_step_size());
                adapt_->reinit_stepsize(new_eps);
            } else {
                arma::vec theta = model.get_vectorized_parameters();
                auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                    return model.logp_and_gradient(params).second;
                };
                auto joint_fn = [&model](const arma::vec& params)
                    -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient(params);
                };
                arma::vec active_inv_mass = model.get_active_inv_mass();
                double new_eps = heuristic_initial_step_size(
                    theta, grad_fn, joint_fn, active_inv_mass, rng,
                    0.625, adapt_->current_step_size());
                adapt_->reinit_stepsize(new_eps);
            }
        }

        // Update step_size_ from controller (may have changed due to mass update)
        step_size_ = adapt_->current_step_size();

        return result;
    }

    double get_step_size() const { return step_size_; }
    double get_final_step_size() const override { return step_size_; }
    double get_averaged_step_size() const {
        return adapt_ ? adapt_->final_step_size() : step_size_;
    }
    const arma::vec& get_inv_mass() const { return adapt_->inv_mass_diag(); }

protected:
    virtual StepResult do_gradient_step(BaseModel& model) = 0;

    /** Whether to use unconstrained (null-space) integration even when
     *  the model reports constraints. Set by "nuts-nullspace" sampler type. */
    bool uses_constrained_integration(const BaseModel& model) const {
        return !force_nullspace_ && model.has_constraints();
    }

    double step_size_;
    double target_acceptance_;

    // Dense mass transform state
    bool use_dense_mass_;
    arma::mat L_active_;       ///< Lower Cholesky of projected covariance
    arma::mat L_full_;         ///< chol(Σ_full) from Stage 2 (computed once)

    /**
     * Recompute the dense mass transform for the current constraint structure.
     *
     * Uses the pre-stored L_full (Cholesky of full-space covariance) and
     * projects through the model's current null-space bases B:
     *   Σ_active = B^T L_full L_full^T B = (L_full^T B)^T (L_full^T B)
     * Then QR of (L_full^T B) gives the Cholesky of Σ_active as R^T.
     */
    void recompute_dense_transform(BaseModel& model) {
        if (L_full_.is_empty()) {
            // First call: compute L_full from the full covariance
            arma::mat full_cov = adapt_->dense_covariance();
            full_cov += 1e-8 * arma::eye(full_cov.n_rows, full_cov.n_cols);
            L_full_ = arma::chol(full_cov, "lower");
        }

        arma::mat B = model.get_projection_matrix();

        if (B.n_cols == B.n_rows) {
            // Square B (typically identity or permutation) — L_active = L_full
            if (B.n_cols == L_full_.n_cols) {
                L_active_ = L_full_;
                return;
            }
        }

        // C = L_full^T * B  (upper-tri × dense, full_dim × active_dim)
        arma::mat C = L_full_.t() * B;

        // QR of C: C = Q R, so C^T C = R^T R → L_active = R^T
        arma::mat Q, R;
        arma::qr_econ(Q, R, C);

        // Ensure positive diagonal (standard Cholesky sign convention)
        for (arma::uword i = 0; i < R.n_rows; ++i) {
            if (R(i, i) < 0) {
                R.row(i) *= -1;
            }
        }
        L_active_ = R.t();
    }

    /**
     * @return true when a valid dense transform is available.
     */
    bool has_dense_transform() const {
        return use_dense_mass_ && !L_active_.is_empty();
    }

public:
    /**
     * Initialize the adaptation controller and run the step-size heuristic.
     * Called by the runner before the MCMC loop.
     */
    void initialize(BaseModel& model) override {
        if (initialized_) return;
        do_initialize(model);
        initialized_ = true;
    }

private:
    void do_initialize(BaseModel& model) {
        int dim = static_cast<int>(model.full_parameter_dimension());
        SafeRNG& rng = model.get_rng();

        // Initialize inverse mass to ones
        arma::vec init_inv_mass = arma::ones<arma::vec>(dim);
        model.set_inv_mass(init_inv_mass);

        double init_eps;

        if (uses_constrained_integration(model)) {
            // Project initial position onto constraint manifold before
            // computing step size. The MLE initialization may violate
            // K_ij = 0 constraints for excluded edges.
            arma::vec x = model.get_full_position();
            arma::vec r_dummy = arma::zeros<arma::vec>(x.n_elem);
            model.project_position(x);
            model.project_momentum(r_dummy, x);
            model.set_full_position(x);

            x = model.get_full_position();
            auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                return model.logp_and_gradient_full(params).second;
            };
            auto joint_fn = [&model](const arma::vec& params)
                -> std::pair<double, arma::vec> {
                return model.logp_and_gradient_full(params);
            };
            init_eps = heuristic_initial_step_size(
                x, grad_fn, joint_fn, rng, target_acceptance_);
        } else {
            arma::vec theta = model.get_vectorized_parameters();
            auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                return model.logp_and_gradient(params).second;
            };
            auto joint_fn = [&model](const arma::vec& params)
                -> std::pair<double, arma::vec> {
                return model.logp_and_gradient(params);
            };
            init_eps = heuristic_initial_step_size(
                theta, grad_fn, joint_fn, rng, target_acceptance_);
        }

        step_size_ = init_eps;

        // Construct the adaptation controller with the shared schedule
        adapt_ = std::make_unique<HMCAdaptationController>(
            dim, init_eps, target_acceptance_, schedule_,
            /*learn_mass_matrix=*/true,
            /*learn_dense_mass=*/use_dense_mass_);
    }

    WarmupSchedule& schedule_;
    bool initialized_;
    bool stage3c_initialized_ = false;
    bool force_nullspace_;
    std::unique_ptr<HMCAdaptationController> adapt_;
};
