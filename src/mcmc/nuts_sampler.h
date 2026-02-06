#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <vector>
#include <memory>
#include "base_sampler.h"
#include "mcmc_utils.h"
#include "mcmc_nuts.h"
#include "mcmc_adaptation.h"
#include "sampler_config.h"
#include "../base_model.h"

class NUTSSampler : public BaseSampler {
public:
    explicit NUTSSampler(const SamplerConfig& config, int n_warmup = 1000)
        : step_size_(config.initial_step_size),
          target_acceptance_(config.target_acceptance),
          max_tree_depth_(config.max_tree_depth),
          no_warmup_(config.no_warmup),
          n_warmup_(n_warmup),
          warmup_iteration_(0),
          initialized_(false),
          step_adapter_(config.initial_step_size)
    {
        build_warmup_schedule(n_warmup);
    }

    SamplerResult warmup_step(BaseModel& model) override {
        if (!initialized_) {
            initialize(model);
            initialized_ = true;
        }

        SamplerResult result = do_nuts_step(model);

        // Adapt step size during all warmup phases
        step_adapter_.update(result.accept_prob, target_acceptance_);
        step_size_ = step_adapter_.current();

        // During Stage 2, accumulate samples for mass matrix estimation
        if (in_stage2()) {
            arma::vec full_params = model.get_full_vectorized_parameters();
            mass_accumulator_->update(full_params);

            if (at_window_end()) {
                // Stan convention: inv_mass = variance (high-variance params move more)
                inv_mass_ = mass_accumulator_->variance();
                mass_accumulator_->reset();

                // Push adapted mass matrix to model
                model.set_inv_mass(inv_mass_);

                arma::vec theta = model.get_vectorized_parameters();
                SafeRNG& rng = model.get_rng();
                auto log_post = [&model](const arma::vec& params) -> double {
                    return model.logp_and_gradient(params).first;
                };
                auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                    return model.logp_and_gradient(params).second;
                };
                arma::vec active_inv_mass = model.get_active_inv_mass();

                double new_eps = heuristic_initial_step_size(
                    theta, log_post, grad_fn, active_inv_mass, rng,
                    0.625, step_size_);
                step_size_ = new_eps;
                step_adapter_.restart(new_eps);
            }
        }

        warmup_iteration_++;
        return result;
    }

    void finalize_warmup() override {
        step_size_ = step_adapter_.averaged();
    }

    SamplerResult sample_step(BaseModel& model) override {
        return do_nuts_step(model);
    }

    bool has_nuts_diagnostics() const override { return true; }
    double get_step_size() const { return step_size_; }
    double get_averaged_step_size() const { return step_adapter_.averaged(); }
    const arma::vec& get_inv_mass() const { return inv_mass_; }

private:
    void build_warmup_schedule(int n_warmup) {
        stage1_end_ = static_cast<int>(0.075 * n_warmup);
        stage3_start_ = n_warmup - static_cast<int>(0.10 * n_warmup);

        window_ends_.clear();
        int cur = stage1_end_;
        int wsize = 25;

        while (cur < stage3_start_) {
            int win = std::min(wsize, stage3_start_ - cur);
            window_ends_.push_back(cur + win);
            cur += win;
            wsize = std::min(wsize * 2, stage3_start_ - cur);
        }
    }

    bool in_stage2() const {
        return warmup_iteration_ >= stage1_end_ && warmup_iteration_ < stage3_start_;
    }

    bool at_window_end() const {
        for (int end : window_ends_) {
            if (warmup_iteration_ + 1 == end) return true;
        }
        return false;
    }

    void initialize(BaseModel& model) {
        arma::vec theta = model.get_vectorized_parameters();
        SafeRNG& rng = model.get_rng();

        inv_mass_ = arma::ones<arma::vec>(model.full_parameter_dimension());
        model.set_inv_mass(inv_mass_);

        mass_accumulator_ = std::make_unique<DiagMassMatrixAccumulator>(
            static_cast<int>(model.full_parameter_dimension()));

        auto log_post = [&model](const arma::vec& params) -> double {
            return model.logp_and_gradient(params).first;
        };
        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };

        step_size_ = heuristic_initial_step_size(
            theta, log_post, grad_fn, rng, target_acceptance_);

        step_adapter_.restart(step_size_);
    }

    SamplerResult do_nuts_step(BaseModel& model) {
        arma::vec theta = model.get_vectorized_parameters();
        SafeRNG& rng = model.get_rng();

        auto joint_fn = [&model](const arma::vec& params)
            -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(params);
        };

        arma::vec active_inv_mass = model.get_active_inv_mass();

        SamplerResult result = nuts_sampler_joint(
            theta, step_size_, joint_fn,
            active_inv_mass, rng, max_tree_depth_
        );

        model.set_vectorized_parameters(result.state);
        return result;
    }

    // Configuration
    double step_size_;
    double target_acceptance_;
    int max_tree_depth_;
    int no_warmup_;
    int n_warmup_;

    // State tracking
    int warmup_iteration_;
    bool initialized_;

    // Step size adaptation
    DualAveraging step_adapter_;

    // Mass matrix adaptation
    arma::vec inv_mass_;
    std::unique_ptr<DiagMassMatrixAccumulator> mass_accumulator_;

    // Warmup schedule
    int stage1_end_;
    int stage3_start_;
    std::vector<int> window_ends_;
};
