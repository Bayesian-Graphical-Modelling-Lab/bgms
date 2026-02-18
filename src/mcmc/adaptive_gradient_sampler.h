#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include <vector>
#include <memory>
#include "mcmc/base_sampler.h"
#include "mcmc/mcmc_utils.h"
#include "mcmc/mcmc_adaptation.h"
#include "mcmc/sampler_config.h"
#include "models/base_model.h"

/**
 * AdaptiveGradientSampler - Base for gradient-based MCMC with warmup adaptation
 *
 * Shared warmup logic for NUTS and HMC:
 *  Stage 1: step-size adaptation only
 *  Stage 2: mass matrix estimation in doubling windows + step-size re-tuning
 *  Stage 3: final step-size averaging
 */
class AdaptiveGradientSampler : public BaseSampler {
public:
    AdaptiveGradientSampler(double step_size, double target_acceptance, int n_warmup)
        : step_size_(step_size),
          target_acceptance_(target_acceptance),
          n_warmup_(n_warmup),
          warmup_iteration_(0),
          initialized_(false),
          step_adapter_(step_size)
    {
        build_warmup_schedule(n_warmup);
    }

    SamplerResult warmup_step(BaseModel& model) override {
        if (!initialized_) {
            initialize(model);
            initialized_ = true;
        }

        SamplerResult result = do_gradient_step(model);

        step_adapter_.update(result.accept_prob, target_acceptance_);
        step_size_ = step_adapter_.current();

        if (in_stage2()) {
            arma::vec full_params = model.get_full_vectorized_parameters();
            mass_accumulator_->update(full_params);

            if (at_window_end()) {
                inv_mass_ = mass_accumulator_->variance();
                mass_accumulator_->reset();
                model.set_inv_mass(inv_mass_);

                arma::vec theta = model.get_vectorized_parameters();
                SafeRNG& rng = model.get_rng();
                auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
                    return model.logp_and_gradient(params).second;
                };
                auto joint_fn = [&model](const arma::vec& params) -> std::pair<double, arma::vec> {
                    return model.logp_and_gradient(params);
                };
                arma::vec active_inv_mass = model.get_active_inv_mass();

                double new_eps = heuristic_initial_step_size(
                    theta, grad_fn, joint_fn, active_inv_mass, rng,
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
        return do_gradient_step(model);
    }

    double get_step_size() const { return step_size_; }
    double get_averaged_step_size() const { return step_adapter_.averaged(); }
    const arma::vec& get_inv_mass() const { return inv_mass_; }

protected:
    virtual SamplerResult do_gradient_step(BaseModel& model) = 0;

    double step_size_;
    double target_acceptance_;

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

        auto grad_fn = [&model](const arma::vec& params) -> arma::vec {
            return model.logp_and_gradient(params).second;
        };
        auto joint_fn = [&model](const arma::vec& params) -> std::pair<double, arma::vec> {
            return model.logp_and_gradient(params);
        };

        step_size_ = heuristic_initial_step_size(
            theta, grad_fn, joint_fn, rng, target_acceptance_);

        step_adapter_.restart(step_size_);
    }

    int n_warmup_;
    int warmup_iteration_;
    bool initialized_;

    DualAveraging step_adapter_;
    arma::vec inv_mass_;
    std::unique_ptr<DiagMassMatrixAccumulator> mass_accumulator_;

    int stage1_end_;
    int stage3_start_;
    std::vector<int> window_ends_;
};
