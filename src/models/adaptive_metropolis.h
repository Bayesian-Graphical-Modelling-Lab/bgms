#pragma once

#include <RcppArmadillo.h>
#include <stdexcept>

class AdaptiveProposal {

public:

    AdaptiveProposal(size_t num_params, size_t adaptation_window = 50, double target_accept = 0.44) {
        proposal_sds_ = arma::vec(num_params, arma::fill::ones) * 0.25; // Initial SD, need to tweak this somehow?
        acceptance_counts_ = arma::ivec(num_params, arma::fill::zeros);
        adaptation_window_ = adaptation_window;
        target_accept_ = target_accept;
    }

    double get_proposal_sd(size_t param_index) const {
        validate_index(param_index);
        return proposal_sds_[param_index];
    }

    void update_proposal_sd(size_t param_index) {

        if (!adapting_) {
            return;
        }

        double current_sd = get_proposal_sd(param_index);
        double observed_acceptance_probability = acceptance_counts_[param_index] / static_cast<double>(iterations_ + 1);
        double rm_weight = std::pow(iterations_, -decay_rate_);

        // Robbins-Monro update step
        double updated_sd = current_sd + (observed_acceptance_probability - target_accept_) * rm_weight;
        updated_sd = std::clamp(updated_sd, rm_lower_bound, rm_upper_bound);

        proposal_sds_(param_index) = updated_sd;
    }

    void increment_accepts(size_t param_index) {
        validate_index(param_index);
        acceptance_counts_[param_index]++;
    }

    void increment_iteration() {
        iterations_++;
        if (iterations_ >= adaptation_window_) {
            adapting_ = false;
        }
    }

private:
    arma::vec   proposal_sds_;
    arma::ivec  acceptance_counts_;
    int         iterations_ = 0,
                adaptation_window_;
    double      target_accept_ = 0.44,
                decay_rate_ = 0.75,
                rm_lower_bound = 0.001,
                rm_upper_bound = 2.0;
    bool        adapting_ = true;

    void validate_index(size_t index) const {
        if (index >= proposal_sds_.n_elem) {
            throw std::out_of_range("Parameter index out of range");
        }
    }

};
