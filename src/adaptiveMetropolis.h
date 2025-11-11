#pragma once

#include <RcppArmadillo.h>
#include <stdexcept>

class AdaptiveProposal {

public:

    AdaptiveProposal(size_t num_params, size_t adaption_window = 50, double target_accept = 0.44) {
        proposal_sds_ = arma::vec(num_params, arma::fill::ones) * 0.25; // Initial SD
        // acceptance_counts_ = arma::ivec(num_params, arma::fill::zeros);
        // total_proposals_ = arma::ivec(num_params, arma::fill::zeros);
        adaptation_window_ = adaption_window;
        target_accept_ = target_accept;
    }

    double get_proposal_sd(size_t param_index) const {
        validate_index(param_index);
        return proposal_sds_[param_index];
    }

    void update_proposal_sd(size_t param_index, double alpha) {

        if (!adapting_) {
            return;
        }

        double current_sd = get_proposal_sd(param_index);
        double updated_sd = current_sd + std::pow(1.0 / iterations_, 0.6) * (alpha - target_accept_);
        // proposal_sds_[param_index] = std::min(20.0, std::max(1.0 / std::sqrt(n), updated_sd));
        proposal_sds_(param_index) = std::min(20.0, updated_sd);
    }

    void increment_iteration() {
        iterations_++;
        if (iterations_ >= adaptation_window_) {
            adapting_ = false;
        }
    }

private:
    arma::vec proposal_sds_;
    int iterations_ = 0,
        adaptation_window_;
    double target_accept_ = 0.44;
    bool adapting_ = true;

    void validate_index(size_t index) const {
        if (index >= proposal_sds_.n_elem) {
            throw std::out_of_range("Parameter index out of range");
        }
    }

};
