#pragma once

#include <string>
#include <RcppArmadillo.h>

class ChainResultNew {

public:
    ChainResultNew() {}

    bool        error = false,
                userInterrupt = false;
    std::string error_msg;
    int         chain_id;

    arma::mat   samples;

    void reserve(const size_t param_dim, const size_t n_iter) {
        samples.set_size(param_dim, n_iter);
    }
    void store_sample(const size_t iter, const arma::vec& sample) {
        samples.col(iter) = sample;
    }

    // arma::imat indicator_samples;

    // other samples
    // arma::ivec treedepth_samples;
    // arma::ivec divergent_samples;
    // arma::vec energy_samples;
    // arma::imat allocation_samples;
};
