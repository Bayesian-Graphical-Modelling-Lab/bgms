#pragma once

#include <string>
#include <RcppArmadillo.h>

class ChainResultNew {

public:
    ChainResultNew() {}

    bool error;
    std::string error_msg;
    int chain_id;
    bool userInterrupt;

    arma::mat samples;

    void reserve(size_t param_dim, size_t n_iter) {
        samples.set_size(param_dim, n_iter);
    }
    void store_sample(size_t iter, const arma::vec& sample) {
        samples.col(iter) = sample;
    }

    // arma::imat indicator_samples;

    // other samples
    // arma::ivec treedepth_samples;
    // arma::ivec divergent_samples;
    // arma::vec energy_samples;
    // arma::imat allocation_samples;
};
