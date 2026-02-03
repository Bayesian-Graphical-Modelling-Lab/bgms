#pragma once

#include <memory>
#include "base_model.h"
#include "adaptiveMetropolis.h"
#include "rng/rng_utils.h"


class SkeletonVariables : public BaseModel {
public:

    // constructor from raw data
    SkeletonVariables(
            const arma::mat& observations,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection = true
    ) :
    {}

    // copy constructor
    SkeletonVariables(const SkeletonVariables& other)
        : BaseModel(other),
    {}

    std::unique_ptr<BaseModel> clone() const override {
        return std::make_unique<SkeletonVariables>(*this); // uses copy constructor
    }

    bool has_gradient()    const          { return false; }
    bool has_adaptive_mh() const override { return true; }

    double logp(const arma::vec& parameters) override {
        // Implement log probability computation
        return 0.0;
    }

    void do_one_mh_step() override;

    size_t parameter_dimension() const override {
        return dim_;
    }

    void set_seed(int seed) override {
        rng_ = SafeRNG(seed);
    }

    // arma::vec get_vectorized_parameters() override {
    //     // upper triangle of precision_matrix_
    //     size_t e = 0;
    //     for (size_t j = 0; j < p_; ++j) {
    //         for (size_t i = 0; i <= j; ++i) {
    //             vectorized_parameters_(e) = precision_matrix_(i, j);
    //             ++e;
    //         }
    //     }
    //     return vectorized_parameters_;
    // }

    // arma::ivec get_vectorized_indicator_parameters() override {
    //     // upper triangle of precision_matrix_
    //     size_t e = 0;
    //     for (size_t j = 0; j < p_; ++j) {
    //         for (size_t i = 0; i <= j; ++i) {
    //             vectorized_indicator_parameters_(e) = edge_indicators_(i, j);
    //             ++e;
    //         }
    //     }
    //     return vectorized_indicator_parameters_;
    // }


private:
    // data
    size_t n_ = 0;
    size_t p_ = 0;
    size_t dim_ = 0;


};