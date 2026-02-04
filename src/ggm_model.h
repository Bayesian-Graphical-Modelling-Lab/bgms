#pragma once

#include <memory>
#include "base_model.h"
#include "adaptiveMetropolis.h"
#include "rng/rng_utils.h"


class GaussianVariables : public BaseModel {
public:

    // constructor from raw data
    GaussianVariables(
            const arma::mat& observations,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection = true
    ) : n_(observations.n_rows),
        p_(observations.n_cols),
        // TODO: need to estimate the means! so + 1
        dim_((p_ * (p_ + 1)) / 2),
        // TODO: need to store sample means!
        suf_stat_(observations.t() * observations),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        proposal_(AdaptiveProposal(dim_, 500)),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),

        edge_indicators_(initial_edge_indicators),

        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none)),
        constants_(6)
    {}

    // constructor from sufficient statistics
    // TODO: needs to implement same TODOs as above constructor
    GaussianVariables(
            const int n,
            const arma::mat& suf_stat,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection = true
    ) : n_(n),
        p_(suf_stat.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(suf_stat),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        proposal_(AdaptiveProposal(dim_, 500)),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none)),
        constants_(6)
    {}

    // copy constructor
    GaussianVariables(const GaussianVariables& other)
        : BaseModel(other),
          dim_(other.dim_),
          suf_stat_(other.suf_stat_),
          n_(other.n_),
          p_(other.p_),
          inclusion_probability_(other.inclusion_probability_),
          edge_selection_(other.edge_selection_),
          precision_matrix_(other.precision_matrix_),
          cholesky_of_precision_(other.cholesky_of_precision_),
          inv_cholesky_of_precision_(other.inv_cholesky_of_precision_),
          covariance_matrix_(other.covariance_matrix_),
          edge_indicators_(other.edge_indicators_),
          vectorized_parameters_(other.vectorized_parameters_),
          vectorized_indicator_parameters_(other.vectorized_indicator_parameters_),
          proposal_(other.proposal_),
          rng_(other.rng_),
          precision_proposal_(other.precision_proposal_),
          constants_(other.constants_)
    {}

    //     // rng_ = SafeRNG(123);

    // }

    void set_adaptive_proposal(AdaptiveProposal proposal) {
        proposal_ = proposal;
    }

    bool has_gradient()    const          { return false; }
    bool has_adaptive_mh() const override { return true; }
    bool has_edge_selection() const override { return edge_selection_; }

    double logp(const arma::vec& parameters) override {
        // Implement log probability computation
        return 0.0;
    }

    // TODO: this can be done more efficiently, no need for the Cholesky!
    double log_likelihood(const arma::mat& omega) const { return log_density_impl(omega,  arma::chol(omega)); };
    double log_likelihood()                       const { return log_density_impl(precision_matrix_, cholesky_of_precision_); }

    void do_one_mh_step() override;

    size_t parameter_dimension() const override {
        return dim_;
    }

    // For GGM, full dimension is the same as parameter dimension (no edge selection filtering)
    size_t full_parameter_dimension() const override {
        return dim_;
    }

    void set_seed(int seed) override {
        rng_ = SafeRNG(seed);
    }

    arma::vec get_vectorized_parameters() const override {
        // upper triangle of precision_matrix_
        arma::vec result(dim_);
        size_t e = 0;
        for (size_t j = 0; j < p_; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                result(e) = precision_matrix_(i, j);
                ++e;
            }
        }
        return result;
    }

    // For GGM, full and active parameter vectors are the same
    arma::vec get_full_vectorized_parameters() const override {
        arma::vec result(dim_);
        size_t e = 0;
        for (size_t j = 0; j < p_; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                result(e) = precision_matrix_(i, j);
                ++e;
            }
        }
        return result;
    }

    arma::ivec get_vectorized_indicator_parameters() override {
        // upper triangle of precision_matrix_
        size_t e = 0;
        for (size_t j = 0; j < p_; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                vectorized_indicator_parameters_(e) = edge_indicators_(i, j);
                ++e;
            }
        }
        return vectorized_indicator_parameters_;
    }

    std::unique_ptr<BaseModel> clone() const override {
        return std::make_unique<GaussianVariables>(*this); // uses copy constructor
    }

private:
    // data
    size_t n_;
    size_t p_;
    size_t dim_;
    arma::mat suf_stat_;
    arma::mat inclusion_probability_;
    bool edge_selection_;

    // parameters
    arma::mat precision_matrix_, cholesky_of_precision_, inv_cholesky_of_precision_, covariance_matrix_;
    arma::imat edge_indicators_;
    arma::vec vectorized_parameters_;
    arma::ivec vectorized_indicator_parameters_;


    AdaptiveProposal proposal_;

    SafeRNG rng_;

    // internal helper variables
    arma::mat precision_proposal_;
    arma::vec constants_; // Phi_q1q, Phi_q1q1, c[1], c[2], c[3], c[4]

    arma::vec v1_ = {0, -1};
    arma::vec v2_ = {0, 0};
    arma::vec vf1_ = arma::zeros<arma::vec>(p_);
    arma::vec vf2_ = arma::zeros<arma::vec>(p_);
    arma::vec u1_ = arma::zeros<arma::vec>(p_);
    arma::vec u2_ = arma::zeros<arma::vec>(p_);

    // Parameter group updates with optimized likelihood evaluations
    void update_edge_parameter(size_t i, size_t j);
    void update_diagonal_parameter(size_t i);
    void update_edge_indicator_parameter_pair(size_t i, size_t j);

    // Helper methods
    void get_constants(size_t i, size_t j);
    double compute_inv_submatrix_i(const arma::mat& A, const size_t i, const size_t ii, const size_t jj) const;
    double R(const double x) const;

    double log_density_impl(const arma::mat& omega, const arma::mat& phi) const;
    double log_density_impl_edge(size_t i, size_t j) const;
    double log_density_impl_diag(size_t j) const;
    double get_log_det(arma::mat triangular_A) const;
    void   cholesky_update_after_edge(double omega_ij_old, double omega_jj_old, size_t i, size_t j);
    void   cholesky_update_after_diag(double omega_ii_old, size_t i);
    // double find_reasonable_step_size_edge(const arma::mat& omega, size_t i, size_t j);
    // double find_reasonable_step_size_diag(const arma::mat& omega, size_t i);
    // double edge_log_ratio(const arma::mat& omega, size_t i, size_t j, double proposal);
    // double diag_log_ratio(const arma::mat& omega, size_t i, double proposal);
};


GaussianVariables createGaussianVariablesFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection = true
);
