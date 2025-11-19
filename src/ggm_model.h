#pragma once

#include <memory>
#include "base_model.h"
#include "adaptiveMetropolis.h"
#include "rng_utils.h"

class GGMModel : public BaseModel {
public:

    GGMModel(
            const arma::mat& X,
            const arma::mat& prior_inclusion_prob,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection = true
    ) : n_(X.n_rows),
        p_(X.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(X.t() * X),
        prior_inclusion_prob_(prior_inclusion_prob),
        edge_selection_(edge_selection),
        proposal_(AdaptiveProposal(dim_, 500)),
        omega_(arma::eye<arma::mat>(p_, p_)),
        phi_(arma::eye<arma::mat>(p_, p_)),
        inv_phi_(arma::eye<arma::mat>(p_, p_)),
        inv_omega_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        omega_prop_(arma::mat(p_, p_, arma::fill::none)),
        constants_(6)
    {}

    GGMModel(const GGMModel& other)
        : BaseModel(other),
          dim_(other.dim_),
          suf_stat_(other.suf_stat_),
          n_(other.n_),
          p_(other.p_),
          prior_inclusion_prob_(other.prior_inclusion_prob_),
          edge_selection_(other.edge_selection_),
          omega_(other.omega_),
          phi_(other.phi_),
          inv_phi_(other.inv_phi_),
          inv_omega_(other.inv_omega_),
          edge_indicators_(other.edge_indicators_),
          vectorized_parameters_(other.vectorized_parameters_),
          vectorized_indicator_parameters_(other.vectorized_indicator_parameters_),
          proposal_(other.proposal_),
          rng_(other.rng_),
          omega_prop_(other.omega_prop_),
          constants_(other.constants_)
    {}

    //     // rng_ = SafeRNG(123);

    // }

    void set_adaptive_proposal(AdaptiveProposal proposal) {
        proposal_ = proposal;
    }

    bool has_gradient()    const          { return false; }
    bool has_adaptive_mh() const override { return true; }

    double logp(const arma::vec& parameters) override {
        // Implement log probability computation
        return 0.0;
    }

    // TODO: this can be done more efficiently, no need for the Cholesky!
    double log_density(const arma::mat& omega) const { return log_density_impl(omega,  arma::chol(omega)); };
    double log_density()                       const { return log_density_impl(omega_, phi_); }

    void do_one_mh_step() override;

    size_t parameter_dimension() const override {
        return dim_;
    }

    void set_seed(int seed) override {
        rng_ = SafeRNG(seed);
    }

    arma::vec get_vectorized_parameters() override {
        // upper triangle of omega_
        size_t e = 0;
        for (size_t j = 0; j < p_; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                vectorized_parameters_(e) = omega_(i, j);
                ++e;
            }
        }
        return vectorized_parameters_;
    }

    arma::ivec get_vectorized_indicator_parameters() override {
        // upper triangle of omega_
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
        return std::make_unique<GGMModel>(*this); // uses copy constructor
    }

private:
    // data
    size_t n_;
    size_t p_;
    size_t dim_;
    arma::mat suf_stat_;
    arma::mat prior_inclusion_prob_;
    bool edge_selection_;

    // parameters
    arma::mat omega_, phi_, inv_phi_, inv_omega_;
    arma::imat edge_indicators_;
    arma::vec vectorized_parameters_;
    arma::ivec vectorized_indicator_parameters_;


    AdaptiveProposal proposal_;
    SafeRNG rng_;

    // internal helper variables
    arma::mat omega_prop_;
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
