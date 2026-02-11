#pragma once

#include <array>
#include <memory>
#include "models/base_model.h"
#include "models/adaptive_metropolis.h"
#include "rng/rng_utils.h"


/**
 * GGMModel - Gaussian Graphical Model
 *
 * Bayesian inference on the precision matrix (inverse covariance) of a
 * multivariate Gaussian via element-wise Metropolis-Hastings. Edge
 * selection uses a spike-and-slab prior with Cauchy slab.
 *
 * The Cholesky factor of the precision matrix is maintained incrementally
 * through rank-1 updates/downdates after each element change.
 */
class GGMModel : public BaseModel {
public:

    // Construct from raw observations
    GGMModel(
            const arma::mat& observations,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection = true,
            const double pairwise_scale = 2.5
    ) : n_(observations.n_rows),
        p_(observations.n_cols),
        // TODO: we need to adjust the algorithm to also sample the means!
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(observations.t() * observations),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        pairwise_scale_(pairwise_scale),
        proposal_(AdaptiveProposal(dim_, 500)),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none))
    {}

    // Construct from sufficient statistics
    GGMModel(
            const int n,
            const arma::mat& suf_stat,
            const arma::mat& inclusion_probability,
            const arma::imat& initial_edge_indicators,
            const bool edge_selection = true,
            const double pairwise_scale = 2.5
    ) : n_(n),
        p_(suf_stat.n_cols),
        dim_((p_ * (p_ + 1)) / 2),
        suf_stat_(suf_stat),
        inclusion_probability_(inclusion_probability),
        edge_selection_(edge_selection),
        pairwise_scale_(pairwise_scale),
        proposal_(AdaptiveProposal(dim_, 500)),
        precision_matrix_(arma::eye<arma::mat>(p_, p_)),
        cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        inv_cholesky_of_precision_(arma::eye<arma::mat>(p_, p_)),
        covariance_matrix_(arma::eye<arma::mat>(p_, p_)),
        edge_indicators_(initial_edge_indicators),
        vectorized_parameters_(dim_),
        vectorized_indicator_parameters_(edge_selection_ ? dim_ : 0),
        precision_proposal_(arma::mat(p_, p_, arma::fill::none))
    {}

    GGMModel(const GGMModel& other)
        : BaseModel(other),
          dim_(other.dim_),
          suf_stat_(other.suf_stat_),
          n_(other.n_),
          p_(other.p_),
          inclusion_probability_(other.inclusion_probability_),
          edge_selection_(other.edge_selection_),
          pairwise_scale_(other.pairwise_scale_),
          precision_matrix_(other.precision_matrix_),
          cholesky_of_precision_(other.cholesky_of_precision_),
          inv_cholesky_of_precision_(other.inv_cholesky_of_precision_),
          covariance_matrix_(other.covariance_matrix_),
          edge_indicators_(other.edge_indicators_),
          vectorized_parameters_(other.vectorized_parameters_),
          vectorized_indicator_parameters_(other.vectorized_indicator_parameters_),
          proposal_(other.proposal_),
          rng_(other.rng_),
          precision_proposal_(other.precision_proposal_)
    {}


    void set_adaptive_proposal(AdaptiveProposal proposal) {
        proposal_ = proposal;
    }

    bool has_gradient()        const override { return false; }
    bool has_adaptive_mh()     const override { return true; }
    bool has_edge_selection()  const override { return edge_selection_; }

    void set_edge_selection_active(bool active) override {
        edge_selection_active_ = active;
    }

    void initialize_graph() override;

    // GGM handles edge indicator updates inside do_one_mh_step()
    void update_edge_indicators() override {}

    // GGM uses component-wise MH; logp is unused.
    double logp(const arma::vec& parameters) override { return 0.0; }

    double log_likelihood(const arma::mat& omega) const { return log_density_impl(omega,  arma::chol(omega)); };
    double log_likelihood()                       const { return log_density_impl(precision_matrix_, cholesky_of_precision_); }

    void do_one_mh_step() override;

    size_t parameter_dimension() const override { return dim_; }
    size_t full_parameter_dimension() const override { return dim_; }

    void set_seed(int seed) override {
        rng_ = SafeRNG(seed);
    }

    arma::vec get_vectorized_parameters() const override {
        return extract_upper_triangle();
    }

    arma::vec get_full_vectorized_parameters() const override {
        return extract_upper_triangle();
    }

    arma::ivec get_vectorized_indicator_parameters() override {
        size_t e = 0;
        for (size_t j = 0; j < p_; ++j) {
            for (size_t i = 0; i <= j; ++i) {
                vectorized_indicator_parameters_(e) = edge_indicators_(i, j);
                ++e;
            }
        }
        return vectorized_indicator_parameters_;
    }

    SafeRNG& get_rng() override { return rng_; }

    const arma::imat& get_edge_indicators() const override {
        return edge_indicators_;
    }

    arma::mat& get_inclusion_probability() override {
        return inclusion_probability_;
    }

    int get_num_variables() const override {
        return static_cast<int>(p_);
    }

    int get_num_pairwise() const override {
        return static_cast<int>(p_ * (p_ - 1) / 2);
    }

    std::unique_ptr<BaseModel> clone() const override {
        return std::make_unique<GGMModel>(*this);
    }

private:

    arma::vec extract_upper_triangle() const {
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

    // Data
    size_t n_;
    size_t p_;
    size_t dim_;
    arma::mat suf_stat_;
    arma::mat inclusion_probability_;
    bool edge_selection_;
    bool edge_selection_active_ = false;
    double pairwise_scale_;

    // Parameters
    arma::mat precision_matrix_, cholesky_of_precision_, inv_cholesky_of_precision_, covariance_matrix_;
    arma::imat edge_indicators_;
    arma::vec vectorized_parameters_;
    arma::ivec vectorized_indicator_parameters_;

    AdaptiveProposal proposal_;
    SafeRNG rng_;

    // Scratch space
    arma::mat precision_proposal_;

    // Workspace for conditional precision reparametrization.
    // [0] Phi_q1q, [1] Phi_q1q1, [2] omega_ij - Phi_q1q*Phi_q1q1,
    // [3] Phi_q1q1, [4] omega_jj - Phi_q1q^2, [5] constrained diagonal at x=0.
    std::array<double, 6> constants_{};

    // Work vectors for rank-2 Cholesky update.
    // A symmetric rank-2 update  A + vf1*vf2' + vf2*vf1'  is decomposed into
    // two rank-1 updates via  u1 = (vf1+vf2)/sqrt(2),  u2 = (vf1-vf2)/sqrt(2).
    arma::vec v1_ = {0, -1};
    arma::vec v2_ = {0, 0};
    arma::vec vf1_ = arma::zeros<arma::vec>(p_);
    arma::vec vf2_ = arma::zeros<arma::vec>(p_);
    arma::vec u1_ = arma::zeros<arma::vec>(p_);
    arma::vec u2_ = arma::zeros<arma::vec>(p_);

    // MH updates
    void update_edge_parameter(size_t i, size_t j);
    void update_diagonal_parameter(size_t i);
    void update_edge_indicator_parameter_pair(size_t i, size_t j);

    // Helpers
    void get_constants(size_t i, size_t j);
    double compute_inv_submatrix_i(const arma::mat& A, const size_t i, const size_t ii, const size_t jj) const;

    // Conditional precision constraint: returns the required diagonal
    // value omega_jj that keeps the precision matrix positive definite
    // after changing the off-diagonal element to x.
    double constrained_diagonal(const double x) const;

    double log_density_impl(const arma::mat& omega, const arma::mat& phi) const;
    double log_density_impl_edge(size_t i, size_t j) const;
    double log_density_impl_diag(size_t j) const;
    double get_log_det(arma::mat triangular_A) const;
    void   cholesky_update_after_edge(double omega_ij_old, double omega_jj_old, size_t i, size_t j);
    void   cholesky_update_after_diag(double omega_ii_old, size_t i);
};


GGMModel createGGMModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection = true,
    const double pairwise_scale = 2.5
);
