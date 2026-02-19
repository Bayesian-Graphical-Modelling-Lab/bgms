#include "models/ggm/ggm_model.h"
#include "models/adaptive_metropolis.h"
#include "rng/rng_utils.h"
#include "models/ggm/cholupdate.h"

double GGMModel::compute_inv_submatrix_i(const arma::mat& A, const size_t i, const size_t ii, const size_t jj) const {
    return(A(ii, jj) - A(ii, i) * A(jj, i) / A(i, i));
}

void GGMModel::get_constants(size_t i, size_t j) {

    double logdet_omega = get_log_det(cholesky_of_precision_);

    double log_adj_omega_ii = logdet_omega + std::log(std::abs(covariance_matrix_(i, i)));
    double log_adj_omega_ij = logdet_omega + std::log(std::abs(covariance_matrix_(i, j)));
    double log_adj_omega_jj = logdet_omega + std::log(std::abs(covariance_matrix_(j, j)));

    double inv_omega_sub_j1j1 = compute_inv_submatrix_i(covariance_matrix_, i, j, j);
    double log_abs_inv_omega_sub_jj = log_adj_omega_ii + std::log(std::abs(inv_omega_sub_j1j1));
    double Phi_q1q  = (2 * std::signbit(covariance_matrix_(i, j)) - 1) * std::exp(
        (log_adj_omega_ij - (log_adj_omega_jj + log_abs_inv_omega_sub_jj) / 2)
    );
    double Phi_q1q1 = std::exp((log_adj_omega_jj - log_abs_inv_omega_sub_jj) / 2);

    constants_[0] = Phi_q1q;
    constants_[1] = Phi_q1q1;
    constants_[2] = precision_matrix_(i, j) - Phi_q1q * Phi_q1q1;
    constants_[3] = Phi_q1q1;
    constants_[4] = precision_matrix_(j, j) - Phi_q1q * Phi_q1q;
    constants_[5] = constants_[4] + constants_[2] * constants_[2] / (constants_[3] * constants_[3]);

}

double GGMModel::constrained_diagonal(const double x) const {
    if (x == 0) {
        return constants_[5];
    } else {
        return constants_[4] + std::pow((x - constants_[2]) / constants_[3], 2);
    }
}

double GGMModel::get_log_det(arma::mat triangular_A) const {
    // log-determinant of A'A where A is upper-triangular Cholesky factor
    return 2 * arma::accu(arma::log(triangular_A.diag()));
}

double GGMModel::log_density_impl(const arma::mat& omega, const arma::mat& phi) const {

    double logdet_omega = get_log_det(phi);
    // TODO: why not just dot(omega, suf_stat_)?
    double trace_prod = arma::accu(omega % suf_stat_);

    double log_likelihood = n_ * (p_ * log(2 * arma::datum::pi) / 2 + logdet_omega / 2) - trace_prod / 2;

    return log_likelihood;
}

double GGMModel::log_density_impl_edge(size_t i, size_t j) const {

    // Log-likelihood ratio (not the full log-likelihood)

    double Ui2 = precision_matrix_(i, j) - precision_proposal_(i, j);
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;

    double cc11 = 0 + covariance_matrix_(j, j);
    double cc12 = 1 - (covariance_matrix_(i, j) * Ui2 + covariance_matrix_(j, j) * Uj2);
    double cc22 = 0 + Ui2 * Ui2 * covariance_matrix_(i, i) + 2 * Ui2 * Uj2 * covariance_matrix_(i, j) + Uj2 * Uj2 * covariance_matrix_(j, j);

    double logdet = std::log(std::abs(cc11 * cc22 - cc12 * cc12));
    // logdet - (logdet(aOmega_prop) - logdet(aOmega))

    double trace_prod = -2 * (suf_stat_(j, j) * Uj2 + suf_stat_(i, j) * Ui2);

    double log_likelihood_ratio = (n_ * logdet - trace_prod) / 2;
    return log_likelihood_ratio;

}

double GGMModel::log_density_impl_diag(size_t j) const {
    // same as above but for i == j, so Ui2 = 0
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;

    double cc11 = 0 + covariance_matrix_(j, j);
    double cc12 = 1 - covariance_matrix_(j, j) * Uj2;
    double cc22 = 0 + Uj2 * Uj2 * covariance_matrix_(j, j);

    double logdet = std::log(std::abs(cc11 * cc22 - cc12 * cc12));
    double trace_prod = -2 * suf_stat_(j, j) * Uj2;

    double log_likelihood_ratio = (n_ * logdet - trace_prod) / 2;
    return log_likelihood_ratio;

}

void GGMModel::update_edge_parameter(size_t i, size_t j) {

    if (edge_indicators_(i, j) == 0) {
        return; // Edge is not included; skip update
    }

    get_constants(i, j);
    double Phi_q1q  = constants_[0];
    (void)constants_[1]; // Phi_q1q1 computed in get_constants but unused here

    size_t e = j * (j + 1) / 2 + i; // parameter index in vectorized form (column-major upper triangle)
    double proposal_sd = proposal_.get_proposal_sd(e);

    double phi_prop       = rnorm(rng_, Phi_q1q, proposal_sd);
    double omega_prop_q1q = constants_[2] + constants_[3] * phi_prop;
    double omega_prop_qq  = constrained_diagonal(omega_prop_q1q);

    // form full proposal matrix for Omega
    precision_proposal_ = precision_matrix_;
    precision_proposal_(i, j) = omega_prop_q1q;
    precision_proposal_(j, i) = omega_prop_q1q;
    precision_proposal_(j, j) = omega_prop_qq;

    double ln_alpha = log_density_impl_edge(i, j);

    ln_alpha += R::dcauchy(precision_proposal_(i, j), 0.0, pairwise_scale_, true);
    ln_alpha -= R::dcauchy(precision_matrix_(i, j), 0.0, pairwise_scale_, true);

    if (std::log(runif(rng_)) < ln_alpha) {
        // accept proposal
        proposal_.increment_accepts(e);

        double omega_ij_old = precision_matrix_(i, j);
        double omega_jj_old = precision_matrix_(j, j);


        precision_matrix_(i, j) = omega_prop_q1q;
        precision_matrix_(j, i) = omega_prop_q1q;
        precision_matrix_(j, j) = omega_prop_qq;

        cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);

    }

    proposal_.update_proposal_sd(e);
}

void GGMModel::cholesky_update_after_edge(double omega_ij_old, double omega_jj_old, size_t i, size_t j)
{

    v2_[0] = omega_ij_old - precision_proposal_(i, j);
    v2_[1] = (omega_jj_old - precision_proposal_(j, j)) / 2;

    vf1_[i] = v1_[0];
    vf1_[j] = v1_[1];
    vf2_[i] = v2_[0];
    vf2_[j] = v2_[1];

    // we now have
    // aOmega_prop - (aOmega + vf1 %*% t(vf2) + vf2 %*% t(vf1))

    u1_ = (vf1_ + vf2_) / sqrt(2);
    u2_ = (vf1_ - vf2_) / sqrt(2);

    // update phi (2x O(p^2))
    cholesky_update(cholesky_of_precision_, u1_);
    cholesky_downdate(cholesky_of_precision_, u2_);

    // update inverse (2x O(p^2))
    arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
    covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();

    // reset for next iteration
    vf1_[i] = 0.0;
    vf1_[j] = 0.0;
    vf2_[i] = 0.0;
    vf2_[j] = 0.0;

}

void GGMModel::update_diagonal_parameter(size_t i) {
    // Implementation of diagonal parameter update
    // 1-3) from before
    double logdet_omega = get_log_det(cholesky_of_precision_);
    double logdet_omega_sub_ii = logdet_omega + std::log(covariance_matrix_(i, i));

    size_t e = i * (i + 3) / 2; // parameter index in vectorized form (column-major upper triangle, i==j)
    double proposal_sd = proposal_.get_proposal_sd(e);

    double theta_curr = (logdet_omega - logdet_omega_sub_ii) / 2;
    double theta_prop = rnorm(rng_, theta_curr, proposal_sd);

    //4) Replace and rebuild omega
    precision_proposal_ = precision_matrix_;
    precision_proposal_(i, i) = precision_matrix_(i, i) - std::exp(theta_curr) * std::exp(theta_curr) + std::exp(theta_prop) * std::exp(theta_prop);

    double ln_alpha = log_density_impl_diag(i);

    ln_alpha += R::dgamma(exp(theta_prop), 1.0, 1.0, true);
    ln_alpha -= R::dgamma(exp(theta_curr), 1.0, 1.0, true);
    ln_alpha += theta_prop - theta_curr; // Jacobian adjustment ?

    if (std::log(runif(rng_)) < ln_alpha) {

        proposal_.increment_accepts(e);

        double omega_ii = precision_matrix_(i, i);
        precision_matrix_(i, i) = precision_proposal_(i, i);

        cholesky_update_after_diag(omega_ii, i);

    }

    proposal_.update_proposal_sd(e);
}

void GGMModel::cholesky_update_after_diag(double omega_ii_old, size_t i)
{

    double delta = omega_ii_old - precision_proposal_(i, i);

    bool s = delta > 0;
    vf1_(i) = std::sqrt(std::abs(delta));

    if (s)
        cholesky_downdate(cholesky_of_precision_, vf1_);
    else
        cholesky_update(cholesky_of_precision_, vf1_);

    // update inverse (2x O(p^2))
    arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
    covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();

    // reset for next iteration
    vf1_(i) = 0.0;
}


void GGMModel::update_edge_indicator_parameter_pair(size_t i, size_t j) {

    size_t e = j * (j + 1) / 2 + i; // parameter index in vectorized form (column-major upper triangle)
    double proposal_sd = proposal_.get_proposal_sd(e);

    if (edge_indicators_(i, j) == 1) {
        // Propose to turn OFF the edge
        precision_proposal_ = precision_matrix_;
        precision_proposal_(i, j) = 0.0;
        precision_proposal_(j, i) = 0.0;

        // Update diagonal to preserve positive-definiteness
        get_constants(i, j);
        precision_proposal_(j, j) = constrained_diagonal(0.0);

        // double ln_alpha = log_likelihood(precision_proposal_) - log_likelihood();
        double ln_alpha = log_density_impl_edge(i, j);
        // {
        //     double ln_alpha_ref = log_likelihood(precision_proposal_) - log_likelihood();
        //     if (std::abs(ln_alpha - ln_alpha_ref) > 1e-6) {
        //         Rcpp::Rcout << "Warning: log density implementations do not match for edge indicator (" << i << ", " << j << ")" << std::endl;
        //         precision_matrix_.print(Rcpp::Rcout, "Current omega:");
        //         precision_proposal_.print(Rcpp::Rcout, "Proposed omega:");
        //         Rcpp::Rcout << "ln_alpha: " << ln_alpha << ", ln_alpha_ref: " << ln_alpha_ref << std::endl;
        //     }
        // }


        ln_alpha += std::log(1.0 - inclusion_probability_(i, j)) - std::log(inclusion_probability_(i, j));

        ln_alpha += R::dnorm(precision_matrix_(i, j) / constants_[3], 0.0, proposal_sd, true) - std::log(constants_[3]);
        ln_alpha -= R::dcauchy(precision_matrix_(i, j), 0.0, pairwise_scale_, true);

        if (std::log(runif(rng_)) < ln_alpha) {

            // Store old values for Cholesky update
            double omega_ij_old = precision_matrix_(i, j);
            double omega_jj_old = precision_matrix_(j, j);

            // Update omega
            precision_matrix_(i, j) = 0.0;
            precision_matrix_(j, i) = 0.0;
            precision_matrix_(j, j) = precision_proposal_(j, j);

            // Update edge indicator
            edge_indicators_(i, j) = 0;
            edge_indicators_(j, i) = 0;

            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);

        }

    } else {
        // Propose to turn ON the edge
        double epsilon = rnorm(rng_, 0.0, proposal_sd);

        // Get constants for current state (with edge OFF)
        get_constants(i, j);
        double omega_prop_ij = constants_[3] * epsilon;
        double omega_prop_jj = constrained_diagonal(omega_prop_ij);

        precision_proposal_ = precision_matrix_;
        precision_proposal_(i, j) = omega_prop_ij;
        precision_proposal_(j, i) = omega_prop_ij;
        precision_proposal_(j, j) = omega_prop_jj;

        // double ln_alpha = log_likelihood(precision_proposal_) - log_likelihood();
        double ln_alpha = log_density_impl_edge(i, j);
        // {
        //     double ln_alpha_ref = log_likelihood(precision_proposal_) - log_likelihood();
        //     if (std::abs(ln_alpha - ln_alpha_ref) > 1e-6) {
        //         Rcpp::Rcout << "Warning: log density implementations do not match for edge indicator (" << i << ", " << j << ")" << std::endl;
        //         precision_matrix_.print(Rcpp::Rcout, "Current omega:");
        //         precision_proposal_.print(Rcpp::Rcout, "Proposed omega:");
        //         Rcpp::Rcout << "ln_alpha: " << ln_alpha << ", ln_alpha_ref: " << ln_alpha_ref << std::endl;
        //     }
        // }
        ln_alpha += std::log(inclusion_probability_(i, j)) - std::log(1.0 - inclusion_probability_(i, j));

        // Prior change: add slab (Cauchy prior)
        ln_alpha += R::dcauchy(omega_prop_ij, 0.0, pairwise_scale_, true);

        // Proposal term: proposed edge value given it was generated from truncated normal
        ln_alpha -= R::dnorm(omega_prop_ij / constants_[3], 0.0, proposal_sd, true) - std::log(constants_[3]);

        if (std::log(runif(rng_)) < ln_alpha) {
            // Accept: turn ON the edge
            proposal_.increment_accepts(e);

            // Store old values for Cholesky update
            double omega_ij_old = precision_matrix_(i, j);
            double omega_jj_old = precision_matrix_(j, j);

            // Update omega
            precision_matrix_(i, j) = omega_prop_ij;
            precision_matrix_(j, i) = omega_prop_ij;
            precision_matrix_(j, j) = omega_prop_jj;

            // Update edge indicator
            edge_indicators_(i, j) = 1;
            edge_indicators_(j, i) = 1;

            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);

        }
    }
}

void GGMModel::do_one_mh_step() {

    // Update off-diagonals (upper triangle)
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            update_edge_parameter(i, j);
        }
    }

    // Update diagonals
    for (size_t i = 0; i < p_; ++i) {
        update_diagonal_parameter(i);
    }

    if (edge_selection_active_) {
        for (size_t i = 0; i < p_ - 1; ++i) {
            for (size_t j = i + 1; j < p_; ++j) {
                update_edge_indicator_parameter_pair(i, j);
            }
        }
    }

    // could also be called in the main MCMC loop
    proposal_.increment_iteration();
}

void GGMModel::initialize_graph() {
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            double p = inclusion_probability_(i, j);
            int draw = (runif(rng_) < p) ? 1 : 0;
            edge_indicators_(i, j) = draw;
            edge_indicators_(j, i) = draw;
            if (!draw) {
                precision_proposal_ = precision_matrix_;
                precision_proposal_(i, j) = 0.0;
                precision_proposal_(j, i) = 0.0;
                get_constants(i, j);
                precision_proposal_(j, j) = constrained_diagonal(0.0);

                double omega_ij_old = precision_matrix_(i, j);
                double omega_jj_old = precision_matrix_(j, j);
                precision_matrix_(j, j) = precision_proposal_(j, j);
                precision_matrix_(i, j) = 0.0;
                precision_matrix_(j, i) = 0.0;
                cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
            }
        }
    }
}


GGMModel createGGMModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection,
    const double pairwise_scale
) {

    if (inputFromR.containsElementNamed("n") && inputFromR.containsElementNamed("suf_stat")) {
        int n = Rcpp::as<int>(inputFromR["n"]);
        arma::mat suf_stat = Rcpp::as<arma::mat>(inputFromR["suf_stat"]);
        return GGMModel(
            n,
            suf_stat,
            prior_inclusion_prob,
            initial_edge_indicators,
            edge_selection,
            pairwise_scale
        );
    } else if (inputFromR.containsElementNamed("X")) {
        arma::mat X = Rcpp::as<arma::mat>(inputFromR["X"]);
        return GGMModel(
            X,
            prior_inclusion_prob,
            initial_edge_indicators,
            edge_selection,
            pairwise_scale
        );
    } else {
        throw std::invalid_argument("Input list must contain either 'X' or both 'n' and 'suf_stat'.");
    }

}
