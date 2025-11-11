#include "ggm_model.h"
#include "adaptiveMetropolis.h"
#include "rng_utils.h"
#include "cholupdate.h"

double GGMModel::compute_inv_submatrix_i(const arma::mat& A, const size_t i, const size_t ii, const size_t jj) const {
    return(A(ii, jj) - A(ii, i) * A(jj, i) / A(i, i));
}

void GGMModel::get_constants(size_t i, size_t j) {

    // TODO: helper function?
    double logdet_omega = 0.0;
    for (size_t i = 0; i < p_; i++) {
        logdet_omega += std::log(phi_(i, i));
    }

    double log_adj_omega_ii = logdet_omega + log(abs(inv_omega_(i, i)));
    double log_adj_omega_ij = logdet_omega + log(abs(inv_omega_(i, j)));
    double log_adj_omega_jj = logdet_omega + log(abs(inv_omega_(j, j)));

    double inv_omega_sub_j1j1 = compute_inv_submatrix_i(inv_omega_, i, j, j);
    double log_abs_inv_omega_sub_jj = log_adj_omega_ii + log(abs(inv_omega_sub_j1j1));

    double Phi_q1q  = (-1 * (2 * std::signbit(inv_omega_(i, j)) - 1)) * std::exp(
        (log_adj_omega_ij - (log_adj_omega_jj + log_abs_inv_omega_sub_jj) / 2)
    );
    double Phi_q1q1 = std::exp((log_adj_omega_jj - log_abs_inv_omega_sub_jj) / 2);

    constants_[1] = Phi_q1q;
    constants_[2] = Phi_q1q1;
    constants_[3] = omega_(i, j) - Phi_q1q * Phi_q1q1;
    constants_[4] = Phi_q1q1;
    constants_[5] = omega_(j, j) - Phi_q1q * Phi_q1q;
    constants_[6] = constants_[5] + constants_[3] * constants_[3] / (constants_[4] * constants_[4]);

}

double GGMModel::R(const double x) const {
      if (x == 0) {
        return constants_[6];
    } else {
        return constants_[3] + std::pow((x - constants_[3]) / constants_[4], 2);
    }
}

void GGMModel::update_edge_parameter(size_t i, size_t j) {

    if (edge_indicators_(i, j) == 0) {
        return; // Edge is not included; skip update
    }

    get_constants(i, j);
    double Phi_q1q  = constants_[1];
    double Phi_q1q1 = constants_[2];

    size_t e = i * (i + 1) / 2 + j; // parameter index in vectorized form
    double proposal_sd = proposal_.get_proposal_sd(e);

    double phi_prop       = rnorm(rng_, Phi_q1q, proposal_sd);
    double omega_prop_q1q = constants_[3] + constants_[4] * phi_prop;
    double omega_prop_qq  = R(omega_prop_q1q);

    // form full proposal matrix for Omega
    omega_prop_ = omega_; // TODO: needs to be a copy!
    omega_prop_(i, j) = omega_prop_q1q;
    omega_prop_(j, i) = omega_prop_q1q;
    omega_prop_(j, j) = omega_prop_qq;

    double ln_alpha = log_density(omega_prop_) - log_density();
    ln_alpha += R::dcauchy(omega_prop_(i, j), 0.0, 2.5, true);
    ln_alpha -= R::dcauchy(omega_(i, j), 0.0, 2.5, true);

    double u = runif(rng_);
    if (ln_alpha > log(u)) {
        // accept proposal

        double omega_ij = omega_(i, j);
        double omega_jj = omega_(j, j);

        omega_(i, j) = omega_prop_q1q;
        omega_(j, i) = omega_prop_q1q;
        omega_(j, j) = omega_prop_qq;

        // TODO: preallocate?
        // find v for low rank update
        arma::vec v1 = {0, -1};
        arma::vec v2 = {omega_ij - omega_prop_(i, j), (omega_jj - omega_prop_(j, j)) / 2};

        arma::vec vf1 = arma::zeros<arma::vec>(p_);
        arma::vec vf2 = arma::zeros<arma::vec>(p_);
        vf1[i] = v1[1];
        vf1[j] = v1[2];
        vf2[i] = v2[1];
        vf2[j] = v2[2];

        // we now have
        // aOmega_prop - (aOmega + vf1 %*% t(vf2) + vf2 %*% t(vf1))

        arma::vec u1 = (vf1 + vf2) / sqrt(2);
        arma::vec u2 = (vf1 - vf2) / sqrt(2);

        // we now have
        // omega_prop_ - (aOmega + u1 %*% t(u1) - u2 %*% t(u2))
        // and also
        // aOmega_prop - (aOmega + cbind(vf1, vf2) %*% matrix(c(0, 1, 1, 0), 2, 2) %*% t(cbind(vf1, vf2)))

        // update phi
        cholesky_update(phi_, u1);
        cholesky_downdate(phi_, u2);

        // update inverse
        inv_omega_ = phi_.t() * phi_;

    }

    double alpha = std::min(1.0, std::exp(ln_alpha));
    proposal_.update_proposal_sd(e, alpha);
}

void GGMModel::do_one_mh_step() {

    // Update off-diagonals (upper triangle)
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            Rcpp::Rcout << "Updating edge parameter (" << i << ", " << j << ")" << std::endl;
            update_edge_parameter(i, j);
        }
    }

    // Update diagonals
    for (size_t i = 0; i < p_; ++i) {
        Rcpp::Rcout << "Updating diagonal parameter " << i << std::endl;
        update_diagonal_parameter(i);
    }

    // if (edge_selection_) {
    //     for (size_t i = 0; i < p_ - 1; ++i) {
    //         for (size_t j = i + 1; j < p_; ++j) {
    //             update_edge_indicator_parameter_pair(i, j);
    //         }
    //     }
    // }
    proposal_.increment_iteration();
}

double GGMModel::log_density_impl(const arma::mat& omega, const arma::mat& phi) const {
    double logdet_omega = 0.0;
    for (size_t i = 0; i < p_; i++) {
        logdet_omega += std::log(phi(i, i));
    }

    // TODO: does this allocate?
    double trace_prod = arma::accu(omega % suf_stat_);

    double log_likelihood = n_ * (p_ * log(2 * arma::datum::pi) / 2 + logdet_omega / 2) - trace_prod / 2;

    return log_likelihood;
}


void GGMModel::update_diagonal_parameter(size_t i) {
    // Implementation of diagonal parameter update
    // 1-3) from before
    double logdet_omega = 0.0;
    for (size_t i = 0; i < p_; i++) {
        logdet_omega += std::log(phi_(i, i));
    }

    double logdet_omega_sub_ii = logdet_omega + std::log(inv_omega_(i, i));

    size_t e = i * (i + 1) / 2 + i; // parameter index in vectorized form
    double proposal_sd = proposal_.get_proposal_sd(e);

    double theta_curr = (logdet_omega - logdet_omega_sub_ii) / 2;
    double theta_prop = rnorm(rng_, theta_curr, proposal_sd);

    //4) Replace and rebuild omega
    omega_prop_ = omega_;
    omega_prop_(i, i) = omega_(i, i) - std::exp(theta_curr) * std::exp(theta_curr) + std::exp(theta_prop) * std::exp(theta_prop);

    // 5) Acceptance ratio
    double ln_alpha = log_density(omega_prop_) - log_density();
    ln_alpha += R::dgamma(exp(theta_prop), 1.0, 1.0, true);
    ln_alpha -= R::dgamma(exp(theta_curr), 1.0, 1.0, true);
    ln_alpha += theta_prop - theta_curr; // Jacobian adjustment ?


    if (log(runif(rng_)) < ln_alpha) {

        double omega_ii = omega_(i, i);

        arma::vec u(p_, arma::fill::zeros);
        double delta = omega_ii - omega_prop_(i, i);
        bool s = delta > 0;
        u(i) = sqrt(abs(delta));

        omega_(i, i) = omega_prop_(i, i);

        if (!s)
            cholesky_update(phi_, u);
        else
            cholesky_downdate(phi_, u);

        inv_omega_ = phi_.t() * phi_;

    }

    double alpha = std::min(1.0, std::exp(ln_alpha));
    proposal_.update_proposal_sd(e, alpha);
}

void GGMModel::update_edge_indicator_parameter_pair(size_t i, size_t j) {
    // Implementation of edge indicator parameter pair update
}
