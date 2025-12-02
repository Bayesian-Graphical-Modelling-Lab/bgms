#include "ggm_model.h"
#include "adaptiveMetropolis.h"
#include "rng_utils.h"
#include "cholupdate.h"

double GGMModel::compute_inv_submatrix_i(const arma::mat& A, const size_t i, const size_t ii, const size_t jj) const {
    return(A(ii, jj) - A(ii, i) * A(jj, i) / A(i, i));
}

void GGMModel::get_constants(size_t i, size_t j) {

    // TODO: helper function?
    double logdet_omega = get_log_det(phi_);

    double log_adj_omega_ii = logdet_omega + std::log(std::abs(inv_omega_(i, i)));
    double log_adj_omega_ij = logdet_omega + std::log(std::abs(inv_omega_(i, j)));
    double log_adj_omega_jj = logdet_omega + std::log(std::abs(inv_omega_(j, j)));

    double inv_omega_sub_j1j1 = compute_inv_submatrix_i(inv_omega_, i, j, j);
    double log_abs_inv_omega_sub_jj = log_adj_omega_ii + std::log(std::abs(inv_omega_sub_j1j1));
    double Phi_q1q  = (2 * std::signbit(inv_omega_(i, j)) - 1) * std::exp(
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
        return constants_[5] + std::pow((x - constants_[3]) / constants_[4], 2);
    }
}

double GGMModel::get_log_det(arma::mat triangular_A) const {
    // assume A is an (upper) triangular cholesky factor
    // returns the log determinant of A'A

    // TODO: should we just do
    // log_det(val, sign, trimatu(A))?
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

    // this is the log likelihood ratio, not the full log likelihood like GGMModel::log_density_impl

    double Ui2 = omega_(i, j) - omega_prop_(i, j);
    // only reached from R
    // if (omega_(j, j) == omega_prop_(j, j)) {
    //     k = i;
    //     i = j;
    //     j = k;
    //   }
    double Uj2 = (omega_(j, j) - omega_prop_(j, j)) / 2;


    // W <- matrix(c(0, 1, 1, 0), 2, 2)
    // U0 <- matrix(c(0, -1, Ui2, Uj2))
    // U <- matrix(0, nrow(aOmega), 2)
    // U[c(i, j), 1] <- c(0, -1)
    // U[c(i, j), 2] <- c(Ui2, Uj2)
    // aOmega_prop - (aOmega + U %*% W %*% t(U))
    // det(aOmega_prop) - det(aOmega + U %*% W %*% t(U))
    // det(aOmega_prop) - det(W + t(U) %*% inv_aOmega %*% U) * det(W) * det(aOmega)
    // below computes logdet(W + t(U) %*% inv_aOmega %*% U) directly (this is a 2x2 matrix)

    double cc11 = 0 + inv_omega_(j, j);
    double cc12 = 1 - (inv_omega_(i, j) * Ui2 + inv_omega_(j, j) * Uj2);
    double cc22 = 0 + Ui2 * Ui2 * inv_omega_(i, i) + 2 * Ui2 * Uj2 * inv_omega_(i, j) + Uj2 * Uj2 * inv_omega_(j, j);

    double logdet = std::log(std::abs(cc11 * cc22 - cc12 * cc12));
    // logdet - (logdet(aOmega_prop) - logdet(aOmega))

    double trace_prod = -2 * (suf_stat_(j, j) * Uj2 + suf_stat_(i, j) * Ui2);

    double log_likelihood_ratio = (n_ * logdet - trace_prod) / 2;
    return log_likelihood_ratio;

}

double GGMModel::log_density_impl_diag(size_t j) const {
    // same as above but for i == j, so Ui2 = 0
    double Uj2 = (omega_(j, j) - omega_prop_(j, j)) / 2;

    double cc11 = 0 + inv_omega_(j, j);
    double cc12 = 1 - inv_omega_(j, j) * Uj2;
    double cc22 = 0 + Uj2 * Uj2 * inv_omega_(j, j);

    double logdet = std::log(std::abs(cc11 * cc22 - cc12 * cc12));
    double trace_prod = -2 * suf_stat_(j, j) * Uj2;

    // This function uses the fact that the determinant doesn't change during edge updates.
    // double trace_prod = 0.0;
    // // TODO: we only need one of the two lines below, but it's not entirely clear which one
    // trace_prod +=     suf_stat_(j, j) * (omega_prop(j, j) - omega(j, j));
    // trace_prod +=     suf_stat_(i, i) * (omega_prop(i, i) - omega(i, i));
    // trace_prod += 2 * suf_stat_(i, j) * (omega_prop(i, j) - omega(i, j));
    // trace_prod - sum((aOmega_prop - aOmega) * SufStat)

    double log_likelihood_ratio = (n_ * logdet - trace_prod) / 2;
    return log_likelihood_ratio;

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

    // Rcpp::Rcout << "i: " << i << ", j: " << j <<
    //                ", proposed phi: " << phi_prop <<
    //                ", proposal_sd omega_ij: " << proposal_sd <<
    //                ", proposed omega_ij: " << omega_prop_q1q <<
    //                ", proposed omega_jj: " << omega_prop_qq << std::endl;
    // constants_.print(Rcpp::Rcout, "Constants:");
    // omega_prop_.print(Rcpp::Rcout, "Proposed omega:");

    // arma::vec eigval = eig_sym(omega_prop_);
    // if (arma::any(eigval <= 0)) {
    //     Rcpp::Rcout << "Warning: omega_prop_ is not positive definite for edge (" << i << ", " << j << ")" << std::endl;

    //     Rcpp::Rcout <<
    //             ", proposed phi: " << phi_prop <<
    //             ", proposal_sd omega_ij: " << proposal_sd <<
    //             ", proposed omega_ij: " << omega_prop_q1q <<
    //             ", proposed omega_jj: " << omega_prop_qq << std::endl;
    //     constants_.print(Rcpp::Rcout, "Constants:");
    //     omega_prop_.print(Rcpp::Rcout, "Proposed omega:");
    //     omega_.print(Rcpp::Rcout, "Current omega:");
    //     phi_.print(Rcpp::Rcout, "Phi:");
    //     inv_omega_.print(Rcpp::Rcout, "Inv(Omega):");

    // }

    // double ln_alpha = log_density(omega_prop_) - log_density();
    double ln_alpha = log_density_impl_edge(i, j);

    // {
    //     double ln_alpha_ref = log_density(omega_prop_) - log_density();
    //     if (std::abs(ln_alpha - ln_alpha_ref) > 1e-6) {
    //         Rcpp::Rcout << "Warning: log density implementations do not match for edge (" << i << ", " << j << ")" << std::endl;
    //         omega_.print(Rcpp::Rcout, "Current omega:");
    //         omega_prop_.print(Rcpp::Rcout, "Proposed omega:");
    //         Rcpp::Rcout << "ln_alpha: " << ln_alpha << ", ln_alpha_ref: " << ln_alpha_ref << std::endl;
    //     }
    // }

    ln_alpha += R::dcauchy(omega_prop_(i, j), 0.0, 2.5, true);
    ln_alpha -= R::dcauchy(omega_(i, j), 0.0, 2.5, true);

    if (std::log(runif(rng_)) < ln_alpha) {
        // accept proposal
        proposal_.increment_accepts(e);

        double omega_ij_old = omega_(i, j);
        double omega_jj_old = omega_(j, j);


        omega_(i, j) = omega_prop_q1q;
        omega_(j, i) = omega_prop_q1q;
        omega_(j, j) = omega_prop_qq;

        cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);

        // // TODO: preallocate?
        // // find v for low rank update
        // arma::vec v1 = {0, -1};
        // arma::vec v2 = {omega_ij - omega_prop_(i, j), (omega_jj - omega_prop_(j, j)) / 2};

        // arma::vec vf1 = arma::zeros<arma::vec>(p_);
        // arma::vec vf2 = arma::zeros<arma::vec>(p_);
        // vf1[i] = v1[0];
        // vf1[j] = v1[1];
        // vf2[i] = v2[0];
        // vf2[j] = v2[1];

        // // we now have
        // // aOmega_prop - (aOmega + vf1 %*% t(vf2) + vf2 %*% t(vf1))

        // arma::vec u1 = (vf1 + vf2) / sqrt(2);
        // arma::vec u2 = (vf1 - vf2) / sqrt(2);

        // // we now have
        // // omega_prop_ - (aOmega + u1 %*% t(u1) - u2 %*% t(u2))
        // // and also
        // // aOmega_prop - (aOmega + cbind(vf1, vf2) %*% matrix(c(0, 1, 1, 0), 2, 2) %*% t(cbind(vf1, vf2)))

        // // update phi (2x O(p^2))
        // cholesky_update(phi_, u1);
        // cholesky_downdate(phi_, u2);

        // // update inverse (2x O(p^2))
        // arma::inv(inv_phi_, arma::trimatu(phi_));
        // inv_omega_ = inv_phi_ * inv_phi_.t();

    }

    proposal_.update_proposal_sd(e);
}

void GGMModel::cholesky_update_after_edge(double omega_ij_old, double omega_jj_old, size_t i, size_t j)
{

    v2_[0] = omega_ij_old - omega_prop_(i, j);
    v2_[1] = (omega_jj_old - omega_prop_(j, j)) / 2;

    vf1_[i] = v1_[0];
    vf1_[j] = v1_[1];
    vf2_[i] = v2_[0];
    vf2_[j] = v2_[1];

    // we now have
    // aOmega_prop - (aOmega + vf1 %*% t(vf2) + vf2 %*% t(vf1))

    u1_ = (vf1_ + vf2_) / sqrt(2);
    u2_ = (vf1_ - vf2_) / sqrt(2);

    // we now have
    // omega_prop_ - (aOmega + u1 %*% t(u1) - u2 %*% t(u2))
    // and also
    // aOmega_prop - (aOmega + cbind(vf1, vf2) %*% matrix(c(0, 1, 1, 0), 2, 2) %*% t(cbind(vf1, vf2)))

    // update phi (2x O(p^2))
    cholesky_update(phi_, u1_);
    cholesky_downdate(phi_, u2_);

    // update inverse (2x O(p^2))
    arma::inv(inv_phi_, arma::trimatu(phi_));
    inv_omega_ = inv_phi_ * inv_phi_.t();

    // reset for next iteration
    vf1_[i] = 0.0;
    vf1_[j] = 0.0;
    vf2_[i] = 0.0;
    vf2_[j] = 0.0;

}

void GGMModel::update_diagonal_parameter(size_t i) {
    // Implementation of diagonal parameter update
    // 1-3) from before
    double logdet_omega = get_log_det(phi_);
    double logdet_omega_sub_ii = logdet_omega + std::log(inv_omega_(i, i));

    size_t e = i * (i + 1) / 2 + i; // parameter index in vectorized form
    double proposal_sd = proposal_.get_proposal_sd(e);

    double theta_curr = (logdet_omega - logdet_omega_sub_ii) / 2;
    double theta_prop = rnorm(rng_, theta_curr, proposal_sd);

    //4) Replace and rebuild omega
    omega_prop_ = omega_;
    omega_prop_(i, i) = omega_(i, i) - std::exp(theta_curr) * std::exp(theta_curr) + std::exp(theta_prop) * std::exp(theta_prop);

    // Rcpp::Rcout << "i: " << i <<
    //                ", current theta: " << theta_curr <<
    //                ", proposed theta: " << theta_prop <<
    //                ", proposal_sd: " << proposal_sd << std::endl;
    // omega_prop_.print(Rcpp::Rcout, "Proposed omega:");

    // 5) Acceptance ratio
    // double ln_alpha = log_density(omega_prop_) - log_density();
    double ln_alpha = log_density_impl_diag(i);
    // {
    //     double ln_alpha_ref = log_density(omega_prop_) - log_density();
    //     if (std::abs(ln_alpha - ln_alpha_ref) > 1e-6) {
    //         Rcpp::Rcout << "Warning: log density implementations do not match for diag (" << i << ", " << i << ")" << std::endl;
    //         // omega_.print(Rcpp::Rcout, "Current omega:");
    //         // omega_prop_.print(Rcpp::Rcout, "Proposed omega:");
    //         Rcpp::Rcout << "ln_alpha: " << ln_alpha << ", ln_alpha_ref: " << ln_alpha_ref << std::endl;
    //         Rcpp::Rcout << "1e4 * diff: " << 10000 * (ln_alpha - ln_alpha_ref) << std::endl;
    //     }
    // }

    ln_alpha += R::dgamma(exp(theta_prop), 1.0, 1.0, true);
    ln_alpha -= R::dgamma(exp(theta_curr), 1.0, 1.0, true);
    ln_alpha += theta_prop - theta_curr; // Jacobian adjustment ?

    if (std::log(runif(rng_)) < ln_alpha) {

        proposal_.increment_accepts(e);

        double omega_ii = omega_(i, i);
        omega_(i, i) = omega_prop_(i, i);

        cholesky_update_after_diag(omega_ii, i);

        // arma::vec u(p_, arma::fill::zeros);
        // double delta = omega_ii - omega_prop_(i, i);
        // bool s = delta > 0;
        // u(i) = std::sqrt(std::abs(delta));


        // if (s)
        //     cholesky_downdate(phi_, u);
        // else
        //     cholesky_update(phi_, u);

        // // update inverse (2x O(p^2))
        // arma::inv(inv_phi_, arma::trimatu(phi_));
        // inv_omega_ = inv_phi_ * inv_phi_.t();


    }

    proposal_.update_proposal_sd(e);
}

void GGMModel::cholesky_update_after_diag(double omega_ii_old, size_t i)
{

    double delta = omega_ii_old - omega_prop_(i, i);

    bool s = delta > 0;
    vf1_(i) = std::sqrt(std::abs(delta));

    if (s)
        cholesky_downdate(phi_, vf1_);
    else
        cholesky_update(phi_, vf1_);

    // update inverse (2x O(p^2))
    arma::inv(inv_phi_, arma::trimatu(phi_));
    inv_omega_ = inv_phi_ * inv_phi_.t();

    // reset for next iteration
    vf1_(i) = 0.0;
}


void GGMModel::update_edge_indicator_parameter_pair(size_t i, size_t j) {

    size_t e = i * (i + 1) / 2 + j; // parameter index in vectorized form
    double proposal_sd = proposal_.get_proposal_sd(e);

    if (edge_indicators_(i, j) == 1) {
        // Propose to turn OFF the edge
        omega_prop_ = omega_;
        omega_prop_(i, j) = 0.0;
        omega_prop_(j, i) = 0.0;

        // Update diagonal using R function with omega_ij = 0
        get_constants(i, j);
        omega_prop_(j, j) = R(0.0);

        // double ln_alpha = log_density(omega_prop_) - log_density();
        double ln_alpha = log_density_impl_edge(i, j);
        // {
        //     double ln_alpha_ref = log_density(omega_prop_) - log_density();
        //     if (std::abs(ln_alpha - ln_alpha_ref) > 1e-6) {
        //         Rcpp::Rcout << "Warning: log density implementations do not match for edge indicator (" << i << ", " << j << ")" << std::endl;
        //         omega_.print(Rcpp::Rcout, "Current omega:");
        //         omega_prop_.print(Rcpp::Rcout, "Proposed omega:");
        //         Rcpp::Rcout << "ln_alpha: " << ln_alpha << ", ln_alpha_ref: " << ln_alpha_ref << std::endl;
        //     }
        // }


        ln_alpha += std::log(1.0 - prior_inclusion_prob_(i, j)) - std::log(prior_inclusion_prob_(i, j));

        ln_alpha += R::dnorm(omega_(i, j) / constants_[4], 0.0, proposal_sd, true) - std::log(constants_[4]);
        ln_alpha -= R::dcauchy(omega_(i, j), 0.0, 2.5, true);

        if (std::log(runif(rng_)) < ln_alpha) {

            // Store old values for Cholesky update
            double omega_ij_old = omega_(i, j);
            double omega_jj_old = omega_(j, j);

            // Update omega
            omega_(i, j) = 0.0;
            omega_(j, i) = 0.0;
            omega_(j, j) = omega_prop_(j, j);

            // Update edge indicator
            edge_indicators_(i, j) = 0;
            edge_indicators_(j, i) = 0;

            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
            // // Cholesky update vectors
            // arma::vec v1 = {0, -1};
            // arma::vec v2 = {omega_ij_old - 0.0, (omega_jj_old - omega_(j, j)) / 2};

            // arma::vec vf1 = arma::zeros<arma::vec>(p_);
            // arma::vec vf2 = arma::zeros<arma::vec>(p_);
            // vf1[i] = v1[0];
            // vf1[j] = v1[1];
            // vf2[i] = v2[0];
            // vf2[j] = v2[1];

            // arma::vec u1 = (vf1 + vf2) / sqrt(2);
            // arma::vec u2 = (vf1 - vf2) / sqrt(2);

            // // Update Cholesky factor
            // cholesky_update(phi_, u1);
            // cholesky_downdate(phi_, u2);

            // // Update inverse
            // arma::inv(inv_phi_, arma::trimatu(phi_));
            // inv_omega_ = inv_phi_ * inv_phi_.t();
        }

    } else {
        // Propose to turn ON the edge
        double epsilon = rnorm(rng_, 0.0, proposal_sd);

        // Get constants for current state (with edge OFF)
        get_constants(i, j);
        double omega_prop_ij = constants_[4] * epsilon;
        double omega_prop_jj = R(omega_prop_ij);

        omega_prop_ = omega_;
        omega_prop_(i, j) = omega_prop_ij;
        omega_prop_(j, i) = omega_prop_ij;
        omega_prop_(j, j) = omega_prop_jj;

        // double ln_alpha = log_density(omega_prop_) - log_density();
        double ln_alpha = log_density_impl_edge(i, j);
        // {
        //     double ln_alpha_ref = log_density(omega_prop_) - log_density();
        //     if (std::abs(ln_alpha - ln_alpha_ref) > 1e-6) {
        //         Rcpp::Rcout << "Warning: log density implementations do not match for edge indicator (" << i << ", " << j << ")" << std::endl;
        //         omega_.print(Rcpp::Rcout, "Current omega:");
        //         omega_prop_.print(Rcpp::Rcout, "Proposed omega:");
        //         Rcpp::Rcout << "ln_alpha: " << ln_alpha << ", ln_alpha_ref: " << ln_alpha_ref << std::endl;
        //     }
        // }
        ln_alpha += std::log(prior_inclusion_prob_(i, j)) - std::log(1.0 - prior_inclusion_prob_(i, j));

        // Prior change: add slab (Cauchy prior)
        ln_alpha += R::dcauchy(omega_prop_ij, 0.0, 2.5, true);

        // Proposal term: proposed edge value given it was generated from truncated normal
        ln_alpha -= R::dnorm(omega_prop_ij / constants_[4], 0.0, proposal_sd, true) - std::log(constants_[4]);

        // TODO: this can be factored out?
        if (std::log(runif(rng_)) < ln_alpha) {
            // Accept: turn ON the edge
            proposal_.increment_accepts(e);

            // Store old values for Cholesky update
            double omega_ij_old = omega_(i, j);
            double omega_jj_old = omega_(j, j);

            // Update omega
            omega_(i, j) = omega_prop_ij;
            omega_(j, i) = omega_prop_ij;
            omega_(j, j) = omega_prop_jj;

            // Update edge indicator
            edge_indicators_(i, j) = 1;
            edge_indicators_(j, i) = 1;

            cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
            // // Cholesky update vectors
            // arma::vec v1 = {0, -1};
            // arma::vec v2 = {omega_ij_old - omega_(i, j), (omega_jj_old - omega_(j, j)) / 2};

            // arma::vec vf1 = arma::zeros<arma::vec>(p_);
            // arma::vec vf2 = arma::zeros<arma::vec>(p_);
            // vf1[i] = v1[0];
            // vf1[j] = v1[1];
            // vf2[i] = v2[0];
            // vf2[j] = v2[1];

            // arma::vec u1 = (vf1 + vf2) / sqrt(2);
            // arma::vec u2 = (vf1 - vf2) / sqrt(2);

            // // Update Cholesky factor
            // cholesky_update(phi_, u1);
            // cholesky_downdate(phi_, u2);

            // // Update inverse
            // arma::inv(inv_phi_, arma::trimatu(phi_));
            // inv_omega_ = inv_phi_ * inv_phi_.t();
        }
    }
}

void GGMModel::do_one_mh_step() {

    // Update off-diagonals (upper triangle)
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            // Rcpp::Rcout << "Updating edge parameter (" << i << ", " << j << ")" << std::endl;
            update_edge_parameter(i, j);
            // if (!arma:: approx_equal(omega_ * inv_omega_, arma::eye<arma::mat>(p_, p_), "absdiff", 1e-6)) {
            //     Rcpp::Rcout << "Warning: Omega * Inv(Omega) not equal to identity after updating edge (" << i << ", " << j << ")" << std::endl;
            //     (omega_ * inv_omega_).print(Rcpp::Rcout, "Omega * Inv(Omega):");
            //     (phi_ * inv_phi_).print(Rcpp::Rcout, "Phi * Inv(Phi):");
            // }
            // if (!arma:: approx_equal(omega_, phi_.t() * phi_, "absdiff", 1e-6)) {
            //     Rcpp::Rcout << "Warning: Omega not equal to Phi.t() * Phi after updating edge (" << i << ", " << j << ")" << std::endl;
            //     omega_.print(Rcpp::Rcout, "Omega:");
            //     phi_.print(Rcpp::Rcout, "Phi:");
            // }
            // if (!arma:: approx_equal(phi_ * inv_phi_, arma::eye<arma::mat>(p_, p_), "absdiff", 1e-6)) {
            //     Rcpp::Rcout << "Warning: Phi * Inv(Phi) not equal to identity after updating edge (" << i << ", " << j << ")" << std::endl;
            //     (omega_ * inv_omega_).print(Rcpp::Rcout, "Omega * Inv(Omega):");
            //     (phi_ * inv_phi_).print(Rcpp::Rcout, "Phi * Inv(Phi):");
            // }
            // if (!arma:: approx_equal(inv_omega_, inv_phi_ * inv_phi_.t(), "absdiff", 1e-6)) {
            //     Rcpp::Rcout << "Warning: Inv(Omega) not equal to Inv(Phi) * Inv(Phi).t() after updating edge (" << i << ", " << j << ")" << std::endl;
            //     omega_.print(Rcpp::Rcout, "Omega:");
            //     phi_.print(Rcpp::Rcout, "Phi:");
            //     inv_omega_.print(Rcpp::Rcout, "Inv(Omega):");
            //     inv_phi_.print(Rcpp::Rcout, "Inv(Phi):");
            // }
        }
    }

    // Update diagonals
    for (size_t i = 0; i < p_; ++i) {
        // Rcpp::Rcout << "Updating diagonal parameter " << i << std::endl;
        update_diagonal_parameter(i);

        // if (!arma:: approx_equal(omega_ * inv_omega_, arma::eye<arma::mat>(p_, p_), "absdiff", 1e-6)) {
        //     Rcpp::Rcout << "Warning: Omega * Inv(Omega) not equal to identity after updating diagonal " << i << std::endl;
        //     (omega_ * inv_omega_).print(Rcpp::Rcout, "Omega * Inv(Omega):");
        //     (phi_ * inv_phi_).print(Rcpp::Rcout, "Phi * Inv(Phi):");
        // }
        // if (!arma:: approx_equal(omega_, phi_.t() * phi_, "absdiff", 1e-6)) {
        //     Rcpp::Rcout << "Warning: Omega not equal to Phi.t() * Phi after updating diagonal " << i << std::endl;
        //     omega_.print(Rcpp::Rcout, "Omega:");
        //     phi_.print(Rcpp::Rcout, "Phi:");
        // }
        // if (!arma:: approx_equal(phi_ * inv_phi_, arma::eye<arma::mat>(p_, p_), "absdiff", 1e-6)) {
        //     Rcpp::Rcout << "Warning: Phi * Inv(Phi) not equal to identity after updating diagonal " << i << std::endl;
        //     (omega_ * inv_omega_).print(Rcpp::Rcout, "Omega * Inv(Omega):");
        //     (phi_ * inv_phi_).print(Rcpp::Rcout, "Phi * Inv(Phi):");
        // }
        // if (!arma:: approx_equal(inv_omega_, inv_phi_ * inv_phi_.t(), "absdiff", 1e-6)) {
        //     Rcpp::Rcout << "Warning: Inv(Omega) not equal to Inv(Phi) * Inv(Phi).t() after updating diagonal " << i << std::endl;
        //     omega_.print(Rcpp::Rcout, "Omega:");
        //     phi_.print(Rcpp::Rcout, "Phi:");
        //     inv_omega_.print(Rcpp::Rcout, "Inv(Omega):");
        //     inv_phi_.print(Rcpp::Rcout, "Inv(Phi):");
        // }
    }

    if (edge_selection_) {
        for (size_t i = 0; i < p_ - 1; ++i) {
            for (size_t j = i + 1; j < p_; ++j) {
                // Rcpp::Rcout << "Between model move for edge (" << i << ", " << j << ")" << std::endl;
                update_edge_indicator_parameter_pair(i, j);
                // if (!arma:: approx_equal(omega_ * inv_omega_, arma::eye<arma::mat>(p_, p_), "absdiff", 1e-6)) {
                //     Rcpp::Rcout << "Warning: Omega * Inv(Omega) not equal to identity after updating edge (" << i << ", " << j << ")" << std::endl;
                //     (omega_ * inv_omega_).print(Rcpp::Rcout, "Omega * Inv(Omega):");
                //     (phi_ * inv_phi_).print(Rcpp::Rcout, "Phi * Inv(Phi):");
                // }
                // if (!arma:: approx_equal(omega_, phi_.t() * phi_, "absdiff", 1e-6)) {
                //     Rcpp::Rcout << "Warning: Omega not equal to Phi.t() * Phi after updating edge (" << i << ", " << j << ")" << std::endl;
                //     omega_.print(Rcpp::Rcout, "Omega:");
                //     phi_.print(Rcpp::Rcout, "Phi:");
                // }
                // if (!arma:: approx_equal(phi_ * inv_phi_, arma::eye<arma::mat>(p_, p_), "absdiff", 1e-6)) {
                //     Rcpp::Rcout << "Warning: Phi * Inv(Phi) not equal to identity after updating edge (" << i << ", " << j << ")" << std::endl;
                //     (omega_ * inv_omega_).print(Rcpp::Rcout, "Omega * Inv(Omega):");
                //     (phi_ * inv_phi_).print(Rcpp::Rcout, "Phi * Inv(Phi):");
                // }
                // if (!arma:: approx_equal(inv_omega_, inv_phi_ * inv_phi_.t(), "absdiff", 1e-6)) {
                //     Rcpp::Rcout << "Warning: Inv(Omega) not equal to Inv(Phi) * Inv(Phi).t() after updating edge (" << i << ", " << j << ")" << std::endl;
                //     omega_.print(Rcpp::Rcout, "Omega:");
                //     phi_.print(Rcpp::Rcout, "Phi:");
                //     inv_omega_.print(Rcpp::Rcout, "Inv(Omega):");
                //     inv_phi_.print(Rcpp::Rcout, "Inv(Phi):");
                // }
            }
        }
    }

    // could also be called in the main MCMC loop
    proposal_.increment_iteration();
}


GGMModel createGGMFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection
) {

    if (inputFromR.containsElementNamed("n") && inputFromR.containsElementNamed("suf_stat")) {
        int n = Rcpp::as<int>(inputFromR["n"]);
        arma::mat suf_stat = Rcpp::as<arma::mat>(inputFromR["suf_stat"]);
        return GGMModel(
            n,
            suf_stat,
            prior_inclusion_prob,
            initial_edge_indicators,
            edge_selection
        );
    } else if (inputFromR.containsElementNamed("X")) {
        arma::mat X = Rcpp::as<arma::mat>(inputFromR["X"]);
        return GGMModel(
            X,
            prior_inclusion_prob,
            initial_edge_indicators,
            edge_selection
        );
    } else {
        throw std::invalid_argument("Input list must contain either 'X' or both 'n' and 'suf_stat'.");
    }

}
