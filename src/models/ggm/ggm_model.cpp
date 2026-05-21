#include "models/ggm/ggm_model.h"
#include "rng/rng_utils.h"
#include "math/explog_macros.h"
#include "math/cholupdate.h"
#include "mcmc/execution/step_result.h"
#include "mcmc/execution/warmup_schedule.h"

// =====================================================================
// NUTS gradient support
// =====================================================================

void GGMModel::ensure_constraint_structure() {
    if (!constraint_dirty_) return;
    constraint_structure_.build(edge_indicators_);
    gradient_engine_.rebuild(constraint_structure_, n_, suf_stat_, *interaction_prior_, *diagonal_prior_, determinant_tilt_);
    constraint_dirty_ = false;
    theta_valid_ = false;
}

void GGMModel::recompute_theta() const {
    if (theta_valid_) return;

    // Build constraint structure (const-safe: structure is already built
    // by ensure_constraint_structure before any gradient call)
    const auto& cs = constraint_structure_;
    theta_.set_size(cs.active_dim);

    arma::mat Aq_buf;

    for (size_t q = 0; q < p_; ++q) {
        const auto& col = cs.columns[q];
        size_t offset = cs.theta_offsets[q];

        // psi_q = log(phi_qq)
        double psi_q = std::log(cholesky_of_precision_(q, q));
        theta_(offset + col.d_q) = psi_q;

        if (q == 0 || col.d_q == 0) continue;

        // Build A_q, compute null-space basis N_q via Givens QR
        arma::mat Q_tmp, R_tmp;
        arma::vec R_diag;
        std::vector<GivensRotation> rots_tmp;
        GGMGradientEngine::build_Aq(cholesky_of_precision_, col, q, Aq_buf);
        GGMGradientEngine::givens_qr(Aq_buf.t(), Q_tmp, R_tmp, R_diag, rots_tmp);
        arma::mat Nq = Q_tmp.cols(col.m_q, q - 1);

        // f_q = N_q^T x_q
        arma::vec x_q = cholesky_of_precision_.col(q).head(q);
        arma::vec f_q = Nq.t() * x_q;

        for (size_t k = 0; k < col.d_q; ++k) {
            theta_(offset + k) = f_q(k);
        }
    }

    theta_valid_ = true;
}

size_t GGMModel::parameter_dimension() const {
    // Lazy: if constraint structure hasn't been built, use full dimension
    if (constraint_dirty_) {
        return p_ + p_ * (p_ - 1) / 2;
    }
    return constraint_structure_.active_dim;
}

size_t GGMModel::full_parameter_dimension() const {
    return p_ + p_ * (p_ - 1) / 2;
}

arma::vec GGMModel::get_vectorized_parameters() const {
    // Ensure the constraint structure is built so we can compute theta
    if (constraint_dirty_) {
        // const_cast is safe: ensure_constraint_structure only modifies
        // the constraint cache, not the model state
        const_cast<GGMModel*>(this)->ensure_constraint_structure();
    }
    recompute_theta();
    return theta_;
}

arma::vec GGMModel::get_full_vectorized_parameters() const {
    if (constraint_dirty_) {
        const_cast<GGMModel*>(this)->ensure_constraint_structure();
    }
    recompute_theta();

    const auto& cs = constraint_structure_;
    arma::vec full(cs.full_dim, arma::fill::zeros);

    for (size_t q = 0; q < p_; ++q) {
        const auto& col = cs.columns[q];
        size_t active_offset = cs.theta_offsets[q];
        size_t full_offset = cs.full_theta_offsets[q];

        // Copy f_q entries into their matching slots in the full vector.
        // In the full vector, column q has q slots for off-diagonal + 1 for diagonal.
        // The included indices map to specific positions.
        for (size_t k = 0; k < col.d_q; ++k) {
            // The k-th included index maps to position included_indices[k] in
            // the column's off-diagonal block
            size_t full_pos = full_offset + col.included_indices[k];
            full(full_pos) = theta_(active_offset + k);
        }

        // psi_q is at the end of the column's block in both layouts
        full(cs.full_psi_offset(q)) = theta_(active_offset + col.d_q);
    }

    return full;
}

void GGMModel::set_vectorized_parameters(const arma::vec& parameters) {
    ensure_constraint_structure();

    // Run forward map: theta -> Phi -> K
    ForwardMapResult fm = gradient_engine_.forward_map(parameters);

    // Update internal state
    precision_matrix_ = fm.K;
    cholesky_of_precision_ = fm.Phi;
    bool ok = arma::solve(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_),
                          arma::eye(p_, p_), arma::solve_opts::fast);
    if (!ok) {
        refresh_cholesky();
    } else {
        covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    }

    // Cache theta
    theta_ = parameters;
    theta_valid_ = true;
}

std::pair<double, arma::vec> GGMModel::logp_and_gradient(
    const arma::vec& parameters)
{
    ensure_constraint_structure();
    return gradient_engine_.logp_and_gradient(parameters);
}


arma::vec GGMModel::get_active_inv_mass() const {
    if (constraint_dirty_) {
        const_cast<GGMModel*>(this)->ensure_constraint_structure();
    }

    const auto& cs = constraint_structure_;

    if (inv_mass_.n_elem == 0) {
        return arma::ones<arma::vec>(cs.active_dim);
    }

    // The theta-space (unconstrained) NUTS path is the only caller of this
    // function. get_full_vectorized_parameters() scatters each active theta
    // entry into a definite slot in the full-dim vector (included off-diag
    // slot for f_q[k]; full_psi_offset for psi_q; excluded slots stay 0).
    // Welford in NUTSAdaptationController therefore tracks the variance of
    // the theta entry at that slot directly, so the active inverse mass is
    // simply the included-slot subset — no N_q rotation needed.
    if (inv_mass_.n_elem == cs.full_dim) {
        arma::vec active(cs.active_dim);
        for (size_t q = 0; q < p_; ++q) {
            const auto& col = cs.columns[q];
            size_t active_offset = cs.theta_offsets[q];
            size_t full_offset = cs.full_theta_offsets[q];
            for (size_t k = 0; k < col.d_q; ++k) {
                active(active_offset + k) =
                    inv_mass_(full_offset + col.included_indices[k]);
            }
            active(cs.psi_offset(q)) = inv_mass_(cs.full_psi_offset(q));
        }
        return active;
    }

    // Fallback: return inv_mass_ as-is (dimensions should match active_dim)
    return inv_mass_;
}


void GGMModel::get_constants(size_t i, size_t j) {

    double logdet_omega = cholesky_helpers::get_log_det(cholesky_of_precision_);

    double log_adj_omega_ii = logdet_omega + MY_LOG(std::abs(covariance_matrix_(i, i)));
    double log_adj_omega_ij = logdet_omega + MY_LOG(std::abs(covariance_matrix_(i, j)));
    double log_adj_omega_jj = logdet_omega + MY_LOG(std::abs(covariance_matrix_(j, j)));

    double inv_omega_sub_j1j1 = cholesky_helpers::compute_inv_submatrix_i(covariance_matrix_, i, j, j);
    double log_abs_inv_omega_sub_jj = log_adj_omega_ii + MY_LOG(std::abs(inv_omega_sub_j1j1));
    double Phi_q1q  = (2 * std::signbit(covariance_matrix_(i, j)) - 1) * MY_EXP(
        (log_adj_omega_ij - (log_adj_omega_jj + log_abs_inv_omega_sub_jj) / 2)
    );
    double Phi_q1q1 = MY_EXP((log_adj_omega_jj - log_abs_inv_omega_sub_jj) / 2);

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

double GGMModel::log_density_impl(const arma::mat& omega, const arma::mat& phi) const {

    double logdet_omega = cholesky_helpers::get_log_det(phi);
    double trace_prod = arma::accu(omega % suf_stat_);

    double log_likelihood = n_ * (p_ * MY_LOG(2 * arma::datum::pi) / 2 + logdet_omega / 2) - trace_prod / 2;

    return log_likelihood;
}

double GGMModel::log_det_ratio_edge(size_t i, size_t j) const {
    // Rank-2 matrix-determinant lemma: log|K_prop| - log|K_curr| where K_prop
    // differs from K_curr at entries (i,j), (j,i), and (j,j). cc11, cc12, cc22
    // are the entries of the 2x2 update Gram matrix I + V^T Sigma U used by
    // the lemma; |det(...)| is its determinant.
    double Ui2 = precision_matrix_(i, j) - precision_proposal_(i, j);
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;

    double cc11 = covariance_matrix_(j, j);
    double cc12 = 1 - (covariance_matrix_(i, j) * Ui2
                       + covariance_matrix_(j, j) * Uj2);
    double cc22 = Ui2 * Ui2 * covariance_matrix_(i, i)
                + 2 * Ui2 * Uj2 * covariance_matrix_(i, j)
                + Uj2 * Uj2 * covariance_matrix_(j, j);

    return MY_LOG(std::abs(cc11 * cc22 - cc12 * cc12));
}

double GGMModel::log_det_ratio_diag(size_t j) const {
    // Rank-1 specialisation of log_det_ratio_edge (Ui2 = 0).
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;

    double cc11 = covariance_matrix_(j, j);
    double cc12 = 1 - covariance_matrix_(j, j) * Uj2;
    double cc22 = Uj2 * Uj2 * covariance_matrix_(j, j);

    return MY_LOG(std::abs(cc11 * cc22 - cc12 * cc12));
}

double GGMModel::log_density_impl_edge(size_t i, size_t j) const {
    // Log-likelihood ratio (not the full log-likelihood).
    double Ui2 = precision_matrix_(i, j) - precision_proposal_(i, j);
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;
    double logdet = log_det_ratio_edge(i, j);
    double trace_prod = -2 * (suf_stat_(j, j) * Uj2 + suf_stat_(i, j) * Ui2);
    return (n_ * logdet - trace_prod) / 2;
}

double GGMModel::log_density_impl_diag(size_t j) const {
    double Uj2 = (precision_matrix_(j, j) - precision_proposal_(j, j)) / 2;
    double logdet = log_det_ratio_diag(j);
    double trace_prod = -2 * suf_stat_(j, j) * Uj2;
    return (n_ * logdet - trace_prod) / 2;
}

double GGMModel::update_edge_parameter(size_t i, size_t j) {

    if (edge_indicators_(i, j) == 0) {
        return 0.0; // Edge is not included; skip update (AR irrelevant, masked out)
    }

    get_constants(i, j);
    double Phi_q1q  = constants_[0];
    (void)constants_[1]; // Phi_q1q1 computed in get_constants but unused here

    size_t e = j * (j + 1) / 2 + i; // parameter index in vectorized form (column-major upper triangle)
    double proposal_sd = proposal_sds_(e);

    double phi_prop       = rnorm(rng_, Phi_q1q, proposal_sd);
    double omega_prop_q1q = constants_[2] + constants_[3] * phi_prop;
    double omega_prop_qq  = constrained_diagonal(omega_prop_q1q);

    // form full proposal matrix for Omega
    precision_proposal_ = precision_matrix_;
    precision_proposal_(i, j) = omega_prop_q1q;
    precision_proposal_(j, i) = omega_prop_q1q;
    precision_proposal_(j, j) = omega_prop_qq;

    double ln_alpha = log_density_impl_edge(i, j);

    // Determinant-tilt prior: |K|^delta contributes
    //   delta * (log|K_prop| - log|K_curr|)
    // to the MH ratio. log_det_ratio_edge uses the rank-2 matrix-determinant
    // lemma in O(p) via the cached covariance, so this is essentially free.
    if (determinant_tilt_ != 0.0) {
        ln_alpha += determinant_tilt_ * log_det_ratio_edge(i, j);
    }

    // Interaction prior on K_yy_{ij} = -0.5 * Omega_{ij}
    ln_alpha += interaction_prior_->logp(-0.5 * precision_proposal_(i, j));
    ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j));

    // Gamma(shape, rate) prior on changed diagonal K_jj. The Roverato move
    // slaves K_jj = c_3 + phi_{q-1,q}^2, so K_jj moves with the off-diagonal
    // and its prior must be re-evaluated.
    ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(j, j));
    ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(j, j));

    if (MY_LOG(runif(rng_)) < ln_alpha) {
        double omega_ij_old = precision_matrix_(i, j);
        double omega_jj_old = precision_matrix_(j, j);

        precision_matrix_(i, j) = omega_prop_q1q;
        precision_matrix_(j, i) = omega_prop_q1q;
        precision_matrix_(j, j) = omega_prop_qq;

        cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
    }

    return std::min(1.0, std::exp(ln_alpha));
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

    // update inverse — fall back to full recomputation if rank-1
    // updates have caused numerical drift
    bool ok = arma::solve(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_),
                          arma::eye(p_, p_), arma::solve_opts::fast);
    if (!ok) {
        refresh_cholesky();
    } else {
        covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    }

    // reset for next iteration
    vf1_[i] = 0.0;
    vf1_[j] = 0.0;
    vf2_[i] = 0.0;
    vf2_[j] = 0.0;

}

double GGMModel::update_diagonal_parameter(size_t i) {
    double logdet_omega = cholesky_helpers::get_log_det(cholesky_of_precision_);
    double logdet_omega_sub_ii = logdet_omega + MY_LOG(covariance_matrix_(i, i));

    size_t e = i * (i + 3) / 2; // parameter index in vectorized form (column-major upper triangle, i==j)
    double proposal_sd = proposal_sds_(e);

    double theta_curr = (logdet_omega - logdet_omega_sub_ii) / 2;
    double theta_prop = rnorm(rng_, theta_curr, proposal_sd);

    precision_proposal_ = precision_matrix_;
    precision_proposal_(i, i) = precision_matrix_(i, i) - MY_EXP(theta_curr) * MY_EXP(theta_curr) + MY_EXP(theta_prop) * MY_EXP(theta_prop);

    double ln_alpha = log_density_impl_diag(i);

    // Determinant-tilt prior: |K|^delta contributes delta * log_det_ratio
    // to the MH ratio. Rank-1 update => O(1) via the cached covariance.
    if (determinant_tilt_ != 0.0) {
        ln_alpha += determinant_tilt_ * log_det_ratio_diag(i);
    }

    ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(i, i));
    ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(i, i));
    ln_alpha += 2.0 * (theta_prop - theta_curr); // Jacobian: dK_ii/dtheta = 2*exp(2*theta)

    if (MY_LOG(runif(rng_)) < ln_alpha) {
        double omega_ii = precision_matrix_(i, i);
        precision_matrix_(i, i) = precision_proposal_(i, i);
        cholesky_update_after_diag(omega_ii, i);
    }

    return std::min(1.0, std::exp(ln_alpha));
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

    // update inverse — fall back to full recomputation if rank-1
    // updates have caused numerical drift
    bool ok = arma::solve(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_),
                          arma::eye(p_, p_), arma::solve_opts::fast);
    if (!ok) {
        refresh_cholesky();
    } else {
        covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    }

    // reset for next iteration
    vf1_(i) = 0.0;
}


void GGMModel::update_edge_indicator_parameter_pair(size_t i, size_t j) {

    size_t e = j * (j + 1) / 2 + i; // parameter index in vectorized form (column-major upper triangle)
    double proposal_sd = proposal_sds_(e);

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

        // Determinant-tilt prior: |K|^delta contributes delta * log_det_ratio
        // to the MH ratio. The rank-2 update at (i,j),(j,j) makes this O(p).
        if (determinant_tilt_ != 0.0) {
            ln_alpha += determinant_tilt_ * log_det_ratio_edge(i, j);
        }

        ln_alpha += MY_LOG(1.0 - inclusion_probability_(i, j)) - MY_LOG(inclusion_probability_(i, j));

        ln_alpha += R::dnorm(precision_matrix_(i, j) / constants_[3], 0.0, proposal_sd, true) - MY_LOG(constants_[3]);
        // Slab in K_yy coords; proposal in K_ij coords. Jacobian |dK_yy/dK_ij| = 1/2.
        ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j)) - MY_LOG(2.0);

        // Gamma(shape, rate) prior on changed diagonal K_jj. The Roverato move
        // slaves K_jj = c_3 + phi_{q-1,q}^2, so K_jj moves with the off-diagonal
        // and its prior must be re-evaluated.
        ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(j, j));
        ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(j, j));

        if (MY_LOG(runif(rng_)) < ln_alpha) {

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

            constraint_dirty_ = true;
            theta_valid_ = false;
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

        // Determinant-tilt prior: |K|^delta contributes delta * log_det_ratio
        // to the MH ratio.
        if (determinant_tilt_ != 0.0) {
            ln_alpha += determinant_tilt_ * log_det_ratio_edge(i, j);
        }

        ln_alpha += MY_LOG(inclusion_probability_(i, j)) - MY_LOG(1.0 - inclusion_probability_(i, j));

        // Slab in K_yy coords; proposal in K_ij coords. Jacobian |dK_yy/dK_ij| = 1/2.
        ln_alpha += interaction_prior_->logp(-0.5 * omega_prop_ij) - MY_LOG(2.0);

        // Gamma(shape, rate) prior on changed diagonal K_jj. The Roverato move
        // slaves K_jj = c_3 + phi_{q-1,q}^2, so K_jj moves with the off-diagonal
        // and its prior must be re-evaluated.
        ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(j, j));
        ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(j, j));

        // Proposal term: proposed edge value given it was generated from truncated normal
        ln_alpha -= R::dnorm(omega_prop_ij / constants_[3], 0.0, proposal_sd, true) - MY_LOG(constants_[3]);

        if (MY_LOG(runif(rng_)) < ln_alpha) {
            // Accept: turn ON the edge
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

            constraint_dirty_ = true;
            theta_valid_ = false;
        }
    }
}

void GGMModel::do_one_metropolis_step(int iteration) {
    // Collect per-slot accept probabilities for the Robbins-Monro adapter.
    // proposal_sds_ is stored as a flat dim_-length vec indexed by the
    // upper-triangle scheme `e = j * (j + 1) / 2 + i`; we mirror that here
    // as a dim_ x 1 matrix.
    arma::mat accept_prob(dim_, 1, arma::fill::zeros);
    arma::umat index_mask(dim_, 1, arma::fill::zeros);

    // Update off-diagonals (upper triangle)
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            double ap = update_edge_parameter(i, j);
            if (edge_indicators_(i, j) == 1) {
                size_t e = j * (j + 1) / 2 + i;
                accept_prob(e, 0) = ap;
                index_mask(e, 0) = 1;
            }
        }
    }

    // Update diagonals
    for (size_t i = 0; i < p_; ++i) {
        double ap = update_diagonal_parameter(i);
        size_t e = i * (i + 3) / 2;
        accept_prob(e, 0) = ap;
        index_mask(e, 0) = 1;
    }

    if (metropolis_adapter_) {
        metropolis_adapter_->update(index_mask, accept_prob, iteration);
    }
}

void GGMModel::init_metropolis_adaptation(const WarmupSchedule& schedule) {
    metropolis_adapter_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sds_, schedule, target_accept_);
}

void GGMModel::prepare_iteration() {
    // Shuffle edge visit order for random-scan edge selection.
    // Called unconditionally to keep RNG state consistent.
    shuffled_edge_order_ = arma_randperm(rng_, num_pairwise_);
}

void GGMModel::update_edge_indicators() {
    for (size_t idx = 0; idx < num_pairwise_; ++idx) {
        size_t flat = shuffled_edge_order_(idx);
        // Convert flat index to (i, j) upper-triangle pair.
        // flat = 0..(num_pairwise_-1), row-major: (0,1),(0,2),...,(0,p-1),(1,2),...
        size_t i = 0, j = 0;
        size_t acc = 0;
        for (size_t row = 0; row < p_ - 1; ++row) {
            size_t cols_in_row = p_ - 1 - row;
            if (flat < acc + cols_in_row) {
                i = row;
                j = row + 1 + (flat - acc);
                break;
            }
            acc += cols_in_row;
        }
        update_edge_indicator_parameter_pair(i, j);
    }
}

void GGMModel::tune_proposal_sd(int iteration, const WarmupSchedule& schedule) {
    auto rm_weight_opt = schedule.rm_weight_for_proposal_sd(iteration);
    if (!rm_weight_opt) return;
    const double rm_weight = *rm_weight_opt;
    const double target_accept = target_accept_;

    // Off-diagonal sweeps
    for (size_t i = 0; i < p_ - 1; ++i) {
        for (size_t j = i + 1; j < p_; ++j) {
            if (edge_indicators_(i, j) == 0) continue;

            get_constants(i, j);
            double Phi_q1q = constants_[0];
            size_t e = j * (j + 1) / 2 + i;
            double proposal_sd = proposal_sds_(e);

            double phi_prop = rnorm(rng_, Phi_q1q, proposal_sd);
            double omega_prop_q1q = constants_[2] + constants_[3] * phi_prop;
            double omega_prop_qq = constrained_diagonal(omega_prop_q1q);

            precision_proposal_ = precision_matrix_;
            precision_proposal_(i, j) = omega_prop_q1q;
            precision_proposal_(j, i) = omega_prop_q1q;
            precision_proposal_(j, j) = omega_prop_qq;

            double ln_alpha = log_density_impl_edge(i, j);
            if (determinant_tilt_ != 0.0) {
                ln_alpha += determinant_tilt_ * log_det_ratio_edge(i, j);
            }
            // Interaction prior on K_yy_{ij} = -0.5 * Omega_{ij}
            ln_alpha += interaction_prior_->logp(-0.5 * precision_proposal_(i, j));
            ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j));

            // Gamma(shape, rate) prior on changed diagonal K_jj. The Roverato move
            // slaves K_jj = c_3 + phi_{q-1,q}^2, so K_jj moves with the off-diagonal
            // and its prior must be re-evaluated.
            ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(j, j));
            ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(j, j));

            if (MY_LOG(runif(rng_)) < ln_alpha) {
                double omega_ij_old = precision_matrix_(i, j);
                double omega_jj_old = precision_matrix_(j, j);
                precision_matrix_(i, j) = omega_prop_q1q;
                precision_matrix_(j, i) = omega_prop_q1q;
                precision_matrix_(j, j) = omega_prop_qq;
                cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);
            }

            proposal_sds_(e) = update_proposal_sd_with_robbins_monro(
                proposal_sds_(e), ln_alpha, rm_weight, target_accept);
        }
    }

    // Diagonal sweeps
    for (size_t i = 0; i < p_; ++i) {
        double logdet_omega = cholesky_helpers::get_log_det(cholesky_of_precision_);
        double logdet_omega_sub_ii = logdet_omega + MY_LOG(covariance_matrix_(i, i));

        size_t e = i * (i + 3) / 2;
        double proposal_sd = proposal_sds_(e);

        double theta_curr = (logdet_omega - logdet_omega_sub_ii) / 2;
        double theta_prop = rnorm(rng_, theta_curr, proposal_sd);

        precision_proposal_ = precision_matrix_;
        precision_proposal_(i, i) = precision_matrix_(i, i)
            - MY_EXP(theta_curr) * MY_EXP(theta_curr)
            + MY_EXP(theta_prop) * MY_EXP(theta_prop);

        double ln_alpha = log_density_impl_diag(i);
        if (determinant_tilt_ != 0.0) {
            ln_alpha += determinant_tilt_ * log_det_ratio_diag(i);
        }
        ln_alpha += diagonal_prior_->logp(0.5 * precision_proposal_(i, i));
        ln_alpha -= diagonal_prior_->logp(0.5 * precision_matrix_(i, i));
        ln_alpha += 2.0 * (theta_prop - theta_curr); // Jacobian: dK_ii/dtheta = 2*exp(2*theta)

        if (MY_LOG(runif(rng_)) < ln_alpha) {
            double omega_ii = precision_matrix_(i, i);
            precision_matrix_(i, i) = precision_proposal_(i, i);
            cholesky_update_after_diag(omega_ii, i);
        }

        proposal_sds_(e) = update_proposal_sd_with_robbins_monro(
            proposal_sds_(e), ln_alpha, rm_weight, target_accept);
    }

    // Invalidate gradient cache after MH updates
    constraint_dirty_ = true;
    theta_valid_ = false;
}

void GGMModel::refresh_cholesky() {
    cholesky_of_precision_ = arma::chol(precision_matrix_, "upper");
    arma::solve(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_),
                arma::eye(p_, p_), arma::solve_opts::fast);
    covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
}


void GGMModel::initialize_precision_from_mle() {
    // With n=0 there is no data; keep the identity initialization.
    if (n_ == 0) return;

    // Regularized MLE: K = n * inv(S + delta * I).
    // delta = trace(S) / (p * n) gives scale-appropriate shrinkage toward I.
    double trace_s = arma::trace(suf_stat_);
    double delta = trace_s / static_cast<double>(p_ * n_);
    arma::mat S_reg = suf_stat_ + delta * arma::eye(p_, p_);
    arma::mat K_init;
    if (arma::inv_sympd(K_init, S_reg)) {
        precision_matrix_ = static_cast<double>(n_) * K_init;

        // For fixed sparse graphs, zero out excluded edges and
        // recompute the diagonal to maintain positive definiteness.
        if (has_sparse_graph_) {
            for (size_t i = 0; i < p_ - 1; ++i) {
                for (size_t j = i + 1; j < p_; ++j) {
                    if (edge_indicators_(i, j) == 0) {
                        precision_matrix_(i, j) = 0.0;
                        precision_matrix_(j, i) = 0.0;
                    }
                }
            }
            // Make diagonally dominant to ensure PD after zeroing.
            for (size_t i = 0; i < p_; ++i) {
                double row_sum = 0.0;
                for (size_t j = 0; j < p_; ++j) {
                    if (j != i) row_sum += std::abs(precision_matrix_(i, j));
                }
                if (precision_matrix_(i, i) <= row_sum) {
                    precision_matrix_(i, i) = row_sum + 0.1;
                }
            }
        }

        refresh_cholesky();
    }
    // If inv_sympd fails, keep the identity initialization.
}


// =============================================================================
// Missing data imputation
// =============================================================================

void GGMModel::update_suf_stat_for_imputation(int variable, int person, double delta) {
    // INVARIANT: observations_(person, variable) must still hold x_old when
    // this function is called. The loop adds 2 * delta * x_old to the (v,v)
    // entry; the delta^2 correction completes the diagonal update.
    for (size_t q = 0; q < p_; q++) {
        suf_stat_(variable, q) += delta * observations_(person, q);
        suf_stat_(q, variable) += delta * observations_(person, q);
    }
    suf_stat_(variable, variable) += delta * delta;
}

void GGMModel::impute_missing() {
    if (!has_missing_) return;

    const int num_missings = missing_index_.n_rows;

    for (int miss = 0; miss < num_missings; miss++) {
        const int person = missing_index_(miss, 0);
        const int variable = missing_index_(miss, 1);

        // Compute conditional mean: mu = -sum_{k != v} omega_{vk} * x_{ik} / omega_{vv}
        double conditional_mean = 0.0;
        for (size_t k = 0; k < p_; k++) {
            if (k != static_cast<size_t>(variable)) {
                conditional_mean += precision_matrix_(variable, k) * observations_(person, k);
            }
        }
        conditional_mean = -conditional_mean / precision_matrix_(variable, variable);

        // Conditional variance: 1 / omega_{vv}
        double conditional_sd = std::sqrt(1.0 / precision_matrix_(variable, variable));

        // Sample new value
        double x_new = rnorm(rng_, conditional_mean, conditional_sd);
        double x_old = observations_(person, variable);
        double delta = x_new - x_old;

        // Incrementally update suf_stat_ (observations_ still holds x_old)
        update_suf_stat_for_imputation(variable, person, delta);

        // Now update the observation
        observations_(person, variable) = x_new;
    }

    // Full recompute at end of sweep to eliminate floating-point drift
    // (matches OMRF pattern; cost is O(np^2), negligible for typical sizes)
    suf_stat_ = observations_.t() * observations_;
}


// =============================================================================
// Factory function
// =============================================================================

GGMModel createGGMModelFromR(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const bool edge_selection,
    std::unique_ptr<BaseParameterPrior> interaction_prior,
    std::unique_ptr<BaseParameterPrior> diagonal_prior,
    const bool na_impute
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
            std::move(interaction_prior),
            std::move(diagonal_prior)
        );
    } else if (inputFromR.containsElementNamed("X")) {
        arma::mat X = Rcpp::as<arma::mat>(inputFromR["X"]);
        return GGMModel(
            X,
            prior_inclusion_prob,
            initial_edge_indicators,
            edge_selection,
            std::move(interaction_prior),
            std::move(diagonal_prior),
            na_impute
        );
    } else {
        throw std::invalid_argument("Input list must contain either 'X' or both 'n' and 'suf_stat'.");
    }

}
