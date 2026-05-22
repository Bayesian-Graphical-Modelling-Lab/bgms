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

    // Interaction prior on K_yy_{ij} = -0.5 * Omega_{ij}. Pass (i, j) so an
    // edgewise prior (GraphicalGPrior) can resolve V_ij; single-scale priors
    // ignore the coords via the BaseParameterPrior default override.
    ln_alpha += interaction_prior_->logp(-0.5 * precision_proposal_(i, j),
                                         static_cast<int>(i), static_cast<int>(j));
    ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j),
                                         static_cast<int>(i), static_cast<int>(j));

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
        ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j),
                                              static_cast<int>(i), static_cast<int>(j))
                    - MY_LOG(2.0);

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
        ln_alpha += interaction_prior_->logp(-0.5 * omega_prop_ij,
                                              static_cast<int>(i), static_cast<int>(j))
                    - MY_LOG(2.0);

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


// =====================================================================
// Graphical G-prior
// =====================================================================

void GGMModel::compute_gg_V_ij_() {
    // V_ij = 1 / (4 n · S̄_ii · S̄_jj)  with  S̄ = suf_stat / n
    //      = n / (4 · suf_stat_(i,i) · suf_stat_(j,j))
    //
    // Uses n_ (= n_rows − 1, the GGM's centered-data degrees of freedom).
    // Diagonal entries are left at zero — they are never indexed.
    V_ij_.set_size(p_, p_);
    V_ij_.zeros();
    for (size_t i = 0; i < p_; ++i) {
        for (size_t j = 0; j < p_; ++j) {
            if (i == j) continue;
            const double s_ii = suf_stat_(i, i);
            const double s_jj = suf_stat_(j, j);
            if (s_ii > 0.0 && s_jj > 0.0) {
                V_ij_(i, j) = static_cast<double>(n_) /
                              (4.0 * s_ii * s_jj);
            } else {
                // Degenerate variable (zero sufficient statistic on the
                // diagonal) — set V_ij to NaN so any subsequent prior eval
                // surfaces the bad-data condition explicitly.
                V_ij_(i, j) = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }
}


void GGMModel::enable_gg_prior(
    GGHyperprior hyperprior,
    double g_init,
    double tcch_a,
    double tcch_b,
    double tcch_r,
    double tcch_s,
    double tcch_u
) {
    if (g_init <= 0.0 || !std::isfinite(g_init)) {
        throw std::invalid_argument(
            "enable_gg_prior: g_init must be a positive finite number.");
    }
    use_gg_prior_  = true;
    gg_hyperprior_ = hyperprior;
    t_             = std::sqrt(g_init);
    gg_tcch_a_     = tcch_a;
    gg_tcch_b_     = tcch_b;
    gg_tcch_r_     = tcch_r;
    gg_tcch_s_     = tcch_s;
    gg_tcch_u_     = tcch_u;
    // Snapshot n for hyper-g/n: under prior-only mode n_ is later zeroed
    // but the hyperprior on g still references the data-defined sample size.
    n_at_gg_setup_ = static_cast<int>(n_);

    compute_gg_V_ij_();

    // Replace user-supplied priors with GG-bound versions. The user-supplied
    // priors are silently dropped; sample_ggm warns at the R level when a
    // user tries to pair graphical_g_prior with a non-default precision
    // scale prior.
    interaction_prior_ = std::make_unique<GraphicalGPrior>(&V_ij_, &t_);
    diagonal_prior_    = std::make_unique<GraphicalGDiag>(&t_);

    // Future gradient calls must rebuild against the new priors.
    constraint_dirty_  = true;
}


void GGMModel::set_prior_only() {
    // Zero the sufficient statistic and sample count so that every data-
    // contribution term in the MH ratios and NUTS gradient evaluates to zero.
    //   log_likelihood = n_ * (...) + log|K|/2 - tr(K·S)/2  → 0
    //   block_log_likelihood: n_·logdet term + suf_stat trace contributions  → 0
    // The GG-prior V_ij was already cached in V_ij_ at enable_gg_prior() time
    // from the original (n_, suf_stat_), so the prior over (K, Gamma) is
    // unchanged. Diagonal Gamma prior, slab Normal prior, edge prior and
    // determinant tilt remain active and continue to drive the chain.
    prior_only_       = true;
    n_                = 0;
    suf_stat_.zeros();
    // Force a gradient-engine rebuild so its cached (n_, suf_stat_) views
    // refresh to the zeroed state.
    gradient_engine_.rebuild(
        constraint_structure_, n_, suf_stat_,
        *interaction_prior_, *diagonal_prior_, determinant_tilt_);
    constraint_dirty_ = false;
}


void GGMModel::gg_rebind_priors_() {
    // Cloned priors may still hold pointers into the source model's V_ij_/t_.
    // dynamic_cast is safe on any other prior — falls through silently.
    if (auto* slab = dynamic_cast<GraphicalGPrior*>(interaction_prior_.get())) {
        slab->bind_to(&V_ij_, &t_);
    }
    if (auto* diag = dynamic_cast<GraphicalGDiag*>(diagonal_prior_.get())) {
        diag->bind_to(&t_);
    }
}


void GGMModel::gg_update_t_() {
    if (!use_gg_prior_) return;
    if (gg_hyperprior_ == GGHyperprior::Fixed) return;

    // ConjugateGamma branch (v1): full conditional on t under
    //   t ~ Gamma(a₀, b₀) prior  ⇒
    //   t | Θ, Γ, Y ~ Gamma(a₀ + n·p/2 + p·δ,  b₀ + ½ tr(η · S))
    // with η = Θ / t. tr(η · S) = tr(Θ · S) / t at current state.
    //
    // After sampling t_new, rescale the chain state Θ ← (t_new/t_old) · Θ
    // and propagate to the Cholesky factor and the covariance cache.
    if (gg_hyperprior_ == GGHyperprior::ConjugateGamma) {
        const double t_old   = t_;
        const double tr_Theta_S = arma::dot(precision_matrix_, suf_stat_);
        const double tr_eta_S   = tr_Theta_S / t_old;
        const double shape = gg_tcch_a_
                             + 0.5 * static_cast<double>(n_) * static_cast<double>(p_)
                             + static_cast<double>(p_) * determinant_tilt_;
        const double rate  = gg_tcch_b_ + 0.5 * tr_eta_S;
        // Guard against pathological rate (numerical).
        if (!(rate > 0.0) || !std::isfinite(rate)) return;
        const double t_new = rgamma(rng_, shape, 1.0 / rate);
        if (!(t_new > 0.0) || !std::isfinite(t_new)) return;
        const double alpha = t_new / t_old;

        // Rescale state: K_new = α · K_old.  Cholesky L_new = √α · L_old.
        // Inverse Cholesky and covariance scale by 1/√α and 1/α respectively.
        precision_matrix_         *= alpha;
        const double sqrt_alpha    = std::sqrt(alpha);
        cholesky_of_precision_    *= sqrt_alpha;
        if (inv_cholesky_of_precision_.n_elem > 0) {
            inv_cholesky_of_precision_ /= sqrt_alpha;
        }
        if (covariance_matrix_.n_elem > 0) {
            covariance_matrix_     /= alpha;
        }
        t_ = t_new;
        // Theta caches and gradient engine state depend on the slab scale,
        // which has changed.
        theta_valid_      = false;
        constraint_dirty_ = true;
        return;
    }

    // ZellnerSiow / HyperG / HyperGOverN: MH-on-log(g) with scale-matched
    // joint proposal (Θ, g_old) → (α·Θ, g_new), α = √(g_new/g_old). The
    // scale-matching collapses the slab contribution to the MH ratio so
    // only the hyperprior, diagonal Gamma, det-tilt, likelihood, and
    // Jacobian terms remain.
    if (gg_hyperprior_ == GGHyperprior::ZellnerSiow ||
        gg_hyperprior_ == GGHyperprior::HyperG     ||
        gg_hyperprior_ == GGHyperprior::HyperGOverN) {
        gg_mh_update_log_g_();
        return;
    }

    // TCCH: not yet implemented (needs the truncated compound confluent
    // hypergeometric pdf evaluation). Fall through with a one-shot
    // warning per chain so the user sees the configuration is unsupported
    // without flooding the console.
    if (!tcch_warned_) {
        Rcpp::warning("GGMModel: tCCH g-hyperprior is not yet implemented; "
                      "g is held fixed.");
        tcch_warned_ = true;
    }
}


double GGMModel::gg_log_hyperprior_(double g) const {
    // Returns log π(g) up to an additive constant (constants drop in MH).
    switch (gg_hyperprior_) {
        case GGHyperprior::ZellnerSiow: {
            // Standard Zellner-Siow: g ~ IG(½, b/2),
            //   π(g) ∝ g^(-3/2) · exp(-b / (2g))
            // The scale b is read from gg_tcch_b_ (default 1.0 for the
            // canonical form). Note this is the "no-n" variant; we leave
            // the dimension-corrected n·b form to the user via tcch_b.
            const double b = gg_tcch_b_;
            return -1.5 * std::log(g) - 0.5 * b / g;
        }
        case GGHyperprior::HyperG: {
            // hyper-g (Liang et al. 2008):
            //   π(g) = (a - 2) / 2 · (1 + g)^(-a/2),  a > 2
            // Default a = 3 gives the recommended uniform-on-shrinkage prior.
            const double a = gg_tcch_a_;
            return -0.5 * a * std::log1p(g);
        }
        case GGHyperprior::HyperGOverN: {
            // hyper-g/n: π(g) = (a-2)/(2n) · (1 + g/n)^(-a/2)
            const double a = gg_tcch_a_;
            const double n = std::max(1.0,
                static_cast<double>(n_at_gg_setup_));
            return -0.5 * a * std::log1p(g / n);
        }
        default:
            return 0.0;
    }
}


void GGMModel::gg_mh_update_log_g_() {
    // Joint move (Θ, g_old) → (α·Θ, g_new) under symmetric Gaussian
    // random walk on log(g). With α = t_new/t_old = √(g_new/g_old) the
    // slab contribution to the MH ratio cancels (the rescaled ω = -K/2
    // ratio against the rescaled scale leaves the quadratic identical).
    // The surviving terms are:
    //   • hyperprior:     log π(g_new) − log π(g_old)
    //   • likelihood:     (n·p/2) log(α) − (α − 1)/2 · tr(Θ·S)
    //   • det-tilt:       p·δ · log(α)
    //   • Jacobian (Θ):   (p + q_active) · log(α)
    //   • Jacobian (g):   2 · log(α)        (log-Normal symmetric proposal)
    //   • diagonal Gamma: −p · log(α)
    // Several terms collapse: −p − q_active (slab norm) + p + q_active
    // (Jacobian) = 0. Net coefficient of log(α):
    //   n·p/2 + p·δ + 2.

    const double g_old     = t_ * t_;
    if (!(g_old > 0.0) || !std::isfinite(g_old)) return;
    const double log_g_old = std::log(g_old);
    const double log_g_new =
        log_g_old + gg_log_g_proposal_sd_ * rnorm(rng_);
    const double g_new = std::exp(log_g_new);
    if (!(g_new > 0.0) || !std::isfinite(g_new)) return;
    const double t_new = std::sqrt(g_new);
    const double alpha = t_new / t_;
    const double log_alpha = std::log(alpha);

    const double log_prior_diff =
        gg_log_hyperprior_(g_new) - gg_log_hyperprior_(g_old);

    const double tr_Theta_S = arma::dot(precision_matrix_, suf_stat_);
    const double n_d = static_cast<double>(n_);
    const double p_d = static_cast<double>(p_);
    const double coef_log_alpha =
        0.5 * n_d * p_d + p_d * determinant_tilt_ + 2.0;
    const double log_ratio = log_prior_diff
                           + coef_log_alpha * log_alpha
                           - 0.5 * (alpha - 1.0) * tr_Theta_S;

    ++gg_log_g_n_total_;
    if (std::log(runif(rng_)) < log_ratio) {
        precision_matrix_      *= alpha;
        const double sqrt_alpha = std::sqrt(alpha);
        cholesky_of_precision_ *= sqrt_alpha;
        if (inv_cholesky_of_precision_.n_elem > 0) {
            inv_cholesky_of_precision_ /= sqrt_alpha;
        }
        if (covariance_matrix_.n_elem > 0) {
            covariance_matrix_ /= alpha;
        }
        t_ = t_new;
        theta_valid_      = false;
        constraint_dirty_ = true;
        ++gg_log_g_n_accept_;
    }
}


void GGMModel::prepare_iteration() {
    // Shuffle edge visit order for random-scan edge selection.
    // Called unconditionally to keep RNG state consistent.
    shuffled_edge_order_ = arma_randperm(rng_, num_pairwise_);
    // Graphical G-prior: end-of-sweep Gibbs/MH on t = √g. Called from
    // prepare_iteration() so the update fires once per sweep regardless
    // of sampler (MH or NUTS). Internal no-op when use_gg_prior_ is false.
    gg_update_t_();
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
            ln_alpha += interaction_prior_->logp(-0.5 * precision_proposal_(i, j),
                                                  static_cast<int>(i), static_cast<int>(j));
            ln_alpha -= interaction_prior_->logp(-0.5 * precision_matrix_(i, j),
                                                  static_cast<int>(i), static_cast<int>(j));

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
