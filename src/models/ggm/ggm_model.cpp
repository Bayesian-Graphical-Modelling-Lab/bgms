#include "models/ggm/ggm_model.h"
#include "rng/rng_utils.h"
#include "math/explog_macros.h"
#include "math/cholupdate.h"
#include "mcmc/execution/step_result.h"
#include "mcmc/execution/warmup_schedule.h"
#include "models/ggm/log_z_nlo.h"
#include "models/ggm/manuscript_nlo.h"
#include "models/ggm/z_ratio_estimator.h"
#include <stdexcept>

// ----------------------------------------------------------------------
// log_Z_NLO closed-form selector. Returns the centring log Z_NLO at G:
//   - manuscript App C NLO (eq:NLO-decomp) when use_manuscript_nlo_ && α=1
//   - bgms's pre-2026-05-21 log_Z_NLO_gamma otherwise (also the fallback
//     at α ≠ 1, where the manuscript form's Na/Nb/Nc terms are not ported).
// File-local helper so it doesn't pollute the GGMModel public API.
// ----------------------------------------------------------------------
static inline double compute_log_Z_NLO_centre(
    const arma::imat& G,
    double alpha, double beta, double sigma, double delta,
    bool use_manuscript
) {
    if (use_manuscript && alpha == 1.0) {
        return ggm_nlo::log_Z_manuscript_NLO_alpha1(G, beta, sigma, delta);
    }
    return log_Z_NLO_gamma(G, alpha, beta, sigma, /*include_F=*/false, delta);
}

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

    // K_new = K_old + vf1 vf2^T + vf2 vf1^T = K_old + u1 u1^T - u2 u2^T,
    // where u1 = (vf1 + vf2) / sqrt(2), u2 = (vf1 - vf2) / sqrt(2). The
    // change-of-basis diagonalises the symmetric rank-2 update so the chol
    // factor advances via one rank-1 update + one rank-1 downdate.
    u1_ = (vf1_ + vf2_) / sqrt(2);
    u2_ = (vf1_ - vf2_) / sqrt(2);

    // L update (2 x O(p^2)). Still required so cholesky_of_precision_ tracks
    // K (used by get_log_det in the next iteration).
    cholesky_update(cholesky_of_precision_, u1_);
    cholesky_downdate(cholesky_of_precision_, u2_);

    // Sherman-Morrison-Woodbury rank-2 update of covariance_matrix_ = inv(K).
    // Replaces the prior arma::solve(trimatu(L), I) (O(p^3)) + inv_L * inv_L.t()
    // (O(p^3)) refresh with 4 x O(p^2) work:
    //   K_new = K_old + M D M^T, M = [u1, u2], D = diag(+1, -1)
    //   inv(K_new) = inv(K_old) - A C^{-1} A^T,
    //     A = inv(K_old) M = [a1, a2],
    //     C = D^{-1} + M^T inv(K_old) M = diag(+1, -1) + symmetric 2x2.
    // Capacitance singularity (|det C| ~ 0) falls back to refresh_cholesky.
    // inv_cholesky_of_precision_ is no longer maintained per accept — only
    // refresh_cholesky() updates it (it's a scratch artefact of the prior
    // path, not read between accepts).
    arma::vec a1 = covariance_matrix_ * u1_;
    arma::vec a2 = covariance_matrix_ * u2_;
    double c11 =  1.0 + arma::dot(u1_, a1);
    double c12 =        arma::dot(u1_, a2);
    double c22 = -1.0 + arma::dot(u2_, a2);
    double det = c11 * c22 - c12 * c12;
    if (!std::isfinite(det) || std::abs(det) < 1e-14) {
        refresh_cholesky();
    } else {
        double inv_c00 =  c22 / det;
        double inv_c11 =  c11 / det;
        double inv_c01 = -c12 / det;
        covariance_matrix_ -= inv_c00 * (a1 * a1.t());
        covariance_matrix_ -= inv_c11 * (a2 * a2.t());
        covariance_matrix_ -= inv_c01 * (a1 * a2.t() + a2 * a1.t());
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

    // L rank-1 update so cholesky_of_precision_ tracks K (used by get_log_det).
    if (s)
        cholesky_downdate(cholesky_of_precision_, vf1_);
    else
        cholesky_update(cholesky_of_precision_, vf1_);

    // SMW rank-1 update of covariance_matrix_ = inv(K). Replaces the prior
    // O(p^3) solve+matmul refresh with one O(p^2) outer-product update.
    //   K_new = K_old + alpha * e_i e_i^T, alpha = K_new(i,i) - K_old(i,i)
    //   inv(K_new) = inv(K_old) - alpha * c_i c_i^T / (1 + alpha * c_i[i]),
    //     c_i = inv(K_old).col(i).
    // alpha = precision_proposal_(i,i) - omega_ii_old (note sign: delta is
    // defined the other way around above so the chol update/downdate branch
    // matches K_new > or < K_old). Refresh-fall-back guards near-singular
    // denom.
    double alpha = precision_proposal_(i, i) - omega_ii_old;
    arma::vec ci = covariance_matrix_.col(i);
    double denom = 1.0 + alpha * ci(i);
    if (!std::isfinite(denom) || std::abs(denom) < 1e-14) {
        refresh_cholesky();
    } else {
        double coeff = alpha / denom;
        covariance_matrix_ -= coeff * (ci * ci.t());
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
        // NOTE: under the Roverato slaving (K_jj <- c_3 + phi^2 with phi chosen
        // so K_ij <- 0), |K| is invariant to machine precision (proven via the
        // 2x2 cofactor identity and verified numerically at q<=10). So this
        // term is identically zero in practice; it's kept here defensively for
        // any future non-Roverato proposal variant.
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

        // Hierarchical-spec correction. The joint-spec MH ratio implicitly
        // targets π(Γ)·Z(Γ) marginally; to convert to the hierarchical
        // target with marginal π(Γ), multiply by Z(Γ_curr)/Z(Γ_star). With
        // V(Γ) ≈ 1/Z(Γ), this is V(Γ_star) / V(Γ_curr). In log form:
        //   ln_alpha += log|V(Γ_star)| - log|V(Γ_curr)|.
        // Lyne (2015) RR debias with the DEGORD-permuted Bartlett-Cholesky
        // inner sampler.
        double log_Z_NLO_star = log_Z_NLO_curr_;  // tentative; set below if hierarchical
        // F2: V(Γ_star) carried to the accept block so we can advance the
        // running V-diagnostic only on accept (set inside the hier_active
        // branch below).
        int    V_star_sign_for_diag   = 1;
        double V_star_log_abs_for_diag = std::numeric_limits<double>::quiet_NaN();
        bool hier_active = (graph_prior_spec_ == GraphPriorSpec::Hierarchical);
        if (hier_active) {
            ++n_hier_del_attempts_;
            ensure_hierarchical_state_();
            // Γ_star: this branch DELETES edge (i, j).
            arma::imat G_star = edge_indicators_;
            G_star(i, j) = 0;
            G_star(j, i) = 0;
            // log_Z_NLO_star via the cheap O(deg²) incremental at α=1 under
            // the default formula; full-recompute when use_manuscript_nlo_ is
            // set (no manuscript-incremental ported yet). At α ≠ 1 the
            // selector falls back to bgms's formula either way.
            if (prior_alpha_ == 1.0 && !use_manuscript_nlo_) {
                double d = log_Z_NLO_gamma_delta_incr_alpha1(
                    edge_indicators_, static_cast<int>(i), static_cast<int>(j),
                    prior_beta_, prior_sigma_, determinant_tilt_, false);
                log_Z_NLO_star = log_Z_NLO_curr_ + d;
            } else {
                log_Z_NLO_star = compute_log_Z_NLO_centre(
                    G_star, prior_alpha_, prior_beta_, prior_sigma_,
                    determinant_tilt_, use_manuscript_nlo_);
            }
            // V evaluated under the DEGORD permutation π that sends (i, j)
            // to (q-2, q-1).
            arma::ivec pi = degord::degord_permutation(
                static_cast<int>(p_), static_cast<int>(i), static_cast<int>(j));
            arma::imat G_pi_curr = degord::permute_graph(edge_indicators_, pi);
            arma::imat G_pi_star = degord::permute_graph(G_star, pi);
            // Log-space V: avoids underflow in c = kappa * exp(log_Z_NLO) at
            // large p (log_Z_NLO is ~ -3500 at p=100, δ=1 → c flushes to 0).
            // log_kappa cancels in the MH ratio, but log_c per Γ is needed
            // to evaluate log|expm1(log_Zhat_m - log_c)| pointwise.
            //
            // Paired call shares the inner Phi-build across Γ_curr / Γ_star
            // by caching (rw_head, S_trail) under a_curr and re-evaluating
            // only row q-2 under a_star — halves per-pool work vs two
            // independent V_log_at_Gamma_pi_degord calls.
            double log_kappa = std::log(v_kappa_);
            double log_c_curr = log_kappa + log_Z_NLO_curr_;
            double log_c_star = log_kappa + log_Z_NLO_star;
            auto V_pair = degord::V_log_pair_at_Gamma_curr_star_degord(
                v_K_depth_, v_pools_t_,
                G_pi_curr, G_pi_star, chain_aux_degord_,
                log_c_curr, log_c_star, v_rho_);
            // F2: initialise running V-diagnostic from V_pair.curr the first
            // time we see a finite value (so the side-car has a meaningful
            // entry from iteration 1 even before any accept). On accept the
            // value is overwritten below from V_pair.star.
            if (!v_diag_initialized_ &&
                std::isfinite(V_pair.curr.first) && V_pair.curr.second != 0) {
                current_sign_V_     = V_pair.curr.second;
                current_log_abs_V_  = V_pair.curr.first;
                v_diag_initialized_ = true;
                last_v_pi_i_ = static_cast<int>(i);
                last_v_pi_j_ = static_cast<int>(j);
            }
            V_star_sign_for_diag    = V_pair.star.second;
            V_star_log_abs_for_diag = V_pair.star.first;
            // Auto-reject on non-finite log|V| (sentinel for V = 0 or
            // non-finite Zhat) or on sign flip across Γ_curr / Γ_star. The
            // sign-flip reject remains until a proper Lyne sign accumulator
            // composes downstream (F3).
            if (!std::isfinite(V_pair.curr.first) || V_pair.curr.second == 0 ||
                !std::isfinite(V_pair.star.first) || V_pair.star.second == 0 ||
                V_pair.curr.second != V_pair.star.second) {
                // Classify the auto-reject for diagnostics.
                if (!std::isfinite(V_pair.curr.first) ||
                    !std::isfinite(V_pair.star.first)) {
                    ++n_hier_del_nonfinite_;
                } else if (V_pair.curr.second == 0 || V_pair.star.second == 0) {
                    ++n_hier_del_signzero_;
                } else {
                    ++n_hier_del_signflip_;
                }
                ln_alpha = -std::numeric_limits<double>::infinity();
            } else {
                ln_alpha += V_pair.star.first - V_pair.curr.first;
            }
        }

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
            if (hier_active) {
                log_Z_NLO_curr_ = log_Z_NLO_star;
                // Γ_star is now Γ_curr; advance the running V state.
                current_sign_V_    = V_star_sign_for_diag;
                current_log_abs_V_ = V_star_log_abs_for_diag;
                v_diag_initialized_ = true;
                last_v_pi_i_ = static_cast<int>(i);
                last_v_pi_j_ = static_cast<int>(j);
            }
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
        // to the MH ratio. See the DELETE branch for the Roverato-invariance
        // note - kept here defensively.
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

        // Hierarchical-spec correction (ADD branch): see DELETE branch for the
        // rationale. Γ_star ADDS edge (i, j) here, so log_Z_NLO_star differs
        // from log_Z_NLO_curr by the +add direction of the incremental.
        double log_Z_NLO_star_add = log_Z_NLO_curr_;
        // F2: V(Γ_star) carried to the accept block (mirrors DELETE branch).
        int    V_star_sign_for_diag_add   = 1;
        double V_star_log_abs_for_diag_add = std::numeric_limits<double>::quiet_NaN();
        bool hier_active_add = (graph_prior_spec_ == GraphPriorSpec::Hierarchical);
        if (hier_active_add) {
            ++n_hier_add_attempts_;
            ensure_hierarchical_state_();
            arma::imat G_star = edge_indicators_;
            G_star(i, j) = 1;
            G_star(j, i) = 1;
            if (prior_alpha_ == 1.0 && !use_manuscript_nlo_) {
                double d = log_Z_NLO_gamma_delta_incr_alpha1(
                    edge_indicators_, static_cast<int>(i), static_cast<int>(j),
                    prior_beta_, prior_sigma_, determinant_tilt_, false);
                log_Z_NLO_star_add = log_Z_NLO_curr_ + d;
            } else {
                log_Z_NLO_star_add = compute_log_Z_NLO_centre(
                    G_star, prior_alpha_, prior_beta_, prior_sigma_,
                    determinant_tilt_, use_manuscript_nlo_);
            }
            arma::ivec pi = degord::degord_permutation(
                static_cast<int>(p_), static_cast<int>(i), static_cast<int>(j));
            arma::imat G_pi_curr = degord::permute_graph(edge_indicators_, pi);
            arma::imat G_pi_star = degord::permute_graph(G_star, pi);
            // Log-space V with within-toggle cache reuse (see DELETE branch
            // for the underflow rationale and the sign-flip auto-reject
            // contract).
            double log_kappa = std::log(v_kappa_);
            double log_c_curr = log_kappa + log_Z_NLO_curr_;
            double log_c_star = log_kappa + log_Z_NLO_star_add;
            auto V_pair = degord::V_log_pair_at_Gamma_curr_star_degord(
                v_K_depth_, v_pools_t_,
                G_pi_curr, G_pi_star, chain_aux_degord_,
                log_c_curr, log_c_star, v_rho_);
            // F2: same diagnostic seeding as in the DELETE branch.
            if (!v_diag_initialized_ &&
                std::isfinite(V_pair.curr.first) && V_pair.curr.second != 0) {
                current_sign_V_     = V_pair.curr.second;
                current_log_abs_V_  = V_pair.curr.first;
                v_diag_initialized_ = true;
                last_v_pi_i_ = static_cast<int>(i);
                last_v_pi_j_ = static_cast<int>(j);
            }
            V_star_sign_for_diag_add    = V_pair.star.second;
            V_star_log_abs_for_diag_add = V_pair.star.first;
            if (!std::isfinite(V_pair.curr.first) || V_pair.curr.second == 0 ||
                !std::isfinite(V_pair.star.first) || V_pair.star.second == 0 ||
                V_pair.curr.second != V_pair.star.second) {
                if (!std::isfinite(V_pair.curr.first) ||
                    !std::isfinite(V_pair.star.first)) {
                    ++n_hier_add_nonfinite_;
                } else if (V_pair.curr.second == 0 || V_pair.star.second == 0) {
                    ++n_hier_add_signzero_;
                } else {
                    ++n_hier_add_signflip_;
                }
                ln_alpha = -std::numeric_limits<double>::infinity();
            } else {
                ln_alpha += V_pair.star.first - V_pair.curr.first;
            }
        }

        if (MY_LOG(runif(rng_)) < ln_alpha) {
            // Accept: turn ON the edge
            // Store old values for Cholesky update
            double omega_ij_old = precision_matrix_(i, j);
            double omega_jj_old = precision_matrix_(j, j);

            if (hier_active_add) {
                log_Z_NLO_curr_ = log_Z_NLO_star_add;
                current_sign_V_    = V_star_sign_for_diag_add;
                current_log_abs_V_ = V_star_log_abs_for_diag_add;
                v_diag_initialized_ = true;
                last_v_pi_i_ = static_cast<int>(i);
                last_v_pi_j_ = static_cast<int>(j);
            }

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

    // Catch SMW-accumulated drift in covariance_matrix_ over a long chain.
    // O(p^2); the refresh path is only taken when drift exceeds tolerance.
    check_and_refresh_if_drift_();
}

void GGMModel::init_metropolis_adaptation(const WarmupSchedule& schedule) {
    metropolis_adapter_ = std::make_unique<MetropolisAdaptationController>(
        proposal_sds_, schedule, target_accept_);
}

void GGMModel::prepare_iteration() {
    // Shuffle edge visit order for random-scan edge selection.
    // Called unconditionally to keep RNG state consistent.
    shuffled_edge_order_ = arma_randperm(rng_, num_pairwise_);
    // Refresh the V/RR U-pool at iteration start.
    //   Legacy (mh_U_ = false): unconditional draw from μ(U). This breaks
    //     PMMH invariance on (Γ, K, U, N) — the conditional under the
    //     augmented target is V·μ, not μ alone — and yields a small Γ-
    //     marginal bias (~−0.001 nats at p=20, p_inc=0.05).
    //   Fixed (mh_U_ = true): V-ratio MH step on U, accepting on
    //     log|V_new| − log|V_old|. Companion-AI delivery 2026-05-21.
    //   On the very first prepare_iteration after lazy init, the state has
    //     just been seeded with a fresh U via ensure_hierarchical_state_();
    //     no comparison is possible, so we skip the MH step and treat the
    //     init draw as the iteration-0 U.
    if (graph_prior_spec_ == GraphPriorSpec::Hierarchical) {
        bool was_built = hierarchical_state_built_;
        ensure_hierarchical_state_();
        if (!was_built) {
            // First iteration: state seeded by ensure_hierarchical_state_;
            // nothing more to do.
            return;
        }
        if (mh_U_) {
            mh_on_U_step_();
        } else {
            refresh_z_ratio_pool_();
        }
    }
}


// ----------------------------------------------------------------------
// Hierarchical-spec (Phase 4) implementation
// ----------------------------------------------------------------------

void GGMModel::set_graph_prior_spec(GraphPriorSpec spec) {
    if (spec == graph_prior_spec_) return;
    graph_prior_spec_ = spec;
    hierarchical_state_built_ = false;  // force rebuild on next use
}


void GGMModel::set_z_ratio_tuning(int M_inner, double kappa, double rho) {
    if (M_inner < 1) throw std::runtime_error("M_inner must be >= 1");
    if (!(rho > 0.0 && rho < 1.0))
        throw std::runtime_error("rho must be in (0, 1)");
    if (!(kappa > 0.0))
        throw std::runtime_error("kappa must be > 0");
    v_M_inner_ = M_inner;
    v_kappa_   = kappa;
    v_rho_     = rho;
    hierarchical_state_built_ = false;
}


void GGMModel::ensure_hierarchical_state_() {
    if (hierarchical_state_built_) return;
    // Validate prior family. The closed-form log_Z_NLO_gamma machinery only
    // covers slab = Normal(0, σ) and diag = Gamma(α, β) on K_ii/2.
    const auto* slab = dynamic_cast<const NormalPrior*>(interaction_prior_.get());
    if (slab == nullptr)
        throw std::runtime_error(
            "Hierarchical graph_prior_spec requires a Normal slab "
            "(NormalPrior). Re-fit with interaction_prior_type = 'normal'.");
    const auto* diag = dynamic_cast<const GammaScalePrior*>(diagonal_prior_.get());
    if (diag == nullptr)
        throw std::runtime_error(
            "Hierarchical graph_prior_spec requires a Gamma diagonal prior "
            "(GammaScalePrior).");

    prior_sigma_ = slab->scale();
    prior_alpha_ = diag->shape();
    prior_beta_  = diag->rate();
    double delta = determinant_tilt_;

    chain_aux_degord_ = degord::make_chain_aux(
        static_cast<int>(p_), prior_alpha_, prior_beta_, prior_sigma_, delta);

    // Analytic centring at the current Γ (full-recompute; the incremental
    // form is only used on accept). Use F = false to match the production
    // convention (NLO without the F-piece — the F overcorrects at α > 1).
    // The selector picks the manuscript App C NLO at α = 1 when the
    // use_manuscript_nlo_ flag is set; otherwise it uses bgms's pre-
    // 2026-05-21 log_Z_NLO_gamma.
    log_Z_NLO_curr_ = compute_log_Z_NLO_centre(
        edge_indicators_, prior_alpha_, prior_beta_, prior_sigma_,
        delta, use_manuscript_nlo_);

    refresh_z_ratio_pool_();
    hierarchical_state_built_ = true;
}


void GGMModel::refresh_z_ratio_pool_() {
    degord::draw_U_degord_rr(
        rng_, v_K_depth_, v_pools_t_, v_M_inner_, static_cast<int>(p_), v_rho_);
}


// MH step on (U, K_depth) at the augmented target. The proposal is a fresh
// draw from μ(U)·P(N), so the proposal density on the forward and reverse
// move cancels and the MH ratio reduces to V_new / V_old. Auto-rejects on
// any non-finite log|V| or sign(V_new) ≠ sign(V_old) — same convention as
// the between-Γ MH (deferred Lyne sign accumulator). Choice of permutation
// π: arbitrary; any π yields an unbiased V estimator. We pick π = (0, 1)
// for simplicity (gives the canonical degord reordering that maps the first
// two vertices to themselves).
void GGMModel::mh_on_U_step_() {
    if (graph_prior_spec_ != GraphPriorSpec::Hierarchical) return;

    ++n_mh_U_attempts_;

    // V_old at (Γ_curr, U_old). Reuse the cached running V state when it's
    // available: the between-Γ MH stores V evaluated at the last-toggled
    // edge's DEGORD permutation π_{last_v_pi_i_, last_v_pi_j_}, and the V
    // estimator is unbiased under any permutation. Reusing the cache skips
    // a full V_log_at_Gamma_pi_degord call — which dominates per-iter cost
    // when K_depth is large.
    //
    // If the cache hasn't been seeded yet (first iteration with no finite
    // V_pair), fall back to recomputing V_old at the canonical π = (0, 1).
    double log_c = std::log(v_kappa_) + log_Z_NLO_curr_;
    int    pi_i, pi_j;
    std::pair<double, int> V_old;
    if (v_diag_initialized_) {
        V_old = { current_log_abs_V_, current_sign_V_ };
        pi_i = last_v_pi_i_;
        pi_j = last_v_pi_j_;
    } else {
        pi_i = 0;
        pi_j = 1;
        arma::ivec pi_canon = degord::degord_permutation(
            static_cast<int>(p_), pi_i, pi_j);
        arma::imat G_pi_canon = degord::permute_graph(edge_indicators_, pi_canon);
        V_old = degord::V_log_at_Gamma_pi_degord(
            v_K_depth_, v_pools_t_, G_pi_canon, chain_aux_degord_,
            log_c, v_rho_);
    }

    // Proposal: fresh U_new ~ μ, K_depth_new ~ P(N).
    int K_depth_new;
    std::vector<arma::mat> pools_new;
    degord::draw_U_degord_rr(
        rng_, K_depth_new, pools_new, v_M_inner_, static_cast<int>(p_), v_rho_);

    // V_new at (Γ_curr, U_new) evaluated under the SAME π as the cached V_old.
    arma::ivec pi_vec = degord::degord_permutation(
        static_cast<int>(p_), pi_i, pi_j);
    arma::imat G_pi_curr = degord::permute_graph(edge_indicators_, pi_vec);
    auto V_new = degord::V_log_at_Gamma_pi_degord(
        K_depth_new, pools_new, G_pi_curr, chain_aux_degord_,
        log_c, v_rho_);

    // Auto-reject paths.
    if (!std::isfinite(V_old.first) || !std::isfinite(V_new.first)) {
        ++n_mh_U_nonfinite_;
        return;
    }
    if (V_old.second == 0 || V_new.second == 0) {
        ++n_mh_U_signzero_;
        return;
    }
    if (V_old.second != V_new.second) {
        ++n_mh_U_signflip_;
        return;
    }

    double log_alpha = V_new.first - V_old.first;
    if (MY_LOG(runif(rng_)) < log_alpha) {
        // Accept: install proposed pool and update running V state to V_new.
        // last_v_pi_i_ / last_v_pi_j_ already point at this π, so they don't
        // need updating.
        v_pools_t_ = std::move(pools_new);
        v_K_depth_ = K_depth_new;
        current_log_abs_V_ = V_new.first;
        current_sign_V_    = V_new.second;
        v_diag_initialized_ = true;
        ++n_mh_U_accepts_;
    }
    // Reject: current_log_abs_V_ / current_sign_V_ already correspond to
    // (Γ_curr, U_old) under last_v_pi_i_,last_v_pi_j_; nothing to update.
}


// NOTE: the on-accept update of log_Z_NLO_curr_ lives inline in
// update_edge_indicator_parameter_pair (both branches set log_Z_NLO_curr_ to
// the pre-computed log_Z_NLO_star{,_add} inside their MH accept blocks).
// The incremental form is the alpha=1 fast path; alpha != 1 full-recomputes.

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
    // SMW drift check; same rationale as the end-of-MH-step path.
    check_and_refresh_if_drift_();
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

void GGMModel::check_and_refresh_if_drift_() {
    // diag(cov * K) should be ones. Compute the max abs deviation in O(p^2)
    // via the elementwise product (K is symmetric so cov(i,:) * K(:,i)
    // equals arma::dot(cov.row(i), K.row(i))).
    arma::vec d = arma::sum(covariance_matrix_ % precision_matrix_, 1) - 1.0;
    double drift = arma::abs(d).max();
    if (!std::isfinite(drift) || drift > kCovDriftTol_) {
        refresh_cholesky();
    }
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
