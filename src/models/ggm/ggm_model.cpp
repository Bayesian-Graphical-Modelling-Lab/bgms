#include "models/ggm/ggm_model.h"
#include "rng/rng_utils.h"
#include "math/explog_macros.h"
#include "math/cholupdate.h"
#include "mcmc/execution/step_result.h"
#include "mcmc/execution/warmup_schedule.h"
#include "models/ggm/log_z_nlo.h"
#include "models/ggm/manuscript_nlo.h"
#include "models/ggm/sd_density_at_zero.h"
#include "models/ggm/sd_density_l_space.h"
#include "models/ggm/sd_density_l_space_quad.h"
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

    // prior_only_ skips the likelihood ratio (chain targets π(Γ, K) instead
    // of π(Γ, K | Y)). PD-ness of the proposal still must be checked: the
    // determinant tilt and Cholesky update both rely on it.
    double ln_alpha = prior_only_ ? 0.0 : log_density_impl_edge(i, j);
    if (prior_only_) {
        arma::mat R_chk;
        if (!arma::chol(R_chk, precision_proposal_)) return 0.0;
    }

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

    // prior_only_: skip likelihood; still validate PD via the implicit
    // positive-K_ii proposal (theta_prop on log-scale → K_ii > 0).
    double ln_alpha = prior_only_ ? 0.0 : log_density_impl_diag(i);

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
            if (plug_in_nlo_) {
                // Plug-in mode: replace V_pair entirely with the deterministic
                // closed-form ratio log Z(Γ_curr) − log Z(Γ_star), approximated
                // by mNLO. No U, no K, no pool work. Flat cost in p.
                ln_alpha += log_Z_NLO_curr_ - log_Z_NLO_star;
            } else {
                // V evaluated under the DEGORD permutation π that sends (i, j)
                // to (q-2, q-1).
                arma::ivec pi = degord::degord_permutation(
                    static_cast<int>(p_), static_cast<int>(i), static_cast<int>(j));
                arma::imat G_pi_curr = degord::permute_graph(edge_indicators_, pi);
                arma::imat G_pi_star = degord::permute_graph(G_star, pi);
                double log_kappa = std::log(v_kappa_);
                double log_c_curr = log_kappa + log_Z_NLO_curr_;
                double log_c_star = log_kappa + log_Z_NLO_star;
                auto V_pair = degord::V_log_pair_at_Gamma_curr_star_degord(
                    v_K_depth_, v_pools_t_,
                    G_pi_curr, G_pi_star, chain_aux_degord_,
                    log_c_curr, log_c_star, v_rho_);
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
                if (!std::isfinite(V_pair.curr.first) || V_pair.curr.second == 0 ||
                    !std::isfinite(V_pair.star.first) || V_pair.star.second == 0) {
                    if (!std::isfinite(V_pair.curr.first) ||
                        !std::isfinite(V_pair.star.first)) {
                        ++n_hier_del_nonfinite_;
                    } else {
                        ++n_hier_del_signzero_;
                    }
                    ln_alpha = -std::numeric_limits<double>::infinity();
                } else {
                    if (V_pair.curr.second != V_pair.star.second) {
                        ++n_hier_del_signflip_;
                    }
                    ln_alpha += V_pair.star.first - V_pair.curr.first;
                }
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
                if (!plug_in_nlo_) {
                    // Γ_star is now Γ_curr; advance the running V state.
                    current_sign_V_    = V_star_sign_for_diag;
                    current_log_abs_V_ = V_star_log_abs_for_diag;
                    v_diag_initialized_ = true;
                    last_v_pi_i_ = static_cast<int>(i);
                    last_v_pi_j_ = static_cast<int>(j);
                }
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
            if (plug_in_nlo_) {
                ln_alpha += log_Z_NLO_curr_ - log_Z_NLO_star_add;
            } else {
                arma::ivec pi = degord::degord_permutation(
                    static_cast<int>(p_), static_cast<int>(i), static_cast<int>(j));
                arma::imat G_pi_curr = degord::permute_graph(edge_indicators_, pi);
                arma::imat G_pi_star = degord::permute_graph(G_star, pi);
                double log_kappa = std::log(v_kappa_);
                double log_c_curr = log_kappa + log_Z_NLO_curr_;
                double log_c_star = log_kappa + log_Z_NLO_star_add;
                auto V_pair = degord::V_log_pair_at_Gamma_curr_star_degord(
                    v_K_depth_, v_pools_t_,
                    G_pi_curr, G_pi_star, chain_aux_degord_,
                    log_c_curr, log_c_star, v_rho_);
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
                    !std::isfinite(V_pair.star.first) || V_pair.star.second == 0) {
                    if (!std::isfinite(V_pair.curr.first) ||
                        !std::isfinite(V_pair.star.first)) {
                        ++n_hier_add_nonfinite_;
                    } else {
                        ++n_hier_add_signzero_;
                    }
                    ln_alpha = -std::numeric_limits<double>::infinity();
                } else {
                    if (V_pair.curr.second != V_pair.star.second) {
                        ++n_hier_add_signflip_;
                    }
                    ln_alpha += V_pair.star.first - V_pair.curr.first;
                }
            }
        }

        if (MY_LOG(runif(rng_)) < ln_alpha) {
            // Accept: turn ON the edge
            // Store old values for Cholesky update
            double omega_ij_old = precision_matrix_(i, j);
            double omega_jj_old = precision_matrix_(j, j);

            if (hier_active_add) {
                log_Z_NLO_curr_ = log_Z_NLO_star_add;
                if (!plug_in_nlo_) {
                    current_sign_V_    = V_star_sign_for_diag_add;
                    current_log_abs_V_ = V_star_log_abs_for_diag_add;
                    v_diag_initialized_ = true;
                    last_v_pi_i_ = static_cast<int>(i);
                    last_v_pi_j_ = static_cast<int>(j);
                }
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

// =====================================================================
// Savage-Dickey between-Γ MH (handoff 2026-05-22 from Z).
//
// Reference: ~/SV/Z/R/src/branchB_chain_SD.cpp::between_step_SD (Z spec
// 2026-05-22_message-to-bgms-companion-SD-pivot.md §5).
//
// The marginalised-γ MH ratio collapses K_ij and uses the 1D conditional
// posterior density at zero (via Savage-Dickey) as the Bayes-factor input:
//
//   log_BF_1:0 = log π_enc(K_ij = 0 | Y, K_{-ij}) - log N(0; 0, σ²)
//
// where the numerator is the standalone primitive in ggm_sd. The cone-
// conditioning factor cancels between numerator and denominator at fixed α,
// so it is never computed.
//
//   ADD (γ_ij: 0 → 1): log α = log(p_inc/(1-p_inc)) - log_BF_1:0
//   DEL (γ_ij: 1 → 0): log α = log((1-p_inc)/p_inc) + log_BF_1:0
//
// Under prior_only_, log_BF_1:0 = 0 (no data) and the MH ratio is the prior
// odds alone — ergodic average of γ_ij should equal inclusion_probability_.
//
// On ADD accept K_ij is drawn from the 1D Laplace proposal N(x*, 1/κ) at
// K_{-ij} (or from the slab N(0, σ²) under prior_only_); on DEL accept K_ij
// is set to zero. K_jj is not touched in this step — the rank-2 (Σ, L)
// update reduces to dK_jj = 0.
// =====================================================================
void GGMModel::update_edge_indicator_parameter_pair_sd(size_t i, size_t j) {
    ensure_prior_params_extracted_();
    const int curr_g = edge_indicators_(i, j);

    // Call the SD primitive for the posterior Laplace (we always need x_-,
    // x_+ for the truncated proposal, and in posterior mode we also need
    // log_density for the BF). apply_pd_truncation=true so the returned
    // log_density is the truncated-Laplace integral.
    ggm_sd::SDResult sd = ggm_sd::density_at_zero_one(
        precision_matrix_,
        static_cast<int>(i), static_cast<int>(j),
        suf_stat_, static_cast<int>(n_),
        determinant_tilt_, prior_sigma_,
        /*nlo=*/true,
        /*apply_pd_truncation=*/true);
    if (sd.status == 1) return;  // K_0 not PD — chain stays put.
    if (!prior_only_ && (sd.status != 0 || !std::isfinite(sd.log_density))) {
        return;  // Laplace invalid in posterior mode → treat as reject.
    }

    // Derive log_BF and proposal parameters.
    // Prior-only: Z's shortcut log_BF = 0 — the marginal BF over K_{-ij} is
    // exactly 1 in expectation (SD identity with no data). With α independent
    // of K_{-ij}, the γ chain converges to Bernoulli(p_inc) regardless of how
    // K_ij is drawn on accept. The TN proposal at slab parameters (0, σ) is
    // pure PD bookkeeping.
    // Posterior: log_BF = sd.log_density - log_slab_at_0; the TN proposal
    // at the Laplace's (x*, 1/√κ) cancels the density's truncation factor
    // in the MH ratio.
    double log_BF_1_to_0;
    double proposal_mean;
    double proposal_sd;
    if (prior_only_) {
        log_BF_1_to_0 = 0.0;
        proposal_mean = 0.0;
        proposal_sd   = prior_sigma_;
    } else {
        const double log_slab_at_0 =
            -0.5 * std::log(2.0 * arma::datum::pi) - std::log(prior_sigma_);
        log_BF_1_to_0 = sd.log_density - log_slab_at_0;
        proposal_mean = sd.x_mode;
        proposal_sd   = 1.0 / std::sqrt(sd.curvature);
    }

    // ---- MH ratio (Z's form; truncation cancels into the density) -------
    const double p_inc = inclusion_probability_(i, j);
    double log_alpha;
    if (curr_g == 0) {
        log_alpha = MY_LOG(p_inc / (1.0 - p_inc)) - log_BF_1_to_0;
    } else {
        log_alpha = MY_LOG((1.0 - p_inc) / p_inc) + log_BF_1_to_0;
    }
    if (!std::isfinite(log_alpha)) return;
    if (MY_LOG(runif(rng_)) >= log_alpha) return;

    // ---- Accept --------------------------------------------------------
    const double omega_ij_old = precision_matrix_(i, j);
    const double omega_jj_old = precision_matrix_(j, j);

    if (curr_g == 0) {
        // ADD: sample K_ij from TN(proposal_mean, proposal_sd, x_-, x_+)
        // via inverse-CDF. By construction the draw is in (x_-, x_+) and
        // K stays PD, so the incremental cholupdate path is safe.
        const double F_lo = R::pnorm5(sd.x_minus, proposal_mean, proposal_sd, 1, 0);
        const double F_hi = R::pnorm5(sd.x_plus,  proposal_mean, proposal_sd, 1, 0);
        const double u    = runif(rng_);
        const double F_x  = F_lo + u * (F_hi - F_lo);
        double new_kij    = R::qnorm5(F_x, proposal_mean, proposal_sd, 1, 0);
        // Defensive clamp on the (rare) F_x → 1 case where qnorm returns inf.
        if (!std::isfinite(new_kij)) {
            new_kij = std::min(std::max(new_kij, sd.x_minus + 1e-12),
                               sd.x_plus - 1e-12);
        }
        precision_proposal_ = precision_matrix_;
        precision_proposal_(i, j) = new_kij;
        precision_proposal_(j, i) = new_kij;
        // K_jj unchanged.
        precision_matrix_(i, j) = new_kij;
        precision_matrix_(j, i) = new_kij;
        edge_indicators_(i, j) = 1;
        edge_indicators_(j, i) = 1;
    } else {
        // DEL: K_ij = 0 deterministically. K stays PD (zeroing an off-diag
        // of a PD matrix preserves PD via the 2x2 cofactor identity).
        precision_proposal_ = precision_matrix_;
        precision_proposal_(i, j) = 0.0;
        precision_proposal_(j, i) = 0.0;
        precision_matrix_(i, j) = 0.0;
        precision_matrix_(j, i) = 0.0;
        edge_indicators_(i, j) = 0;
        edge_indicators_(j, i) = 0;
    }

    // K stays PD by construction (TN proposal honours the PD cone), so the
    // existing rank-2 cholupdate + SMW Σ-update path is safe.
    cholesky_update_after_edge(omega_ij_old, omega_jj_old, i, j);

    constraint_dirty_ = true;
    theta_valid_ = false;
}


// =====================================================================
// L-space Savage-Dickey between-Γ MH (α = 1 closed-form Gibbs variant).
//
// Derivation: notes/2026-05-22_message-to-Z-companion-L-space-SD-derivation.md
//
// Per edge (i, j):
//   1. Σ_corner = Σ_{(i,j),(i,j)}; invert (2×2) to get Schur S.
//   2. l_ii = √S_11, l_ji = S_12/l_ii, l_jj = √(S_22 - S_12²/S_11).
//   3. a^T b = K_ij - l_ii · l_ji  ⟹  m_ij = -a^T b/l_ii = l_ji - K_ij/l_ii.
//   4. s_jj = K_jj - l_ji² - l_jj² = b^T b.
//   5. τ_post = l_ii²/σ² + 2β + S_jj   (S_jj entry of suf_stat_, not Schur S)
//   6. μ_post = (l_ii² m_ij/σ² - S_ij l_ii) / τ_post
//   7. log_BF_local = ½ log(τ_post/2π) - ½ τ_post (m_ij - μ_post)²
//                     - log l_ii - log N(0; 0, σ²)
//   8. MH:  log α(ADD) = log(p/(1-p)) - log_BF_local
//           log α(DEL) = log((1-p)/p) + log_BF_local
//   9. On accept (ADD): l_ji_new ~ N(μ_post, 1/τ_post)
//      On accept (DEL): l_ji_new = m_ij
//  10. Δl_ji = l_ji_new - l_ji
//      ΔK_ij = l_ii · Δl_ji
//      ΔK_jj = l_ji_new² - l_ji² = Δl_ji · (l_ji_new + l_ji)
//
// Update K (in place) and Σ via the existing rank-2 cholupdate+SMW path.
//
// Prior-only mode: zero out the S_ij and S_jj terms (set log_BF to come
// from the prior-only Laplace, which at α=1 is also closed-form). The
// marginal γ is Bernoulli(p_inc) under the encompassing prior with the
// |K|^δ tilt accounted for in the Gaussian conditional on l_ji.
// =====================================================================
void GGMModel::update_edge_indicator_parameter_pair_sd_lspace(size_t i, size_t j) {
    ensure_prior_params_extracted_();
    const int curr_g = edge_indicators_(i, j);
    const double p_inc = inclusion_probability_(i, j);
    const double sigma = prior_sigma_;
    const double beta  = prior_beta_;
    const double sigma2     = sigma * sigma;
    const double inv_sigma2 = 1.0 / sigma2;

    // Extract the trailing 2×2 (l_ii, l_ji) of the DEGORD-permuted Cholesky
    // factor via the Σ_BB block, read directly from the maintained
    // covariance_matrix_ (no triangular solves — the accept block below
    // refreshes Σ from K via arma::chol on every accept, so Σ is exactly in
    // sync with K).
    const double s_ii_cov = covariance_matrix_(i, i);
    const double s_jj_cov = covariance_matrix_(j, j);
    const double s_ij_cov = covariance_matrix_(i, j);
    const double det_cov  = s_ii_cov * s_jj_cov - s_ij_cov * s_ij_cov;
    if (!(det_cov > 0.0)) return;
    const double S_11 =  s_jj_cov / det_cov;
    const double S_12 = -s_ij_cov / det_cov;
    if (!(S_11 > 0.0)) return;
    const double l_ii = std::sqrt(S_11);
    const double l_ji = S_12 / l_ii;

    // ---- m_ij (Roverato slave in L-space) ----
    const double K_ij = precision_matrix_(i, j);
    const double m_ij = l_ji - K_ij / l_ii;

    // ---- Data terms ----
    // S_ij and S_jj here are entries of suf_stat_ (data); not the Schur matrix.
    const double S_ij_data = prior_only_ ? 0.0 : suf_stat_(i, j);
    const double S_jj_data = prior_only_ ? 0.0 : suf_stat_(j, j);
    const double n_eff     = prior_only_ ? 0.0 : static_cast<double>(n_);
    // Note: log_det_factor for the data piece in the conditional log-posterior
    // is (n/2) on log|K|, which is l_ji-independent. So the n_eff enters only
    // through S_ij_data * K_ij and S_jj_data * K_jj in tr(SK). In prior-only
    // mode we drop these.
    (void) n_eff;

    // ---- Branch by α: α=1 has closed-form Gaussian conditional (exact
    //      Gibbs proposal); α>1 uses Gauss-Hermite quadrature for the
    //      SD log-BF and Laplace-Gaussian as proposal with an explicit MH
    //      correction term log[π_target / q_proposal] evaluated at the
    //      sample point (ADD) or current state (DEL).
    double log_BF_1_to_0;
    double proposal_mu, proposal_sd;
    double A_post_save = 0.0, B_post_save = 0.0, s_jj_save = 0.0;  // for α>1 correction
    if (prior_alpha_ == 1.0) {
        const double tau_post  = l_ii * l_ii * inv_sigma2 + 2.0 * beta + S_jj_data;
        const double tau_prior = l_ii * l_ii * inv_sigma2 + 2.0 * beta;
        if (!(tau_post > 0.0) || !(tau_prior > 0.0)) return;
        const double mu_post   = (l_ii * l_ii * m_ij * inv_sigma2
                                  - S_ij_data * l_ii) / tau_post;
        const double mu_prior  = (l_ii * l_ii * m_ij * inv_sigma2) / tau_prior;
        const double log_pi_post  = 0.5 * std::log(tau_post  / (2.0 * arma::datum::pi))
                                  - 0.5 * tau_post  * (m_ij - mu_post)  * (m_ij - mu_post);
        const double log_pi_prior = 0.5 * std::log(tau_prior / (2.0 * arma::datum::pi))
                                  - 0.5 * tau_prior * (m_ij - mu_prior) * (m_ij - mu_prior);
        log_BF_1_to_0  = log_pi_post - log_pi_prior;
        proposal_mu    = mu_post;
        proposal_sd    = 1.0 / std::sqrt(tau_post);
    } else {
        // α > 1: GH quadrature for log_BF (reliable across cells), Laplace
        // for the proposal mode/curvature (cheap; if Laplace fails, use the
        // pure-Gaussian fallback B/(2A), 1/√(2A)).
        const double A_post  = 0.5 * (l_ii * l_ii * inv_sigma2
                                      + 2.0 * beta + S_jj_data);
        const double B_post  = l_ii * l_ii * m_ij * inv_sigma2
                               - S_ij_data * l_ii;
        const double A_prior = 0.5 * (l_ii * l_ii * inv_sigma2 + 2.0 * beta);
        const double B_prior = l_ii * l_ii * m_ij * inv_sigma2;
        const double s_jj    = precision_matrix_(j, j) - l_ji * l_ji;
        if (!(s_jj > 0.0) || !(A_post > 0.0) || !(A_prior > 0.0)) return;
        const auto gh_post  = ggm_sd::density_at_l_ji_gh(
            m_ij, A_post,  B_post,  s_jj, prior_alpha_);
        const auto gh_prior = ggm_sd::density_at_l_ji_gh(
            m_ij, A_prior, B_prior, s_jj, prior_alpha_);
        if (gh_post.status != 0 || gh_prior.status != 0
            || !std::isfinite(gh_post.log_density)
            || !std::isfinite(gh_prior.log_density)) {
            ++n_pd_reverts_;
            return;
        }
        log_BF_1_to_0 = gh_post.log_density - gh_prior.log_density;

        const auto lp_post = ggm_sd::density_at_l_ji_one(
            0.0, A_post, B_post, s_jj, prior_alpha_, /*nlo=*/false);
        if (lp_post.status == 0
            && std::isfinite(lp_post.curvature) && lp_post.curvature > 0.0) {
            proposal_mu = lp_post.x_mode;
            proposal_sd = 1.0 / std::sqrt(lp_post.curvature);
        } else {
            proposal_mu = B_post / (2.0 * A_post);
            proposal_sd = 1.0 / std::sqrt(2.0 * A_post);
        }
        A_post_save = A_post;
        B_post_save = B_post;
        s_jj_save   = s_jj;
    }

    // ---- MH ratio (proposal MH correction added below at α>1) ----
    double log_alpha;
    if (curr_g == 0) {
        log_alpha = MY_LOG(p_inc / (1.0 - p_inc)) - log_BF_1_to_0;
    } else {
        log_alpha = MY_LOG((1.0 - p_inc) / p_inc) + log_BF_1_to_0;
    }

    // ---- Sample proposal ----
    double l_ji_new;
    if (curr_g == 0) {
        // ADD: draw from proposal Gaussian. At α=1 this is exact Gibbs; at
        // α>1 it's the Laplace approximation with an explicit MH correction.
        l_ji_new = rnorm(rng_, proposal_mu, proposal_sd);
    } else {
        // DEL: deterministic slave.
        l_ji_new = m_ij;
    }

    // ---- MH proposal correction at α>1 ----
    // ADD: log_α += log π_target(l_ji_new) − log q_proposal(l_ji_new)
    // DEL: log_α −= log π_target(l_ji_curr) − log q_proposal(l_ji_curr)
    // (At α=1 the proposal IS the conditional, correction ≡ 0.)
    if (prior_alpha_ != 1.0) {
        const double x_corr  = (curr_g == 0) ? l_ji_new : l_ji;
        const auto gh_corr   = ggm_sd::density_at_l_ji_gh(
            x_corr, A_post_save, B_post_save, s_jj_save, prior_alpha_);
        if (gh_corr.status != 0 || !std::isfinite(gh_corr.log_density)) {
            ++n_pd_reverts_;
            return;
        }
        const double dx     = (x_corr - proposal_mu) / proposal_sd;
        const double log_q  = -0.5 * std::log(2.0 * arma::datum::pi)
                              - std::log(proposal_sd)
                              - 0.5 * dx * dx;
        const double correction = gh_corr.log_density - log_q;
        log_alpha += (curr_g == 0 ? +correction : -correction);
    }

    if (!std::isfinite(log_alpha)) return;
    if (MY_LOG(runif(rng_)) >= log_alpha) return;

    // Translate Δl_ji into ΔK_ij, ΔK_jj.
    const double d_l_ji = l_ji_new - l_ji;
    const double K_ij_new = K_ij + l_ii * d_l_ji;
    const double K_jj_new = precision_matrix_(j, j) + (l_ji_new + l_ji) * d_l_ji;

    // Accept: install the new K and rebuild L + Σ from scratch via
    // arma::chol. O(p³) per accept, but the only path proven correct.
    // Both architectural alternatives we tried — SMW-only on Σ, and the
    // drift-bounded rank-2 chol-up on L — produce K that disagrees with
    // L L^T after enough accepts, and the next periodic refresh throws on
    // a numerically non-PD K. The disagreement is unavoidable as long as
    // K is updated algebraically and L (or Σ) is updated by a non-exact
    // rank-r transform.
    const double omega_ij_old = precision_matrix_(i, j);
    const double omega_jj_old = precision_matrix_(j, j);
    precision_matrix_(i, j) = K_ij_new;
    precision_matrix_(j, i) = K_ij_new;
    precision_matrix_(j, j) = K_jj_new;

    if (curr_g == 0) {
        edge_indicators_(i, j) = 1;
        edge_indicators_(j, i) = 1;
    } else {
        edge_indicators_(i, j) = 0;
        edge_indicators_(j, i) = 0;
    }

    // PD canary via two-arg arma::chol (returns bool, no throw).
    arma::mat U_check;
    if (!arma::chol(U_check, precision_matrix_, "upper")) {
        ++n_pd_reverts_;
        precision_matrix_(i, j) = omega_ij_old;
        precision_matrix_(j, i) = omega_ij_old;
        precision_matrix_(j, j) = omega_jj_old;
        if (curr_g == 0) {
            edge_indicators_(i, j) = 0;
            edge_indicators_(j, i) = 0;
        } else {
            edge_indicators_(i, j) = 1;
            edge_indicators_(j, i) = 1;
        }
        return;
    }
    cholesky_of_precision_ = U_check;
    arma::solve(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_),
                arma::eye(p_, p_), arma::solve_opts::fast);
    covariance_matrix_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    constraint_dirty_ = true;
    theta_valid_ = false;
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
}


// V/RR U-pool refresh — only relevant when between-Γ moves are active.
// The chain runner calls this iff WarmupSchedule::u_refresh_enabled(iter)
// is true (i.e., stage 3c + sampling), so this method itself doesn't
// re-check the schedule — it just executes whichever refresh strategy is
// configured (mh_U or legacy fresh-from-prior). Running this earlier
// would let U/K_depth drift via PMMH dynamics with no Γ moves consuming
// it, polluting the sampling phase with inflated K_depth.
void GGMModel::refresh_auxiliary_u() {
    if (graph_prior_spec_ != GraphPriorSpec::Hierarchical) return;
    // Plug-in mode skips the U machinery entirely — the closed-form mNLO
    // ratio carries the hierarchical correction deterministically.
    if (plug_in_nlo_) return;

    bool was_built = hierarchical_state_built_;
    ensure_hierarchical_state_();
    if (!was_built) {
        // First refresh after lazy init: state was just seeded with a
        // fresh U; no MH comparison is possible yet.
        return;
    }
    if (mh_U_) {
        mh_on_U_step_();
    } else {
        refresh_z_ratio_pool_();
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


void GGMModel::ensure_prior_params_extracted_() {
    if (prior_params_extracted_) return;
    // Validate prior family. The closed-form log_Z_NLO_gamma and SD
    // density-at-zero machinery only covers slab = Normal(0, σ) and
    // diag = Gamma(α, β) on K_ii/2.
    const auto* slab = dynamic_cast<const NormalPrior*>(interaction_prior_.get());
    if (slab == nullptr)
        throw std::runtime_error(
            "This code path requires a Normal slab (NormalPrior). "
            "Re-fit with interaction_prior_type = 'normal'.");
    const auto* diag = dynamic_cast<const GammaScalePrior*>(diagonal_prior_.get());
    if (diag == nullptr)
        throw std::runtime_error(
            "This code path requires a Gamma diagonal prior "
            "(GammaScalePrior).");
    prior_sigma_ = slab->scale();
    prior_alpha_ = diag->shape();
    prior_beta_  = diag->rate();
    prior_params_extracted_ = true;
}


void GGMModel::ensure_hierarchical_state_() {
    if (hierarchical_state_built_) return;
    ensure_prior_params_extracted_();
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
// Dispatcher: routes to local-K or global proposal depending on flag and
// mixture frequency. Bit-equality with the pre-flag binary is preserved
// when mh_U_local_K_ == false (no RNG draws on the new branch).
void GGMModel::mh_on_U_step_() {
    if (graph_prior_spec_ != GraphPriorSpec::Hierarchical) return;

    if (mh_U_local_K_) {
        // Local + global mixture. The global slice (fresh-from-prior K)
        // keeps the long-jump escape route alive at low frequency.
        bool do_global = (runif(rng_) < mh_U_local_K_global_freq_);
        if (do_global) {
            ++n_mh_U_local_global_steps_;
            mh_on_U_step_global_();
        } else {
            mh_on_U_step_local_K_();
        }
        return;
    }
    mh_on_U_step_global_();
}


void GGMModel::mh_on_U_step_global_() {
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

    // Auto-reject only on degenerate sentinels (Lyne 2015: chain targets
    // π·|V|·μ·P(N), sign flips are tracked not rejected).
    if (!std::isfinite(V_old.first) || !std::isfinite(V_new.first)) {
        ++n_mh_U_nonfinite_;
        return;
    }
    if (V_old.second == 0 || V_new.second == 0) {
        ++n_mh_U_signzero_;
        return;
    }
    if (V_old.second != V_new.second) {
        ++n_mh_U_signflip_;  // diagnostic only; sign flips no longer reject
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


// Local random-walk on K_depth with reflection at 0. The geometric prior
// P(K=k) = (1−ρ)·ρ^k anchors the chain near small K via the prior ratio
// ρ^(K_new − K_old) in the MH numerator; this attacks the K-dwell trap
// that fresh-from-prior K proposals cannot escape efficiently.
//
// Pool reconstruction is product-form: shrink drops the last pool, grow
// appends one fresh draw. The auxiliary draw's μ density cancels with the
// corresponding μ factor in the augmented target — no Jacobian.
//
// Boundary corrections (q proposal asymmetry under reflection at 0):
//   K_old = 0 → K_new = 1 (forced):       log_q_ratio = log(1/2)
//   K_old = 1 → K_new = 0 (down, p=1/2):  log_q_ratio = log(2)
//   interior:                              log_q_ratio = 0
void GGMModel::mh_on_U_step_local_K_() {
    ++n_mh_U_attempts_;

    // V_old fetch — same cache logic as the global path.
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

    // Propose K_new
    const int K_old = v_K_depth_;
    int    K_new;
    double log_q_ratio = 0.0;
    bool   moves_up;
    if (K_old == 0) {
        K_new   = 1;
        moves_up = true;
        log_q_ratio = -MY_LOG(2.0);  // forced up; reverse picks down with p=1/2
    } else {
        moves_up = (runif(rng_) < 0.5);
        K_new   = moves_up ? K_old + 1 : K_old - 1;
        if (K_old == 1 && K_new == 0) {
            log_q_ratio = MY_LOG(2.0);  // reverse from K=0 is forced up
        }
    }
    if (moves_up) ++n_mh_U_local_up_attempts_;
    else          ++n_mh_U_local_down_attempts_;

    // Construct pools_new: shrink drops last; grow appends fresh.
    std::vector<arma::mat> pools_new;
    if (moves_up) {
        pools_new = v_pools_t_;
        pools_new.push_back(degord::draw_bartlett_pool(
            rng_, static_cast<int>(p_), v_M_inner_));
    } else {
        pools_new.assign(
            v_pools_t_.begin(),
            v_pools_t_.begin() + static_cast<std::ptrdiff_t>(K_new));
    }

    arma::ivec pi_vec = degord::degord_permutation(
        static_cast<int>(p_), pi_i, pi_j);
    arma::imat G_pi_curr = degord::permute_graph(edge_indicators_, pi_vec);
    auto V_new = degord::V_log_at_Gamma_pi_degord(
        K_new, pools_new, G_pi_curr, chain_aux_degord_, log_c, v_rho_);

    if (!std::isfinite(V_old.first) || !std::isfinite(V_new.first)) {
        ++n_mh_U_nonfinite_; return;
    }
    if (V_old.second == 0 || V_new.second == 0) {
        ++n_mh_U_signzero_; return;
    }
    if (V_old.second != V_new.second) {
        ++n_mh_U_signflip_;
    }

    // MH ratio: log|V_new/V_old| + log(P(K_new)/P(K_old)) + log(q_rev/q_fwd)
    double log_alpha = V_new.first - V_old.first
                     + static_cast<double>(K_new - K_old) * std::log(v_rho_)
                     + log_q_ratio;

    if (MY_LOG(runif(rng_)) < log_alpha) {
        v_pools_t_ = std::move(pools_new);
        v_K_depth_ = K_new;
        current_log_abs_V_ = V_new.first;
        current_sign_V_    = V_new.second;
        v_diag_initialized_ = true;
        ++n_mh_U_accepts_;
        if (moves_up) ++n_mh_U_local_up_accepts_;
        else          ++n_mh_U_local_down_accepts_;
    }
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
        if (use_sd_between_step_) {
            if (use_sd_lspace_) {
                update_edge_indicator_parameter_pair_sd_lspace(i, j);
            } else {
                update_edge_indicator_parameter_pair_sd(i, j);
            }
        } else {
            update_edge_indicator_parameter_pair(i, j);
        }
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
