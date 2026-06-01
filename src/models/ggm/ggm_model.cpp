#include "models/ggm/ggm_model.h"
#include "rng/rng_utils.h"
#include "math/explog_macros.h"
#include "math/cholupdate.h"
#include "mcmc/execution/step_result.h"
#include "mcmc/execution/warmup_schedule.h"
#include "math/savage_dickey/sinh_midpoint.h"
#include "math/savage_dickey/cubic_mode.h"
#include "math/savage_dickey/cauchy_omega.h"
#include <stdexcept>

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
        // givens_qr now consumes A_q directly and stores R in transposed
        // orientation; Q is still (q × q), unaffected by R's layout.
        GGMGradientEngine::givens_qr(Aq_buf, Q_tmp, R_tmp, R_diag, rots_tmp);
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
    // the lemma. The formula's `cc11*cc22 - cc12²` equals MINUS the true
    // determinant ratio (sign-flipped by the lemma's convention here), so
    // std::abs recovers |K*|/|K|. PD-violation protection is provided by the
    // explicit arma::chol canaries in update_edge_parameter and
    // update_diagonal_parameter (prior_only mode), not by detecting the sign
    // here.
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
    // Rank-1 specialisation of log_det_ratio_edge (Ui2 = 0). Same sign-flip
    // convention; std::abs recovers the true |K*|/|K|.
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

    apply_rank2_chol_smw_update_();

    // reset for next iteration
    vf1_[i] = 0.0;
    vf1_[j] = 0.0;
    vf2_[i] = 0.0;
    vf2_[j] = 0.0;

}

void GGMModel::apply_rank2_chol_smw_update_()
{
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
        const double inv_c00 =  c22 / det;
        const double inv_c11 =  c11 / det;
        const double inv_c01 = -c12 / det;
        // ΔΣ = -inv_c00 · a1 a1ᵀ - inv_c11 · a2 a2ᵀ - inv_c01 · (a1 a2ᵀ + a2 a1ᵀ)
        // Refactor to two rank-1 outer products with bundled scalar weights:
        //   ΔΣ = -(a1 b1ᵀ + a2 b2ᵀ),  with
        //     b1 = inv_c00 · a1 + inv_c01 · a2
        //     b2 = inv_c01 · a1 + inv_c11 · a2
        // Saves 2 of the 4 arma p×p temporaries the 3-line form materialised.
        const arma::vec b1 = inv_c00 * a1 + inv_c01 * a2;
        const arma::vec b2 = inv_c01 * a1 + inv_c11 * a2;
        covariance_matrix_ -= a1 * b1.t();
        covariance_matrix_ -= a2 * b2.t();
    }
}

void GGMModel::set_within_step_kind(WithinStepKind kind) {
    if (kind == WithinStepKind::RowBlockGibbs && !row_block_gibbs_eligible()) {
        // Warn-and-downgrade: opt-in scripts under unusual prior combinations
        // (Cauchy slab, Gamma α ≠ 1, δ ≠ 0) keep running on the AM path
        // until the planned PR-4..PR-6 coverage lands.
        Rcpp::warning(
            "within_step_kind = 'row_block_gibbs' requested but the prior "
            "configuration is not yet supported (need Normal slab, Gamma "
            "alpha=1, delta=0). Falling back to 'adaptive_metropolis'.");
        within_step_kind_ = WithinStepKind::AdaptiveMetropolis;
        return;
    }
    within_step_kind_ = kind;
}

bool GGMModel::row_block_gibbs_eligible() {
    // Normal slab + Gamma(alpha=1, .) on K_ii/2 + zero determinant tilt is
    // the exact-conjugate scope of PR-2. Other prior families and alpha/delta
    // values land in PR-4..PR-6 as post-step corrections.
    if (dynamic_cast<const NormalPrior*>(interaction_prior_.get()) == nullptr)
        return false;
    const auto* diag = dynamic_cast<const GammaScalePrior*>(diagonal_prior_.get());
    if (diag == nullptr) return false;
    if (std::abs(diag->shape() - 1.0) > 1e-12) return false;
    if (determinant_tilt_ != 0.0) return false;
    return true;
}

void GGMModel::update_row_block_gibbs(size_t i) {
    // Conjugate Gaussian–Gamma draw of (β = K_{N_i, i}, kii = K_{i,i}) given
    // A = K_{-i, -i}, S, and the Normal-slab × Gamma(α=1) prior. δ = 0 here;
    // PR-4 will add the determinant-tilt shape shift.
    ensure_prior_params_extracted_();
    const double sigma = prior_sigma_;     // Normal slab std on K_yy_ij = -K_ij/2
    const double beta0 = prior_beta_;      // Gamma rate on K_ii/2
    const double s_ii  = suf_stat_(i, i);

    // Active neighbour set N_i (in row order).
    std::vector<size_t> Ni;
    Ni.reserve(p_ - 1);
    for (size_t k = 0; k < p_; ++k) {
        if (k != i && edge_indicators_(i, k) == 1) Ni.push_back(k);
    }
    const size_t q = Ni.size();

    // Stash old K column entries so the rank-2 update can encode the delta.
    const double kii_old = precision_matrix_(i, i);
    arma::vec beta_old(q);
    for (size_t k = 0; k < q; ++k) beta_old(k) = precision_matrix_(i, Ni[k]);

    // ξ shape α = 1, δ = 0 → n/2 + 1.
    const double xi_shape = static_cast<double>(n_) / 2.0 + 1.0;
    const double xi_rate  = (beta0 + s_ii) / 2.0;

    arma::vec beta_new(q, arma::fill::zeros);
    double kii_new;

    if (q == 0) {
        // No active neighbours: K_{i,i} = ξ, no β draw.
        kii_new = rgamma(rng_, xi_shape, xi_rate);
    } else {
        // C = (A⁻¹)_{N_i, N_i} via Schur on Σ:
        //   C_{kl} = Σ_{N_i[k], N_i[l]} − Σ_{N_i[k], i} Σ_{i, N_i[l]} / Σ_{ii}
        const double sigma_ii = covariance_matrix_(i, i);
        arma::vec sigma_iNi(q);
        for (size_t k = 0; k < q; ++k) sigma_iNi(k) = covariance_matrix_(i, Ni[k]);
        arma::mat C(q, q);
        for (size_t k = 0; k < q; ++k) {
            for (size_t l = 0; l < q; ++l) {
                C(k, l) = covariance_matrix_(Ni[k], Ni[l])
                          - sigma_iNi(k) * sigma_iNi(l) / sigma_ii;
            }
        }

        // M = (β₀ + S_ii) C + (1/(4σ²)) I — symmetric positive-definite.
        const double inv_4sig2 = 1.0 / (4.0 * sigma * sigma);
        arma::mat M = (beta0 + s_ii) * C;
        for (size_t k = 0; k < q; ++k) M(k, k) += inv_4sig2;

        arma::mat L_M;
        if (!arma::chol(L_M, M, "lower")) {
            // M is PD by construction (C PD as submatrix of A⁻¹, prior precision > 0).
            // A failure here means numerical trouble; skip this row's update
            // rather than corrupt K.
            return;
        }

        // S_{N_i, i} vector.
        arma::vec s_Ni_i(q);
        for (size_t k = 0; k < q; ++k) s_Ni_i(k) = suf_stat_(Ni[k], i);

        // Mean μ = −M⁻¹ S_{N_i, i}. Two triangular solves: L y = −S, Lᵀ μ = y.
        arma::vec y  = arma::solve(arma::trimatl(L_M), -s_Ni_i);
        arma::vec mu = arma::solve(arma::trimatu(L_M.t()), y);

        // β = μ + Lᵀ⁻¹ z, z ~ N(0, I_q): one triangular solve.
        arma::vec z = arma_rnorm_vec(rng_, q);
        arma::vec w = arma::solve(arma::trimatu(L_M.t()), z);
        beta_new = mu + w;

        // ξ ~ Gamma; K_{i,i} = ξ + βᵀ C β.
        const double xi = rgamma(rng_, xi_shape, xi_rate);
        kii_new = xi + arma::as_scalar(beta_new.t() * C * beta_new);
    }

    // Write the new K column / row, then apply the symmetric rank-2 update
    // to chol(K) and Σ. The rank-2 decomposition is
    //   ΔK = e_i vf2ᵀ + vf2 e_iᵀ,
    //   (vf2)_i      = (kii_new − kii_old) / 2,
    //   (vf2)_{N_i} = β_new − β_old,
    //   (vf2)_k      = 0  otherwise
    // which reproduces the sparse column-change at row/col i without touching
    // the (k, l) entries off row i.
    precision_matrix_(i, i) = kii_new;
    for (size_t k = 0; k < q; ++k) {
        precision_matrix_(i, Ni[k]) = beta_new(k);
        precision_matrix_(Ni[k], i) = beta_new(k);
    }

    vf1_[i] = 1.0;
    vf2_[i] = (kii_new - kii_old) / 2.0;
    for (size_t k = 0; k < q; ++k) vf2_[Ni[k]] = beta_new(k) - beta_old(k);

    apply_rank2_chol_smw_update_();

    vf1_[i] = 0.0;
    vf2_[i] = 0.0;
    for (size_t k = 0; k < q; ++k) vf2_[Ni[k]] = 0.0;
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

    // prior_only_: skip likelihood, but explicitly PD-check the proposal —
    // K_ii > 0 by construction is not sufficient (a small K_ii relative to
    // off-diagonals can still violate PD). Without this canary the chain
    // accumulates non-PD K accepts at shallow det-tilt (cf. 2026-05-23).
    double ln_alpha = prior_only_ ? 0.0 : log_density_impl_diag(i);
    if (prior_only_) {
        arma::mat R_chk;
        if (!arma::chol(R_chk, precision_proposal_)) return 0.0;
    }

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
    // bgms convention: NormalPrior(scale = σ_user) is on K_yy = -K/2, so the
    // slab on K_ij has sd σ_K = 2·σ_user (variance 4σ_user²). The SD primitives
    // in math/savage_dickey/ take the slab as N(0, σ²) on K_ij directly, so
    // we pass σ_K = 2·prior_sigma_ to match the within-K parameterization.
    //
    // Cauchy slab: K_ij | gamma=1, omega ~ N(0, sigma^2 * omega_ij) with
    // omega_ij ~ IG(1/2, 1/2) marginalised to Cauchy(0, sigma). The kernel
    // structure is unchanged; the per-edge omega_(i,j) just rescales the
    // effective slab variance. Normal slab has omega_ij == 1.0 from the
    // constructor so this multiplication is a no-op there.
    const double sigma = 2.0 * prior_sigma_;
    const double beta  = prior_beta_;
    const double omega_ij_slab = slab_is_cauchy_ ? omega_(i, j) : 1.0;
    const double sigma2     = sigma * sigma * omega_ij_slab;
    const double inv_sigma2 = 1.0 / sigma2;

    // Extract the trailing 2×2 (l_ii, l_ji) of the DEGORD-permuted Cholesky
    // factor via the Σ_BB block. Σ_corner is computed on demand as
    // Σ_kl = inv(L)[k,:] · inv(L)[l,:]^T (one dot product per entry). The
    // accept block keeps inv(L) fresh per accept (chol + triangular solve)
    // and skips the full Σ matmul; Σ is refreshed once at the end of
    // update_edge_indicators() for the within-step's reads.
    const arma::rowvec ri = inv_cholesky_of_precision_.row(i);
    const arma::rowvec rj = inv_cholesky_of_precision_.row(j);
    const double s_ii_cov = arma::dot(ri, ri);
    const double s_jj_cov = arma::dot(rj, rj);
    const double s_ij_cov = arma::dot(ri, rj);
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
    //
    // Diag prior convention: bgms evaluates diag at K_jj/2 ~ Gamma(α, β),
    // implying K_jj ~ Gamma(α, β/2). Under Roverato slaving K_jj = c_3 + l²,
    // the resulting contribution to the kernel -A l² + B l is A = β/2, so
    // tau = 2A picks up β (not 2β). Match the within-K convention by using
    // `beta` rather than `2β` below. (Companion to the slab σ_K = 2σ_user
    // fix at line 977.)
    double log_BF_1_to_0;
    double proposal_mu, proposal_sd;
    double A_post_save = 0.0, B_post_save = 0.0, s_jj_save = 0.0;  // for α>1 correction
    if (prior_alpha_ == 1.0) {
        const double tau_post  = l_ii * l_ii * inv_sigma2 + beta + S_jj_data;
        const double tau_prior = l_ii * l_ii * inv_sigma2 + beta;
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
        // α > 1: sinh-substitution + midpoint-rule quadrature
        // (math/savage_dickey/sinh_midpoint.h) for the local Bayes factor. The
        // substitution phi = sqrt(s) sinh(t) pins the integrand's
        // branch points at t = +- i pi/2 (constant distance from the
        // real axis, independent of s), giving uniform exponential
        // convergence across the (A, B, s, alpha) plane. At N = 128
        // the worst-case error in the chain-realistic regime is
        // ~1e-10 (machine roundoff floor).
        //
        // The proposal Gaussian uses the cubic critical-point solver
        // (math/savage_dickey/cubic_mode.h) for closed-form mode and
        // curvature in O(1) -- no Newton iteration, no saddle / triple-root
        // failure modes.
        const double A_post  = 0.5 * (l_ii * l_ii * inv_sigma2
                                      + beta + S_jj_data);
        const double B_post  = l_ii * l_ii * m_ij * inv_sigma2
                               - S_ij_data * l_ii;
        const double A_prior = 0.5 * (l_ii * l_ii * inv_sigma2 + beta);
        const double B_prior = l_ii * l_ii * m_ij * inv_sigma2;
        const double s_jj    = precision_matrix_(j, j) - l_ji * l_ji;
        if (!(s_jj > 0.0) || !(A_post > 0.0) || !(A_prior > 0.0)) return;
        const auto sinh_post  = savage_dickey::density_at_l_ji_sinh(
            m_ij, A_post,  B_post,  s_jj, prior_alpha_);
        const auto sinh_prior = savage_dickey::density_at_l_ji_sinh(
            m_ij, A_prior, B_prior, s_jj, prior_alpha_);
        if (sinh_post.status != 0 || sinh_prior.status != 0
            || !std::isfinite(sinh_post.log_density)
            || !std::isfinite(sinh_prior.log_density)) {
            ++n_pd_reverts_;
            return;
        }
        log_BF_1_to_0 = sinh_post.log_density - sinh_prior.log_density;

        // Proposal Gaussian: closed-form mode + curvature via the
        // critical-point cubic. The cubic returns every real root in
        // O(1) and labels modes by ell_pp sign; we use the global mode.
        // If the cubic flags no usable mode (PD revert, near-zero
        // curvature), fall back to the alpha = 1 reference Gaussian
        // (B_post/(2 A_post), 2 A_post).
        const auto cubic_post = savage_dickey::solve_sd_cubic(
            A_post, B_post, s_jj, prior_alpha_);
        const double kappa_floor = 1e-10 * 2.0 * A_post;
        if (cubic_post.status == 0
            && cubic_post.n_modes >= 1
            && cubic_post.global_mode_index >= 0
            && cubic_post.roots[cubic_post.global_mode_index].curvature
                > kappa_floor) {
            const auto& m = cubic_post.roots[cubic_post.global_mode_index];
            proposal_mu = m.phi;
            proposal_sd = 1.0 / std::sqrt(m.curvature);
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
    //
    // Uses the same sinh-substitution + midpoint-rule primitive as the
    // BF computation; log_Z is identical between this call and the
    // earlier sinh_post (deterministic in A, B, s_jj, alpha), so this
    // is effectively a kernel re-evaluation at x_corr.
    if (prior_alpha_ != 1.0) {
        const double x_corr = (curr_g == 0) ? l_ji_new : l_ji;
        const auto sinh_corr = savage_dickey::density_at_l_ji_sinh(
            x_corr, A_post_save, B_post_save, s_jj_save, prior_alpha_);
        if (sinh_corr.status != 0 || !std::isfinite(sinh_corr.log_density)) {
            ++n_pd_reverts_;
            return;
        }
        const double dx     = (x_corr - proposal_mu) / proposal_sd;
        const double log_q  = -0.5 * std::log(2.0 * arma::datum::pi)
                              - std::log(proposal_sd)
                              - 0.5 * dx * dx;
        const double correction = sinh_corr.log_density - log_q;
        log_alpha += (curr_g == 0 ? +correction : -correction);
    }

    // Translate Δl_ji into ΔK_ij, ΔK_jj. Compute before the MH check so we
    // can evaluate the log|K| ratio for the MH correction below.
    const double d_l_ji = l_ji_new - l_ji;
    const double K_ij_new = K_ij + l_ii * d_l_ji;
    const double K_jj_new = precision_matrix_(j, j) + (l_ji_new + l_ji) * d_l_ji;

    // ---- log|K| MH correction (likelihood log-det + determinant tilt) ----
    // The L-space SD primitive's kernel f(l_ji) = -A l_ji² + B l_ji
    //   + (α-1) log(s_jj + l_ji²) targets π_FALSE = π_TRUE / |K(l_ji)|^{n/2+δ}.
    // (The (n/2) log|K| comes from the data likelihood; the δ log|K| from the
    // tilt. Both depend on l_ji via K_ij(l_ji) and K_jj(l_ji), but neither is
    // in f.) Shifting the chain to target π_TRUE via standard MH theory:
    //   log α += [log π_TRUE(after) − log π_TRUE(before)]
    //          − [log π_FALSE(after) − log π_FALSE(before)]
    //          = (n/2 + δ) [log|K_after| − log|K_before|].
    //
    // Implementation: inline the rank-2 matrix-determinant lemma (the same
    // formula as log_det_ratio_edge but with PD detection via the sign of
    // the 2×2 determinant — negative = PD under the lemma's sign-flip
    // convention, positive = K_proposed non-PD). This avoids an arma::chol
    // per SD attempt and short-circuits on PD violations.
    const double log_det_factor =
        0.5 * static_cast<double>(n_) + determinant_tilt_;
    if (log_det_factor != 0.0) {
        const double Ui2 = K_ij - K_ij_new;
        const double Uj2 = (precision_matrix_(j, j) - K_jj_new) * 0.5;
        // Reuse the Σ corner already computed above (s_ii_cov, s_jj_cov,
        // s_ij_cov) so we never read covariance_matrix_ here. log_det MH
        // agrees with the log_BF moments computed at the start of the attempt.
        const double sii = s_ii_cov;
        const double sjj = s_jj_cov;
        const double sij = s_ij_cov;
        const double cc11 = sjj;
        const double cc12 = 1.0 - (sij * Ui2 + sjj * Uj2);
        const double cc22 = Ui2 * Ui2 * sii + 2.0 * Ui2 * Uj2 * sij
                          + Uj2 * Uj2 * sjj;
        const double det2x2 = cc11 * cc22 - cc12 * cc12;
        if (det2x2 >= 0.0) {
            // Under the lemma's sign-flip convention a non-negative value
            // signals K_proposed is non-PD; reject without further work.
            ++n_pd_reverts_;
            return;
        }
        log_alpha += log_det_factor * MY_LOG(-det2x2);
    }

    if (!std::isfinite(log_alpha)) return;
    if (MY_LOG(runif(rng_)) >= log_alpha) return;

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

    // PD canary via two-arg arma::chol + min-pivot guard. arma::chol returns
    // false only for strict non-PD; barely-PD K (one pivot ~ 1e-10) still
    // passes but its κ(K) is huge, and the next within-K NUTS step computes
    // K^{-1} entries on the order of 1/min_pivot — gradient overflow, theta
    // becomes inf, K becomes garbage (we saw K_jj samples ~ 1e13 at
    // extreme priors). Treat min_diag(U_check) < kSDMinCholDiag as a PD
    // violation and revert. The threshold corresponds to a minimum K
    // eigenvalue of ~kSDMinCholDiag² (= 1e-12 at 1e-6), well clear of
    // double-precision noise but loose enough not to perturb well-conditioned
    // chains.
    constexpr double kSDMinCholDiag = 1.0e-6;
    arma::mat U_check;
    if (!arma::chol(U_check, precision_matrix_, "upper")
        || U_check.diag().min() < kSDMinCholDiag) {
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
    // Σ is intentionally NOT refreshed here — Σ_corner reads inside the
    // sweep go directly off inv(L) row inner products, and update_edge_
    // indicators() refreshes covariance_matrix_ once at the end of the
    // between sweep for the within-step's reads.
    constraint_dirty_ = true;
    theta_valid_ = false;
}


void GGMModel::do_one_metropolis_step(int iteration) {
    if (within_step_kind_ == WithinStepKind::RowBlockGibbs) {
        // Conjugate row sweep — no acceptance, no per-slot adaptation. The
        // within-slot entries of proposal_sds_ become inert (no mask bits
        // are ever set under this branch); between-step adaptation still
        // fires from its own driver.
        for (size_t i = 0; i < p_; ++i) {
            update_row_block_gibbs(i);
        }
        check_and_refresh_if_drift_();
        return;
    }

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


// ----------------------------------------------------------------------
// Hierarchical-spec configuration
// ----------------------------------------------------------------------

void GGMModel::set_graph_prior_spec(GraphPriorSpec spec) {
    graph_prior_spec_ = spec;
}


void GGMModel::ensure_prior_params_extracted_() {
    if (prior_params_extracted_) return;
    // Validate prior family. The SD density-at-zero machinery covers the
    // Normal slab directly and the Cauchy slab via the scale-mixture
    // representation (K_ij | omega ~ N(0, sigma^2 * omega); omega ~ IG(1/2, 1/2)).
    // Diag must be Gamma(alpha, beta) on K_ii/2.
    if (const auto* slab_n = dynamic_cast<const NormalPrior*>(interaction_prior_.get())) {
        prior_sigma_     = slab_n->scale();
        slab_is_cauchy_  = false;
    } else if (const auto* slab_c = dynamic_cast<const CauchyPrior*>(interaction_prior_.get())) {
        prior_sigma_     = slab_c->scale();
        slab_is_cauchy_  = true;
    } else {
        throw std::runtime_error(
            "This code path requires a Normal or Cauchy slab. "
            "Re-fit with interaction_prior_type = 'normal' or 'cauchy'.");
    }
    const auto* diag = dynamic_cast<const GammaScalePrior*>(diagonal_prior_.get());
    if (diag == nullptr)
        throw std::runtime_error(
            "This code path requires a Gamma diagonal prior "
            "(GammaScalePrior).");
    prior_alpha_ = diag->shape();
    prior_beta_  = diag->rate();
    prior_params_extracted_ = true;
}


// Thin GGM-side wrapper around the model-agnostic Cauchy omega primitive
// in math/savage_dickey/cauchy_omega.h. Extracts (l_ii, m_ij) from the
// maintained Cholesky/covariance state, computes the GGM-specific
// diag-prior contributions
//
//     A_diag_K =  beta / (2 * l_ii^2),
//     B_K      = -beta * m_ij / l_ii,
//
// and forwards to savage_dickey::slice_sample_cauchy_omega_active. At
// gamma = 0 the slab carries no information about omega; we sample from
// the IG(1/2, 1/2) prior via sample_cauchy_omega_prior.
void GGMModel::slice_sample_omega_ij_(size_t i, size_t j) {
    if (!slab_is_cauchy_) return;

    if (edge_indicators_(i, j) == 0) {
        const double w = savage_dickey::sample_cauchy_omega_prior(rng_);
        omega_(i, j) = w;
        omega_(j, i) = w;
        return;
    }

    const arma::rowvec ri_cov = inv_cholesky_of_precision_.row(i);
    const arma::rowvec rj_cov = inv_cholesky_of_precision_.row(j);
    const double s_ii_cov = arma::dot(ri_cov, ri_cov);
    const double s_jj_cov = arma::dot(rj_cov, rj_cov);
    const double s_ij_cov = arma::dot(ri_cov, rj_cov);
    const double det_cov  = s_ii_cov * s_jj_cov - s_ij_cov * s_ij_cov;
    if (!(det_cov > 0.0)) return;
    const double S_11 = s_jj_cov / det_cov;
    if (!(S_11 > 0.0)) return;
    const double l_ii = std::sqrt(S_11);
    const double l_ji = (-s_ij_cov / det_cov) / l_ii;
    const double K_ij = precision_matrix_(i, j);
    const double m_ij = l_ji - K_ij / l_ii;

    const double sigma    = 2.0 * prior_sigma_;
    const double sigma2   = sigma * sigma;
    const double beta     = prior_beta_;
    const double A_diag_K = 0.5 * beta / (l_ii * l_ii);
    const double B_K      = -beta * m_ij / l_ii;

    const double w_new = savage_dickey::slice_sample_cauchy_omega_active(
        K_ij, sigma2, A_diag_K, B_K, omega_(i, j), rng_);
    omega_(i, j) = w_new;
    omega_(j, i) = w_new;
}


void GGMModel::update_edge_indicators() {
    // Hierarchical SD between-step reads Σ_corner on demand via inv(L) row
    // inner products. Within-step's SMW updates have made inv_L stale wrt
    // the L modified by cholesky_update_after_edge/diag; refresh inv_L once
    // at the top of the sweep so the corner reads are consistent.
    if (graph_prior_spec_ == GraphPriorSpec::Hierarchical) {
        arma::solve(inv_cholesky_of_precision_,
                    arma::trimatu(cholesky_of_precision_),
                    arma::eye(p_, p_), arma::solve_opts::fast);
    }

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
        if (graph_prior_spec_ == GraphPriorSpec::Hierarchical) {
            update_edge_indicator_parameter_pair_sd_lspace(i, j);
            // Cauchy slab: refresh omega_(i,j) from its L-space-consistent
            // full conditional after the (gamma, K_ij) MH step.  No-op for
            // Normal slab (omega_ij is fixed at 1.0).
            if (slab_is_cauchy_) slice_sample_omega_ij_(i, j);
        } else {
            update_edge_indicator_parameter_pair(i, j);
        }
    }
    // Hierarchical: Σ has been bypassed throughout the between sweep
    // (corner reads come from inv_L row inner products). Refresh Σ once here
    // so the within-step that follows reads a fresh covariance_matrix_. One
    // matmul per iter, vs one matmul per SD accept under the old refresh-
    // per-accept scheme.
    if (graph_prior_spec_ == GraphPriorSpec::Hierarchical) {
        covariance_matrix_ = inv_cholesky_of_precision_
                           * inv_cholesky_of_precision_.t();
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
    arma::mat U_new;
    if (!arma::chol(U_new, precision_matrix_, "upper")) {
        // Diagnostic dump on failure.
        arma::vec eigs;
        bool eig_ok = arma::eig_sym(eigs, precision_matrix_);
        arma::vec diag_K = precision_matrix_.diag();
        arma::mat off = precision_matrix_;
        off.diag().zeros();
        double max_abs_off = arma::abs(off).max();
        double drift_max  = arma::abs(arma::sum(covariance_matrix_ % precision_matrix_, 1) - 1.0).max();
        int n_edges = (static_cast<int>(arma::accu(edge_indicators_)) - static_cast<int>(p_)) / 2;
        Rcpp::Rcout << "[refresh_cholesky] non-PD K detected:"
                    << "  p="          << p_
                    << "  n_edges="    << n_edges
                    << "  min(diag K)="  << diag_K.min()
                    << "  max(diag K)="  << diag_K.max()
                    << "  max|offdiag|=" << max_abs_off
                    << "  cov*K drift=" << drift_max;
        if (eig_ok) {
            Rcpp::Rcout << "  min_eig=" << eigs.min() << "  max_eig=" << eigs.max();
        } else {
            Rcpp::Rcout << "  (eig_sym failed too)";
        }
        Rcpp::Rcout << std::endl;
        Rcpp::stop("refresh_cholesky: precision matrix is not positive definite");
    }
    cholesky_of_precision_ = U_new;
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
