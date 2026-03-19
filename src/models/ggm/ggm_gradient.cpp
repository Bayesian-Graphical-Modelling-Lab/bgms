#include "models/ggm/ggm_gradient.h"

#include <cmath>
#include <stdexcept>

// =====================================================================
// rebuild
// =====================================================================

void GGMGradientEngine::rebuild(
    const GraphConstraintStructure& structure,
    size_t n,
    const arma::mat& suf_stat,
    double pairwise_scale)
{
    structure_ = &structure;
    n_ = n;
    p_ = structure.p;
    suf_stat_ = &suf_stat;
    pairwise_scale_ = pairwise_scale;
}

// =====================================================================
// build_Aq
// =====================================================================
// Build constraint matrix A_q for column q.
//
// For each excluded index i (rows where edge (i,q) is absent),
// row r of A_q is: A_q[r, 0:i] = Phi[0:i, i].
// This encodes the zero-constraint on K[i,q] = sum_l Phi[l,i]*Phi[l,q].

void GGMGradientEngine::build_Aq(
    const arma::mat& Phi,
    const ColumnConstraints& col,
    size_t q,
    arma::mat& Aq)
{
    Aq.zeros(col.m_q, q);
    for (size_t r = 0; r < col.m_q; ++r) {
        size_t i = col.excluded_indices[r];
        for (size_t l = 0; l <= i; ++l) {
            Aq(r, l) = Phi(l, i);
        }
    }
}

// =====================================================================
// compute_null_basis_and_rdiag
// =====================================================================
// Compute null-space basis N_q and QR R-diagonal from A_q.
//
// Uses QR of A_q^T (which is q x m_q). The last (q - rank) columns
// of Q form the null-space basis. The diagonal of R gives the
// Jacobian contribution.

void GGMGradientEngine::compute_null_basis_and_rdiag(
    const arma::mat& Aq,
    arma::mat& Nq,
    arma::vec& R_diag)
{
    size_t m_q = Aq.n_rows;
    size_t q_minus_1 = Aq.n_cols;  // q (0-based), so this is the number of rows above diagonal

    if (m_q == 0) {
        // No constraints: null space is the full space
        Nq = arma::eye(q_minus_1, q_minus_1);
        R_diag.reset();
        return;
    }

    if (m_q >= q_minus_1) {
        // Fully constrained: no free parameters
        Nq.reset();
        // Still need R_diag for the Jacobian
        arma::mat Q, R;
        arma::qr_econ(Q, R, Aq.t());
        R_diag.set_size(R.n_rows);
        for (size_t j = 0; j < R.n_rows; ++j) {
            R_diag(j) = std::abs(R(j, j));
        }
        return;
    }

    // General case: QR of A_q^T to get null space
    arma::mat AqT = Aq.t();  // (q-1) x m_q
    arma::mat Q, R;
    arma::qr(Q, R, AqT);    // Full QR: Q is (q-1) x (q-1)

    size_t rank = m_q;  // A_q should have full row rank

    // R_diag: absolute diagonal of R (length = rank)
    R_diag.set_size(rank);
    for (size_t j = 0; j < rank; ++j) {
        R_diag(j) = std::abs(R(j, j));
    }

    // N_q: last (q-1 - rank) columns of Q
    Nq = Q.cols(rank, q_minus_1 - 1);
}

// =====================================================================
// forward_map
// =====================================================================

ForwardMapResult GGMGradientEngine::forward_map(const arma::vec& theta) const {
    ForwardMapResult result;
    result.Phi.zeros(p_, p_);
    result.psi.set_size(p_);
    result.Nq.resize(p_);
    result.R_diag.resize(p_);

    arma::mat Aq_buf;  // reusable buffer for A_q

    for (size_t q = 0; q < p_; ++q) {
        const auto& col = structure_->columns[q];
        size_t offset = structure_->theta_offsets[q];

        if (q == 0) {
            // Column 0: just the diagonal
            double psi_q = theta(offset);
            result.psi(q) = psi_q;
            result.Phi(0, 0) = std::exp(psi_q);
            continue;
        }

        // Extract f_q from theta
        size_t d_q = col.d_q;
        arma::vec f_q;
        if (d_q > 0) {
            f_q = theta.subvec(offset, offset + d_q - 1);
        }

        // psi_q is after f_q
        double psi_q = theta(offset + d_q);
        result.psi(q) = psi_q;
        result.Phi(q, q) = std::exp(psi_q);

        // Build constraint matrix A_q
        size_t m_q = col.m_q;

        if (m_q == 0 && d_q == q) {
            // No constraints: x_q = f_q directly (N_q = I)
            result.Nq[q] = arma::eye(q, q);
            result.R_diag[q].reset();
            for (size_t k = 0; k < d_q; ++k) {
                result.Phi(k, q) = f_q(k);
            }
        } else if (d_q == 0) {
            // Fully constrained: x_q = 0
            // Still compute R_diag for Jacobian
            if (m_q > 0) {
                build_Aq(result.Phi, col, q, Aq_buf);
                arma::mat Nq_dummy;
                compute_null_basis_and_rdiag(Aq_buf, Nq_dummy, result.R_diag[q]);
            }
            result.Nq[q].reset();
            // x_q stays zero (already zeroed)
        } else {
            // General case: build A_q, QR, null space
            build_Aq(result.Phi, col, q, Aq_buf);
            compute_null_basis_and_rdiag(Aq_buf, result.Nq[q], result.R_diag[q]);

            // x_q = N_q f_q
            arma::vec x_q = result.Nq[q] * f_q;
            for (size_t k = 0; k < q; ++k) {
                result.Phi(k, q) = x_q(k);
            }
        }
    }

    // Jacobian: log|det J| = p*log(2) + 2*sum(psi) + sum_{i<p}(p-i)*psi_i
    //                       - sum_q sum_j log|R_{q,jj}|
    double ldj = static_cast<double>(p_) * std::log(2.0);
    for (size_t q = 0; q < p_; ++q) {
        ldj += 2.0 * result.psi(q);
    }
    for (size_t i = 0; i + 1 < p_; ++i) {
        ldj += static_cast<double>(p_ - 1 - i) * result.psi(i);
    }
    for (size_t q = 1; q < p_; ++q) {
        const auto& rd = result.R_diag[q];
        for (size_t j = 0; j < rd.n_elem; ++j) {
            ldj -= std::log(rd(j));
        }
    }
    result.log_det_jacobian = ldj;

    // K = Phi^T Phi
    result.K = result.Phi.t() * result.Phi;

    return result;
}

// =====================================================================
// logp_and_gradient
// =====================================================================

std::pair<double, arma::vec> GGMGradientEngine::logp_and_gradient(
    const arma::vec& theta) const
{
    // --- Forward pass ---
    ForwardMapResult fm = forward_map(theta);
    const arma::mat& Phi = fm.Phi;
    const arma::mat& K = fm.K;
    const arma::mat& S = *suf_stat_;
    double n = static_cast<double>(n_);
    double scale2 = pairwise_scale_ * pairwise_scale_;

    // --- Log-posterior value ---
    double log_det_K = 2.0 * arma::accu(fm.psi);
    double log_lik = (n / 2.0) * log_det_K - 0.5 * arma::accu(K % S);

    // Cauchy slab prior on included off-diagonal K entries
    double log_slab = 0.0;
    for (size_t q = 1; q < p_; ++q) {
        for (size_t i : structure_->columns[q].included_indices) {
            double kij = K(i, q);
            log_slab += R::dcauchy(kij, 0.0, pairwise_scale_, 1);
        }
    }

    // Gamma(1,1) prior on diagonal K entries
    double log_diag_prior = 0.0;
    for (size_t i = 0; i < p_; ++i) {
        log_diag_prior += R::dgamma(K(i, i), 1.0, 1.0, 1);
    }

    double lp = log_lik + log_slab + log_diag_prior + fm.log_det_jacobian;

    // --- Backward pass ---
    arma::vec gradient(theta.n_elem, arma::fill::zeros);

    // Step 1: K_bar from likelihood + priors
    // Likelihood: (n/2) K^{-1} - (1/2) S
    // Use Phi^{-1} to compute K^{-1} = Phi^{-1} Phi^{-T} (more robust
    // than inv_sympd(K) when the leapfrog integrator pushes theta to
    // extreme values)
    arma::mat Phi_inv = arma::inv(arma::trimatu(Phi));
    arma::mat K_inv = Phi_inv * Phi_inv.t();
    arma::mat K_bar = (n / 2.0) * K_inv - 0.5 * S;

    // Cauchy prior on off-diagonal included edges (upper triangle only)
    // d/dk log(dcauchy(k; 0, s)) = -2k / (s^2 + k^2)
    for (size_t q = 1; q < p_; ++q) {
        for (size_t i : structure_->columns[q].included_indices) {
            double kij = K(i, q);
            K_bar(i, q) += -2.0 * kij / (scale2 + kij * kij);
        }
    }

    // Gamma(1,1) prior on diagonals: d/dk [log dgamma(k;1,1)] = -1
    for (size_t i = 0; i < p_; ++i) {
        K_bar(i, i) -= 1.0;
    }

    // Step 2: Phi_bar from K_bar
    // K = Phi^T Phi, adjoint: Phi_bar = Phi * (K_bar + K_bar^T)
    arma::mat K_bar_sym = K_bar + K_bar.t();
    arma::mat Phi_bar = Phi * K_bar_sym;

    // Step 3: Process columns in reverse order -> theta gradient
    for (size_t q = p_; q-- > 0; ) {
        const auto& col = structure_->columns[q];
        size_t offset = structure_->theta_offsets[q];
        size_t d_q = col.d_q;

        // --- psi_q gradient ---
        // Chain rule through phi_qq = exp(psi_q)
        double psi_bar = Phi_bar(q, q) * Phi(q, q);

        // Jacobian contribution: 2 + (p - 1 - q) for q < p-1, or 2 for q = p-1
        psi_bar += 2.0;
        if (q + 1 < p_) {
            psi_bar += static_cast<double>(p_ - 1 - q);
        }

        gradient(offset + d_q) = psi_bar;

        if (q == 0) continue;

        // --- f_q gradient ---
        arma::vec x_bar = Phi_bar.col(q).head(q);

        if (d_q > 0) {
            const arma::mat& Nq = fm.Nq[q];
            arma::vec f_bar = Nq.t() * x_bar;
            for (size_t k = 0; k < d_q; ++k) {
                gradient(offset + k) = f_bar(k);
            }
        }

        // --- Cross-column adjoint (Phi_{<q} via N_q and R_diag) ---
        size_t m_q = col.m_q;
        if (m_q == 0) continue;

        // Extract f_q for FD perturbations
        arma::vec f_q;
        if (d_q > 0) {
            f_q = theta.subvec(offset, offset + d_q - 1);
        }

        // Finite differences: perturb Phi[l, i] for each excluded
        // neighbor i of column q, and each l <= i.
        double eps = 1e-7;
        arma::mat Phi_pert = Phi;
        arma::mat Aq_plus, Aq_minus;
        arma::mat Nq_plus, Nq_minus;
        arma::vec Rd_plus, Rd_minus;

        for (size_t r = 0; r < m_q; ++r) {
            size_t i = col.excluded_indices[r];
            for (size_t l = 0; l <= i; ++l) {
                double orig = Phi(l, i);

                // Plus perturbation
                Phi_pert(l, i) = orig + eps;
                build_Aq(Phi_pert, col, q, Aq_plus);
                compute_null_basis_and_rdiag(Aq_plus, Nq_plus, Rd_plus);

                // Minus perturbation
                Phi_pert(l, i) = orig - eps;
                build_Aq(Phi_pert, col, q, Aq_minus);
                compute_null_basis_and_rdiag(Aq_minus, Nq_minus, Rd_minus);

                // Restore
                Phi_pert(l, i) = orig;

                // x_q contribution (when d_q > 0)
                if (d_q > 0) {
                    arma::vec x_q_plus = Nq_plus * f_q;
                    arma::vec x_q_minus = Nq_minus * f_q;
                    arma::vec dx_q = (x_q_plus - x_q_minus) / (2.0 * eps);
                    Phi_bar(l, i) += arma::dot(x_bar, dx_q);
                }

                // Jacobian R-diagonal contribution
                double jac_plus = 0.0, jac_minus = 0.0;
                for (size_t j = 0; j < Rd_plus.n_elem; ++j) {
                    jac_plus -= std::log(Rd_plus(j));
                }
                for (size_t j = 0; j < Rd_minus.n_elem; ++j) {
                    jac_minus -= std::log(Rd_minus(j));
                }
                double d_jac = (jac_plus - jac_minus) / (2.0 * eps);
                Phi_bar(l, i) += d_jac;
            }
        }
    }

    return {lp, gradient};
}
