#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "rng/rng_utils.h"
#include "mcmc/execution/step_result.h"
#include "math/explog_macros.h"


// =============================================================================
// Beta-type prior used for all main effects (ordinal thresholds and BC α/β).
// Matches OMRFModel::log_pseudoposterior_main_component.
// =============================================================================

static double log_beta_prior(double x, double alpha, double beta) {
    return x * alpha - std::log1p(MY_EXP(x)) * (alpha + beta);
}


// =============================================================================
// update_main_effect
// =============================================================================
// MH update for one main-effect parameter.
//   Ordinal: mux_(s, c) = threshold for category c+1  (c in [0, C_s-1])
//   Blume-Capel: mux_(s, 0) = linear α, mux_(s, 1) = quadratic β
//                (c indexes 0 or 1 for BC)
//
// The accept/reject uses log_conditional_omrf(s) + beta-type prior.
// =============================================================================

void MixedMRFModel::update_main_effect(int s, int c) {
    double& current = mux_(s, c);
    double proposal_sd = prop_sd_mux_(s, c);

    double current_val = current;
    double proposed = rnorm(rng_, current_val, proposal_sd);

    // Current log-posterior
    double ll_curr = (use_marginal_pl_ ? log_marginal_omrf(s) : log_conditional_omrf(s))
                   + log_beta_prior(current_val, main_alpha_, main_beta_);

    // Proposed log-posterior
    current = proposed;
    double ll_prop = (use_marginal_pl_ ? log_marginal_omrf(s) : log_conditional_omrf(s))
                   + log_beta_prior(proposed, main_alpha_, main_beta_);

    double ln_alpha = ll_prop - ll_curr;

    if(std::log(runif(rng_)) >= ln_alpha) {
        current = current_val;  // reject
    }
}


// =============================================================================
// update_continuous_mean
// =============================================================================
// MH update for one continuous mean parameter muy_(j).
// The accept/reject uses log_conditional_ggm() + Normal(0, 1) prior.
// Must save/restore conditional_mean_ around the proposal.
// =============================================================================

void MixedMRFModel::update_continuous_mean(int j) {
    double current_val = muy_(j);
    double proposed = rnorm(rng_, current_val, prop_sd_muy_(j));

    // Current log-posterior (Normal(0,1) prior: -x^2/2 up to constant)
    double ll_curr = log_conditional_ggm() + R::dnorm(current_val, 0.0, 1.0, true);
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ll_curr += log_marginal_omrf(s);
    }

    // Set proposed value and refresh conditional_mean_
    arma::mat cond_mean_saved = conditional_mean_;
    muy_(j) = proposed;
    recompute_conditional_mean();

    double ll_prop = log_conditional_ggm() + R::dnorm(proposed, 0.0, 1.0, true);
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ll_prop += log_marginal_omrf(s);
    }

    double ln_alpha = ll_prop - ll_curr;

    if(std::log(runif(rng_)) >= ln_alpha) {
        muy_(j) = current_val;  // reject
        conditional_mean_ = std::move(cond_mean_saved);
    }
}


// =============================================================================
// update_Kxx
// =============================================================================
// MH update for one discrete-discrete interaction Kxx_(i, j).
// Symmetric: sets both (i,j) and (j,i).
// Acceptance: log_conditional_omrf(i) + log_conditional_omrf(j) + Cauchy prior.
// =============================================================================

void MixedMRFModel::update_Kxx(int i, int j) {
    double current_val = Kxx_(i, j);
    double proposed = rnorm(rng_, current_val, prop_sd_Kxx_(i, j));

    // Current log-posterior
    double ll_curr, ll_prop;
    if(use_marginal_pl_) {
        ll_curr = log_marginal_omrf(i) + log_marginal_omrf(j)
                + R::dcauchy(current_val, 0.0, pairwise_scale_, true);

        Kxx_(i, j) = proposed;
        Kxx_(j, i) = proposed;
        recompute_Theta();

        ll_prop = log_marginal_omrf(i) + log_marginal_omrf(j)
                + R::dcauchy(proposed, 0.0, pairwise_scale_, true);
    } else {
        ll_curr = log_conditional_omrf(i) + log_conditional_omrf(j)
                + R::dcauchy(current_val, 0.0, pairwise_scale_, true);

        Kxx_(i, j) = proposed;
        Kxx_(j, i) = proposed;

        ll_prop = log_conditional_omrf(i) + log_conditional_omrf(j)
                + R::dcauchy(proposed, 0.0, pairwise_scale_, true);
    }

    double ln_alpha = ll_prop - ll_curr;

    if(std::log(runif(rng_)) >= ln_alpha) {
        Kxx_(i, j) = current_val;  // reject
        Kxx_(j, i) = current_val;
        if(use_marginal_pl_) recompute_Theta();
    }
}


// =============================================================================
// Rank-1 Kyy proposal helpers (permutation-free)
// =============================================================================
// Direct analogs of GGMModel::get_constants / constrained_diagonal,
// operating on Kyy_, Kyy_chol_, and covariance_yy_.
// =============================================================================

void MixedMRFModel::get_kyy_constants(int i, int j) {
    double logdet = cholesky_helpers::get_log_det(Kyy_chol_);

    double log_adj_ii = logdet + std::log(std::abs(covariance_yy_(i, i)));
    double log_adj_ij = logdet + std::log(std::abs(covariance_yy_(i, j)));
    double log_adj_jj = logdet + std::log(std::abs(covariance_yy_(j, j)));

    double inv_sub_jj = cholesky_helpers::compute_inv_submatrix_i(covariance_yy_, i, j, j);
    double log_abs_inv_sub_jj = log_adj_ii + std::log(std::abs(inv_sub_jj));

    double Phi_q1q  = (2 * std::signbit(covariance_yy_(i, j)) - 1) * std::exp(
        (log_adj_ij - (log_adj_jj + log_abs_inv_sub_jj) / 2)
    );
    double Phi_q1q1 = std::exp((log_adj_jj - log_abs_inv_sub_jj) / 2);

    kyy_constants_[0] = Phi_q1q;
    kyy_constants_[1] = Phi_q1q1;
    kyy_constants_[2] = Kyy_(i, j) - Phi_q1q * Phi_q1q1;
    kyy_constants_[3] = Phi_q1q1;
    kyy_constants_[4] = Kyy_(j, j) - Phi_q1q * Phi_q1q;
    kyy_constants_[5] = kyy_constants_[4] +
        kyy_constants_[2] * kyy_constants_[2] / (kyy_constants_[3] * kyy_constants_[3]);
}

double MixedMRFModel::kyy_constrained_diagonal(double x) const {
    if(x == 0.0) {
        return kyy_constants_[5];
    } else {
        double t = (x - kyy_constants_[2]) / kyy_constants_[3];
        return kyy_constants_[4] + t * t;
    }
}


// =============================================================================
// log_ggm_ratio_edge
// =============================================================================
// Log-likelihood ratio for a rank-2 off-diagonal Kyy change using the
// matrix determinant lemma for the log-det part and Woodbury for the
// quadratic-form part.  Assumes precision_yy_proposal_ is filled.
//
// B+.10 will replace the O(npq + nq²) quadratic-form computation with
// an O(nq) rank-2 shortcut.
// =============================================================================

double MixedMRFModel::log_ggm_ratio_edge(int i, int j) const {
    size_t ui = static_cast<size_t>(i);
    size_t uj = static_cast<size_t>(j);

    // --- Log-determinant ratio via matrix determinant lemma ---
    // ΔΩ has 3 nonzero entries: (i,j), (j,i), (j,j).
    // Ui = old - new off-diag, Uj = (old - new diag) / 2
    double Ui = Kyy_(ui, uj) - precision_yy_proposal_(ui, uj);
    double Uj = (Kyy_(uj, uj) - precision_yy_proposal_(uj, uj)) / 2.0;

    double cc11 = covariance_yy_(uj, uj);
    double cc12 = 1.0 - (covariance_yy_(ui, uj) * Ui +
                          covariance_yy_(uj, uj) * Uj);
    double cc22 = Ui * Ui * covariance_yy_(ui, ui) +
                  2.0 * Ui * Uj * covariance_yy_(ui, uj) +
                  Uj * Uj * covariance_yy_(uj, uj);

    double logdet_ratio = std::log(std::abs(cc11 * cc22 - cc12 * cc12));

    // --- Proposed covariance via Woodbury ---
    // ΔΩ = vf1 vf2' + vf2 vf1' where vf1 = [0,...,-1,...] (j-th),
    //   vf2 = [0,...,Ui,...,Uj,...] (i-th and j-th).
    // s1 = Σ vf1 = -Σ[:,j], s2 = Σ vf2 = Ui*Σ[:,i] + Uj*Σ[:,j]
    arma::vec s1 = -covariance_yy_.col(uj);
    arma::vec s2 = Ui * covariance_yy_.col(ui) + Uj * covariance_yy_.col(uj);

    // 2×2 core matrix T = I + [vf2,vf1]' [s1,s2]
    // T = [1 + vf2's1,  vf2's2;  vf1's1,  1 + vf1's2]
    double t11 = 1.0 + Ui * s1(ui) + Uj * s1(uj);     // 1 + vf2' s1
    double t12 = Ui * s2(ui) + Uj * s2(uj);            // vf2' s2
    double t21 = -s1(uj);                              // vf1' s1 = Σ(j,j)
    double t22 = 1.0 - s2(uj);                         // 1 + vf1' s2

    double det_T = t11 * t22 - t12 * t21;

    // T^{-1}
    double inv_t11 =  t22 / det_T;
    double inv_t12 = -t12 / det_T;
    double inv_t21 = -t21 / det_T;
    double inv_t22 =  t11 / det_T;

    // Σ' = Σ - [s1,s2] T^{-1} [s2',s1']
    //     = Σ - (inv_t11*s1 + inv_t21*s2)*s2' - (inv_t12*s1 + inv_t22*s2)*s1'
    arma::vec w1 = inv_t11 * s1 + inv_t21 * s2;  // coefficient for s2' row
    arma::vec w2 = inv_t12 * s1 + inv_t22 * s2;  // coefficient for s1' row
    arma::mat cov_prop = covariance_yy_ - w1 * s2.t() - w2 * s1.t();

    // --- Proposed conditional mean ---
    // M' = 1*μ_y' + 2 * X * Kxy * Σ'
    arma::mat cond_mean_prop = arma::repmat(muy_.t(), n_, 1) +
                               2.0 * discrete_observations_dbl_ * Kxy_ * cov_prop;

    // --- Quadratic form difference ---
    arma::mat D_curr = continuous_observations_ - conditional_mean_;
    arma::mat D_prop = continuous_observations_ - cond_mean_prop;

    double quad_curr = arma::accu((D_curr * Kyy_) % D_curr);
    double quad_prop = arma::accu((D_prop * precision_yy_proposal_) % D_prop);

    double n = static_cast<double>(n_);
    return n / 2.0 * logdet_ratio - (quad_prop - quad_curr) / 2.0;
}


// =============================================================================
// log_ggm_ratio_diag
// =============================================================================
// Log-likelihood ratio for a rank-1 diagonal Kyy change.
// Same structure as log_ggm_ratio_edge but simpler (Ui = 0).
// =============================================================================

double MixedMRFModel::log_ggm_ratio_diag(int i) const {
    size_t ui = static_cast<size_t>(i);

    // --- Log-determinant ratio (rank-1) ---
    double Uj = (Kyy_(ui, ui) - precision_yy_proposal_(ui, ui)) / 2.0;

    double cc11 = covariance_yy_(ui, ui);
    double cc12 = 1.0 - covariance_yy_(ui, ui) * Uj;
    double cc22 = Uj * Uj * covariance_yy_(ui, ui);

    double logdet_ratio = std::log(std::abs(cc11 * cc22 - cc12 * cc12));

    // --- Proposed covariance via Sherman-Morrison (rank-1 special case) ---
    // ΔΩ = -2Uj * e_i e_i', so Σ' = Σ + 2Uj * Σ[:,i] Σ[i,:]' / (1 - 2Uj * Σ(i,i))
    arma::vec s = covariance_yy_.col(ui);
    double denom = 1.0 - 2.0 * Uj * covariance_yy_(ui, ui);
    arma::mat cov_prop = covariance_yy_ + (2.0 * Uj / denom) * s * s.t();

    // --- Proposed conditional mean ---
    arma::mat cond_mean_prop = arma::repmat(muy_.t(), n_, 1) +
                               2.0 * discrete_observations_dbl_ * Kxy_ * cov_prop;

    // --- Quadratic form difference ---
    arma::mat D_curr = continuous_observations_ - conditional_mean_;
    arma::mat D_prop = continuous_observations_ - cond_mean_prop;

    double quad_curr = arma::accu((D_curr * Kyy_) % D_curr);
    double quad_prop = arma::accu((D_prop * precision_yy_proposal_) % D_prop);

    double n = static_cast<double>(n_);
    return n / 2.0 * logdet_ratio - (quad_prop - quad_curr) / 2.0;
}


// =============================================================================
// cholesky_update_after_kyy_edge
// =============================================================================
// Rank-2 Cholesky update after accepting an off-diagonal Kyy change.
// Decomposes ΔΩ = vf1*vf2' + vf2*vf1' into two rank-1 ops.
// Then recomputes inv_cholesky_yy_ and covariance_yy_.
// =============================================================================

void MixedMRFModel::cholesky_update_after_kyy_edge(
    double old_ij, double old_jj, int i, int j)
{
    kyy_v2_[0] = old_ij - precision_yy_proposal_(i, j);
    kyy_v2_[1] = (old_jj - precision_yy_proposal_(j, j)) / 2.0;

    kyy_vf1_[i] = kyy_v1_[0];   // 0
    kyy_vf1_[j] = kyy_v1_[1];   // -1
    kyy_vf2_[i] = kyy_v2_[0];
    kyy_vf2_[j] = kyy_v2_[1];

    kyy_u1_ = (kyy_vf1_ + kyy_vf2_) / std::sqrt(2.0);
    kyy_u2_ = (kyy_vf1_ - kyy_vf2_) / std::sqrt(2.0);

    cholesky_update(Kyy_chol_, kyy_u1_);
    cholesky_downdate(Kyy_chol_, kyy_u2_);

    arma::inv(inv_cholesky_yy_, arma::trimatu(Kyy_chol_));
    covariance_yy_ = inv_cholesky_yy_ * inv_cholesky_yy_.t();
    Kyy_log_det_ = cholesky_helpers::get_log_det(Kyy_chol_);

    kyy_vf1_[i] = 0.0;
    kyy_vf1_[j] = 0.0;
    kyy_vf2_[i] = 0.0;
    kyy_vf2_[j] = 0.0;
}


// =============================================================================
// cholesky_update_after_kyy_diag
// =============================================================================
// Rank-1 Cholesky update after accepting a diagonal Kyy change.
// =============================================================================

void MixedMRFModel::cholesky_update_after_kyy_diag(double old_ii, int i) {
    double delta = old_ii - precision_yy_proposal_(i, i);
    bool downdate = delta > 0.0;

    kyy_vf1_[i] = std::sqrt(std::abs(delta));

    if(downdate)
        cholesky_downdate(Kyy_chol_, kyy_vf1_);
    else
        cholesky_update(Kyy_chol_, kyy_vf1_);

    arma::inv(inv_cholesky_yy_, arma::trimatu(Kyy_chol_));
    covariance_yy_ = inv_cholesky_yy_ * inv_cholesky_yy_.t();
    Kyy_log_det_ = cholesky_helpers::get_log_det(Kyy_chol_);

    kyy_vf1_[i] = 0.0;
}


// =============================================================================
// update_Kyy_offdiag
// =============================================================================
// MH update for one off-diagonal element of the precision matrix Kyy_(i, j).
// Uses rank-1 Cholesky infrastructure (GGM-style, no permutation):
//   1. Extract constants from covariance_yy_ and Kyy_chol_
//   2. Propose on the unconstrained Cholesky scale
//   3. Map to Kyy space with constrained diagonal
//   4. Evaluate rank-2 log-likelihood ratio
//   5. On accept: rank-1 Cholesky update
//
// Prior: Cauchy(0, pairwise_scale_) on off-diag, Gamma(1, 1) on diagonal.
// =============================================================================

void MixedMRFModel::update_Kyy_offdiag(int i, int j) {
    get_kyy_constants(i, j);

    double phi_curr = kyy_constants_[0];  // Phi_q1q
    double phi_prop = rnorm(rng_, phi_curr, prop_sd_Kyy_(i, j));

    double omega_prop_ij = kyy_constants_[2] + kyy_constants_[3] * phi_prop;
    double omega_prop_jj = kyy_constrained_diagonal(omega_prop_ij);
    double diag_curr = Kyy_(j, j);

    // Fill proposal matrix (only the 3 changed entries matter)
    precision_yy_proposal_ = Kyy_;
    precision_yy_proposal_(i, j) = omega_prop_ij;
    precision_yy_proposal_(j, i) = omega_prop_ij;
    precision_yy_proposal_(j, j) = omega_prop_jj;

    double ln_alpha = log_ggm_ratio_edge(i, j);

    // Marginal mode: add OMRF ratio with proposed Theta
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ln_alpha -= log_marginal_omrf(s);

        arma::mat Theta_saved = Theta_;
        arma::mat Kyy_saved = Kyy_;
        Kyy_ = precision_yy_proposal_;
        recompute_Theta();
        for(size_t s = 0; s < p_; ++s)
            ln_alpha += log_marginal_omrf(s);
        Kyy_ = Kyy_saved;
        Theta_ = std::move(Theta_saved);
    }

    // Prior ratio: Cauchy on off-diag + Gamma(1,1) on diagonal
    ln_alpha += R::dcauchy(omega_prop_ij, 0.0, pairwise_scale_, true);
    ln_alpha -= R::dcauchy(Kyy_(i, j), 0.0, pairwise_scale_, true);
    ln_alpha += R::dgamma(omega_prop_jj, 1.0, 1.0, true);
    ln_alpha -= R::dgamma(diag_curr, 1.0, 1.0, true);

    if(std::log(runif(rng_)) < ln_alpha) {
        double old_ij = Kyy_(i, j);
        double old_jj = Kyy_(j, j);

        Kyy_(i, j) = omega_prop_ij;
        Kyy_(j, i) = omega_prop_ij;
        Kyy_(j, j) = omega_prop_jj;

        cholesky_update_after_kyy_edge(old_ij, old_jj, i, j);
        recompute_conditional_mean();
        if(use_marginal_pl_) recompute_Theta();
    }
}


// =============================================================================
// update_Kyy_diag
// =============================================================================
// MH update for one diagonal element of Kyy.
// Proposes on the log-Cholesky scale to ensure positivity.
// Uses rank-1 Cholesky update on accept.
// Prior: Gamma(1, 1) on the diagonal element + Jacobian for log-scale proposal.
// =============================================================================

void MixedMRFModel::update_Kyy_diag(int i) {
    double logdet = cholesky_helpers::get_log_det(Kyy_chol_);
    double logdet_sub_ii = logdet + std::log(covariance_yy_(i, i));

    double theta_curr = (logdet - logdet_sub_ii) / 2.0;
    double theta_prop = rnorm(rng_, theta_curr, prop_sd_Kyy_(i, i));

    precision_yy_proposal_ = Kyy_;
    precision_yy_proposal_(i, i) = Kyy_(i, i)
        - std::exp(theta_curr) * std::exp(theta_curr)
        + std::exp(theta_prop) * std::exp(theta_prop);

    double ln_alpha = log_ggm_ratio_diag(i);

    // Marginal mode: add OMRF ratio with proposed Theta
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ln_alpha -= log_marginal_omrf(s);

        arma::mat Theta_saved = Theta_;
        arma::mat Kyy_saved = Kyy_;
        Kyy_ = precision_yy_proposal_;
        recompute_Theta();
        for(size_t s = 0; s < p_; ++s)
            ln_alpha += log_marginal_omrf(s);
        Kyy_ = Kyy_saved;
        Theta_ = std::move(Theta_saved);
    }

    // Prior ratio: Gamma(1,1) on diagonal
    ln_alpha += R::dgamma(precision_yy_proposal_(i, i), 1.0, 1.0, true);
    ln_alpha -= R::dgamma(Kyy_(i, i), 1.0, 1.0, true);

    // Jacobian for log-scale proposal
    ln_alpha += theta_prop - theta_curr;

    if(std::log(runif(rng_)) < ln_alpha) {
        double old_ii = Kyy_(i, i);
        Kyy_(i, i) = precision_yy_proposal_(i, i);

        cholesky_update_after_kyy_diag(old_ii, i);
        recompute_conditional_mean();
        if(use_marginal_pl_) recompute_Theta();
    }
}


// =============================================================================
// update_Kxy
// =============================================================================
// MH update for one cross-type interaction Kxy_(i, j).
// Acceptance: log_conditional_omrf(i) + log_conditional_ggm() + Cauchy prior.
// Must save/restore conditional_mean_ around the proposal.
// =============================================================================

void MixedMRFModel::update_Kxy(int i, int j) {
    double current_val = Kxy_(i, j);
    double proposed = rnorm(rng_, current_val, prop_sd_Kxy_(i, j));

    // Current log-posterior
    double ll_curr = log_conditional_ggm()
                   + R::dcauchy(current_val, 0.0, pairwise_scale_, true);
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ll_curr += log_marginal_omrf(s);
    } else {
        ll_curr += log_conditional_omrf(i);
    }

    // Set proposed value and refresh caches
    arma::mat cond_mean_saved = conditional_mean_;
    arma::mat Theta_saved;
    if(use_marginal_pl_) Theta_saved = Theta_;
    Kxy_(i, j) = proposed;
    recompute_conditional_mean();
    if(use_marginal_pl_) recompute_Theta();

    double ll_prop = log_conditional_ggm()
                   + R::dcauchy(proposed, 0.0, pairwise_scale_, true);
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ll_prop += log_marginal_omrf(s);
    } else {
        ll_prop += log_conditional_omrf(i);
    }

    double ln_alpha = ll_prop - ll_curr;

    if(std::log(runif(rng_)) >= ln_alpha) {
        Kxy_(i, j) = current_val;  // reject
        conditional_mean_ = std::move(cond_mean_saved);
        if(use_marginal_pl_) Theta_ = std::move(Theta_saved);
    }
}
