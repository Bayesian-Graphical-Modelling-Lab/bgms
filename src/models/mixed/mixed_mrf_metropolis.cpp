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
//   Ordinal: main_effects_discrete_(s, c) = threshold for category c+1  (c in [0, C_s-1])
//   Blume-Capel: main_effects_discrete_(s, 0) = linear α, main_effects_discrete_(s, 1) = quadratic β
//                (c indexes 0 or 1 for BC)
//
// The accept/reject uses log_conditional_omrf(s) + beta-type prior.
// =============================================================================

void MixedMRFModel::update_main_effect(int s, int c, int iteration) {
    double& current = main_effects_discrete_(s, c);
    double proposal_sd = proposal_sd_main_discrete_(s, c);

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

    if(MY_LOG(runif(rng_)) >= ln_alpha) {
        current = current_val;  // reject
    }

    if(iteration >= 1 && iteration < total_warmup_) {
        double rm_weight = std::pow(iteration, -0.75);
        proposal_sd_main_discrete_(s, c) = update_proposal_sd_with_robbins_monro(
            proposal_sd_main_discrete_(s, c), ln_alpha, rm_weight, 0.44);
    }
}


// =============================================================================
// update_continuous_mean
// =============================================================================
// MH update for one continuous mean parameter main_effects_continuous_(j).
// The accept/reject uses log_conditional_ggm() + Normal(0, 1) prior.
// Must save/restore conditional_mean_ around the proposal.
// =============================================================================

void MixedMRFModel::update_continuous_mean(int j, int iteration) {
    double current_val = main_effects_continuous_(j);
    double proposed = rnorm(rng_, current_val, proposal_sd_main_continuous_(j));

    // Current log-posterior (Normal(0,1) prior: -x^2/2 up to constant)
    double ll_curr = log_conditional_ggm() + R::dnorm(current_val, 0.0, 1.0, true);
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ll_curr += log_marginal_omrf(s);
    }

    // Set proposed value and refresh conditional_mean_
    arma::mat cond_mean_saved = conditional_mean_;
    main_effects_continuous_(j) = proposed;
    recompute_conditional_mean();

    double ll_prop = log_conditional_ggm() + R::dnorm(proposed, 0.0, 1.0, true);
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ll_prop += log_marginal_omrf(s);
    }

    double ln_alpha = ll_prop - ll_curr;

    if(MY_LOG(runif(rng_)) >= ln_alpha) {
        main_effects_continuous_(j) = current_val;  // reject
        conditional_mean_ = std::move(cond_mean_saved);
    }

    if(iteration >= 1 && iteration < total_warmup_) {
        double rm_weight = std::pow(iteration, -0.75);
        proposal_sd_main_continuous_(j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_main_continuous_(j), ln_alpha, rm_weight, 0.44);
    }
}


// =============================================================================
// update_pairwise_discrete
// =============================================================================
// MH update for one discrete-discrete interaction pairwise_effects_discrete_(i, j).
// Symmetric: sets both (i,j) and (j,i).
// Acceptance: log_conditional_omrf(i) + log_conditional_omrf(j) + Cauchy prior.
// =============================================================================

void MixedMRFModel::update_pairwise_discrete(int i, int j, int iteration) {
    double current_val = pairwise_effects_discrete_(i, j);
    double proposed = rnorm(rng_, current_val, proposal_sd_pairwise_discrete_(i, j));

    // Current log-posterior
    // Kxx prior: Cauchy(0, pairwise_scale/2) on K-scale
    const double kxx_scale = 0.5 * pairwise_scale_;
    double ll_curr, ll_prop;
    if(use_marginal_pl_) {
        ll_curr = log_marginal_omrf(i) + log_marginal_omrf(j)
                + R::dcauchy(current_val, 0.0, kxx_scale, true);

        pairwise_effects_discrete_(i, j) = proposed;
        pairwise_effects_discrete_(j, i) = proposed;
        recompute_Theta();

        ll_prop = log_marginal_omrf(i) + log_marginal_omrf(j)
                + R::dcauchy(proposed, 0.0, kxx_scale, true);
    } else {
        ll_curr = log_conditional_omrf(i) + log_conditional_omrf(j)
                + R::dcauchy(current_val, 0.0, kxx_scale, true);

        pairwise_effects_discrete_(i, j) = proposed;
        pairwise_effects_discrete_(j, i) = proposed;

        ll_prop = log_conditional_omrf(i) + log_conditional_omrf(j)
                + R::dcauchy(proposed, 0.0, kxx_scale, true);
    }

    double ln_alpha = ll_prop - ll_curr;

    if(MY_LOG(runif(rng_)) >= ln_alpha) {
        pairwise_effects_discrete_(i, j) = current_val;  // reject
        pairwise_effects_discrete_(j, i) = current_val;
        if(use_marginal_pl_) recompute_Theta();
    }

    if(iteration >= 1 && iteration < total_warmup_) {
        double rm_weight = std::pow(iteration, -0.75);
        proposal_sd_pairwise_discrete_(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_discrete_(i, j), ln_alpha, rm_weight, 0.44);
    }
}


// =============================================================================
// Rank-1 precision proposal helpers (permutation-free)
// =============================================================================
// Direct analogs of GGMModel::get_constants / constrained_diagonal,
// operating on Theta = -2 Kyy (positive-definite precision),
// cholesky_of_precision_, and covariance_continuous_.
//
// All constants and proposals live in Theta space.  Conversion to/from
// K-scale (pairwise_effects_continuous_) happens at the outer call sites.
// =============================================================================

void MixedMRFModel::get_precision_constants(int i, int j) {
    double logdet = cholesky_helpers::get_log_det(cholesky_of_precision_);

    double log_adj_ii = logdet + MY_LOG(std::abs(covariance_continuous_(i, i)));
    double log_adj_ij = logdet + MY_LOG(std::abs(covariance_continuous_(i, j)));
    double log_adj_jj = logdet + MY_LOG(std::abs(covariance_continuous_(j, j)));

    double inv_sub_jj = cholesky_helpers::compute_inv_submatrix_i(covariance_continuous_, i, j, j);
    double log_abs_inv_sub_jj = log_adj_ii + MY_LOG(std::abs(inv_sub_jj));

    double Phi_q1q  = (2 * std::signbit(covariance_continuous_(i, j)) - 1) * MY_EXP(
        (log_adj_ij - (log_adj_jj + log_abs_inv_sub_jj) / 2)
    );
    double Phi_q1q1 = MY_EXP((log_adj_jj - log_abs_inv_sub_jj) / 2);

    // Read Theta = -2 Kyy for constants extraction
    double theta_ij = -2.0 * pairwise_effects_continuous_(i, j);
    double theta_jj = -2.0 * pairwise_effects_continuous_(j, j);

    kyy_constants_[0] = Phi_q1q;
    kyy_constants_[1] = Phi_q1q1;
    kyy_constants_[2] = theta_ij - Phi_q1q * Phi_q1q1;
    kyy_constants_[3] = Phi_q1q1;
    kyy_constants_[4] = theta_jj - Phi_q1q * Phi_q1q;
    kyy_constants_[5] = kyy_constants_[4] +
        kyy_constants_[2] * kyy_constants_[2] / (kyy_constants_[3] * kyy_constants_[3]);
}

double MixedMRFModel::precision_constrained_diagonal(double x) const {
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
// Log-likelihood ratio for a rank-2 off-diagonal precision change using the
// matrix determinant lemma for the log-det part and Woodbury for the
// quadratic-form part.  Assumes precision_proposal_ is filled.
//
// TODO: replace the O(npq + nq²) quadratic-form computation with
// an O(nq) rank-2 shortcut.
// =============================================================================

double MixedMRFModel::log_ggm_ratio_edge(int i, int j) const {
    size_t ui = static_cast<size_t>(i);
    size_t uj = static_cast<size_t>(j);

    // Current Theta = -2 Kyy (positive-definite precision)
    arma::mat Theta_curr = -2.0 * pairwise_effects_continuous_;

    // --- Log-determinant ratio via matrix determinant lemma ---
    // ΔΩ has 3 nonzero entries: (i,j), (j,i), (j,j).
    // Ui = old - new off-diag, Uj = (old - new diag) / 2
    double Ui = Theta_curr(ui, uj) - precision_proposal_(ui, uj);
    double Uj = (Theta_curr(uj, uj) - precision_proposal_(uj, uj)) / 2.0;

    double cc11 = covariance_continuous_(uj, uj);
    double cc12 = 1.0 - (covariance_continuous_(ui, uj) * Ui +
                          covariance_continuous_(uj, uj) * Uj);
    double cc22 = Ui * Ui * covariance_continuous_(ui, ui) +
                  2.0 * Ui * Uj * covariance_continuous_(ui, uj) +
                  Uj * Uj * covariance_continuous_(uj, uj);

    double logdet_ratio = MY_LOG(std::abs(cc11 * cc22 - cc12 * cc12));

    // --- Proposed covariance via Woodbury ---
    // ΔΩ = vf1 vf2' + vf2 vf1' where vf1 = [0,...,-1,...] (j-th),
    //   vf2 = [0,...,Ui,...,Uj,...] (i-th and j-th).
    // s1 = Σ vf1 = -Σ[:,j], s2 = Σ vf2 = Ui*Σ[:,i] + Uj*Σ[:,j]
    arma::vec s1 = -covariance_continuous_.col(uj);
    arma::vec s2 = Ui * covariance_continuous_.col(ui) + Uj * covariance_continuous_.col(uj);

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
    arma::mat cov_prop = covariance_continuous_ - w1 * s2.t() - w2 * s1.t();

    // --- Proposed conditional mean ---
    // M' = μ_y' + 2 X K_xy Σ'
    arma::mat cond_mean_prop = arma::repmat(main_effects_continuous_.t(), n_, 1) +
                               2.0 * discrete_observations_dbl_ * pairwise_effects_cross_ * cov_prop;

    // --- Quadratic form difference ---
    arma::mat D_curr = continuous_observations_ - conditional_mean_;
    arma::mat D_prop = continuous_observations_ - cond_mean_prop;

    double quad_curr = arma::accu((D_curr * Theta_curr) % D_curr);
    double quad_prop = arma::accu((D_prop * precision_proposal_) % D_prop);

    double n = static_cast<double>(n_);
    return n / 2.0 * logdet_ratio - (quad_prop - quad_curr) / 2.0;
}


// =============================================================================
// log_ggm_ratio_diag
// =============================================================================
// Log-likelihood ratio for a rank-1 diagonal precision change.
// Same structure as log_ggm_ratio_edge but simpler (Ui = 0).
// =============================================================================

double MixedMRFModel::log_ggm_ratio_diag(int i) const {
    size_t ui = static_cast<size_t>(i);

    // Current Theta = -2 Kyy
    double theta_ii = -2.0 * pairwise_effects_continuous_(ui, ui);

    // --- Log-determinant ratio (rank-1) ---
    double Uj = (theta_ii - precision_proposal_(ui, ui)) / 2.0;

    double cc11 = covariance_continuous_(ui, ui);
    double cc12 = 1.0 - covariance_continuous_(ui, ui) * Uj;
    double cc22 = Uj * Uj * covariance_continuous_(ui, ui);

    double logdet_ratio = MY_LOG(std::abs(cc11 * cc22 - cc12 * cc12));

    // --- Proposed covariance via Sherman-Morrison (rank-1 special case) ---
    // ΔΩ = -2Uj * e_i e_i', so Σ' = Σ + 2Uj * Σ[:,i] Σ[i,:]' / (1 - 2Uj * Σ(i,i))
    arma::vec s = covariance_continuous_.col(ui);
    double denom = 1.0 - 2.0 * Uj * covariance_continuous_(ui, ui);
    arma::mat cov_prop = covariance_continuous_ + (2.0 * Uj / denom) * s * s.t();

    // --- Proposed conditional mean ---
    arma::mat cond_mean_prop = arma::repmat(main_effects_continuous_.t(), n_, 1) +
                               2.0 * discrete_observations_dbl_ * pairwise_effects_cross_ * cov_prop;

    // --- Quadratic form difference ---
    arma::mat D_curr = continuous_observations_ - conditional_mean_;
    arma::mat D_prop = continuous_observations_ - cond_mean_prop;

    // Theta_curr for quad: only diagonal changed, use full -2 Kyy
    arma::mat Theta_curr = -2.0 * pairwise_effects_continuous_;
    double quad_curr = arma::accu((D_curr * Theta_curr) % D_curr);
    double quad_prop = arma::accu((D_prop * precision_proposal_) % D_prop);

    double n = static_cast<double>(n_);
    return n / 2.0 * logdet_ratio - (quad_prop - quad_curr) / 2.0;
}


// =============================================================================
// cholesky_update_after_precision_edge
// =============================================================================
// Rank-2 Cholesky update after accepting an off-diagonal precision change.
// Decomposes ΔΩ = vf1*vf2' + vf2*vf1' into two rank-1 ops.
// Then recomputes inv_cholesky_of_precision_ and covariance_continuous_.
// =============================================================================

void MixedMRFModel::cholesky_update_after_precision_edge(
    double old_ij, double old_jj, int i, int j)
{
    kyy_v2_[0] = old_ij - precision_proposal_(i, j);
    kyy_v2_[1] = (old_jj - precision_proposal_(j, j)) / 2.0;

    kyy_vf1_[i] = kyy_v1_[0];   // 0
    kyy_vf1_[j] = kyy_v1_[1];   // -1
    kyy_vf2_[i] = kyy_v2_[0];
    kyy_vf2_[j] = kyy_v2_[1];

    kyy_u1_ = (kyy_vf1_ + kyy_vf2_) / std::sqrt(2.0);
    kyy_u2_ = (kyy_vf1_ - kyy_vf2_) / std::sqrt(2.0);

    cholesky_update(cholesky_of_precision_, kyy_u1_);
    cholesky_downdate(cholesky_of_precision_, kyy_u2_);

    arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
    covariance_continuous_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    log_det_precision_ = cholesky_helpers::get_log_det(cholesky_of_precision_);

    kyy_vf1_[i] = 0.0;
    kyy_vf1_[j] = 0.0;
    kyy_vf2_[i] = 0.0;
    kyy_vf2_[j] = 0.0;
}


// =============================================================================
// cholesky_update_after_precision_diag
// =============================================================================
// Rank-1 Cholesky update after accepting a diagonal precision change.
// =============================================================================

void MixedMRFModel::cholesky_update_after_precision_diag(double old_ii, int i) {
    double delta = old_ii - precision_proposal_(i, i);
    bool downdate = delta > 0.0;

    kyy_vf1_[i] = std::sqrt(std::abs(delta));

    if(downdate)
        cholesky_downdate(cholesky_of_precision_, kyy_vf1_);
    else
        cholesky_update(cholesky_of_precision_, kyy_vf1_);

    arma::inv(inv_cholesky_of_precision_, arma::trimatu(cholesky_of_precision_));
    covariance_continuous_ = inv_cholesky_of_precision_ * inv_cholesky_of_precision_.t();
    log_det_precision_ = cholesky_helpers::get_log_det(cholesky_of_precision_);

    kyy_vf1_[i] = 0.0;
}


// =============================================================================
// update_pairwise_effects_continuous_offdiag
// =============================================================================
// MH update for one off-diagonal element of the precision matrix pairwise_effects_continuous_(i, j).
// Uses rank-1 Cholesky infrastructure (GGM-style, no permutation):
//   1. Extract constants from covariance_continuous_ and cholesky_of_precision_
//   2. Propose on the unconstrained Cholesky scale
//   3. Map to precision space with constrained diagonal
//   4. Evaluate rank-2 log-likelihood ratio
//   5. On accept: rank-1 Cholesky update
//
// Constants and proposals live in Theta space (Theta = -2 Kyy).
// Prior: Cauchy(0, pairwise_scale_) on Kyy off-diag,
//        Gamma(1, 1) on -Kyy diagonal.
// Storage: pairwise_effects_continuous_ = Kyy = -1/2 Theta.
// =============================================================================

void MixedMRFModel::update_pairwise_effects_continuous_offdiag(int i, int j, int iteration) {
    get_precision_constants(i, j);

    double phi_curr = kyy_constants_[0];  // Phi_q1q
    double phi_prop = rnorm(rng_, phi_curr, proposal_sd_pairwise_continuous_(i, j));

    // Propose in Theta space
    double theta_prop_ij = kyy_constants_[2] + kyy_constants_[3] * phi_prop;
    double theta_prop_jj = precision_constrained_diagonal(theta_prop_ij);

    // Current Theta values
    double theta_curr_ij = -2.0 * pairwise_effects_continuous_(i, j);
    double theta_curr_jj = -2.0 * pairwise_effects_continuous_(j, j);

    // Fill proposal matrix in Theta space
    precision_proposal_ = -2.0 * pairwise_effects_continuous_;
    precision_proposal_(i, j) = theta_prop_ij;
    precision_proposal_(j, i) = theta_prop_ij;
    precision_proposal_(j, j) = theta_prop_jj;

    double ln_alpha = log_ggm_ratio_edge(i, j);

    // Marginal mode: add OMRF ratio with proposed Kyy
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ln_alpha -= log_marginal_omrf(s);

        arma::mat Theta_saved = Theta_;
        arma::mat pairwise_effects_continuous_saved = pairwise_effects_continuous_;
        // Temporarily set Kyy to proposed value
        pairwise_effects_continuous_(i, j) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, i) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, j) = -0.5 * theta_prop_jj;
        recompute_Theta();
        for(size_t s = 0; s < p_; ++s)
            ln_alpha += log_marginal_omrf(s);
        pairwise_effects_continuous_ = pairwise_effects_continuous_saved;
        Theta_ = std::move(Theta_saved);
    }

    // Prior ratio on K-scale:
    //   Cauchy(0, scale) on Kyy[i,j] = -1/2 Theta[i,j]
    //   Gamma(1, 1) on -Kyy[j,j] = 1/2 Theta[j,j]
    double kyy_prop_ij = -0.5 * theta_prop_ij;
    double kyy_curr_ij = pairwise_effects_continuous_(i, j);
    double neg_kyy_prop_jj = 0.5 * theta_prop_jj;
    double neg_kyy_curr_jj = -pairwise_effects_continuous_(j, j);

    ln_alpha += R::dcauchy(kyy_prop_ij, 0.0, pairwise_scale_, true);
    ln_alpha -= R::dcauchy(kyy_curr_ij, 0.0, pairwise_scale_, true);
    ln_alpha += R::dgamma(neg_kyy_prop_jj, 1.0, 1.0, true);
    ln_alpha -= R::dgamma(neg_kyy_curr_jj, 1.0, 1.0, true);

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        // Pass old Theta values to Cholesky update
        double old_theta_ij = theta_curr_ij;
        double old_theta_jj = theta_curr_jj;

        // Store Kyy = -1/2 Theta
        pairwise_effects_continuous_(i, j) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, i) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, j) = -0.5 * theta_prop_jj;

        cholesky_update_after_precision_edge(old_theta_ij, old_theta_jj, i, j);
        recompute_conditional_mean();
        if(use_marginal_pl_) recompute_Theta();
    }

    if(iteration >= 1 && iteration < total_warmup_) {
        double rm_weight = std::pow(iteration, -0.75);
        proposal_sd_pairwise_continuous_(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_continuous_(i, j), ln_alpha, rm_weight, 0.44);
    }
}


// =============================================================================
// update_pairwise_effects_continuous_diag
// =============================================================================
// MH update for one diagonal element of the precision matrix.
// Proposes on the log-Cholesky scale to ensure positivity of Theta = -2 Kyy.
// Uses rank-1 Cholesky update on accept.
// Prior: Gamma(1, 1) on -Kyy[i,i] = 1/2 Theta[i,i] + Jacobian for log-scale proposal.
// =============================================================================

void MixedMRFModel::update_pairwise_effects_continuous_diag(int i, int iteration) {
    double logdet = cholesky_helpers::get_log_det(cholesky_of_precision_);
    double logdet_sub_ii = logdet + MY_LOG(covariance_continuous_(i, i));

    double theta_curr = (logdet - logdet_sub_ii) / 2.0;
    double theta_prop = rnorm(rng_, theta_curr, proposal_sd_pairwise_continuous_(i, i));

    // Current Theta diagonal
    double theta_ii_curr = -2.0 * pairwise_effects_continuous_(i, i);
    double theta_ii_prop = theta_ii_curr
        - MY_EXP(theta_curr) * MY_EXP(theta_curr)
        + MY_EXP(theta_prop) * MY_EXP(theta_prop);

    // Fill proposal in Theta space
    precision_proposal_ = -2.0 * pairwise_effects_continuous_;
    precision_proposal_(i, i) = theta_ii_prop;

    double ln_alpha = log_ggm_ratio_diag(i);

    // Marginal mode: add OMRF ratio with proposed Kyy
    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ln_alpha -= log_marginal_omrf(s);

        arma::mat Theta_saved = Theta_;
        double kyy_saved = pairwise_effects_continuous_(i, i);
        pairwise_effects_continuous_(i, i) = -0.5 * theta_ii_prop;
        recompute_Theta();
        for(size_t s = 0; s < p_; ++s)
            ln_alpha += log_marginal_omrf(s);
        pairwise_effects_continuous_(i, i) = kyy_saved;
        Theta_ = std::move(Theta_saved);
    }

    // Prior ratio: Gamma(1,1) on -Kyy[i,i] = 1/2 Theta[i,i]
    double neg_kyy_prop = 0.5 * theta_ii_prop;
    double neg_kyy_curr = -pairwise_effects_continuous_(i, i);
    ln_alpha += R::dgamma(neg_kyy_prop, 1.0, 1.0, true);
    ln_alpha -= R::dgamma(neg_kyy_curr, 1.0, 1.0, true);

    // Jacobian for log-scale proposal
    ln_alpha += theta_prop - theta_curr;

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        // Pass old Theta value to Cholesky update
        double old_theta_ii = theta_ii_curr;

        // Store Kyy = -1/2 Theta
        pairwise_effects_continuous_(i, i) = -0.5 * theta_ii_prop;

        cholesky_update_after_precision_diag(old_theta_ii, i);
        recompute_conditional_mean();
        if(use_marginal_pl_) recompute_Theta();
    }

    if(iteration >= 1 && iteration < total_warmup_) {
        double rm_weight = std::pow(iteration, -0.75);
        proposal_sd_pairwise_continuous_(i, i) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_continuous_(i, i), ln_alpha, rm_weight, 0.44);
    }
}


// =============================================================================
// update_pairwise_cross
// =============================================================================
// MH update for one cross-type interaction pairwise_effects_cross_(i, j).
// Acceptance: log_conditional_omrf(i) + log_conditional_ggm() + Cauchy prior.
// Must save/restore conditional_mean_ around the proposal.
// =============================================================================

void MixedMRFModel::update_pairwise_cross(int i, int j, int iteration) {
    double current_val = pairwise_effects_cross_(i, j);
    double proposed = rnorm(rng_, current_val, proposal_sd_pairwise_cross_(i, j));

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
    pairwise_effects_cross_(i, j) = proposed;
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

    if(MY_LOG(runif(rng_)) >= ln_alpha) {
        pairwise_effects_cross_(i, j) = current_val;  // reject
        conditional_mean_ = std::move(cond_mean_saved);
        if(use_marginal_pl_) Theta_ = std::move(Theta_saved);
    }

    if(iteration >= 1 && iteration < total_warmup_) {
        double rm_weight = std::pow(iteration, -0.75);
        proposal_sd_pairwise_cross_(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd_pairwise_cross_(i, j), ln_alpha, rm_weight, 0.44);
    }
}


// =============================================================================
// update_edge_indicator_discrete
// =============================================================================
// Metropolis-Hastings add-delete move for a discrete-discrete edge (i, j).
//   Add (G=0→1): propose k ~ N(0, σ), accept with slab + Hastings.
//   Delete (G=1→0): set k = 0, accept with reverse terms.
// =============================================================================

void MixedMRFModel::update_edge_indicator_discrete(int i, int j) {
    double k_curr = pairwise_effects_discrete_(i, j);
    double prop_sd = proposal_sd_pairwise_discrete_(i, j);

    int g_curr = gxx(i, j);
    int g_prop = 1 - g_curr;

    double k_prop;
    if(g_prop == 1) {
        k_prop = rnorm(rng_, k_curr, prop_sd);  // k_curr = 0 on a true add
    } else {
        k_prop = 0.0;
    }

    // --- Likelihood ratio ---
    double ll_curr, ll_prop;
    if(use_marginal_pl_) {
        ll_curr = log_marginal_omrf(i) + log_marginal_omrf(j);

        pairwise_effects_discrete_(i, j) = k_prop;
        pairwise_effects_discrete_(j, i) = k_prop;
        recompute_Theta();

        ll_prop = log_marginal_omrf(i) + log_marginal_omrf(j);

        // Restore
        pairwise_effects_discrete_(i, j) = k_curr;
        pairwise_effects_discrete_(j, i) = k_curr;
        recompute_Theta();
    } else {
        ll_curr = log_conditional_omrf(i) + log_conditional_omrf(j);

        pairwise_effects_discrete_(i, j) = k_prop;
        pairwise_effects_discrete_(j, i) = k_prop;

        ll_prop = log_conditional_omrf(i) + log_conditional_omrf(j);

        // Restore
        pairwise_effects_discrete_(i, j) = k_curr;
        pairwise_effects_discrete_(j, i) = k_curr;
    }

    double ln_alpha = ll_prop - ll_curr;

    // Kxx slab prior: Cauchy(0, pairwise_scale/2) on K-scale
    const double kxx_slab_scale = 0.5 * pairwise_scale_;
    if(g_prop == 1) {
        // Add: slab prior, subtract proposal density, inclusion prior
        ln_alpha += R::dcauchy(k_prop, 0.0, kxx_slab_scale, true);
        ln_alpha -= R::dnorm(k_prop, k_curr, prop_sd, true);
        ln_alpha += MY_LOG(inclusion_probability_(i, j))
                  - MY_LOG(1.0 - inclusion_probability_(i, j));
    } else {
        // Delete: subtract slab prior, add reverse proposal density, inclusion prior
        ln_alpha -= R::dcauchy(k_curr, 0.0, kxx_slab_scale, true);
        ln_alpha += R::dnorm(k_curr, k_prop, prop_sd, true);
        ln_alpha -= MY_LOG(inclusion_probability_(i, j))
                  - MY_LOG(1.0 - inclusion_probability_(i, j));
    }

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        pairwise_effects_discrete_(i, j) = k_prop;
        pairwise_effects_discrete_(j, i) = k_prop;
        set_gxx(i, j, g_prop);
        if(use_marginal_pl_) recompute_Theta();
    }
}


// =============================================================================
// update_edge_indicator_continuous
// =============================================================================
// Metropolis-Hastings add-delete move for a continuous-continuous edge (i, j).
// Uses Cholesky reparameterization (permute-free constants extraction).
// All proposals and constants live in Theta space (Theta = -2 Kyy).
//   Add (G=0→1): propose ε ~ N(0, σ), theta_ij = C[3]*ε, constrain diagonal.
//   Delete (G=1→0): set theta_ij = 0, constrain diagonal.
// =============================================================================

void MixedMRFModel::update_edge_indicator_continuous(int i, int j) {
    get_precision_constants(i, j);

    int g_curr = gyy(i, j);
    int g_prop = 1 - g_curr;

    double theta_prop_ij, theta_prop_jj;

    if(g_prop == 1) {
        // Add: propose from N(0, σ) on reparameterized scale
        double epsilon = rnorm(rng_, 0.0, proposal_sd_pairwise_continuous_(i, j));
        theta_prop_ij = kyy_constants_[3] * epsilon;
        theta_prop_jj = precision_constrained_diagonal(theta_prop_ij);
    } else {
        // Delete: set off-diagonal to 0 in Theta space
        theta_prop_ij = 0.0;
        theta_prop_jj = precision_constrained_diagonal(0.0);
    }

    // Fill proposal in Theta space
    precision_proposal_ = -2.0 * pairwise_effects_continuous_;
    precision_proposal_(i, j) = theta_prop_ij;
    precision_proposal_(j, i) = theta_prop_ij;
    precision_proposal_(j, j) = theta_prop_jj;

    // --- Likelihood ratio ---
    double ln_alpha = log_ggm_ratio_edge(i, j);

    if(use_marginal_pl_) {
        for(size_t s = 0; s < p_; ++s)
            ln_alpha -= log_marginal_omrf(s);

        arma::mat Theta_saved = Theta_;
        arma::mat pairwise_effects_continuous_saved = pairwise_effects_continuous_;
        // Temporarily set Kyy to proposed value
        pairwise_effects_continuous_(i, j) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, i) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, j) = -0.5 * theta_prop_jj;
        recompute_Theta();
        for(size_t s = 0; s < p_; ++s)
            ln_alpha += log_marginal_omrf(s);
        pairwise_effects_continuous_ = pairwise_effects_continuous_saved;
        Theta_ = std::move(Theta_saved);
    }

    // --- Spike-and-slab terms on K-scale ---
    // Kyy[i,j] = -1/2 Theta[i,j]
    double kyy_prop_ij = -0.5 * theta_prop_ij;
    double kyy_curr_ij = pairwise_effects_continuous_(i, j);
    // The proposal density is on the Theta reparameterized scale:
    // theta_prop_ij = C[3] * epsilon, so epsilon = theta_prop_ij / C[3]
    if(g_prop == 1) {
        // Add: slab prior on proposed Kyy off-diag
        ln_alpha += R::dcauchy(kyy_prop_ij, 0.0, pairwise_scale_, true);
        // Subtract proposal density (in Theta space, Jacobian -1/2 cancels symmetrically)
        ln_alpha -= R::dnorm(theta_prop_ij / kyy_constants_[3], 0.0,
                             proposal_sd_pairwise_continuous_(i, j), true)
                  - MY_LOG(kyy_constants_[3]);
        // Inclusion prior: log(π / (1-π))
        ln_alpha += MY_LOG(inclusion_probability_(p_ + i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(p_ + i, p_ + j));
    } else {
        // Delete: subtract slab prior on current Kyy off-diag
        ln_alpha -= R::dcauchy(kyy_curr_ij, 0.0, pairwise_scale_, true);
        // Add reverse proposal density
        double theta_curr_ij = -2.0 * kyy_curr_ij;
        ln_alpha += R::dnorm(theta_curr_ij / kyy_constants_[3], 0.0,
                             proposal_sd_pairwise_continuous_(i, j), true)
                  - MY_LOG(kyy_constants_[3]);
        // Inclusion prior: log((1-π) / π)
        ln_alpha -= MY_LOG(inclusion_probability_(p_ + i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(p_ + i, p_ + j));
    }

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        // Pass old Theta values to Cholesky update
        double old_theta_ij = -2.0 * pairwise_effects_continuous_(i, j);
        double old_theta_jj = -2.0 * pairwise_effects_continuous_(j, j);

        // Store Kyy = -1/2 Theta
        pairwise_effects_continuous_(i, j) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, i) = -0.5 * theta_prop_ij;
        pairwise_effects_continuous_(j, j) = -0.5 * theta_prop_jj;

        set_gyy(i, j, g_prop);
        cholesky_update_after_precision_edge(old_theta_ij, old_theta_jj, i, j);
        recompute_conditional_mean();
        if(use_marginal_pl_) recompute_Theta();
    }
}


// =============================================================================
// update_edge_indicator_cross
// =============================================================================
// Metropolis-Hastings add-delete move for a cross-type edge (i, j).
//   Add (G=0→1): propose k ~ N(0, σ).
//   Delete (G=1→0): set k = 0.
// =============================================================================

void MixedMRFModel::update_edge_indicator_cross(int i, int j) {
    double k_curr = pairwise_effects_cross_(i, j);
    double prop_sd = proposal_sd_pairwise_cross_(i, j);

    int g_curr = gxy(i, j);
    int g_prop = 1 - g_curr;

    double k_prop;
    if(g_prop == 1) {
        k_prop = rnorm(rng_, k_curr, prop_sd);  // k_curr = 0 on a true add
    } else {
        k_prop = 0.0;
    }

    // --- Likelihood ratio ---
    double ll_curr, ll_prop;
    if(use_marginal_pl_) {
        ll_curr = log_conditional_ggm();
        for(size_t s = 0; s < p_; ++s)
            ll_curr += log_marginal_omrf(s);

        arma::mat cond_mean_saved = conditional_mean_;
        arma::mat Theta_saved = Theta_;
        pairwise_effects_cross_(i, j) = k_prop;
        recompute_conditional_mean();
        recompute_Theta();

        ll_prop = log_conditional_ggm();
        for(size_t s = 0; s < p_; ++s)
            ll_prop += log_marginal_omrf(s);

        // Restore
        pairwise_effects_cross_(i, j) = k_curr;
        conditional_mean_ = std::move(cond_mean_saved);
        Theta_ = std::move(Theta_saved);
    } else {
        ll_curr = log_conditional_omrf(i) + log_conditional_ggm();

        arma::mat cond_mean_saved = conditional_mean_;
        pairwise_effects_cross_(i, j) = k_prop;
        recompute_conditional_mean();

        ll_prop = log_conditional_omrf(i) + log_conditional_ggm();

        // Restore
        pairwise_effects_cross_(i, j) = k_curr;
        conditional_mean_ = std::move(cond_mean_saved);
    }

    double ln_alpha = ll_prop - ll_curr;

    if(g_prop == 1) {
        // Add
        ln_alpha += R::dcauchy(k_prop, 0.0, pairwise_scale_, true);
        ln_alpha -= R::dnorm(k_prop, k_curr, prop_sd, true);
        ln_alpha += MY_LOG(inclusion_probability_(i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(i, p_ + j));
    } else {
        // Delete
        ln_alpha -= R::dcauchy(k_curr, 0.0, pairwise_scale_, true);
        ln_alpha += R::dnorm(k_curr, k_prop, prop_sd, true);
        ln_alpha -= MY_LOG(inclusion_probability_(i, p_ + j))
                  - MY_LOG(1.0 - inclusion_probability_(i, p_ + j));
    }

    if(MY_LOG(runif(rng_)) < ln_alpha) {
        pairwise_effects_cross_(i, j) = k_prop;
        set_gxy(i, j, g_prop);
        recompute_conditional_mean();
        if(use_marginal_pl_) recompute_Theta();
    }
}
