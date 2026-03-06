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
    double ll_curr = log_conditional_omrf(s) + log_beta_prior(current_val, main_alpha_, main_beta_);

    // Proposed log-posterior
    current = proposed;
    double ll_prop = log_conditional_omrf(s) + log_beta_prior(proposed, main_alpha_, main_beta_);

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

    // Set proposed value and refresh conditional_mean_
    arma::mat cond_mean_saved = conditional_mean_;
    muy_(j) = proposed;
    recompute_conditional_mean();

    double ll_prop = log_conditional_ggm() + R::dnorm(proposed, 0.0, 1.0, true);

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
    double ll_curr = log_conditional_omrf(i) + log_conditional_omrf(j)
                   + R::dcauchy(current_val, 0.0, pairwise_scale_, true);

    // Set proposed value (symmetric)
    Kxx_(i, j) = proposed;
    Kxx_(j, i) = proposed;

    double ll_prop = log_conditional_omrf(i) + log_conditional_omrf(j)
                   + R::dcauchy(proposed, 0.0, pairwise_scale_, true);

    double ln_alpha = ll_prop - ll_curr;

    if(std::log(runif(rng_)) >= ln_alpha) {
        Kxx_(i, j) = current_val;  // reject
        Kxx_(j, i) = current_val;
    }
}


// =============================================================================
// Kyy helpers: permutation-based Cholesky updates
// =============================================================================
// Ported from mixedGM/R/continuous_variable_helper.R.
// These operate on the (q-1, q) block after permuting target indices to
// the last two positions.

// Build a permutation vector that moves indices (i, j) to positions (q-1, q).
// Same involution property as the R version.
static arma::uvec make_perm_offdiag(int i, int j, int q) {
    int A = q - 2;  // 0-based target for first index
    int B = q - 1;  // 0-based target for second index
    arma::uvec perm = arma::regspace<arma::uvec>(0, q - 1);

    if((i == A && j == B) || (i == B && j == A)) {
        // already in position
    } else if(i == A || j == A) {
        int other = (i == A) ? j : i;
        if(other != B) std::swap(perm(other), perm(B));
    } else if(i == B || j == B) {
        int other = (i == B) ? j : i;
        if(other != A) std::swap(perm(other), perm(A));
    } else {
        std::swap(perm(i), perm(A));
        std::swap(perm(j), perm(B));
    }
    return perm;
}

// Build permutation vector that moves index i to position q-1 (last).
static arma::uvec make_perm_diag(int i, int q) {
    arma::uvec perm = arma::regspace<arma::uvec>(0, q - 1);
    if(i != q - 1) std::swap(perm(i), perm(q - 1));
    return perm;
}

// Apply symmetric permutation: M[perm, perm]
static arma::mat permute_matrix(const arma::mat& M, const arma::uvec& perm) {
    return M.submat(perm, perm);
}

// Extract constants from upper Cholesky factor (q×q) for (q-2, q-1) block update.
// Returns: {c1, c2, c3, c4}
struct CholConstants {
    double c1, c2, c3, c4;
};

static CholConstants get_constants(const arma::mat& L, int q) {
    // L is upper-triangular Cholesky factor
    // phi_q and phi_q_1 are columns q-1 and q-2, excluding rows (q-2, q-1)
    int qm1 = q - 1;  // last column (0-based)
    int qm2 = q - 2;  // second-to-last column

    double c1 = 0.0;
    for(int k = 0; k < qm2; ++k) {
        c1 += L(k, qm1) * L(k, qm2);
    }
    double c2 = L(qm2, qm2);
    // c3 excludes only row qm2 (the row being modified), not qm1.
    // In R: sum(cholm[-(q-1), q]^2)  — excludes row q-1 (1-based) = row qm2 (0-based).
    // Must include L(qm1, qm1) (the diagonal entry of the last row).
    double c3 = 0.0;
    for(int k = 0; k < qm2; ++k) {
        c3 += L(k, qm1) * L(k, qm1);
    }
    c3 += L(qm1, qm1) * L(qm1, qm1);  // diagonal entry
    double c4 = c3 + c1 * c1 / (c2 * c2);

    return {c1, c2, c3, c4};
}

// Constrained diagonal: R(omega, constants) from mixedGM
static double constrained_diagonal(double omega, const CholConstants& C) {
    if(omega == 0.0) {
        return C.c4;
    } else {
        double t = (omega - C.c1) / C.c2;
        return C.c3 + t * t;
    }
}


// =============================================================================
// update_Kyy_offdiag
// =============================================================================
// MH update for one off-diagonal element of the precision matrix Kyy_(i, j).
// Uses the mixedGM permute-based Cholesky approach:
//   1. Permute target pair to last two positions
//   2. Cholesky of permuted Kyy
//   3. Propose on the Cholesky scale: phi' ~ N(Phi[q-2, q-1], sd)
//   4. Map back to Kyy space with constrained diagonal
//   5. Unpermute
//   6. Full recompute of GGM likelihood
//
// Prior: Cauchy(0, pairwise_scale_) on off-diag, Gamma(1, 1) on diagonal.
// =============================================================================

void MixedMRFModel::update_Kyy_offdiag(int i, int j) {
    int q = static_cast<int>(q_);

    arma::uvec perm = make_perm_offdiag(i, j, q);
    arma::mat Kyy_perm = permute_matrix(Kyy_, perm);
    arma::mat Phi = arma::chol(Kyy_perm);  // upper Cholesky
    CholConstants C = get_constants(Phi, q);

    // Current Cholesky-scale value
    double phi_curr = Phi(q - 2, q - 1);
    double phi_prop = rnorm(rng_, phi_curr, prop_sd_Kyy_(i, j));

    // Map to Kyy space
    double omega_prop = C.c1 + C.c2 * phi_prop;
    double diag_curr = Kyy_perm(q - 1, q - 1);
    double diag_prop = constrained_diagonal(omega_prop, C);

    // Build proposed permuted Kyy
    arma::mat Kyy_prop_perm = Kyy_perm;
    Kyy_prop_perm(q - 2, q - 1) = omega_prop;
    Kyy_prop_perm(q - 1, q - 2) = omega_prop;
    Kyy_prop_perm(q - 1, q - 1) = diag_prop;

    // Unpermute
    arma::mat Kyy_prop = permute_matrix(Kyy_prop_perm, perm);

    // Save current state
    arma::mat Kyy_saved = Kyy_;
    arma::mat Kyy_inv_saved = Kyy_inv_;
    arma::mat Kyy_chol_saved = Kyy_chol_;
    double Kyy_log_det_saved = Kyy_log_det_;
    arma::mat cond_mean_saved = conditional_mean_;

    // Evaluate proposed likelihood
    Kyy_ = Kyy_prop;
    recompute_Kyy_decomposition();
    recompute_conditional_mean();
    double ll_prop = log_conditional_ggm();

    // Evaluate current likelihood (restore)
    Kyy_ = Kyy_saved;
    Kyy_inv_ = Kyy_inv_saved;
    Kyy_chol_ = Kyy_chol_saved;
    Kyy_log_det_ = Kyy_log_det_saved;
    conditional_mean_ = cond_mean_saved;
    double ll_curr = log_conditional_ggm();

    double ln_alpha = ll_prop - ll_curr;

    // Prior ratio: Cauchy on off-diag + Gamma(1,1) on diagonal that changed
    ln_alpha += R::dcauchy(Kyy_prop(i, j), 0.0, pairwise_scale_, true);
    ln_alpha -= R::dcauchy(Kyy_(i, j), 0.0, pairwise_scale_, true);
    ln_alpha += R::dgamma(diag_prop, 1.0, 1.0, true);
    ln_alpha -= R::dgamma(diag_curr, 1.0, 1.0, true);

    if(std::log(runif(rng_)) < ln_alpha) {
        // Accept: install proposed Kyy and refresh caches
        Kyy_ = Kyy_prop;
        recompute_Kyy_decomposition();
        recompute_conditional_mean();
        if(use_marginal_pl_) recompute_Theta();
    }
}


// =============================================================================
// update_Kyy_diag
// =============================================================================
// MH update for one diagonal element of Kyy.
// Proposes on the log-Cholesky scale to ensure positivity.
// Prior: Gamma(1, 1) on the diagonal element + Jacobian for log-scale proposal.
// =============================================================================

void MixedMRFModel::update_Kyy_diag(int i) {
    int q = static_cast<int>(q_);

    arma::uvec perm = make_perm_diag(i, q);
    arma::mat Kyy_perm = permute_matrix(Kyy_, perm);
    arma::mat L = arma::chol(Kyy_perm);  // upper Cholesky

    // Propose on log(L[q-1, q-1]) scale
    double theta_curr = std::log(L(q - 1, q - 1));
    double theta_prop = rnorm(rng_, theta_curr, prop_sd_Kyy_(i, i));

    // Rebuild Kyy_prop from modified Cholesky
    L(q - 1, q - 1) = std::exp(theta_prop);
    arma::mat Kyy_prop_perm = L.t() * L;
    arma::mat Kyy_prop = permute_matrix(Kyy_prop_perm, perm);

    // Save current state
    arma::mat Kyy_saved = Kyy_;
    arma::mat Kyy_inv_saved = Kyy_inv_;
    arma::mat Kyy_chol_saved = Kyy_chol_;
    double Kyy_log_det_saved = Kyy_log_det_;
    arma::mat cond_mean_saved = conditional_mean_;

    // Evaluate proposed likelihood
    Kyy_ = Kyy_prop;
    recompute_Kyy_decomposition();
    recompute_conditional_mean();
    double ll_prop = log_conditional_ggm();

    // Evaluate current likelihood (restore)
    Kyy_ = Kyy_saved;
    Kyy_inv_ = Kyy_inv_saved;
    Kyy_chol_ = Kyy_chol_saved;
    Kyy_log_det_ = Kyy_log_det_saved;
    conditional_mean_ = cond_mean_saved;
    double ll_curr = log_conditional_ggm();

    double ln_alpha = ll_prop - ll_curr;

    // Prior ratio: Gamma(1, 1) on diagonal
    ln_alpha += R::dgamma(Kyy_prop(i, i), 1.0, 1.0, true);
    ln_alpha -= R::dgamma(Kyy_(i, i), 1.0, 1.0, true);

    // Jacobian for log-scale proposal
    ln_alpha += theta_prop - theta_curr;

    if(std::log(runif(rng_)) < ln_alpha) {
        Kyy_ = Kyy_prop;
        recompute_Kyy_decomposition();
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
    double ll_curr = log_conditional_omrf(i) + log_conditional_ggm()
                   + R::dcauchy(current_val, 0.0, pairwise_scale_, true);

    // Set proposed value and refresh conditional_mean_
    arma::mat cond_mean_saved = conditional_mean_;
    Kxy_(i, j) = proposed;
    recompute_conditional_mean();

    double ll_prop = log_conditional_omrf(i) + log_conditional_ggm()
                   + R::dcauchy(proposed, 0.0, pairwise_scale_, true);

    double ln_alpha = ll_prop - ll_curr;

    if(std::log(runif(rng_)) >= ln_alpha) {
        Kxy_(i, j) = current_val;  // reject
        conditional_mean_ = std::move(cond_mean_saved);
    } else if(use_marginal_pl_) {
        recompute_Theta();
    }
}
