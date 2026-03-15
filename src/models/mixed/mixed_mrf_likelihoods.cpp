#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "utils/variable_helpers.h"
#include "math/explog_macros.h"


// =============================================================================
// log_conditional_omrf
// =============================================================================
// Conditional OMRF pseudolikelihood for discrete variable s:
//   log f(x_s | x_{-s}, y) = numerator - sum_v log(Z_v)
//
// Branches on is_ordinal_variable_(s) to select ordinal thresholds or
// Blume-Capel (linear + quadratic) main effects.
// =============================================================================

double MixedMRFModel::log_conditional_omrf(int s) const {
    int C_s = num_categories_(s);

    // Rest score: contribution from other discrete vars + continuous vars.
    // Factor 2 comes from the joint density x'Kxx x = 2Σ_{i<j} K_ij x_i x_j,
    // so the full conditional coefficient for x_j is 2 K_{sj}.
    arma::vec rest = 2.0 * (discrete_observations_dbl_ * pairwise_effects_discrete_.col(s)
                   - discrete_observations_dbl_.col(s) * pairwise_effects_discrete_(s, s))
                   + 2.0 * continuous_observations_ * pairwise_effects_cross_.row(s).t();

    // Numerator (sufficient-statistic form): dot(x_s, rest) + main-effect sums
    double numer = arma::dot(discrete_observations_dbl_.col(s), rest);

    if(is_ordinal_variable_(s)) {
        // Ordinal: add threshold contributions  sum_{c=1}^{C_s} count_c * main_effects_discrete_(s, c-1)
        for(int c = 1; c <= C_s; ++c) {
            numer += static_cast<double>(counts_per_category_(c, s)) * main_effects_discrete_(s, c - 1);
        }

        // Denominator via compute_denom_ordinal (FAST/SAFE block-split)
        arma::vec main_param = main_effects_discrete_.row(s).cols(0, C_s - 1).t();
        arma::vec bound = static_cast<double>(C_s) * rest;
        arma::vec denom = compute_denom_ordinal(rest, main_param, bound);

        return numer - arma::accu(bound + ARMA_MY_LOG(denom));
    } else {
        // Blume-Capel: alpha * sum(x) + beta * sum(x^2)
        double alpha = main_effects_discrete_(s, 0);
        double beta = main_effects_discrete_(s, 1);
        numer += alpha * static_cast<double>(blume_capel_stats_(0, s))
               + beta * static_cast<double>(blume_capel_stats_(1, s));

        // Denominator via compute_denom_blume_capel (computes bound internally)
        arma::vec bound;
        arma::vec denom = compute_denom_blume_capel(
            rest, alpha, beta, baseline_category_(s), C_s, bound
        );

        return numer - arma::accu(bound + ARMA_MY_LOG(denom));
    }
}


// =============================================================================
// log_marginal_omrf
// =============================================================================
// Marginal OMRF pseudolikelihood for discrete variable s:
//   log f(x_s | x_{-s}) using marginal_interactions = A_xx + 2 A_xy Σ_yy A_xy'
//
// Differs from conditional form in three ways:
//   1. rest score uses marginal_interactions_ instead of pairwise_effects_discrete_, minus self-interaction
//   2. scalar bias 2 A_xy(s,:) μ_y added to rest
//   3. numerator includes marginal(s,s) * sum(x_s^2)
//   4. denominator offsets include c^2 * marginal(s,s)
// =============================================================================

double MixedMRFModel::log_marginal_omrf(int s) const {
    int C_s = num_categories_(s);

    // Rest score: marginal-interaction-based + cross * μ_y bias
    double precision_ss = marginal_interactions_(s, s);
    arma::vec rest = discrete_observations_dbl_ * marginal_interactions_.col(s)
                   - discrete_observations_dbl_.col(s) * precision_ss
                   + 2.0 * arma::dot(pairwise_effects_cross_.row(s), main_effects_continuous_);

    // Numerator: dot(x_s, rest) + precision_ss * dot(x_s, x_s) + main effects
    double numer = arma::dot(discrete_observations_dbl_.col(s), rest)
                 + precision_ss * arma::dot(discrete_observations_dbl_.col(s),
                                        discrete_observations_dbl_.col(s));

    if(is_ordinal_variable_(s)) {
        for(int c = 1; c <= C_s; ++c) {
            numer += static_cast<double>(counts_per_category_(c, s)) * main_effects_discrete_(s, c - 1);
        }

        // Denominator: main_param(c) = μ_x(s,c) + (c+1)^2 Θ_ss
        arma::vec main_param(C_s);
        for(int c = 0; c < C_s; ++c) {
            main_param(c) = main_effects_discrete_(s, c) + static_cast<double>((c + 1) * (c + 1)) * precision_ss;
        }

        arma::vec bound = static_cast<double>(C_s) * rest;
        arma::vec denom = compute_denom_ordinal(rest, main_param, bound);

        return numer - arma::accu(bound + ARMA_MY_LOG(denom));
    } else {
        // Blume-Capel: alpha * sum(x) + beta * sum(x^2)
        double alpha = main_effects_discrete_(s, 0);
        double beta = main_effects_discrete_(s, 1);
        numer += alpha * static_cast<double>(blume_capel_stats_(0, s))
               + beta * static_cast<double>(blume_capel_stats_(1, s));

        // Denominator: theta_c includes marginal_interactions_(s,s) * (c - ref)^2
        int ref = baseline_category_(s);
        double effective_beta = beta + precision_ss;

        arma::vec bound;
        arma::vec denom = compute_denom_blume_capel(
            rest, alpha, effective_beta, ref, C_s, bound
        );

        return numer - arma::accu(bound + ARMA_MY_LOG(denom));
    }
}


// =============================================================================
// log_conditional_ggm
// =============================================================================
// Conditional GGM log-likelihood: log f(y | x)
//   y | x ~ N(conditional_mean_, covariance_continuous_)
//
// Uses cached covariance_continuous_, log_det_precision_, and conditional_mean_.
// The quadratic form uses precision = -2 * pairwise_effects_continuous_.
// =============================================================================

double MixedMRFModel::log_conditional_ggm() const {
    arma::mat D = continuous_observations_ - conditional_mean_;

    // Quadratic form: trace(Precision D'D)
    double quad_sum = arma::accu((D * (-2.0 * pairwise_effects_continuous_)) % D);

    return static_cast<double>(n_) / 2.0 *
           (-static_cast<double>(q_) * MY_LOG(2.0 * arma::datum::pi)
            + log_det_precision_)
         - quad_sum / 2.0;
}
