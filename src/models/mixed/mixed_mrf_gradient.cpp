#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"
#include "utils/variable_helpers.h"
#include "math/explog_macros.h"


// =============================================================================
// Gradient cache
// =============================================================================
// The gradient cache stores precomputed index mappings and observed-statistic
// contributions that do not change during a leapfrog trajectory.  It is
// invalidated whenever edge indicators change (same pattern as the OMRF).
// =============================================================================

void MixedMRFModel::ensure_gradient_cache() {
    if(gradient_cache_valid_) return;

    // --- Build index matrix for pairwise_effects_discrete_ upper-triangular entries ---
    // Maps (i, j) to a position in the flat gradient vector (offset from
    // the start of pairwise_discrete entries, which sits at num_main_).
    kxx_index_cache_.set_size(p_, p_);
    kxx_index_cache_.zeros();

    int num_active_kxx = 0;
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 1) {
                kxx_index_cache_(i, j) = num_main_ + num_active_kxx;
                kxx_index_cache_(j, i) = kxx_index_cache_(i, j);
                num_active_kxx++;
            }
        }
    }

    // --- Build index matrix for pairwise_effects_cross_ entries ---
    // Maps (i, j) to a position in the flat gradient vector (offset from
    // the start of pairwise_cross entries, which sits at num_main_ + active_kxx + q).
    kxy_index_cache_.set_size(p_, q_);
    kxy_index_cache_.zeros();

    int kxy_offset = num_main_ + num_active_kxx + static_cast<int>(q_);
    int num_active_kxy = 0;
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 1) {
                kxy_index_cache_(i, j) = kxy_offset + num_active_kxy;
                num_active_kxy++;
            }
        }
    }

    // main_effects_continuous_ offset in gradient vector
    main_effects_continuous_grad_offset_ = num_main_ + num_active_kxx;

    // --- Precompute observed statistics portion of the gradient ---
    size_t active_dim = num_main_ + num_active_kxx + q_ + num_active_kxy;
    grad_obs_cache_.set_size(active_dim);
    grad_obs_cache_.zeros();

    // Observed statistics for discrete main effects
    int offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            int C_s = num_categories_(s);
            for(int c = 0; c < C_s; ++c) {
                grad_obs_cache_(offset + c) = counts_per_category_(c + 1, s);
            }
            offset += C_s;
        } else {
            grad_obs_cache_(offset)     = blume_capel_stats_(0, s);
            grad_obs_cache_(offset + 1) = blume_capel_stats_(1, s);
            offset += 2;
        }
    }

    // Observed statistics for pairwise_effects_discrete_ edges
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            int loc = kxx_index_cache_(i, j);
            // Factor 2 from the symmetric double-count in the pseudo-likelihood
            grad_obs_cache_(loc) = 2.0 * arma::dot(
                discrete_observations_dbl_.col(i),
                discrete_observations_dbl_.col(j)
            );
        }
    }

    // No precomputed observed stats for means or cross effects — those depend on
    // continuous_observations_ combined with current parameters, so they
    // are computed fresh each logp_and_gradient call.

    // Cache transpose of discrete observations for vectorized pairwise gradient
    discrete_observations_dbl_t_ = discrete_observations_dbl_.t();

    gradient_cache_valid_ = true;
}


void MixedMRFModel::invalidate_gradient_cache() {
    gradient_cache_valid_ = false;
}


// =============================================================================
// Unvectorize NUTS parameters into temporaries
// =============================================================================
// Unpacks a NUTS-dimension parameter vector into temporary matrices without
// mutating model state.  Used during leapfrog trajectory evaluation.
// =============================================================================

void MixedMRFModel::unvectorize_nuts_to_temps(
    const arma::vec& params,
    arma::mat& temp_main_discrete,
    arma::mat& temp_pairwise_discrete,
    arma::vec& temp_main_continuous,
    arma::mat& temp_pairwise_cross
) const {
    size_t idx = 0;

    // 1. main_effects_discrete_
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                temp_main_discrete(s, c) = params(idx++);
            }
        } else {
            temp_main_discrete(s, 0) = params(idx++);
            temp_main_discrete(s, 1) = params(idx++);
        }
    }

    // 2. pairwise_effects_discrete_ upper-triangular (active only)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 1) {
                temp_pairwise_discrete(i, j) = params(idx++);
                temp_pairwise_discrete(j, i) = temp_pairwise_discrete(i, j);
            }
        }
    }

    // 3. main_effects_continuous_
    for(size_t j = 0; j < q_; ++j) {
        temp_main_continuous(j) = params(idx++);
    }

    // 4. pairwise_effects_cross_ row-major (active only)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 1) {
                temp_pairwise_cross(i, j) = params(idx++);
            }
        }
    }
}


// =============================================================================
// gradient
// =============================================================================

arma::vec MixedMRFModel::gradient(const arma::vec& parameters) {
    auto [logp, grad] = logp_and_gradient(parameters);
    return grad;
}


// =============================================================================
// logp_and_gradient — conditional pseudo-likelihood
// =============================================================================
// Computes the log pseudo-posterior and its gradient with respect to the
// NUTS parameters (μ_x, K_xx, μ_y, K_xy).  K_yy is treated as fixed.
//
// The pseudo-log-posterior is:
//   l(θ) = sum_s log p(x_s | x_{-s}, y)   [OMRF conditionals]
//            + log p(y | x)                     [GGM conditional]
//            + log π(θ)                    [priors]
//
// For marginal PL, the OMRF conditionals use Θ = K_xx + 2 K_xy Σ_yy K_xy'
// instead of K_xx directly, and derive rest scores and denominator offsets
// from Θ.  The GGM conditional is the same in both modes.
// =============================================================================

std::pair<double, arma::vec> MixedMRFModel::logp_and_gradient(
    const arma::vec& parameters)
{
    ensure_gradient_cache();

    // --- Unvectorize into temporaries ---
    arma::mat temp_main_discrete = main_effects_discrete_;
    arma::mat temp_pairwise_discrete = pairwise_effects_discrete_;
    arma::vec temp_main_continuous = main_effects_continuous_;
    arma::mat temp_pairwise_cross = pairwise_effects_cross_;
    unvectorize_nuts_to_temps(parameters, temp_main_discrete, temp_pairwise_discrete, temp_main_continuous, temp_pairwise_cross);

    // --- Derived quantities ---
    // Conditional mean: M_i = μ_y' + 2 x_i' K_xy Σ_yy  (n x q)
    arma::mat temp_cond_mean = arma::repmat(temp_main_continuous.t(), n_, 1)
                             + 2.0 * discrete_observations_dbl_ * temp_pairwise_cross * covariance_continuous_;

    // Residual: D = Y - M  (n x q)
    arma::mat D = continuous_observations_ - temp_cond_mean;

    // Theta for marginal PL
    arma::mat temp_Theta;
    if(use_marginal_pl_) {
        temp_Theta = temp_pairwise_discrete + 2.0 * temp_pairwise_cross * covariance_continuous_ * temp_pairwise_cross.t();
    }

    // Start gradient from observed-statistics cache
    arma::vec grad = grad_obs_cache_;

    double logp = 0.0;

    // For marginal PL: precompute K_xy Σ_yy (used in cross-contributions)
    arma::mat cross_times_cov;  // p x q
    if(use_marginal_pl_) {
        cross_times_cov = temp_pairwise_cross * covariance_continuous_;
    }

    // =========================================================================
    // Part 1: OMRF conditionals
    // =========================================================================

    int main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);

        // --- Rest score for variable s ---
        arma::vec rest;
        if(use_marginal_pl_) {
            // Marginal: Θ-based rest + K_xy μ_y bias
            double theta_ss = temp_Theta(s, s);
            rest = discrete_observations_dbl_ * temp_Theta.col(s)
                 - discrete_observations_dbl_.col(s) * theta_ss
                 + 2.0 * arma::dot(temp_pairwise_cross.row(s), temp_main_continuous);
        } else {
            // Conditional: K_xx-based rest + 2 K_xy y
            rest = discrete_observations_dbl_ * temp_pairwise_discrete.col(s)
                 - discrete_observations_dbl_.col(s) * temp_pairwise_discrete(s, s)
                 + 2.0 * continuous_observations_ * temp_pairwise_cross.row(s).t();
        }

        if(is_ordinal_variable_(s)) {
            arma::vec main_param = temp_main_discrete.row(s).cols(0, C_s - 1).t();

            // Marginal PL: absorb Theta_ss into main_param
            if(use_marginal_pl_) {
                double theta_ss = temp_Theta(s, s);
                for(int c = 0; c < C_s; ++c) {
                    main_param(c) += static_cast<double>((c + 1) * (c + 1)) * theta_ss;
                }
            }

            // bound = per-observation upper bound on log-scores for numerical
            // stability. Must cover max_c(main_param(c) + (c+1)*rest(i)).
            // The highest-category term main_param(C_s-1) + C_s*rest dominates
            // when rest > 0; category 0 (score = 0) dominates when rest << 0.
            arma::vec bound = main_param(C_s - 1) + static_cast<double>(C_s) * rest;
            bound = arma::max(bound, arma::zeros<arma::vec>(bound.n_elem));

            LogZAndProbs result = compute_logZ_and_probs_ordinal(
                main_param, rest, bound, C_s
            );

            // log pseudo-posterior contribution
            logp -= arma::accu(result.log_Z);

            // Main-effect gradient: ∂/∂main_effects_discrete_{s,c} = count_c - sum_i prob(c)
            for(int c = 0; c < C_s; ++c) {
                grad(main_effects_discrete_offset + c) -= arma::accu(result.probs.col(c + 1));
            }

            // Expected value E_s[c+1|rest] per observation
            arma::vec weights = arma::regspace<arma::vec>(1, C_s);
            arma::vec E = result.probs.cols(1, C_s) * weights;

            // Pairwise discrete gradient: sum_i x_{i,t} * (x_{i,s}+1 - E_s)
            // (uses pre-transposed discrete observations for BLAS efficiency)
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? kxx_index_cache_(s, t) : kxx_index_cache_(t, s);
                grad(loc) -= pw_grad(t);
            }

            if(use_marginal_pl_) {
                // Additional pairwise_discrete gradient from Θ_ss in denominator:
                // ∂/∂pairwise_effects_discrete_{st} through Θ_ss: zero (∂Θ_ss/∂pairwise_effects_discrete_st = δ_{st})
                // So pairwise_discrete gradient from Θ rest scores is already handled above.

                // Pairwise_cross gradient from marginal OMRF (through Θ):
                // ∂Theta_{st}/∂pairwise_effects_cross_{a,j} has two terms:
                //   = 2 [Σyy pairwise_effects_cross_t']_j δ_{as} + 2 [pairwise_effects_cross_s Σyy]_j δ_{at}
                // Self-contribution (a=s): from rest_s → pairwise_effects_cross_s
                // Cross-contribution (a=t): from rest_s → pairwise_effects_cross_t for each t≠s

                arma::vec weights_sq = arma::square(weights);
                arma::vec E_sq = result.probs.cols(1, C_s) * weights_sq;

                arma::vec diff_pw = discrete_observations_dbl_t_ *
                    (discrete_observations_dbl_.col(s) - E);
                diff_pw(s) = 0.0;

                double diff_diag = arma::dot(
                    discrete_observations_dbl_.col(s),
                    discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

                double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

                // Self-contribution: a = s
                // Off-diagonal Theta: ∂Θ_{st}/∂pairwise_effects_cross_{s,j} = 2 [Σyy pairwise_effects_cross_t']_j
                // Diagonal Theta: ∂Θ_{ss}/∂pairwise_effects_cross_{s,j} = 4 [Σyy pairwise_effects_cross_s']_j
                // Rest-score bias: ∂(2 pairwise_effects_cross_s μy)/∂pairwise_effects_cross_{s,j} = 2 μy_j
                arma::rowvec kxy_self = 2.0 * (diff_pw.t() * temp_pairwise_cross) * covariance_continuous_
                                      + 4.0 * diff_diag * cross_times_cov.row(s)
                                      + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_self(j);
                }

                // Cross-contribution: a = t, for each t ≠ s
                // ∂l_s/∂pairwise_effects_cross_{t,:} = diff_pw(t) * 2 * pairwise_effects_cross_s * Σyy
                arma::rowvec V_s = 2.0 * cross_times_cov.row(s);  // 2 K_xy_s Σ_yy
                for(size_t t = 0; t < p_; ++t) {
                    if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                    for(size_t j = 0; j < q_; ++j) {
                        if(edge_indicators_(t, p_ + j) == 0) continue;
                        int loc = kxy_index_cache_(t, j);
                        grad(loc) += diff_pw(t) * V_s(j);
                    }
                }

                // Continuous mean gradient from marginal OMRF:
                // ∂l_s/∂main_effects_continuous_j = 2 pairwise_effects_cross_{sj} * sum_i (x_{is} - E_s)
                for(size_t j = 0; j < q_; ++j) {
                    grad(main_effects_continuous_grad_offset_ + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
                }
            } else {
                // Conditional PL: pairwise_cross gradient from OMRF rest score
                // ∂/∂pairwise_effects_cross_{s,j} = 2 sum_i y_{ij} (x_{is}+1 - E_s)
                arma::rowvec kxy_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_grad_s(j);
                }
            }

            main_effects_discrete_offset += C_s;
        } else {
            // --- Blume-Capel variable ---
            int ref = baseline_category_(s);
            double lin_eff = temp_main_discrete(s, 0);
            double quad_eff = temp_main_discrete(s, 1);

            // Marginal PL: absorb Theta_ss into quadratic effect
            double effective_quad = quad_eff;
            if(use_marginal_pl_) {
                effective_quad += temp_Theta(s, s);
            }

            arma::vec bc_bound;
            LogZAndProbs result = compute_logZ_and_probs_blume_capel(
                rest, lin_eff, effective_quad, ref, C_s, bc_bound
            );

            logp -= arma::accu(result.log_Z);

            arma::vec score = arma::regspace<arma::vec>(0, C_s) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            // Main-effect gradient
            grad(main_effects_discrete_offset)     -= arma::accu(result.probs * score);
            grad(main_effects_discrete_offset + 1) -= arma::accu(result.probs * sq_score);

            // Expected score per person
            arma::vec E = result.probs * score;

            // Pairwise discrete gradient
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? kxx_index_cache_(s, t) : kxx_index_cache_(t, s);
                grad(loc) -= pw_grad(t);
            }

            if(use_marginal_pl_) {
                // Pairwise_cross gradient from marginal OMRF (same structure as ordinal)
                arma::vec E_sq = result.probs * sq_score;

                arma::vec diff_pw = discrete_observations_dbl_t_ *
                    (discrete_observations_dbl_.col(s) - E);
                diff_pw(s) = 0.0;

                double diff_diag = arma::dot(
                    discrete_observations_dbl_.col(s),
                    discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

                double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

                // Self-contribution: a = s
                arma::rowvec kxy_self = 2.0 * (diff_pw.t() * temp_pairwise_cross) * covariance_continuous_
                                      + 4.0 * diff_diag * cross_times_cov.row(s)
                                      + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_self(j);
                }

                // Cross-contribution: a = t, for each t ≠ s
                arma::rowvec V_s = 2.0 * cross_times_cov.row(s);
                for(size_t t = 0; t < p_; ++t) {
                    if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                    for(size_t j = 0; j < q_; ++j) {
                        if(edge_indicators_(t, p_ + j) == 0) continue;
                        int loc = kxy_index_cache_(t, j);
                        grad(loc) += diff_pw(t) * V_s(j);
                    }
                }

                // Continuous mean gradient from marginal OMRF
                for(size_t j = 0; j < q_; ++j) {
                    grad(main_effects_continuous_grad_offset_ + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
                }
            } else {
                // Conditional PL: pairwise_cross gradient from OMRF rest score
                arma::rowvec kxy_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_grad_s(j);
                }
            }

            main_effects_discrete_offset += 2;
        }
    }

    // Add numerator contribution to logp from discrete sufficient statistics
    // (already in grad_obs_cache_ as counts, but logp needs the actual dot-products)
    main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);
        arma::vec rest;
        if(use_marginal_pl_) {
            double theta_ss = temp_Theta(s, s);
            rest = discrete_observations_dbl_ * temp_Theta.col(s)
                 - discrete_observations_dbl_.col(s) * theta_ss
                 + 2.0 * arma::dot(temp_pairwise_cross.row(s), temp_main_continuous);
            // Theta_ss quadratic contribution
            logp += temp_Theta(s, s) * arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s));
        } else {
            rest = discrete_observations_dbl_ * temp_pairwise_discrete.col(s)
                 - discrete_observations_dbl_.col(s) * temp_pairwise_discrete(s, s)
                 + 2.0 * continuous_observations_ * temp_pairwise_cross.row(s).t();
        }
        // Numerator: dot(x_s, rest) + main-effect sums
        logp += arma::dot(discrete_observations_dbl_.col(s), rest);

        if(is_ordinal_variable_(s)) {
            for(int c = 1; c <= C_s; ++c) {
                logp += static_cast<double>(counts_per_category_(c, s)) * temp_main_discrete(s, c - 1);
            }
        } else {
            logp += temp_main_discrete(s, 0) * static_cast<double>(blume_capel_stats_(0, s))
                  + temp_main_discrete(s, 1) * static_cast<double>(blume_capel_stats_(1, s));
        }
    }

    // =========================================================================
    // Part 2: GGM conditional log-likelihood and gradients
    // =========================================================================
    // log p(y | x) = n/2 (log|K_yy| - q log(2π)) - ½ trace(K_yy D'D)
    // where D = Y - M, M_i = μ_y' + 2 x_i' K_xy Σ_yy
    //
    // K_yy is fixed, so log|K_yy| contributes to logp but not gradient.

    double quad_sum = arma::accu((D * pairwise_effects_continuous_) % D);
    logp += static_cast<double>(n_) / 2.0 *
            (-static_cast<double>(q_) * MY_LOG(2.0 * arma::datum::pi)
             + log_det_precision_)
          - quad_sum / 2.0;

    // ∂/∂μ_y: K_yy D' 1_n = K_yy sum_over_rows(D)
    arma::vec D_colsums = arma::sum(D, 0).t();  // q-vector
    arma::vec grad_main_effects_continuous_ggm = pairwise_effects_continuous_ * D_colsums;

    for(size_t j = 0; j < q_; ++j) {
        grad(main_effects_continuous_grad_offset_ + j) += grad_main_effects_continuous_ggm(j);
    }

    // ∂/∂K_xy: The GGM conditional depends on K_xy through M.
    // ∂M/∂pairwise_effects_cross_{s,j} = 2 x_s [Σ_yy]_{j,:}
    // ∂logp_ggm/∂K_xy = 2 X' D  (shortcut: K_yy Σ_yy = I eliminates K_yy)
    //
    // Correctly: ∂(−½ trace(K_yy D'D))/∂pairwise_effects_cross_{s,j}
    //   = trace(K_yy D' ∂M/∂pairwise_effects_cross_{s,j})
    //   = trace(K_yy D' · 2 x_s [Σ_yy]_{j,:})
    //   = 2 [x_s' D K_yy Σ_yy]_j
    //   = 2 [x_s' D]_j    (since K_yy Σ_yy = I)
    arma::mat grad_pairwise_effects_cross_ggm = 2.0 * discrete_observations_dbl_t_ * D;  // p x q

    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            int loc = kxy_index_cache_(i, j);
            grad(loc) += grad_pairwise_effects_cross_ggm(i, j);
        }
    }

    // =========================================================================
    // Part 3: Prior log-densities and gradient contributions
    // =========================================================================

    // --- main_effects_discrete_ priors: Beta(alpha, beta) on sigmoid scale ---
    main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            int C_s = num_categories_(s);
            for(int c = 0; c < C_s; ++c) {
                double val = temp_main_discrete(s, c);
                logp += val * main_alpha_ -
                        std::log1p(MY_EXP(val)) * (main_alpha_ + main_beta_);
                double p = 1.0 / (1.0 + MY_EXP(-val));
                grad(main_effects_discrete_offset + c) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            main_effects_discrete_offset += C_s;
        } else {
            for(int k = 0; k < 2; ++k) {
                double val = temp_main_discrete(s, k);
                logp += val * main_alpha_ -
                        std::log1p(MY_EXP(val)) * (main_alpha_ + main_beta_);
                double p = 1.0 / (1.0 + MY_EXP(-val));
                grad(main_effects_discrete_offset + k) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            main_effects_discrete_offset += 2;
        }
    }

    // --- pairwise_effects_discrete_ priors: Cauchy(0, pairwise_scale_) ---
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            int loc = kxx_index_cache_(i, j);
            double val = temp_pairwise_discrete(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            grad(loc) -= 2.0 * val / (val * val + pairwise_scale_ * pairwise_scale_);
        }
    }

    // --- main_effects_continuous_ priors: Normal(0, 1) ---
    for(size_t j = 0; j < q_; ++j) {
        double val = temp_main_continuous(j);
        logp += R::dnorm(val, 0.0, 1.0, true);
        grad(main_effects_continuous_grad_offset_ + j) -= val;  // -val from -val^2/2
    }

    // --- pairwise_effects_cross_ priors: Cauchy(0, pairwise_scale_) ---
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            int loc = kxy_index_cache_(i, j);
            double val = temp_pairwise_cross(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            grad(loc) -= 2.0 * val / (val * val + pairwise_scale_ * pairwise_scale_);
        }
    }

    return {logp, grad};
}
