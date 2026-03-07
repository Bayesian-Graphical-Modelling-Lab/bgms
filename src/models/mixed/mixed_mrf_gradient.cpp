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

    // --- Build index matrix for Kxx upper-triangular entries ---
    // Maps (i, j) to a position in the flat gradient vector (offset from
    // the start of Kxx entries, which sits at num_main_).
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

    // --- Build index matrix for Kxy entries ---
    // Maps (i, j) to a position in the flat gradient vector (offset from
    // the start of Kxy entries, which sits at num_main_ + active_kxx + q).
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

    // muy offset in gradient vector
    muy_grad_offset_ = num_main_ + num_active_kxx;

    // --- Precompute observed statistics portion of the gradient ---
    size_t active_dim = num_main_ + num_active_kxx + q_ + num_active_kxy;
    grad_obs_cache_.set_size(active_dim);
    grad_obs_cache_.zeros();

    // Observed statistics for discrete main effects (mux)
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

    // Observed statistics for Kxx edges
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

    // No precomputed observed stats for muy or Kxy — those depend on
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
    arma::mat& temp_mux,
    arma::mat& temp_Kxx,
    arma::vec& temp_muy,
    arma::mat& temp_Kxy
) const {
    size_t idx = 0;

    // 1. mux
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                temp_mux(s, c) = params(idx++);
            }
        } else {
            temp_mux(s, 0) = params(idx++);
            temp_mux(s, 1) = params(idx++);
        }
    }

    // 2. Kxx upper-triangular (active only)
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 1) {
                temp_Kxx(i, j) = params(idx++);
                temp_Kxx(j, i) = temp_Kxx(i, j);
            }
        }
    }

    // 3. muy
    for(size_t j = 0; j < q_; ++j) {
        temp_muy(j) = params(idx++);
    }

    // 4. Kxy row-major (active only)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 1) {
                temp_Kxy(i, j) = params(idx++);
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
// NUTS parameters (mux, Kxx, muy, Kxy).  Kyy is treated as fixed.
//
// The pseudo-log-posterior is:
//   l(theta) = sum_s log p(x_s | x_{-s}, y)   [OMRF conditionals]
//            + log p(y | x)                     [GGM conditional]
//            + log pi(theta)                    [priors]
//
// For marginal PL, the OMRF conditionals use Theta = Kxx + 2 Kxy Σyy Kxy'
// instead of Kxx directly, and derive rest scores and denominator offsets
// from Theta.  The GGM conditional is the same in both modes.
// =============================================================================

std::pair<double, arma::vec> MixedMRFModel::logp_and_gradient(
    const arma::vec& parameters)
{
    ensure_gradient_cache();

    // --- Unvectorize into temporaries ---
    arma::mat temp_mux = mux_;
    arma::mat temp_Kxx = Kxx_;
    arma::vec temp_muy = muy_;
    arma::mat temp_Kxy = Kxy_;
    unvectorize_nuts_to_temps(parameters, temp_mux, temp_Kxx, temp_muy, temp_Kxy);

    // --- Derived quantities ---
    // Conditional mean: M_i = muy' + 2 x_i' Kxy Σyy  (n x q)
    arma::mat temp_cond_mean = arma::repmat(temp_muy.t(), n_, 1)
                             + 2.0 * discrete_observations_dbl_ * temp_Kxy * covariance_yy_;

    // Residual: D = Y - M  (n x q)
    arma::mat D = continuous_observations_ - temp_cond_mean;

    // Theta for marginal PL
    arma::mat temp_Theta;
    if(use_marginal_pl_) {
        temp_Theta = temp_Kxx + 2.0 * temp_Kxy * covariance_yy_ * temp_Kxy.t();
    }

    // Start gradient from observed-statistics cache
    arma::vec grad = grad_obs_cache_;

    double logp = 0.0;

    // =========================================================================
    // Part 1: OMRF conditionals
    // =========================================================================

    int mux_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);

        // --- Rest score for variable s ---
        arma::vec rest;
        if(use_marginal_pl_) {
            // Marginal: Theta-based rest + Kxy*muy bias
            double theta_ss = temp_Theta(s, s);
            rest = discrete_observations_dbl_ * temp_Theta.col(s)
                 - discrete_observations_dbl_.col(s) * theta_ss
                 + 2.0 * arma::dot(temp_Kxy.row(s), temp_muy);
        } else {
            // Conditional: Kxx-based rest + 2 Kxy y
            rest = discrete_observations_dbl_ * temp_Kxx.col(s)
                 - discrete_observations_dbl_.col(s) * temp_Kxx(s, s)
                 + 2.0 * continuous_observations_ * temp_Kxy.row(s).t();
        }

        arma::vec bound = static_cast<double>(C_s) * rest;

        if(is_ordinal_variable_(s)) {
            arma::vec main_param = temp_mux.row(s).cols(0, C_s - 1).t();

            // Marginal PL: absorb Theta_ss into main_param
            if(use_marginal_pl_) {
                double theta_ss = temp_Theta(s, s);
                for(int c = 0; c < C_s; ++c) {
                    main_param(c) += static_cast<double>((c + 1) * (c + 1)) * theta_ss;
                }
            }

            LogZAndProbs result = compute_logZ_and_probs_ordinal(
                main_param, rest, bound, C_s
            );

            // log pseudo-posterior contribution
            logp -= arma::accu(result.log_Z);

            // Main-effect gradient: ∂/∂mux_{s,c} = count_c - sum_i prob(c)
            for(int c = 0; c < C_s; ++c) {
                grad(mux_offset + c) -= arma::accu(result.probs.col(c + 1));
            }

            // Expected value E_s[c+1|rest] per observation
            arma::vec weights = arma::regspace<arma::vec>(1, C_s);
            arma::vec E = result.probs.cols(1, C_s) * weights;

            // Kxx pairwise gradient: sum_i x_{i,t} * (x_{i,s}+1 - E_s)
            // (uses pre-transposed discrete observations for BLAS efficiency)
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? kxx_index_cache_(s, t) : kxx_index_cache_(t, s);
                grad(loc) -= pw_grad(t);
            }

            if(use_marginal_pl_) {
                // Additional Kxx gradient from Theta_ss in denominator:
                // ∂/∂Kxx_{st} through Theta_ss: zero (∂Theta_ss/∂Kxx_st = δ_{st})
                // So Kxx gradient from Theta rest scores is already handled above.

                // Kxy gradient from marginal OMRF (through Theta):
                // E_s[c^2] per observation
                arma::vec weights_sq = arma::square(weights);
                arma::vec E_sq = result.probs.cols(1, C_s) * weights_sq;

                // Through Theta rest scores: for each (s, t!=s):
                //   ∂Theta_{st}/∂Kxy_{a,j} = 2 Σyy_{j,:} Kxy_{t,:}' δ_{as}
                //                           + 2 Kxy_{s,:} Σyy_{:,j} δ_{at}  -- but a=s here
                // The contribution is: sum_i x_{it} * (x_{is}+1 - E_s)
                //   times ∂Theta_{st}/∂Kxy_{s,j} = 2 [Σyy Kxy_t']_j
                // Summed over t!=s: this equals 2 * sum_t!=s pw_grad(t) * [Σyy Kxy_t']
                // Plus the diagonal term: Theta_ss gives
                //   ∂/∂Kxy_{s,j} of Theta_ss = 4 [Σyy Kxy_s']_j
                //   contribution: -sum_i E_sq * 4 [Σyy Kxy_s']_j
                // Actually let's use the vector form:

                // Off-diagonal Theta contribution to Kxy_s gradient:
                // For variable s, the rest score R_s uses Theta_{st} for t!=s
                // ∂R_{is}/∂Kxy_{s,j} = sum_{t!=s} x_{it} ∂Theta_{st}/∂Kxy_{s,j}
                //                     = 2 sum_{t!=s} x_{it} [Σyy Kxy_t']_j
                //                     = 2 [X_{-s} Kxy_{-s} Σyy]_{i,j}
                // Gradient = sum_i (x_{is}+1 - E_s(i)) * ∂R_{is}/∂Kxy_{s,j}
                //          = sum_i (obs_s - E_s) * 2 [X_{-s} Kxy_{-s} Σyy]_{i,j}
                // But X_{-s} means "all columns except s" — complex indexing.
                //
                // Simpler: ∂/∂Kxy_{s,:} of the Theta-based omrf for variable s:
                //   ∂Θ_{st}/∂Kxy_{s,:} = 2 Σyy Kxy_t'  for each t (including t=s)
                //   Full rest-score chain rule:
                //     ∂l_s/∂Kxy_{s,:} = sum_{t!=s} (obs_stat_{st} - exp_stat_{st})
                //                       * 2 Σyy Kxy_t'
                //                     + (obs_sq_s - exp_sq_s) * 2 Σyy Kxy_s'
                //   where obs_stat_{st} = sum_i x_{it} x_{is}
                //                exp_stat_{st} = sum_i x_{it} E_s
                //                obs_sq_s = sum_i x_{is}^2
                //                exp_sq_s = sum_i E_sq

                // obs_minus_E for the pairwise statistics: already pw_grad = X' E
                // obs_stat_{st} = sum_i x_{it}*x_{is} (cached in grad_obs for Kxx edges)
                // So: obs - exp for the pairwise is exactly pw_grad subtracted from obs.
                // But we already computed that as a vector over all t.

                // Vector form for Kxy_s gradient from Theta:
                // diff_t = X'(x_s - E)  (p-vector, where element t = sum_i x_it(x_is - E_s))
                arma::vec diff_pw = discrete_observations_dbl_t_ *
                    (discrete_observations_dbl_.col(s) - E);
                diff_pw(s) = 0.0;  // exclude self-interaction from rest score

                // diff_diag = sum_i (x_is^2 - E_sq)
                double diff_diag = arma::dot(
                    discrete_observations_dbl_.col(s),
                    discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

                // ∂l_s/∂Kxy_{s,:} from Theta = 2 (sum_{t!=s} diff_pw(t) * Kxy_t + diff_diag * Kxy_s) Σyy
                arma::rowvec kxy_contrib = 2.0 * (diff_pw.t() * temp_Kxy + diff_diag * temp_Kxy.row(s)) * covariance_yy_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_contrib(j);
                }

                // muy gradient from marginal OMRF:
                // ∂l_s/∂muy_j = 2 Kxy_{sj} * sum_i (x_{is} - E_s)
                double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);
                for(size_t j = 0; j < q_; ++j) {
                    grad(muy_grad_offset_ + j) += 2.0 * temp_Kxy(s, j) * sum_obs_minus_E;
                }
            } else {
                // Conditional PL: Kxy gradient from OMRF rest score
                // ∂/∂Kxy_{s,j} = 2 sum_i y_{ij} (x_{is}+1 - E_s)
                arma::rowvec kxy_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_grad_s(j);
                }
            }

            mux_offset += C_s;
        } else {
            // --- Blume-Capel variable ---
            int ref = baseline_category_(s);
            double lin_eff = temp_mux(s, 0);
            double quad_eff = temp_mux(s, 1);

            // Marginal PL: absorb Theta_ss into quadratic effect
            double effective_quad = quad_eff;
            if(use_marginal_pl_) {
                effective_quad += temp_Theta(s, s);
            }

            LogZAndProbs result = compute_logZ_and_probs_blume_capel(
                rest, lin_eff, effective_quad, ref, C_s, bound
            );

            logp -= arma::accu(result.log_Z);

            arma::vec score = arma::regspace<arma::vec>(0, C_s) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            // Main-effect gradient
            grad(mux_offset)     -= arma::accu(result.probs * score);
            grad(mux_offset + 1) -= arma::accu(result.probs * sq_score);

            // Expected score per person
            arma::vec E = result.probs * score;

            // Kxx pairwise gradient
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? kxx_index_cache_(s, t) : kxx_index_cache_(t, s);
                grad(loc) -= pw_grad(t);
            }

            if(use_marginal_pl_) {
                // Kxy gradient from marginal OMRF (same structure as ordinal)
                arma::vec E_sq = result.probs * sq_score;

                arma::vec diff_pw = discrete_observations_dbl_t_ *
                    (discrete_observations_dbl_.col(s) - E);
                diff_pw(s) = 0.0;

                double diff_diag = arma::dot(
                    discrete_observations_dbl_.col(s),
                    discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

                arma::rowvec kxy_contrib = 2.0 * (diff_pw.t() * temp_Kxy + diff_diag * temp_Kxy.row(s)) * covariance_yy_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_contrib(j);
                }

                // muy gradient from marginal OMRF
                double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);
                for(size_t j = 0; j < q_; ++j) {
                    grad(muy_grad_offset_ + j) += 2.0 * temp_Kxy(s, j) * sum_obs_minus_E;
                }
            } else {
                // Conditional PL: Kxy gradient from OMRF rest score
                arma::rowvec kxy_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = kxy_index_cache_(s, j);
                    grad(loc) += kxy_grad_s(j);
                }
            }

            mux_offset += 2;
        }
    }

    // Add numerator contribution to logp from discrete sufficient statistics
    // (already in grad_obs_cache_ as counts, but logp needs the actual dot-products)
    mux_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);
        arma::vec rest;
        if(use_marginal_pl_) {
            double theta_ss = temp_Theta(s, s);
            rest = discrete_observations_dbl_ * temp_Theta.col(s)
                 - discrete_observations_dbl_.col(s) * theta_ss
                 + 2.0 * arma::dot(temp_Kxy.row(s), temp_muy);
            // Theta_ss quadratic contribution
            logp += temp_Theta(s, s) * arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s));
        } else {
            rest = discrete_observations_dbl_ * temp_Kxx.col(s)
                 - discrete_observations_dbl_.col(s) * temp_Kxx(s, s)
                 + 2.0 * continuous_observations_ * temp_Kxy.row(s).t();
        }
        // Numerator: dot(x_s, rest) + main-effect sums
        logp += arma::dot(discrete_observations_dbl_.col(s), rest);

        if(is_ordinal_variable_(s)) {
            for(int c = 1; c <= C_s; ++c) {
                logp += static_cast<double>(counts_per_category_(c, s)) * temp_mux(s, c - 1);
            }
        } else {
            logp += temp_mux(s, 0) * static_cast<double>(blume_capel_stats_(0, s))
                  + temp_mux(s, 1) * static_cast<double>(blume_capel_stats_(1, s));
        }
    }

    // =========================================================================
    // Part 2: GGM conditional log-likelihood and gradients
    // =========================================================================
    // log p(y | x) = n/2 * (log|Kyy| - q log(2pi)) - 1/2 trace(Kyy D'D)
    // where D = Y - M, M_i = muy' + 2 x_i' Kxy Σyy
    //
    // Kyy is fixed, so log|Kyy| contributes to logp but not gradient.

    double quad_sum = arma::accu((D * Kyy_) % D);
    logp += static_cast<double>(n_) / 2.0 *
            (-static_cast<double>(q_) * std::log(2.0 * arma::datum::pi)
             + Kyy_log_det_)
          - quad_sum / 2.0;

    // ∂/∂muy: Kyy * D' * 1_n = Kyy * sum_over_rows(D)
    arma::vec D_colsums = arma::sum(D, 0).t();  // q-vector
    arma::vec grad_muy_ggm = Kyy_ * D_colsums;

    for(size_t j = 0; j < q_; ++j) {
        grad(muy_grad_offset_ + j) += grad_muy_ggm(j);
    }

    // ∂/∂Kxy: The GGM conditional depends on Kxy through M.
    // ∂M/∂Kxy_{s,j} = 2 x_s [Σyy]_{j,:}
    // ∂logp_ggm/∂Kxy = 2 X' D  (shortcut: Kyy Σyy = I eliminates Kyy)
    //
    // Correctly: ∂(−½ trace(Kyy D'D))/∂Kxy_{s,j}
    //   = trace(Kyy D' ∂M/∂Kxy_{s,j})
    //   = trace(Kyy D' * 2 x_s [Σyy]_{j,:})
    //   = 2 [x_s' D Kyy Σyy]_j
    //   = 2 [x_s' D]_j    (since Kyy Σyy = I)
    arma::mat grad_Kxy_ggm = 2.0 * discrete_observations_dbl_t_ * D;  // p x q

    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            int loc = kxy_index_cache_(i, j);
            grad(loc) += grad_Kxy_ggm(i, j);
        }
    }

    // =========================================================================
    // Part 3: Prior log-densities and gradient contributions
    // =========================================================================

    // --- mux priors: Beta(alpha, beta) on sigmoid scale ---
    mux_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            int C_s = num_categories_(s);
            for(int c = 0; c < C_s; ++c) {
                double val = temp_mux(s, c);
                logp += val * main_alpha_ -
                        std::log1p(MY_EXP(val)) * (main_alpha_ + main_beta_);
                double p = 1.0 / (1.0 + MY_EXP(-val));
                grad(mux_offset + c) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            mux_offset += C_s;
        } else {
            for(int k = 0; k < 2; ++k) {
                double val = temp_mux(s, k);
                logp += val * main_alpha_ -
                        std::log1p(MY_EXP(val)) * (main_alpha_ + main_beta_);
                double p = 1.0 / (1.0 + MY_EXP(-val));
                grad(mux_offset + k) += main_alpha_ - (main_alpha_ + main_beta_) * p;
            }
            mux_offset += 2;
        }
    }

    // --- Kxx priors: Cauchy(0, pairwise_scale_) ---
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            int loc = kxx_index_cache_(i, j);
            double val = temp_Kxx(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            grad(loc) -= 2.0 * val / (val * val + pairwise_scale_ * pairwise_scale_);
        }
    }

    // --- muy priors: Normal(0, 1) ---
    for(size_t j = 0; j < q_; ++j) {
        double val = temp_muy(j);
        logp += R::dnorm(val, 0.0, 1.0, true);
        grad(muy_grad_offset_ + j) -= val;  // ∂/∂muy: -muy (from -muy^2/2)
    }

    // --- Kxy priors: Cauchy(0, pairwise_scale_) ---
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            int loc = kxy_index_cache_(i, j);
            double val = temp_Kxy(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            grad(loc) -= 2.0 * val / (val * val + pairwise_scale_ * pairwise_scale_);
        }
    }

    return {logp, grad};
}
