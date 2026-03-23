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
    disc_index_cache_.set_size(p_, p_);
    disc_index_cache_.zeros();

    int num_active_disc = 0;
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 1) {
                disc_index_cache_(i, j) = num_main_ + num_active_disc;
                disc_index_cache_(j, i) = disc_index_cache_(i, j);
                num_active_disc++;
            }
        }
    }

    // --- Build index matrix for pairwise_effects_cross_ entries ---
    // Maps (i, j) to a position in the flat gradient vector (offset from
    // the start of pairwise_cross entries, which sits at num_main_ + active_kxx + q).
    cross_index_cache_.set_size(p_, q_);
    cross_index_cache_.zeros();

    int cross_offset = num_main_ + num_active_disc + static_cast<int>(q_);
    int num_active_cross = 0;
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 1) {
                cross_index_cache_(i, j) = cross_offset + num_active_cross;
                num_active_cross++;
            }
        }
    }

    // main_effects_continuous_ offset in gradient vector
    main_effects_continuous_grad_offset_ = num_main_ + num_active_disc;

    // --- Precompute observed statistics portion of the gradient ---
    size_t active_dim = num_main_ + num_active_disc + q_ + num_active_cross
                      + num_cholesky_;
    grad_obs_cache_.set_size(active_dim);
    grad_obs_cache_.zeros();

    // Cholesky block offset in gradient vector
    chol_grad_offset_ = num_main_ + num_active_disc + q_ + num_active_cross;

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
            int loc = disc_index_cache_(i, j);
            // Factor 4: K = σ, and the log-PL has edge (i,j) in two conditionals,
            // giving d/dK [4K·(x^Tx)] = 4·(x^Tx)
            grad_obs_cache_(loc) = 4.0 * arma::dot(
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
// NUTS parameters (μ_x, A_xx, μ_y, A_xy, R) where R is the Cholesky
// factor of the continuous precision matrix Ω = R^T R.
//
// The pseudo-log-posterior is:
//   l(θ) = sum_s log p(x_s | x_{-s}, y)   [OMRF conditionals]
//            + log p(y | x)                [GGM conditional]
//            + log π(θ)                    [priors on all params]
//            + log |det J|                 [Cholesky Jacobian]
//
// For marginal PL, the OMRF conditionals use Θ = A_xx + 2 A_xy Σ_yy A_xy'
// instead of A_xx directly.  The GGM conditional is the same in both modes.
// =============================================================================

std::pair<double, arma::vec> MixedMRFModel::logp_and_gradient(
    const arma::vec& parameters)
{
    ensure_gradient_cache();

    // --- Unvectorize into temporaries (blocks 1–4) ---
    arma::mat temp_main_discrete = main_effects_discrete_;
    arma::mat temp_pairwise_discrete = pairwise_effects_discrete_;
    arma::vec temp_main_continuous = main_effects_continuous_;
    arma::mat temp_pairwise_cross = pairwise_effects_cross_;
    unvectorize_nuts_to_temps(parameters, temp_main_discrete, temp_pairwise_discrete, temp_main_continuous, temp_pairwise_cross);

    // --- Unpack block 5: Cholesky of precision ---
    arma::mat temp_cholesky(q_, q_, arma::fill::zeros);
    size_t chol_idx = static_cast<size_t>(chol_grad_offset_);
    for(size_t j = 0; j < q_; ++j) {
        for(size_t i = 0; i < j; ++i) {
            temp_cholesky(i, j) = parameters(chol_idx++);
        }
        temp_cholesky(j, j) = std::exp(parameters(chol_idx++));
    }

    // Guard against degenerate Cholesky (extreme theta pushed by leapfrog)
    double min_diag = temp_cholesky.diag().min();
    if(!std::isfinite(min_diag) || min_diag < 1e-15) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(parameters.n_elem, arma::fill::zeros)};
    }

    arma::mat temp_precision = temp_cholesky.t() * temp_cholesky;
    arma::mat temp_inv_chol;
    bool solve_ok = arma::solve(temp_inv_chol, arma::trimatu(temp_cholesky),
                                arma::eye(q_, q_), arma::solve_opts::fast);
    if(!solve_ok) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(parameters.n_elem, arma::fill::zeros)};
    }
    arma::mat temp_covariance = temp_inv_chol * temp_inv_chol.t();
    double temp_log_det = 2.0 * arma::sum(arma::log(temp_cholesky.diag()));

    // --- Derived quantities ---
    // Conditional mean: M_i = μ_y' + 2 x_i' A_xy Σ_yy  (n x q)
    arma::mat temp_cond_mean = arma::repmat(temp_main_continuous.t(), n_, 1)
                             + 2.0 * discrete_observations_dbl_ * temp_pairwise_cross * temp_covariance;

    // Residual: D = Y - M  (n x q)
    arma::mat D = continuous_observations_ - temp_cond_mean;

    // Marginal PL effective discrete interaction matrix
    arma::mat temp_marginal;
    if(use_marginal_pl_) {
        temp_marginal = 2.0 * temp_pairwise_discrete + 2.0 * temp_pairwise_cross * temp_covariance * temp_pairwise_cross.t();
    }

    // Start gradient from observed-statistics cache
    arma::vec grad = grad_obs_cache_;

    double logp = 0.0;

    // For marginal PL: precompute A_xy Σ_yy (used in cross-contributions)
    arma::mat cross_times_cov;  // p x q
    arma::mat Theta_bar;        // p x p marginal-PL coupling for precision gradient
    if(use_marginal_pl_) {
        cross_times_cov = temp_pairwise_cross * temp_covariance;
        Theta_bar = arma::zeros<arma::mat>(p_, p_);
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
            // Marginal: Θ-based rest + A_xy μ_y bias
            double precision_ss = temp_marginal(s, s);
            rest = discrete_observations_dbl_ * temp_marginal.col(s)
                 - discrete_observations_dbl_.col(s) * precision_ss
                 + 2.0 * arma::dot(temp_pairwise_cross.row(s), temp_main_continuous);
        } else {
            // Conditional: 2 * discrete_int rest + 2 * cross_int * y
            rest = 2.0 * (discrete_observations_dbl_ * temp_pairwise_discrete.col(s)
                 - discrete_observations_dbl_.col(s) * temp_pairwise_discrete(s, s))
                 + 2.0 * continuous_observations_ * temp_pairwise_cross.row(s).t();
        }

        if(is_ordinal_variable_(s)) {
            arma::vec main_param = temp_main_discrete.row(s).cols(0, C_s - 1).t();

            // Marginal PL: absorb marginal self-interaction into main_param
            if(use_marginal_pl_) {
                double precision_ss = temp_marginal(s, s);
                for(int c = 0; c < C_s; ++c) {
                    main_param(c) += static_cast<double>((c + 1) * (c + 1)) * precision_ss;
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
            // Factor 2: chain rule d/dK = 2 × d/dσ
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? disc_index_cache_(s, t) : disc_index_cache_(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
            }

            if(use_marginal_pl_) {
                // Additional pairwise_discrete gradient from Θ_ss in denominator:
                // ∂/∂pairwise_effects_discrete_{st} through Θ_ss: zero (∂Θ_ss/∂pairwise_effects_discrete_st = δ_{st})
                // So pairwise_discrete gradient from Θ rest scores is already handled above.

                // Pairwise_cross gradient from marginal OMRF (through Θ):
                // ∂marginal_{st}/∂pairwise_effects_cross_{a,j} has two terms:
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

                // Accumulate Θ̄ for precision gradient coupling
                for(size_t t = 0; t < p_; ++t) {
                    if(t != s) Theta_bar(s, t) += diff_pw(t);
                }
                Theta_bar(s, s) += diff_diag;

                // Self-contribution: a = s
                // Off-diagonal effective interaction: ∂Θ_{st}/∂pairwise_effects_cross_{s,j} = 2 [Σyy pairwise_effects_cross_t']_j
                // Diagonal effective interaction: ∂Θ_{ss}/∂pairwise_effects_cross_{s,j} = 4 [Σyy pairwise_effects_cross_s']_j
                // Rest-score bias: ∂(2 pairwise_effects_cross_s μy)/∂pairwise_effects_cross_{s,j} = 2 μy_j
                arma::rowvec cross_self = 2.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                      + 4.0 * diff_diag * cross_times_cov.row(s)
                                      + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = cross_index_cache_(s, j);
                    grad(loc) += cross_self(j);
                }

                // Cross-contribution: a = t, for each t ≠ s
                // ∂l_s/∂pairwise_effects_cross_{t,:} = diff_pw(t) * 2 * pairwise_effects_cross_s * Σyy
                arma::rowvec V_s = 2.0 * cross_times_cov.row(s);  // 2 A_xy_s Σ_yy
                for(size_t t = 0; t < p_; ++t) {
                    if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                    for(size_t j = 0; j < q_; ++j) {
                        if(edge_indicators_(t, p_ + j) == 0) continue;
                        int loc = cross_index_cache_(t, j);
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
                arma::rowvec cross_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = cross_index_cache_(s, j);
                    grad(loc) += cross_grad_s(j);
                }
            }

            main_effects_discrete_offset += C_s;
        } else {
            // --- Blume-Capel variable ---
            int ref = baseline_category_(s);
            double lin_eff = temp_main_discrete(s, 0);
            double quad_eff = temp_main_discrete(s, 1);

            // Marginal PL: absorb marginal self-interaction into quadratic effect
            double effective_quad = quad_eff;
            if(use_marginal_pl_) {
                effective_quad += temp_marginal(s, s);
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
            // Factor 2: chain rule d/dK = 2 × d/dσ
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(edge_indicators_(s, t) == 0 || s == t) continue;
                int loc = (s < t) ? disc_index_cache_(s, t) : disc_index_cache_(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
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

                // Accumulate Θ̄ for precision gradient coupling
                for(size_t t = 0; t < p_; ++t) {
                    if(t != s) Theta_bar(s, t) += diff_pw(t);
                }
                Theta_bar(s, s) += diff_diag;

                // Self-contribution: a = s
                arma::rowvec cross_self = 2.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                      + 4.0 * diff_diag * cross_times_cov.row(s)
                                      + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = cross_index_cache_(s, j);
                    grad(loc) += cross_self(j);
                }

                // Cross-contribution: a = t, for each t ≠ s
                arma::rowvec V_s = 2.0 * cross_times_cov.row(s);
                for(size_t t = 0; t < p_; ++t) {
                    if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                    for(size_t j = 0; j < q_; ++j) {
                        if(edge_indicators_(t, p_ + j) == 0) continue;
                        int loc = cross_index_cache_(t, j);
                        grad(loc) += diff_pw(t) * V_s(j);
                    }
                }

                // Continuous mean gradient from marginal OMRF
                for(size_t j = 0; j < q_; ++j) {
                    grad(main_effects_continuous_grad_offset_ + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
                }
            } else {
                // Conditional PL: pairwise_cross gradient from OMRF rest score
                arma::rowvec cross_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    if(edge_indicators_(s, p_ + j) == 0) continue;
                    int loc = cross_index_cache_(s, j);
                    grad(loc) += cross_grad_s(j);
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
            double precision_ss = temp_marginal(s, s);
            rest = discrete_observations_dbl_ * temp_marginal.col(s)
                 - discrete_observations_dbl_.col(s) * precision_ss
                 + 2.0 * arma::dot(temp_pairwise_cross.row(s), temp_main_continuous);
            // Marginal self-interaction quadratic contribution
            logp += temp_marginal(s, s) * arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s));
        } else {
            rest = 2.0 * (discrete_observations_dbl_ * temp_pairwise_discrete.col(s)
                 - discrete_observations_dbl_.col(s) * temp_pairwise_discrete(s, s))
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
    // log p(y | x) = n/2 (log|Ω| - q log(2π)) - ½ trace(Ω D'D)
    // where Ω = R'R (precision), D = Y - M

    double quad_sum = arma::accu((D * temp_precision) % D);
    logp += static_cast<double>(n_) / 2.0 *
            (-static_cast<double>(q_) * MY_LOG(2.0 * arma::datum::pi)
             + temp_log_det)
          - quad_sum / 2.0;

    // ∂/∂μ_y: Ω * sum_over_rows(D)
    arma::vec D_colsums = arma::sum(D, 0).t();  // q-vector
    arma::vec grad_main_effects_continuous_ggm = temp_precision * D_colsums;

    for(size_t j = 0; j < q_; ++j) {
        grad(main_effects_continuous_grad_offset_ + j) += grad_main_effects_continuous_ggm(j);
    }

    // ∂/∂A_xy: The GGM conditional depends on A_xy through M.
    // ∂M/∂pairwise_effects_cross_{s,j} = 2 x_s [Σ_yy]_{j,:}
    // ∂logp_ggm/∂A_xy = 2 X' D  (shortcut: Θ Σ_yy = I eliminates Θ)
    //
    // Correctly: ∂(−½ trace(Θ D'D))/∂pairwise_effects_cross_{s,j}
    //   = trace(Θ D' ∂M/∂pairwise_effects_cross_{s,j})
    //   = trace(Θ D' · 2 x_s [Σ_yy]_{j,:})
    //   = 2 [x_s' D Θ Σ_yy]_j
    //   = 2 [x_s' D]_j    (since Θ Σ_yy = I)
    arma::mat grad_pairwise_effects_cross_ggm = 2.0 * discrete_observations_dbl_t_ * D;  // p x q

    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            if(edge_indicators_(i, p_ + j) == 0) continue;
            int loc = cross_index_cache_(i, j);
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
    const double disc_scale = pairwise_scale_;
    const double disc_scale_sq = disc_scale * disc_scale;
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            if(edge_indicators_(i, j) == 0) continue;
            int loc = disc_index_cache_(i, j);
            double val = temp_pairwise_discrete(i, j);
            logp += R::dcauchy(val, 0.0, disc_scale, true);
            grad(loc) -= 2.0 * val / (val * val + disc_scale_sq);
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
            int loc = cross_index_cache_(i, j);
            double val = temp_pairwise_cross(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            grad(loc) -= 2.0 * val / (val * val + pairwise_scale_ * pairwise_scale_);
        }
    }

    // =========================================================================
    // Part 4: Precision gradient via Cholesky parameterization
    // =========================================================================
    // Compute Ω̄ = ∂ℓ/∂Ω, then map to R̄ = ∂ℓ/∂R, then to position gradient.
    //
    // Ω̄ = (n/2) Σ − ½ D^T D − 2 Σ A_xy^T X^T D + priors on Ω
    //     + [marginal PL coupling through Θ]

    // --- Phase 1: GGM conditional contribution ---
    arma::mat Omega_bar = 0.5 * static_cast<double>(n_) * temp_covariance
                        - 0.5 * D.t() * D;

    // --- Phase 2: Conditional-mean coupling ---
    // M_i = μ_y + 2 Σ A_xy^T x_i depends on Σ = Ω^{-1}
    // ∂ℓ/∂Σ_{ab} from GGM conditional = 2 [A_xy^T X^T D Ω]_{ab}
    // Mapping: ∂ℓ/∂Ω += −Σ (∂ℓ/∂Σ) Σ = −2 Σ A_xy^T X^T D
    Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
               * discrete_observations_dbl_t_ * D;

    // --- Phase 2b: Marginal PL coupling through Θ ---
    // Θ = 2 A_xx + 2 A_xy Σ A_xy^T depends on Σ
    // ∂Θ/∂Σ = 2 A_xy ⊗ A_xy  →  ∂ℓ/∂Σ = 2 A_xy^T Θ̄ A_xy
    // ∂ℓ/∂Ω += −Σ (∂ℓ/∂Σ) Σ = −2 Σ A_xy^T Θ̄ A_xy Σ
    if(use_marginal_pl_) {
        Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
                   * Theta_bar * temp_pairwise_cross * temp_covariance;
    }

    // --- Phase 3: Priors on precision entries ---
    // Gamma(1, 1) on diagonal Ω_{jj}: log π(Ω_{jj}) = −Ω_{jj} + const
    for(size_t j = 0; j < q_; ++j) {
        Omega_bar(j, j) -= 1.0;
        logp -= temp_precision(j, j);  // Gamma(1,1) log-density
    }
    // Cauchy(0, scale) on off-diagonal Ω_{ij} (upper triangle only)
    // Only add to Omega_bar(i,j), not (j,i): the symmetrization
    // Ω̄ + Ω̄ᵀ in Phase 4 handles the lower triangle automatically.
    const double cont_scale_sq = pairwise_scale_ * pairwise_scale_;
    for(size_t i = 0; i < q_ - 1; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            double val = temp_precision(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            double cauchy_grad = -2.0 * val / (val * val + cont_scale_sq);
            Omega_bar(i, j) += cauchy_grad;
        }
    }

    // --- Phase 4: Map Ω̄ → R̄ → position gradient ---
    // R̄ = R (Ω̄ + Ω̄^T)
    arma::mat Omega_bar_sym = Omega_bar + Omega_bar.t();
    arma::mat R_bar = temp_cholesky * Omega_bar_sym;

    // Cholesky Jacobian: log|det J| = q log 2 + Σ_j (q − j + 1) ψ_j
    // where (q − j + 1) = (q − j) from Bartlett + 1 from exp(ψ_j).
    size_t gidx = static_cast<size_t>(chol_grad_offset_);
    for(size_t j = 0; j < q_; ++j) {
        double psi_j = std::log(temp_cholesky(j, j));
        double jac_weight = static_cast<double>(q_ - j + 1);
        logp += jac_weight * psi_j;

        // Off-diagonal Cholesky entries: ∂ℓ/∂R_{ij} = R̄_{ij}
        for(size_t i = 0; i < j; ++i) {
            grad(gidx++) = R_bar(i, j);
        }
        // Diagonal (log-scale): ∂ℓ/∂ψ_j = R̄_{jj} R_{jj} + (q − j + 1)
        grad(gidx++) = R_bar(j, j) * temp_cholesky(j, j) + jac_weight;
    }
    // Add constant Jacobian term to logp
    logp += static_cast<double>(q_) * std::log(2.0);

    return {logp, grad};
}


// =============================================================================
// logp_and_gradient_full — full-space version for RATTLE
// =============================================================================
// Same computation as logp_and_gradient() but with fixed full-dimension
// indexing: all Kxx, Kxy slots present (zeros for excluded edges).
// No edge gating — gradient is computed for every parameter.
// The Cholesky Jacobian is included (parameterization, not constraint).
// =============================================================================

std::pair<double, arma::vec> MixedMRFModel::logp_and_gradient_full(
    const arma::vec& x)
{
    const size_t full_dim = full_parameter_dimension();

    // --- Unpack all 5 blocks from full-space vector ---
    size_t idx = 0;

    // Block 1: main_effects_discrete_
    arma::mat temp_main_discrete(p_, max_cats_, arma::fill::zeros);
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            for(int c = 0; c < num_categories_(s); ++c) {
                temp_main_discrete(s, c) = x(idx++);
            }
        } else {
            temp_main_discrete(s, 0) = x(idx++);
            temp_main_discrete(s, 1) = x(idx++);
        }
    }

    // Block 2: ALL pairwise_effects_discrete_ upper-triangular
    arma::mat temp_pairwise_discrete(p_, p_, arma::fill::zeros);
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            temp_pairwise_discrete(i, j) = x(idx);
            temp_pairwise_discrete(j, i) = x(idx);
            idx++;
        }
    }

    // Block 3: main_effects_continuous_
    arma::vec temp_main_continuous(q_);
    for(size_t j = 0; j < q_; ++j) {
        temp_main_continuous(j) = x(idx++);
    }

    // Block 4: ALL pairwise_effects_cross_ row-major
    arma::mat temp_pairwise_cross(p_, q_, arma::fill::zeros);
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            temp_pairwise_cross(i, j) = x(idx++);
        }
    }

    // Block 5: Cholesky of precision
    const size_t chol_offset = idx;
    arma::mat temp_cholesky(q_, q_, arma::fill::zeros);
    for(size_t j = 0; j < q_; ++j) {
        for(size_t i = 0; i < j; ++i) {
            temp_cholesky(i, j) = x(idx++);
        }
        temp_cholesky(j, j) = std::exp(x(idx++));
    }

    // Guard against degenerate Cholesky
    double min_diag = temp_cholesky.diag().min();
    if(!std::isfinite(min_diag) || min_diag < 1e-15) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(full_dim, arma::fill::zeros)};
    }

    arma::mat temp_precision = temp_cholesky.t() * temp_cholesky;
    arma::mat temp_inv_chol;
    bool solve_ok = arma::solve(temp_inv_chol, arma::trimatu(temp_cholesky),
                                arma::eye(q_, q_), arma::solve_opts::fast);
    if(!solve_ok) {
        return {-std::numeric_limits<double>::infinity(),
                arma::vec(full_dim, arma::fill::zeros)};
    }
    arma::mat temp_covariance = temp_inv_chol * temp_inv_chol.t();
    double temp_log_det = 2.0 * arma::sum(arma::log(temp_cholesky.diag()));

    // --- Derived quantities ---
    arma::mat temp_cond_mean = arma::repmat(temp_main_continuous.t(), n_, 1)
                             + 2.0 * discrete_observations_dbl_ * temp_pairwise_cross * temp_covariance;
    arma::mat D = continuous_observations_ - temp_cond_mean;

    arma::mat temp_marginal;
    if(use_marginal_pl_) {
        temp_marginal = 2.0 * temp_pairwise_discrete + 2.0 * temp_pairwise_cross * temp_covariance * temp_pairwise_cross.t();
    }

    // Initialize gradient (full dimension)
    arma::vec grad(full_dim, arma::fill::zeros);

    double logp = 0.0;

    // Full-space index offsets (fixed layout)
    const size_t kxx_offset = num_main_;
    const size_t mean_offset = num_main_ + num_pairwise_xx_;
    const size_t kxy_offset = num_main_ + num_pairwise_xx_ + q_;

    // For marginal PL: precompute helpers
    arma::mat cross_times_cov;
    arma::mat Theta_bar;
    if(use_marginal_pl_) {
        cross_times_cov = temp_pairwise_cross * temp_covariance;
        Theta_bar = arma::zeros<arma::mat>(p_, p_);
    }

    // =========================================================================
    // Part 1: OMRF conditionals (same as active-space but full indexing)
    // =========================================================================

    // Helper: flat index for Kxx(i,j) in full vector (i < j)
    auto kxx_idx = [&](size_t i, size_t j) -> size_t {
        // Row-major upper triangle: offset = sum_{r=0}^{i-1} (p-1-r) + (j-i-1)
        return kxx_offset + i * (2 * p_ - 3 - i) / 2 + (j - 1);
    };

    // Helper: flat index for Kxy(i,j) in full vector
    auto kxy_idx = [&](size_t i, size_t j) -> size_t {
        return kxy_offset + i * q_ + j;
    };

    // Precompute observed statistics for Kxx gradient (same as grad_obs_cache_)
    // Main effect observed stats
    int main_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        if(is_ordinal_variable_(s)) {
            int C_s = num_categories_(s);
            for(int c = 0; c < C_s; ++c) {
                grad(main_offset + c) = counts_per_category_(c + 1, s);
            }
            main_offset += C_s;
        } else {
            grad(main_offset)     = blume_capel_stats_(0, s);
            grad(main_offset + 1) = blume_capel_stats_(1, s);
            main_offset += 2;
        }
    }

    // Kxx observed stats
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            grad(kxx_idx(i, j)) = 4.0 * arma::dot(
                discrete_observations_dbl_.col(i),
                discrete_observations_dbl_.col(j));
        }
    }

    // OMRF conditionals loop
    int main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);

        // Rest score
        arma::vec rest;
        if(use_marginal_pl_) {
            double precision_ss = temp_marginal(s, s);
            rest = discrete_observations_dbl_ * temp_marginal.col(s)
                 - discrete_observations_dbl_.col(s) * precision_ss
                 + 2.0 * arma::dot(temp_pairwise_cross.row(s), temp_main_continuous);
        } else {
            rest = 2.0 * (discrete_observations_dbl_ * temp_pairwise_discrete.col(s)
                 - discrete_observations_dbl_.col(s) * temp_pairwise_discrete(s, s))
                 + 2.0 * continuous_observations_ * temp_pairwise_cross.row(s).t();
        }

        if(is_ordinal_variable_(s)) {
            arma::vec main_param = temp_main_discrete.row(s).cols(0, C_s - 1).t();

            if(use_marginal_pl_) {
                double precision_ss = temp_marginal(s, s);
                for(int c = 0; c < C_s; ++c) {
                    main_param(c) += static_cast<double>((c + 1) * (c + 1)) * precision_ss;
                }
            }

            arma::vec bound = main_param(C_s - 1) + static_cast<double>(C_s) * rest;
            bound = arma::max(bound, arma::zeros<arma::vec>(bound.n_elem));

            LogZAndProbs result = compute_logZ_and_probs_ordinal(
                main_param, rest, bound, C_s
            );

            logp -= arma::accu(result.log_Z);

            // Main-effect gradient
            for(int c = 0; c < C_s; ++c) {
                grad(main_effects_discrete_offset + c) -= arma::accu(result.probs.col(c + 1));
            }

            arma::vec weights = arma::regspace<arma::vec>(1, C_s);
            arma::vec E = result.probs.cols(1, C_s) * weights;

            // Pairwise discrete gradient (ALL edges, no gating)
            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(s == t) continue;
                size_t loc = (s < t) ? kxx_idx(s, t) : kxx_idx(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
            }

            if(use_marginal_pl_) {
                arma::vec weights_sq = arma::square(weights);
                arma::vec E_sq = result.probs.cols(1, C_s) * weights_sq;

                arma::vec diff_pw = discrete_observations_dbl_t_ *
                    (discrete_observations_dbl_.col(s) - E);
                diff_pw(s) = 0.0;

                double diff_diag = arma::dot(
                    discrete_observations_dbl_.col(s),
                    discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

                double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

                for(size_t t = 0; t < p_; ++t) {
                    if(t != s) Theta_bar(s, t) += diff_pw(t);
                }
                Theta_bar(s, s) += diff_diag;

                arma::rowvec cross_self = 2.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                      + 4.0 * diff_diag * cross_times_cov.row(s)
                                      + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

                for(size_t j = 0; j < q_; ++j) {
                    grad(kxy_idx(s, j)) += cross_self(j);
                }

                arma::rowvec V_s = 2.0 * cross_times_cov.row(s);
                for(size_t t = 0; t < p_; ++t) {
                    if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                    for(size_t j = 0; j < q_; ++j) {
                        grad(kxy_idx(t, j)) += diff_pw(t) * V_s(j);
                    }
                }

                for(size_t j = 0; j < q_; ++j) {
                    grad(mean_offset + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
                }
            } else {
                arma::rowvec cross_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    grad(kxy_idx(s, j)) += cross_grad_s(j);
                }
            }

            main_effects_discrete_offset += C_s;
        } else {
            // Blume-Capel variable
            int ref = baseline_category_(s);
            double lin_eff = temp_main_discrete(s, 0);
            double quad_eff = temp_main_discrete(s, 1);

            double effective_quad = quad_eff;
            if(use_marginal_pl_) {
                effective_quad += temp_marginal(s, s);
            }

            arma::vec bc_bound;
            LogZAndProbs result = compute_logZ_and_probs_blume_capel(
                rest, lin_eff, effective_quad, ref, C_s, bc_bound
            );

            logp -= arma::accu(result.log_Z);

            arma::vec score = arma::regspace<arma::vec>(0, C_s) - static_cast<double>(ref);
            arma::vec sq_score = arma::square(score);

            grad(main_effects_discrete_offset)     -= arma::accu(result.probs * score);
            grad(main_effects_discrete_offset + 1) -= arma::accu(result.probs * sq_score);

            arma::vec E = result.probs * score;

            arma::vec pw_grad = discrete_observations_dbl_t_ * E;
            for(size_t t = 0; t < p_; ++t) {
                if(s == t) continue;
                size_t loc = (s < t) ? kxx_idx(s, t) : kxx_idx(t, s);
                grad(loc) -= 2.0 * pw_grad(t);
            }

            if(use_marginal_pl_) {
                arma::vec E_sq = result.probs * sq_score;

                arma::vec diff_pw = discrete_observations_dbl_t_ *
                    (discrete_observations_dbl_.col(s) - E);
                diff_pw(s) = 0.0;

                double diff_diag = arma::dot(
                    discrete_observations_dbl_.col(s),
                    discrete_observations_dbl_.col(s)) - arma::accu(E_sq);

                double sum_obs_minus_E = arma::accu(discrete_observations_dbl_.col(s)) - arma::accu(E);

                for(size_t t = 0; t < p_; ++t) {
                    if(t != s) Theta_bar(s, t) += diff_pw(t);
                }
                Theta_bar(s, s) += diff_diag;

                arma::rowvec cross_self = 2.0 * (diff_pw.t() * temp_pairwise_cross) * temp_covariance
                                      + 4.0 * diff_diag * cross_times_cov.row(s)
                                      + 2.0 * sum_obs_minus_E * temp_main_continuous.t();

                for(size_t j = 0; j < q_; ++j) {
                    grad(kxy_idx(s, j)) += cross_self(j);
                }

                arma::rowvec V_s = 2.0 * cross_times_cov.row(s);
                for(size_t t = 0; t < p_; ++t) {
                    if(t == s || std::abs(diff_pw(t)) < 1e-300) continue;
                    for(size_t j = 0; j < q_; ++j) {
                        grad(kxy_idx(t, j)) += diff_pw(t) * V_s(j);
                    }
                }

                for(size_t j = 0; j < q_; ++j) {
                    grad(mean_offset + j) += 2.0 * temp_pairwise_cross(s, j) * sum_obs_minus_E;
                }
            } else {
                arma::rowvec cross_grad_s = 2.0 * (
                    discrete_observations_dbl_.col(s) - E
                ).t() * continuous_observations_;

                for(size_t j = 0; j < q_; ++j) {
                    grad(kxy_idx(s, j)) += cross_grad_s(j);
                }
            }

            main_effects_discrete_offset += 2;
        }
    }

    // Numerator contribution to logp
    main_effects_discrete_offset = 0;
    for(size_t s = 0; s < p_; ++s) {
        int C_s = num_categories_(s);
        arma::vec rest;
        if(use_marginal_pl_) {
            double precision_ss = temp_marginal(s, s);
            rest = discrete_observations_dbl_ * temp_marginal.col(s)
                 - discrete_observations_dbl_.col(s) * precision_ss
                 + 2.0 * arma::dot(temp_pairwise_cross.row(s), temp_main_continuous);
            logp += temp_marginal(s, s) * arma::dot(
                discrete_observations_dbl_.col(s),
                discrete_observations_dbl_.col(s));
        } else {
            rest = 2.0 * (discrete_observations_dbl_ * temp_pairwise_discrete.col(s)
                 - discrete_observations_dbl_.col(s) * temp_pairwise_discrete(s, s))
                 + 2.0 * continuous_observations_ * temp_pairwise_cross.row(s).t();
        }
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
    // Part 2: GGM conditional
    // =========================================================================

    double quad_sum = arma::accu((D * temp_precision) % D);
    logp += static_cast<double>(n_) / 2.0 *
            (-static_cast<double>(q_) * MY_LOG(2.0 * arma::datum::pi)
             + temp_log_det)
          - quad_sum / 2.0;

    arma::vec D_colsums = arma::sum(D, 0).t();
    arma::vec grad_mean_ggm = temp_precision * D_colsums;
    for(size_t j = 0; j < q_; ++j) {
        grad(mean_offset + j) += grad_mean_ggm(j);
    }

    arma::mat grad_cross_ggm = 2.0 * discrete_observations_dbl_t_ * D;
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            grad(kxy_idx(i, j)) += grad_cross_ggm(i, j);
        }
    }

    // =========================================================================
    // Part 3: Priors
    // =========================================================================

    // Main effects priors: Beta(alpha, beta) on sigmoid scale
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

    // Kxx priors: Cauchy(0, pairwise_scale_)
    const double disc_scale = pairwise_scale_;
    const double disc_scale_sq = disc_scale * disc_scale;
    for(size_t i = 0; i < p_ - 1; ++i) {
        for(size_t j = i + 1; j < p_; ++j) {
            double val = temp_pairwise_discrete(i, j);
            logp += R::dcauchy(val, 0.0, disc_scale, true);
            grad(kxx_idx(i, j)) -= 2.0 * val / (val * val + disc_scale_sq);
        }
    }

    // Continuous mean priors: Normal(0, 1)
    for(size_t j = 0; j < q_; ++j) {
        double val = temp_main_continuous(j);
        logp += R::dnorm(val, 0.0, 1.0, true);
        grad(mean_offset + j) -= val;
    }

    // Kxy priors: Cauchy(0, pairwise_scale_)
    for(size_t i = 0; i < p_; ++i) {
        for(size_t j = 0; j < q_; ++j) {
            double val = temp_pairwise_cross(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            grad(kxy_idx(i, j)) -= 2.0 * val / (val * val + pairwise_scale_ * pairwise_scale_);
        }
    }

    // =========================================================================
    // Part 4: Precision gradient via Cholesky parameterization
    // =========================================================================

    arma::mat Omega_bar = 0.5 * static_cast<double>(n_) * temp_covariance
                        - 0.5 * D.t() * D;

    Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
               * discrete_observations_dbl_t_ * D;

    if(use_marginal_pl_) {
        Omega_bar -= 2.0 * temp_covariance * temp_pairwise_cross.t()
                   * Theta_bar * temp_pairwise_cross * temp_covariance;
    }

    // Gamma(1,1) on diagonal
    for(size_t j = 0; j < q_; ++j) {
        Omega_bar(j, j) -= 1.0;
        logp -= temp_precision(j, j);
    }
    // Cauchy on off-diagonal
    const double cont_scale_sq = pairwise_scale_ * pairwise_scale_;
    for(size_t i = 0; i < q_ - 1; ++i) {
        for(size_t j = i + 1; j < q_; ++j) {
            double val = temp_precision(i, j);
            logp += R::dcauchy(val, 0.0, pairwise_scale_, true);
            double cauchy_grad = -2.0 * val / (val * val + cont_scale_sq);
            Omega_bar(i, j) += cauchy_grad;
        }
    }

    // R̄ = R (Ω̄ + Ω̄ᵀ)
    arma::mat Omega_bar_sym = Omega_bar + Omega_bar.t();
    arma::mat R_bar = temp_cholesky * Omega_bar_sym;

    // Cholesky Jacobian and gradient extraction
    size_t gidx = chol_offset;
    for(size_t j = 0; j < q_; ++j) {
        double psi_j = std::log(temp_cholesky(j, j));
        double jac_weight = static_cast<double>(q_ - j + 1);
        logp += jac_weight * psi_j;

        for(size_t i = 0; i < j; ++i) {
            grad(gidx++) = R_bar(i, j);
        }
        grad(gidx++) = R_bar(j, j) * temp_cholesky(j, j) + jac_weight;
    }
    logp += static_cast<double>(q_) * std::log(2.0);

    return {logp, grad};
}
