// test_mixed_mrf.cpp — Rcpp test helpers for MixedMRFModel skeleton (Phase A.3)
//
// These are lightweight exported functions that construct a MixedMRFModel from
// R inputs and return diagnostic information for testthat assertions.

#include <RcppArmadillo.h>
#include "models/mixed/mixed_mrf_model.h"

// [[Rcpp::export]]
Rcpp::List test_mixed_mrf_skeleton(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection,
    const std::string& pseudolikelihood,
    int seed
) {
    MixedMRFModel model(
        discrete_observations,
        continuous_observations,
        num_categories,
        is_ordinal_variable,
        baseline_category,
        inclusion_probability,
        initial_edge_indicators,
        edge_selection,
        pseudolikelihood,
        1.0, 1.0, 2.5,
        seed
    );

    // --- Dimension checks ---
    size_t param_dim = model.parameter_dimension();
    size_t full_param_dim = model.full_parameter_dimension();
    int num_variables = model.get_num_variables();
    int num_pairwise = model.get_num_pairwise();

    // --- Vectorized parameter round-trip ---
    arma::vec params_before = model.get_vectorized_parameters();
    model.set_vectorized_parameters(params_before);
    arma::vec params_after = model.get_vectorized_parameters();
    double roundtrip_max_diff = arma::max(arma::abs(params_before - params_after));

    // --- Full-parameter vector ---
    arma::vec full_params = model.get_full_vectorized_parameters();

    // --- Edge indicators ---
    arma::ivec indicator_vec = model.get_vectorized_indicator_parameters();
    const arma::imat& edge_mat = model.get_edge_indicators();

    // --- Non-trivial parameter round-trip ---
    // Set some non-zero parameters and verify round-trip
    arma::vec test_params(full_param_dim);
    for(size_t i = 0; i < full_param_dim; ++i) {
        test_params(i) = static_cast<double>(i + 1) * 0.01;
    }
    // Kyy must be SPD: overwrite the Kyy block with a diagonal-dominant SPD matrix.
    // Kyy occupies indices: num_main + num_pairwise_xx + q ... num_main + num_pairwise_xx + q + q(q+1)/2 - 1
    // For simplicity, just make it diagonal = large positive
    int p = discrete_observations.n_cols;
    int q = continuous_observations.n_cols;
    size_t num_main = 0;
    for(int s = 0; s < p; ++s) {
        if(is_ordinal_variable(s)) {
            num_main += num_categories(s);
        } else {
            num_main += 2;
        }
    }
    size_t kyy_start = num_main + p * (p - 1) / 2 + q;
    // Write identity-like SPD into Kyy block
    size_t kyy_idx = kyy_start;
    for(int i = 0; i < q; ++i) {
        for(int j = i; j < q; ++j) {
            if(i == j) {
                test_params(kyy_idx) = 2.0 + i * 0.1;  // positive diagonal
            } else {
                test_params(kyy_idx) = 0.01;  // small off-diagonal
            }
            kyy_idx++;
        }
    }

    model.set_vectorized_parameters(test_params);
    arma::vec recovered = model.get_vectorized_parameters();
    double nontrivial_max_diff = arma::max(arma::abs(test_params - recovered));

    // --- Clone round-trip ---
    auto cloned = model.clone();
    arma::vec cloned_params = cloned->get_vectorized_parameters();
    double clone_max_diff = arma::max(arma::abs(recovered - cloned_params));

    return Rcpp::List::create(
        Rcpp::Named("parameter_dimension") = param_dim,
        Rcpp::Named("full_parameter_dimension") = full_param_dim,
        Rcpp::Named("num_variables") = num_variables,
        Rcpp::Named("num_pairwise") = num_pairwise,
        Rcpp::Named("params_length") = params_before.n_elem,
        Rcpp::Named("roundtrip_max_diff") = roundtrip_max_diff,
        Rcpp::Named("nontrivial_roundtrip_max_diff") = nontrivial_max_diff,
        Rcpp::Named("clone_max_diff") = clone_max_diff,
        Rcpp::Named("indicator_length") = indicator_vec.n_elem,
        Rcpp::Named("edge_indicators_rows") = edge_mat.n_rows,
        Rcpp::Named("edge_indicators_cols") = edge_mat.n_cols,
        Rcpp::Named("has_edge_selection") = model.has_edge_selection(),
        Rcpp::Named("has_adaptive_metropolis") = model.has_adaptive_metropolis(),
        Rcpp::Named("has_gradient") = model.has_gradient(),
        Rcpp::Named("has_missing_data") = model.has_missing_data()
    );
}


// [[Rcpp::export]]
Rcpp::List test_mixed_mrf_likelihoods(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection,
    const std::string& pseudolikelihood,
    const arma::vec& params,
    int seed
) {
    MixedMRFModel model(
        discrete_observations,
        continuous_observations,
        num_categories,
        is_ordinal_variable,
        baseline_category,
        inclusion_probability,
        initial_edge_indicators,
        edge_selection,
        pseudolikelihood,
        1.0, 1.0, 2.5,
        seed
    );

    // Set parameters (caller must ensure Kyy block is SPD)
    model.set_vectorized_parameters(params);

    int p = discrete_observations.n_cols;

    // Evaluate log_conditional_omrf for each discrete variable
    Rcpp::NumericVector omrf_ll(p);
    for(int s = 0; s < p; ++s) {
        omrf_ll[s] = model.log_conditional_omrf(s);
    }

    // Evaluate log_marginal_omrf for each discrete variable (marginal mode only)
    Rcpp::NumericVector marg_omrf_ll(p);
    if(pseudolikelihood == "marginal") {
        for(int s = 0; s < p; ++s) {
            marg_omrf_ll[s] = model.log_marginal_omrf(s);
        }
    }

    // Evaluate log_conditional_ggm
    double ggm_ll = model.log_conditional_ggm();

    return Rcpp::List::create(
        Rcpp::Named("omrf_ll") = omrf_ll,
        Rcpp::Named("marg_omrf_ll") = marg_omrf_ll,
        Rcpp::Named("ggm_ll") = ggm_ll
    );
}


// [[Rcpp::export]]
Rcpp::List test_mixed_mrf_sampler(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    bool edge_selection,
    const std::string& pseudolikelihood,
    int n_warmup,
    int n_samples,
    int seed
) {
    MixedMRFModel model(
        discrete_observations,
        continuous_observations,
        num_categories,
        is_ordinal_variable,
        baseline_category,
        inclusion_probability,
        initial_edge_indicators,
        edge_selection,
        pseudolikelihood,
        1.0, 1.0, 2.5,
        seed
    );

    size_t full_dim = model.full_parameter_dimension();

    // Warmup (discard samples, no edge selection yet)
    for(int iter = 0; iter < n_warmup; ++iter) {
        model.prepare_iteration();
        model.do_one_metropolis_step(iter);
    }

    // Activate edge selection after warmup
    if(edge_selection) {
        model.set_edge_selection_active(true);
    }

    // Sampling
    arma::mat samples(n_samples, full_dim);
    int num_indicators = model.get_num_pairwise();
    arma::imat indicator_samples(n_samples, num_indicators);

    for(int iter = 0; iter < n_samples; ++iter) {
        model.prepare_iteration();
        model.do_one_metropolis_step(n_warmup + iter);
        samples.row(iter) = model.get_full_vectorized_parameters().t();
        if(edge_selection) {
            indicator_samples.row(iter) = model.get_vectorized_indicator_parameters().t();
        }
    }

    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("samples") = samples,
        Rcpp::Named("full_parameter_dimension") = full_dim
    );

    if(edge_selection) {
        result["indicator_samples"] = indicator_samples;
    }

    return result;
}


// =============================================================================
// test_mixed_mrf_cholesky — T28 (log-ratio agreement) and T29 (Cholesky fidelity)
// =============================================================================

// [[Rcpp::export]]
Rcpp::List test_mixed_mrf_cholesky(
    const arma::imat& discrete_observations,
    const arma::mat& continuous_observations,
    const arma::ivec& num_categories,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    const arma::mat& inclusion_probability,
    const arma::imat& initial_edge_indicators,
    const arma::vec& params,
    int seed,
    int target_i,
    int target_j
) {
    MixedMRFModel model(
        discrete_observations,
        continuous_observations,
        num_categories,
        is_ordinal_variable,
        baseline_category,
        inclusion_probability,
        initial_edge_indicators,
        false,       // no edge selection
        "conditional",
        1.0, 1.0, 2.5,
        seed
    );

    model.set_vectorized_parameters(params);

    int q = static_cast<int>(model.q_);

    // =========================================================================
    // T28: log_ggm_ratio_edge vs brute-force log_conditional_ggm difference
    // =========================================================================

    // Current log-likelihood
    double ggm_ll_curr = model.log_conditional_ggm();

    // Extract reparameterization constants and build a deterministic proposal
    model.get_kyy_constants(target_i, target_j);
    double phi_curr = model.kyy_constants_[0];
    double phi_prop = phi_curr + 0.05;  // small deterministic shift

    double omega_prop_ij = model.kyy_constants_[2] + model.kyy_constants_[3] * phi_prop;
    double omega_prop_jj = model.kyy_constrained_diagonal(omega_prop_ij);

    model.precision_yy_proposal_ = model.Kyy_;
    model.precision_yy_proposal_(target_i, target_j) = omega_prop_ij;
    model.precision_yy_proposal_(target_j, target_i) = omega_prop_ij;
    model.precision_yy_proposal_(target_j, target_j) = omega_prop_jj;

    // Rank-2 log-ratio via new function
    double ratio_rank2 = model.log_ggm_ratio_edge(target_i, target_j);

    // Brute-force: install proposed Kyy, full recompute, evaluate
    arma::mat Kyy_saved = model.Kyy_;
    arma::mat cov_saved = model.covariance_yy_;
    arma::mat chol_saved = model.Kyy_chol_;
    arma::mat inv_chol_saved = model.inv_cholesky_yy_;
    double logdet_saved = model.Kyy_log_det_;
    arma::mat cmean_saved = model.conditional_mean_;

    model.Kyy_ = model.precision_yy_proposal_;
    model.recompute_Kyy_decomposition();
    model.recompute_conditional_mean();
    double ggm_ll_prop = model.log_conditional_ggm();

    double ratio_bruteforce = ggm_ll_prop - ggm_ll_curr;

    // =========================================================================
    // T28b: log_ggm_ratio_diag vs brute-force
    // =========================================================================

    // Restore to current state
    model.Kyy_ = Kyy_saved;
    model.covariance_yy_ = cov_saved;
    model.Kyy_chol_ = chol_saved;
    model.inv_cholesky_yy_ = inv_chol_saved;
    model.Kyy_log_det_ = logdet_saved;
    model.conditional_mean_ = cmean_saved;

    // Diagonal proposal: small shift on log-Cholesky scale
    double logdet_curr = cholesky_helpers::get_log_det(model.Kyy_chol_);
    double logdet_sub = logdet_curr + std::log(model.covariance_yy_(target_i, target_i));
    double theta_curr = (logdet_curr - logdet_sub) / 2.0;
    double theta_prop = theta_curr + 0.05;

    model.precision_yy_proposal_ = model.Kyy_;
    model.precision_yy_proposal_(target_i, target_i) = model.Kyy_(target_i, target_i)
        - std::exp(theta_curr) * std::exp(theta_curr)
        + std::exp(theta_prop) * std::exp(theta_prop);

    double ratio_diag_rank1 = model.log_ggm_ratio_diag(target_i);

    // Brute-force
    model.Kyy_ = model.precision_yy_proposal_;
    model.recompute_Kyy_decomposition();
    model.recompute_conditional_mean();
    double ggm_ll_diag_prop = model.log_conditional_ggm();

    double ratio_diag_brute = ggm_ll_diag_prop - ggm_ll_curr;

    // =========================================================================
    // T29: Cholesky update fidelity — rank-2 edge update
    // =========================================================================

    // Restore to current state
    model.Kyy_ = Kyy_saved;
    model.covariance_yy_ = cov_saved;
    model.Kyy_chol_ = chol_saved;
    model.inv_cholesky_yy_ = inv_chol_saved;
    model.Kyy_log_det_ = logdet_saved;
    model.conditional_mean_ = cmean_saved;

    // Re-fill off-diagonal proposal
    model.precision_yy_proposal_ = model.Kyy_;
    model.precision_yy_proposal_(target_i, target_j) = omega_prop_ij;
    model.precision_yy_proposal_(target_j, target_i) = omega_prop_ij;
    model.precision_yy_proposal_(target_j, target_j) = omega_prop_jj;

    // Apply rank-1 Cholesky update
    double old_ij = model.Kyy_(target_i, target_j);
    double old_jj = model.Kyy_(target_j, target_j);
    model.Kyy_(target_i, target_j) = omega_prop_ij;
    model.Kyy_(target_j, target_i) = omega_prop_ij;
    model.Kyy_(target_j, target_j) = omega_prop_jj;
    model.cholesky_update_after_kyy_edge(old_ij, old_jj, target_i, target_j);

    arma::mat chol_rank1 = model.Kyy_chol_;
    arma::mat cov_rank1 = model.covariance_yy_;
    double logdet_rank1 = model.Kyy_log_det_;

    // Full recompute for ground truth
    model.recompute_Kyy_decomposition();
    arma::mat chol_full = model.Kyy_chol_;
    arma::mat cov_full = model.covariance_yy_;
    double logdet_full = model.Kyy_log_det_;

    double chol_max_diff = arma::max(arma::max(arma::abs(chol_rank1 - chol_full)));
    double cov_max_diff = arma::max(arma::max(arma::abs(cov_rank1 - cov_full)));
    double logdet_diff = std::abs(logdet_rank1 - logdet_full);

    // =========================================================================
    // T29b: Cholesky update fidelity — rank-1 diagonal update
    // =========================================================================

    // Restore
    model.Kyy_ = Kyy_saved;
    model.covariance_yy_ = cov_saved;
    model.Kyy_chol_ = chol_saved;
    model.inv_cholesky_yy_ = inv_chol_saved;
    model.Kyy_log_det_ = logdet_saved;

    // Diagonal proposal (reuse from T28b)
    model.precision_yy_proposal_ = model.Kyy_;
    model.precision_yy_proposal_(target_i, target_i) = model.Kyy_(target_i, target_i)
        - std::exp(theta_curr) * std::exp(theta_curr)
        + std::exp(theta_prop) * std::exp(theta_prop);

    double old_ii = model.Kyy_(target_i, target_i);
    model.Kyy_(target_i, target_i) = model.precision_yy_proposal_(target_i, target_i);
    model.cholesky_update_after_kyy_diag(old_ii, target_i);

    arma::mat chol_diag_rank1 = model.Kyy_chol_;
    arma::mat cov_diag_rank1 = model.covariance_yy_;
    double logdet_diag_rank1 = model.Kyy_log_det_;

    model.recompute_Kyy_decomposition();
    arma::mat chol_diag_full = model.Kyy_chol_;
    arma::mat cov_diag_full = model.covariance_yy_;
    double logdet_diag_full = model.Kyy_log_det_;

    double chol_diag_max_diff = arma::max(arma::max(arma::abs(chol_diag_rank1 - chol_diag_full)));
    double cov_diag_max_diff = arma::max(arma::max(arma::abs(cov_diag_rank1 - cov_diag_full)));
    double logdet_diag_diff = std::abs(logdet_diag_rank1 - logdet_diag_full);

    return Rcpp::List::create(
        // T28: off-diagonal log-ratio agreement
        Rcpp::Named("ratio_rank2") = ratio_rank2,
        Rcpp::Named("ratio_bruteforce") = ratio_bruteforce,
        // T28b: diagonal log-ratio agreement
        Rcpp::Named("ratio_diag_rank1") = ratio_diag_rank1,
        Rcpp::Named("ratio_diag_brute") = ratio_diag_brute,
        // T29: off-diagonal Cholesky fidelity
        Rcpp::Named("chol_max_diff") = chol_max_diff,
        Rcpp::Named("cov_max_diff") = cov_max_diff,
        Rcpp::Named("logdet_diff") = logdet_diff,
        // T29b: diagonal Cholesky fidelity
        Rcpp::Named("chol_diag_max_diff") = chol_diag_max_diff,
        Rcpp::Named("cov_diag_max_diff") = cov_diag_max_diff,
        Rcpp::Named("logdet_diag_diff") = logdet_diag_diff
    );
}
