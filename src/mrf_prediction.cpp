// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils/variable_helpers.h"
using namespace Rcpp;


// ============================================================================
//   GGM Conditional Prediction
// ============================================================================

/**
 * Compute conditional Gaussian parameters for a GGM.
 *
 * For a GGM with precision matrix Omega, the conditional distribution of
 * X_j given X_{-j} = x_{-j} is:
 *
 *   X_j | X_{-j} ~ N( -omega_{jj}^{-1} sum_{k != j} omega_{jk} x_k,
 *                       omega_{jj}^{-1} )
 *
 * @param observations  n x p matrix of observed continuous data
 * @param predict_vars  0-based indices of variables to predict
 * @param precision     p x p precision matrix (Omega)
 *
 * @return List of n x 2 matrices (one per predicted variable), where
 *         column 0 = conditional mean, column 1 = conditional SD.
 */
// [[Rcpp::export]]
Rcpp::List compute_conditional_ggm(
    const arma::mat& observations,
    const arma::ivec& predict_vars,
    const arma::mat& precision
) {
  int n = observations.n_rows;
  int num_predict_vars = predict_vars.n_elem;

  Rcpp::List result(num_predict_vars);

  for (int pv = 0; pv < num_predict_vars; pv++) {
    int j = predict_vars[pv];

    double omega_jj = precision(j, j);
    double cond_var = 1.0 / omega_jj;
    double cond_sd  = std::sqrt(cond_var);

    // Compute conditional means for all observations at once:
    // cond_mean_i = -(1/omega_jj) * sum_{k != j} omega_{jk} * x_{ik}
    //             = -(1/omega_jj) * (observations * omega_{.,j} - x_{ij} * omega_{jj})
    arma::vec omega_col_j = precision.col(j);
    arma::vec linear_pred = observations * omega_col_j;   // n x 1
    // Subtract the self-contribution x_{ij} * omega_{jj}
    linear_pred -= observations.col(j) * omega_jj;
    arma::vec cond_means = -cond_var * linear_pred;

    // Build n x 2 output matrix: [mean, sd]
    Rcpp::NumericMatrix out(n, 2);
    for (int i = 0; i < n; i++) {
      out(i, 0) = cond_means[i];
      out(i, 1) = cond_sd;
    }

    result[pv] = out;
  }

  return result;
}


// ============================================================================
//   OMRF Conditional Prediction
// ============================================================================

// Compute conditional probabilities P(X_j = c | X_{-j}) for ordinal MRF
// Uses numerically stable vectorized computation from variable_helpers.h

// [[Rcpp::export]]
Rcpp::List compute_conditional_probs(
    arma::imat observations,          // n x p matrix of observed data
    arma::ivec predict_vars,          // which variables to predict (0-based)
    arma::mat pairwise,           // p x p pairwise matrix
    arma::mat main,             // p x max_cat main effect matrix
    arma::ivec num_categories,         // number of categories per variable
    Rcpp::StringVector variable_type, // "ordinal" or "blume-capel" per variable
    arma::ivec baseline_category      // baseline for blume-capel variables
) {
  int no_persons = observations.n_rows;
  int num_variables = observations.n_cols;
  int num_predict_vars = predict_vars.n_elem;

  // Output is a list of probability matrices, one per predict variable
  Rcpp::List result(num_predict_vars);

  for(int pv = 0; pv < num_predict_vars; pv++) {
    int variable = predict_vars[pv];
    int n_cats = num_categories[variable] + 1;  // Include category 0

    // Compute rest scores for all persons at once (vectorized)
    arma::vec rest_scores(no_persons, arma::fill::zeros);

    for(int vertex = 0; vertex < num_variables; vertex++) {
      if(vertex == variable) continue;  // Skip the variable we're predicting

      arma::vec obs_col = arma::conv_to<arma::vec>::from(observations.col(vertex));

      if(std::string(variable_type[vertex]) != "blume-capel") {
        rest_scores += obs_col * pairwise(vertex, variable);
      } else {
        int ref = baseline_category[vertex];
        rest_scores += (obs_col - double(ref)) * pairwise(vertex, variable);
      }
    }

    // Use numerically stable probability computation
    arma::mat probs;

    if(std::string(variable_type[variable]) == "blume-capel") {
      int ref = baseline_category[variable];
      double lin_eff = main(variable, 0);
      double quad_eff = main(variable, 1);
      arma::vec bound;  // Will be computed inside

      probs = compute_probs_blume_capel(
        rest_scores,
        lin_eff,
        quad_eff,
        ref,
        num_categories[variable],
        bound
      );
    } else {
      // Regular ordinal variable
      // Extract main effect parameters for this variable
      arma::vec main_param = main.row(variable).head(num_categories[variable]).t();

      // Compute bounds for numerical stability: max exponent per person
      arma::vec bound(no_persons, arma::fill::zeros);
      for(int c = 0; c < num_categories[variable]; c++) {
        arma::vec exps = main_param[c] + (c + 1) * rest_scores;
        bound = arma::max(bound, exps);
      }

      probs = compute_probs_ordinal(
        main_param,
        rest_scores,
        bound,
        num_categories[variable]
      );
    }

    // Convert arma::mat to Rcpp::NumericMatrix
    Rcpp::NumericMatrix prob_mat(no_persons, n_cats);
    for(int i = 0; i < no_persons; i++) {
      for(int c = 0; c < n_cats; c++) {
        prob_mat(i, c) = probs(i, c);
      }
    }

    result[pv] = prob_mat;
  }

  return result;
}
