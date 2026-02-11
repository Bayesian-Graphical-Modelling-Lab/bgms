// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils/variable_helpers.h"
using namespace Rcpp;

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
