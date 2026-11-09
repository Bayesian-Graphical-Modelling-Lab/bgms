#include <RcppArmadillo.h>
#include "math/explog_macros.h"

// ----------------------------------------------------------------------------|
// The c++ code below is based on the R code accompanying the paper:
//  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
//  Detection With Unknown Number of Communities, Journal of the American
//  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
// ----------------------------------------------------------------------------|

// ----------------------------------------------------------------------------|
// Compute partition coefficient for the MFM - SBM
// ----------------------------------------------------------------------------|
// [[Rcpp::export]]
arma::vec compute_Vn_mfm_sbm(arma::uword num_variables,
                             double dirichlet_alpha,
                             arma::uword t_max,
                             double lambda) {
  arma::vec log_Vn(t_max);
  double r;

  // Sum over B = k + 1 components for T = t + 1 occupied clusters, so B >= T.
  // Past B = T the terms decay at least as fast as a Poisson(lambda) in B - T,
  // up to a polynomial factor, so twenty standard deviations past lambda the
  // omitted mass is below double precision. The cut-off at k = 500 only moves
  // when lambda or T is large enough to need it.
  const arma::uword k_tail = static_cast<arma::uword>(
    std::ceil(lambda + 20.0 * std::sqrt(lambda) + 20.0));

  for(arma::uword t = 0; t < t_max; t++) {
    r = -INFINITY; // initialize log-coefficient at -Inf
    const arma::uword k_max = std::max<arma::uword>(500, t + k_tail);
    for(arma::uword k = t; k <= k_max; k++){
      arma::vec b_linspace_1 = arma::linspace(k-t+1,k+1,t+1); // numerator = b*(b-1)*...*(b-|C|+1)
      arma::vec b_linspace_2 = arma::linspace((k+1)*dirichlet_alpha,(k+1)*dirichlet_alpha+num_variables-1, num_variables); // denominator b*e*(b*e+1)*...*(b*e+p-1)
      double b = arma::accu(ARMA_MY_LOG(b_linspace_1))-arma::accu(ARMA_MY_LOG(b_linspace_2)) + R::dpois((k+1)-1, lambda, true); // sum(log(numerator)) - sum(log(denominator)) + log(P=(k+1|lambda))
      double m = std::max(b,r);  // scaling factor for log-sum-exp formula
      r = MY_LOG(MY_EXP(r-m) +  MY_EXP(b-m)) + m; // update r using log-sum-exp formula to ensure numerical stability and avoid underflow
    }
    log_Vn(t) = r;
  }
  return log_Vn;
}