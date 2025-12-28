#include "math/explog_switch.h"
#include "Rcpp.h"

Rcpp::String get_explog_switch() {
#if USE_CUSTOM_LOG
      return "custom";
  #else
      return "standard";
  #endif
}

Rcpp::NumericVector rcpp_ieee754_exp(Rcpp::NumericVector x) {
  Rcpp::NumericVector y(x.size());
  for (int i = 0; i < x.size(); i++) {
    y[i] = MY_EXP(x[i]);
  }
  return y;
}

Rcpp::NumericVector rcpp_ieee754_log(Rcpp::NumericVector x) {
  Rcpp::NumericVector y(x.size());
  for (int i = 0; i < x.size(); i++) {
    y[i] = MY_LOG(x[i]);
  }
  return y;
}