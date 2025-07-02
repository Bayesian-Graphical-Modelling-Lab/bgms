#include "Rcpp.h"
#include "explog_switch.h"

// [[Rcpp::export]]
Rcpp::String get_explog_switch() {
#if USE_CUSTOM_LOG
      return "custom";
  #else
      return "standard";
  #endif
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ieee754_exp(Rcpp::NumericVector x) {
  Rcpp::NumericVector y(x.size());
  for (int i = 0; i < x.size(); i++) {
    y[i] = MY_EXP(x[i]);
  }
  return y;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ieee754_log(Rcpp::NumericVector x) {
  Rcpp::NumericVector y(x.size());
  for (int i = 0; i < x.size(); i++) {
    y[i] = MY_LOG(x[i]);
  }
  return y;
}