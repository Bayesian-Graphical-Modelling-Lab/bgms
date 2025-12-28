#ifndef BGMS_EXPLOG_SWITCH_H
#define BGMS_EXPLOG_SWITCH_H

#include "RcppArmadillo.h"

double __ieee754_exp(double x); // forward declaration
double __ieee754_log(double x); // forward declaration

// elementwise exp
template<typename T1>
arma::Mat<double> custom_arma_exp(const arma::Base<double, T1>& X)
{
  arma::Mat<double> Xin = X.get_ref();
  arma::Mat<double> out(Xin.n_rows, Xin.n_cols, arma::fill::none);

  const double* in_mem  = Xin.memptr();
  double* out_mem       = out.memptr();
  const arma::uword N   = Xin.n_elem;

  for (arma::uword i = 0; i < N; ++i)
    out_mem[i] = __ieee754_exp(in_mem[i]);

  return out;
}

// elementwise log
template<typename T1>
arma::Mat<double> custom_arma_log(const arma::Base<double, T1>& X)
{
    arma::Mat<double> Xin = X.get_ref();
    arma::Mat<double> out(Xin.n_rows, Xin.n_cols, arma::fill::none);

    double* out_mem       = out.memptr();
    const double* in_mem  = Xin.memptr();
    const arma::uword N   = Xin.n_elem;

    for (arma::uword i = 0; i < N; ++i)
        out_mem[i] = __ieee754_log(in_mem[i]);

    return out;
}


#endif
