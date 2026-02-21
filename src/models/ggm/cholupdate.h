#pragma once

#include <RcppArmadillo.h>

void cholesky_update(  arma::mat& R, arma::vec& u, double eps = 1e-12);
void cholesky_downdate(arma::mat& R, arma::vec& u, double eps = 1e-12);
