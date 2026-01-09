// mcmc_leapfrog.h
#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc/mcmc_memoization.h"



/**
 * Function: leapfrog
 *
 * Performs num_leapfrogs standard leapfrog integration steps using the given gradient function.
 */
std::pair<arma::vec, arma::vec> leapfrog(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag
);



/**
 * Function: leapfrog_memo
 *
 * Performs a leapfrog step using a memoization wrapper to avoid redundant gradient evaluations.
 */
std::pair<arma::vec, arma::vec> leapfrog_memo(
    const arma::vec& theta,
    const arma::vec& r,
    double eps,
    Memoizer& memo,
    const arma::vec& inv_mass_diag
);