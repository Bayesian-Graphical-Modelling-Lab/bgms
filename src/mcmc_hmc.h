#pragma once

#include <RcppArmadillo.h>
#include <functional>
#include "mcmc_utils.h"
struct SafeRNG;

SamplerResult hmc_sampler(
    const arma::vec& init_theta,
    double step_size,
    const std::function<double(const arma::vec&)>& log_post,
    const std::function<arma::vec(const arma::vec&)>& grad,
    const int num_leapfrogs,
    const arma::vec& inv_mass_diag,
    SafeRNG& rng
);