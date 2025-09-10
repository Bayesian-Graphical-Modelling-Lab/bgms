#pragma once

#include <RcppArmadillo.h>
#include <dqrng.h>
#include <dqrng_generator.h>
#include <dqrng_distribution.h>
#include <xoshiro.h>
#include <random>
#include <boost/random/beta_distribution.hpp>

struct SafeRNG {
  dqrng::xoshiro256plusplus eng;

  // Default constructor // TODO: perhaps delete this to require a seed
  SafeRNG() : eng(1) {}

  SafeRNG(const int seed) : eng(seed) {}

  // Required interface for std::distributions
  using result_type = uint64_t;
  static constexpr result_type min() { return dqrng::xoshiro256plusplus::min(); }
  static constexpr result_type max() { return dqrng::xoshiro256plusplus::max(); }
  result_type operator()() { return eng(); }
};


// ============================================================
// Scalar RNG helpers
// ============================================================

// Uniform(0,1)
inline double runif(SafeRNG& rng) {
  return dqrng::uniform_distribution(0.0, 1.0)(rng.eng);
}

// Normal(mu, sigma)
inline double rnorm(SafeRNG& rng, double mu = 0.0, double sigma = 1.0) {
  return dqrng::normal_distribution(mu, sigma)(rng.eng);
}

// Bernoulli(p)
inline int rbern(SafeRNG& rng, double p) {
  return (runif(rng) < p) ? 1 : 0;
}

// Beta(a, b)
inline double rbeta(SafeRNG& rng, double a, double b) {
  return boost::random::beta_distribution<double>(a, b)(rng.eng);
}

// Exponential(lambda)
inline double rexp(SafeRNG& rng, double lambda) {
  return dqrng::exponential_distribution(lambda)(rng.eng);
}

// ============================================================
// Armadillo RNG helpers
// ============================================================

// arma::vec of N(mu, sigma)
inline arma::vec arma_rnorm_vec(SafeRNG& rng,
                                arma::uword n,
                                double mu = 0.0, double sigma = 1.0) {
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i)
    out[i] = rnorm(rng, mu, sigma);
  return out;
}

// arma::mat of N(mu, sigma)
inline arma::mat arma_rnorm_mat(SafeRNG& rng,
                                arma::uword nrow, arma::uword ncol,
                                double mu = 0.0, double sigma = 1.0) {
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = rnorm(rng, mu, sigma);
    }
  }
  return out;
}

// arma::vec of U(0,1)
inline arma::vec arma_runif_vec(SafeRNG& rng, arma::uword n) {
  arma::vec out(n);
  for (arma::uword i = 0; i < n; ++i)
    out[i] = runif(rng);
  return out;
}

// arma::mat of U(0,1)
inline arma::mat arma_runif_mat(SafeRNG& rng,
                                arma::uword nrow, arma::uword ncol) {
  arma::mat out(nrow, ncol);
  for (arma::uword j = 0; j < ncol; ++j) {
    for (arma::uword i = 0; i < nrow; ++i) {
      out(i, j) = runif(rng);
    }
  }
  return out;
}

// Random permutation (like arma::randperm)
inline arma::uvec arma_randperm(SafeRNG& rng, arma::uword n) {
  arma::uvec out(n);
  std::iota(out.begin(), out.end(), 0);
  std::shuffle(out.begin(), out.end(), rng.eng);
  return out;
}
