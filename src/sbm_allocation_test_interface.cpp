// -----------------------------------------------------------------------------
// sbm_allocation_test_interface.cpp
//
// R-facing test entry for the plain (hierarchical) MFM-SBM update. Runs the
// production StochasticBlockEdgePrior::update() for n_iter sweeps, so tests
// can compare the sampled partitions with exact enumeration over all
// partitions of a small graph.
//
//  - resample_graph = false: the graph is held fixed, so the chain targets
//    p(z | G) exactly.
//  - resample_graph = true: each sweep is followed by an exact draw of
//    G | z, theta ~ Bernoulli(theta_{z_i z_j}), so the chain is a prior-only
//    Gibbs sampler on the joint p(G, z, theta) of the hierarchical
//    specification.
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>

#include "priors/edge_prior.h"
#include "rng/rng_utils.h"

// [[Rcpp::export]]
Rcpp::List test_sbm_allocation_chain(
    const arma::imat& edge_indicators,
    int n_iter,
    double alpha,
    double beta,
    double alpha_between,
    double beta_between,
    double dirichlet_alpha,
    double lambda,
    int seed,
    bool resample_graph = false
) {
  const int p = static_cast<int>(edge_indicators.n_rows);
  const int num_pairs = p * (p - 1) / 2;

  SafeRNG rng(seed);
  StochasticBlockEdgePrior prior(
    alpha, beta, alpha_between, beta_between, dirichlet_alpha, lambda
  );

  arma::imat graph = edge_indicators;
  arma::mat inclusion_probability(p, p, arma::fill::zeros);
  arma::imat allocations(n_iter, p);
  arma::imat graphs;
  if (resample_graph) graphs.set_size(n_iter, num_pairs);

  for (int s = 0; s < n_iter; s++) {
    prior.update(graph, inclusion_probability, p, num_pairs, rng);

    if (resample_graph) {
      int e = 0;
      for (int i = 0; i < p - 1; i++) {
        for (int j = i + 1; j < p; j++) {
          int g = runif(rng) < inclusion_probability(i, j) ? 1 : 0;
          graph(i, j) = g;
          graph(j, i) = g;
          graphs(s, e++) = g;
        }
      }
    }

    allocations.row(s) = prior.get_allocations().t();
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("allocations") = allocations
  );
  if (resample_graph) out["graphs"] = graphs;
  return out;
}
