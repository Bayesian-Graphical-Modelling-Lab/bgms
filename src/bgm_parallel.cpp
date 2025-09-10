// [[Rcpp::depends(RcppParallel, RcppArmadillo, dqrng)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "bgm_sampler.h"
#include <tbb/global_control.h>
#include <vector>
#include <string>
#include "rng_utils.h"

using namespace Rcpp;
using namespace RcppParallel;

// -----------------------------------------------------------------------------
// Wrapper to silence Clang warning
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Result struct
// -----------------------------------------------------------------------------
struct ChainResult {
  bool error;
  std::string error_msg;
  int chain_id;
  Rcpp::List result;
};

// -----------------------------------------------------------------------------
// Worker struct
// -----------------------------------------------------------------------------
struct GibbsChainRunner : public Worker {
  // inputs
  const arma::imat& observations;
  const arma::ivec& num_categories;
  double interaction_scale;
  const std::string& edge_prior;
  const arma::mat& inclusion_probability;
  double beta_bernoulli_alpha;
  double beta_bernoulli_beta;
  double dirichlet_alpha;
  double lambda;
  const arma::imat& interaction_index_matrix;
  int iter;
  int burnin;
  const arma::imat& num_obs_categories;
  const arma::imat& sufficient_blume_capel;
  double threshold_alpha;
  double threshold_beta;
  bool na_impute;
  const arma::imat& missing_index;
  const arma::uvec& is_ordinal_variable;
  const arma::ivec& reference_category;
  bool edge_selection;
  const std::string& update_method;
  const arma::imat& pairwise_effect_indices;
  double target_accept;
  const arma::imat& sufficient_pairwise;
  int hmc_num_leapfrogs;
  int nuts_max_depth;
  bool learn_mass_matrix;

  // Wrapped RNG engines
  const std::vector<SafeRNG>& chain_rngs;

  // output buffer
  std::vector<ChainResult>& results;

  GibbsChainRunner(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    double interaction_scale,
    const std::string& edge_prior,
    const arma::mat& inclusion_probability,
    double beta_bernoulli_alpha,
    double beta_bernoulli_beta,
    double dirichlet_alpha,
    double lambda,
    const arma::imat& interaction_index_matrix,
    int iter,
    int burnin,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    double threshold_alpha,
    double threshold_beta,
    bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    const arma::imat& sufficient_pairwise,
    int hmc_num_leapfrogs,
    int nuts_max_depth,
    bool learn_mass_matrix,
    const std::vector<SafeRNG>& chain_rngs,
    std::vector<ChainResult>& results
  ) :
    observations(observations),
    num_categories(num_categories),
    interaction_scale(interaction_scale),
    edge_prior(edge_prior),
    inclusion_probability(inclusion_probability),
    beta_bernoulli_alpha(beta_bernoulli_alpha),
    beta_bernoulli_beta(beta_bernoulli_beta),
    dirichlet_alpha(dirichlet_alpha),
    lambda(lambda),
    interaction_index_matrix(interaction_index_matrix),
    iter(iter),
    burnin(burnin),
    num_obs_categories(num_obs_categories),
    sufficient_blume_capel(sufficient_blume_capel),
    threshold_alpha(threshold_alpha),
    threshold_beta(threshold_beta),
    na_impute(na_impute),
    missing_index(missing_index),
    is_ordinal_variable(is_ordinal_variable),
    reference_category(reference_category),
    edge_selection(edge_selection),
    update_method(update_method),
    pairwise_effect_indices(pairwise_effect_indices),
    target_accept(target_accept),
    sufficient_pairwise(sufficient_pairwise),
    hmc_num_leapfrogs(hmc_num_leapfrogs),
    nuts_max_depth(nuts_max_depth),
    learn_mass_matrix(learn_mass_matrix),
    chain_rngs(chain_rngs),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      ChainResult out;
      out.chain_id = static_cast<int>(i + 1);
      out.error = false;

      try {
        SafeRNG rng = chain_rngs[i];

        Rcpp::List result = run_gibbs_sampler_for_bgm(
          out.chain_id,
          observations,
          num_categories,
          interaction_scale,
          edge_prior,
          inclusion_probability,
          beta_bernoulli_alpha,
          beta_bernoulli_beta,
          dirichlet_alpha,
          lambda,
          interaction_index_matrix,
          iter,
          burnin,
          num_obs_categories,
          sufficient_blume_capel,
          threshold_alpha,
          threshold_beta,
          na_impute,
          missing_index,
          is_ordinal_variable,
          reference_category,
          edge_selection,
          update_method,
          pairwise_effect_indices,
          target_accept,
          sufficient_pairwise,
          hmc_num_leapfrogs,
          nuts_max_depth,
          learn_mass_matrix,
          rng
        );
        out.result = result;

      } catch (std::exception& e) {
        out.error = true;
        out.error_msg = e.what();
      } catch (...) {
        out.error = true;
        out.error_msg = "Unknown error";
      }

      results[i] = out;
    }
  }
};

// -----------------------------------------------------------------------------
// Main entry point
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List run_bgm_parallel(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    double interaction_scale,
    const std::string& edge_prior,
    const arma::mat& inclusion_probability,
    double beta_bernoulli_alpha,
    double beta_bernoulli_beta,
    double dirichlet_alpha,
    double lambda,
    const arma::imat& interaction_index_matrix,
    int iter,
    int burnin,
    const arma::imat& num_obs_categories,
    const arma::imat& sufficient_blume_capel,
    double threshold_alpha,
    double threshold_beta,
    bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& reference_category,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    const arma::imat& sufficient_pairwise,
    int hmc_num_leapfrogs,
    int nuts_max_depth,
    bool learn_mass_matrix,
    int num_chains,
    int nThreads,
    uint64_t seed   // 👈 new
) {
  std::vector<ChainResult> results(num_chains);

  // Prepare one independent RNG per chain via jump()
  std::vector<SafeRNG> chain_rngs(num_chains);
  for (int c = 0; c < num_chains; ++c) {
    chain_rngs[c] = SafeRNG(seed + c);
  }

  GibbsChainRunner worker(
      observations, num_categories, interaction_scale, edge_prior,
      inclusion_probability, beta_bernoulli_alpha, beta_bernoulli_beta,
      dirichlet_alpha, lambda, interaction_index_matrix, iter, burnin,
      num_obs_categories, sufficient_blume_capel, threshold_alpha, threshold_beta,
      na_impute, missing_index, is_ordinal_variable, reference_category,
      edge_selection, update_method, pairwise_effect_indices, target_accept,
      sufficient_pairwise, hmc_num_leapfrogs, nuts_max_depth, learn_mass_matrix,
      chain_rngs, results
  );

  {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, num_chains, worker);
  }

  // gather results
  Rcpp::List output(num_chains);
  for (int i = 0; i < num_chains; ++i) {
    if (results[i].error) {
      output[i] = Rcpp::List::create(
        Rcpp::Named("error") = results[i].error_msg,
        Rcpp::Named("chain_id") = results[i].chain_id
      );
    } else {
      output[i] = results[i].result;
    }
  }

  return output;
}
