// [[Rcpp::depends(RcppParallel, RcppArmadillo, dqrng)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include "bgm_sampler.h"
#include <tbb/global_control.h>
#include <vector>
#include <string>
#include "rng_utils.h"
#include "progress_manager.h"
#include "mcmc_adaptation.h"

using namespace Rcpp;
using namespace RcppParallel;



/**
 * Container for the result of a single MCMC chain.
 *
 * Fields:
 *  - error: True if the chain terminated with an error, false otherwise.
 *  - error_msg: Error message in case of failure (empty if no error).
 *  - chain_id: Integer identifier for the chain (1-based).
 *  - result: Rcpp::List containing the chain’s outputs (samples, diagnostics, etc.).
 *
 * Usage:
 *  - Used in parallel samplers to collect per-chain results.
 *  - Checked after execution to propagate errors or assemble outputs into R.
 */
struct ChainResult {
  bool error;
  std::string error_msg;
  int chain_id;
  Rcpp::List result;
};



/**
 * Worker struct for running a single Gibbs sampling chain in parallel (bgm model).
 *
 * This class wraps all inputs needed to run one chain, executes
 * `run_gibbs_sampler_bgm()` inside the `operator()`, and writes results
 * into the shared output buffer.
 *
 * Fields:
 *  - inputs: model data, priors, sampler settings, adaptation controls.
 *  - chain_rngs: vector of RNG engines (one per chain).
 *  - results: reference to output buffer for storing per-chain results.
 *
 * Usage:
 *  - Constructed once per parallel run with common inputs and result buffer.
 *  - `operator()(begin, end)` executes chains in [begin, end).
 *  - Each chain’s output is written into `results` at the matching index.
 *
 * Error handling:
 *  - If an exception is thrown inside the sampler, the error flag and message
 *    are recorded in the corresponding `ChainResult`.
 */
struct GibbsChainRunner : public Worker {
  const arma::imat& observations;
  const arma::ivec& num_categories;
  double  pairwise_scale;
  const std::string& edge_prior;
  const arma::mat& inclusion_probability;
  double beta_bernoulli_alpha;
  double beta_bernoulli_beta;
  double dirichlet_alpha;
  double lambda;
  const arma::imat& interaction_index_matrix;
  int iter;
  int burnin;
  const arma::imat& counts_per_category;
  const arma::imat& blume_capel_stats;
  double main_alpha;
  double main_beta;
  bool na_impute;
  const arma::imat& missing_index;
  const arma::uvec& is_ordinal_variable;
  const arma::ivec& baseline_category;
  bool edge_selection;
  const std::string& update_method;
  const arma::imat& pairwise_effect_indices;
  double target_accept;
  const arma::imat& pairwise_stats;
  int hmc_num_leapfrogs;
  int nuts_max_depth;
  bool learn_mass_matrix;

  // Wrapped RNG engines
  const std::vector<SafeRNG>& chain_rngs;
  ProgressManager& pm;

  // output buffer
  std::vector<ChainResult>& results;

  GibbsChainRunner(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    double  pairwise_scale,
    const std::string& edge_prior,
    const arma::mat& inclusion_probability,
    double beta_bernoulli_alpha,
    double beta_bernoulli_beta,
    double dirichlet_alpha,
    double lambda,
    const arma::imat& interaction_index_matrix,
    int iter,
    int burnin,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    double main_alpha,
    double main_beta,
    bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    const arma::imat& pairwise_stats,
    int hmc_num_leapfrogs,
    int nuts_max_depth,
    bool learn_mass_matrix,
    const std::vector<SafeRNG>& chain_rngs,
    ProgressManager& pm,
    std::vector<ChainResult>& results
  ) :
    observations(observations),
    num_categories(num_categories),
     pairwise_scale( pairwise_scale),
    edge_prior(edge_prior),
    inclusion_probability(inclusion_probability),
    beta_bernoulli_alpha(beta_bernoulli_alpha),
    beta_bernoulli_beta(beta_bernoulli_beta),
    dirichlet_alpha(dirichlet_alpha),
    lambda(lambda),
    interaction_index_matrix(interaction_index_matrix),
    iter(iter),
    burnin(burnin),
    counts_per_category(counts_per_category),
    blume_capel_stats(blume_capel_stats),
    main_alpha(main_alpha),
    main_beta(main_beta),
    na_impute(na_impute),
    missing_index(missing_index),
    is_ordinal_variable(is_ordinal_variable),
    baseline_category(baseline_category),
    edge_selection(edge_selection),
    update_method(update_method),
    pairwise_effect_indices(pairwise_effect_indices),
    target_accept(target_accept),
    pairwise_stats(pairwise_stats),
    hmc_num_leapfrogs(hmc_num_leapfrogs),
    nuts_max_depth(nuts_max_depth),
    learn_mass_matrix(learn_mass_matrix),
    chain_rngs(chain_rngs),
    pm(pm),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      ChainResult out;
      out.chain_id = static_cast<int>(i + 1);
      out.error = false;

      try {
        SafeRNG rng = chain_rngs[i];

        Rcpp::List result = run_gibbs_sampler_bgm(
          out.chain_id,
          observations,
          num_categories,
           pairwise_scale,
          edge_prior,
          inclusion_probability,
          beta_bernoulli_alpha,
          beta_bernoulli_beta,
          dirichlet_alpha,
          lambda,
          interaction_index_matrix,
          iter,
          burnin,
          counts_per_category,
          blume_capel_stats,
          main_alpha,
          main_beta,
          na_impute,
          missing_index,
          is_ordinal_variable,
          baseline_category,
          edge_selection,
          update_method,
          pairwise_effect_indices,
          target_accept,
          pairwise_stats,
          hmc_num_leapfrogs,
          nuts_max_depth,
          learn_mass_matrix,
          rng,
          pm
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



/**
 * Runs multiple Gibbs sampling chains for the bgm model in parallel.
 *
 * Each chain is executed independently with its own RNG stream.
 * Results are collected into a list of per-chain outputs.
 *
 * Inputs:
 *  - observations: Matrix of categorical observations (persons × variables).
 *  - num_categories: Number of categories per variable.
 *  -  pairwise_scale: Scale parameter of the Cauchy prior on pairwise effects.
 *  - edge_prior: Prior specification for edge inclusion.
 *  - inclusion_probability: Matrix of prior inclusion probabilities.
 *  - beta_bernoulli_alpha, beta_bernoulli_beta: Hyperparameters for Bernoulli priors.
 *  - dirichlet_alpha, lambda: Hyperparameters for other prior components.
 *  - interaction_index_matrix: Indexing matrix for candidate interactions.
 *  - iter: Total number of iterations per chain.
 *  - burnin: Number of burn-in iterations.
 *  - counts_per_category: Category counts per variable.
 *  - blume_capel_stats: Sufficient statistics for Blume–Capel variables.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - na_impute: If true, impute missing data during sampling.
 *  - missing_index: Indices of missing observations.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - edge_selection: If true, update inclusion indicators during sampling.
 *  - update_method: Sampler type ("adaptive-metropolis", "hamiltonian-mc", "nuts").
 *  - pairwise_effect_indices: Indexing matrix for pairwise effects.
 *  - target_accept: Target acceptance probability for MH/HMC updates.
 *  - pairwise_stats: Sufficient statistics for pairwise effects.
 *  - hmc_num_leapfrogs: Number of leapfrog steps for HMC.
 *  - nuts_max_depth: Maximum tree depth for NUTS.
 *  - learn_mass_matrix: If true, adapt the mass matrix during warmup.
 *  - num_chains: Number of chains to run in parallel.
 *  - nThreads: Number of parallel threads.
 *  - seed: RNG seed for chain initialization.
 *
 * Returns:
 *  - An Rcpp::List of length num_chains.
 *    * Each element contains the chain’s samples and diagnostics,
 *      or an error message if the chain failed.
 *
 * Notes:
 *  - Each chain receives an independent RNG stream (`SafeRNG(seed + chain_id)`).
 *  - Errors are caught per chain and returned instead of stopping all chains.
 *  - Parallelism is controlled via Intel TBB (`nThreads` argument).
 */
// [[Rcpp::export]]
Rcpp::List run_bgm_parallel(
    const arma::imat& observations,
    const arma::ivec& num_categories,
    double  pairwise_scale,
    const std::string& edge_prior,
    const arma::mat& inclusion_probability,
    double beta_bernoulli_alpha,
    double beta_bernoulli_beta,
    double dirichlet_alpha,
    double lambda,
    const arma::imat& interaction_index_matrix,
    int iter,
    int burnin,
    const arma::imat& counts_per_category,
    const arma::imat& blume_capel_stats,
    double main_alpha,
    double main_beta,
    bool na_impute,
    const arma::imat& missing_index,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool edge_selection,
    const std::string& update_method,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    const arma::imat& pairwise_stats,
    int hmc_num_leapfrogs,
    int nuts_max_depth,
    bool learn_mass_matrix,
    int num_chains,
    int nThreads,
    uint64_t seed,
    int progress_type
) {
  std::vector<ChainResult> results(num_chains);

  // Prepare one independent RNG per chain via jump()
  std::vector<SafeRNG> chain_rngs(num_chains);
  for (int c = 0; c < num_chains; ++c) {
    chain_rngs[c] = SafeRNG(seed + c);
  }

  // only used to determine the total no. burnin iterations, a bit hacky
  WarmupSchedule warmup_schedule_temp(burnin, edge_selection, (update_method != "adaptive-metropolis"));
  int total_burnin = warmup_schedule_temp.total_burnin;
  ProgressManager pm(num_chains, iter, total_burnin, 50, progress_type);

  GibbsChainRunner worker(
      observations, num_categories,  pairwise_scale, edge_prior,
      inclusion_probability, beta_bernoulli_alpha, beta_bernoulli_beta,
      dirichlet_alpha, lambda, interaction_index_matrix, iter, burnin,
      counts_per_category, blume_capel_stats, main_alpha, main_beta,
      na_impute, missing_index, is_ordinal_variable, baseline_category,
      edge_selection, update_method, pairwise_effect_indices, target_accept,
      pairwise_stats, hmc_num_leapfrogs, nuts_max_depth, learn_mass_matrix,
      chain_rngs, pm, results
  );

  {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, num_chains, worker);
  }

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

  pm.finish();

  return output;
}
