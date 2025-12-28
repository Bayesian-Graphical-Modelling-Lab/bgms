// [[Rcpp::depends(RcppParallel, RcppArmadillo, dqrng)]]
#include <RcppArmadillo.h>
#include "bgmCompare_sampler.h"
#include "rng_utils.h" // must be included before RcppParallel
#include <RcppParallel.h>
#include <tbb/global_control.h>
#include <vector>
#include <string>
#include "progress_manager.h"
#include "bgmCompare_output.h"
#include "mcmc_adaptation.h"
#include "common_helpers.h"

using namespace RcppParallel;



/**
 * Container for the result of a single MCMC chain (bgmCompare model).
 *
 * Fields:
 *  - error: True if the chain terminated with an error, false otherwise.
 *  - error_msg: Error message if an error occurred (empty if none).
 *  - chain_id: Integer identifier for the chain (1-based).
 *  - result: bgmCompareOutput object containing chain results
 *    (samples, diagnostics, metadata).
 *
 * Usage:
 *  - Used in parallel execution to collect results from each chain.
 *  - Checked after execution to propagate errors or assemble outputs
 *    into an R-accessible list.
 */
struct bgmCompareChainResult {
  bool error;
  std::string error_msg;
  int chain_id;
  bgmCompareOutput result;
};



/**
 * Parallel worker for running a single Gibbs sampling chain (bgmCompare model).
 *
 * This struct wraps all inputs needed for one chain and provides an
 * `operator()` so that multiple chains can be launched in parallel with TBB.
 *
 * Workflow per chain:
 *  - Construct a chain-specific RNG from `chain_rngs`.
 *  - Copy master statistics and observation data into per-chain buffers.
 *  - Call `run_gibbs_sampler_bgmCompare()` to execute the full chain.
 *  - Catch and record any errors (sets `error = true` and stores `error_msg`).
 *  - Store results into the shared `results` vector at the chain index.
 *
 * Inputs (stored as const references or values):
 *  - observations_master: Input observation matrix (persons × variables).
 *  - num_groups: Number of groups.
 *  - counts_per_category_master: Group-level category counts.
 *  - blume_capel_stats_master: Group-level Blume–Capel sufficient statistics.
 *  - pairwise_stats_master: Group-level pairwise sufficient statistics.
 *  - num_categories: Number of categories per variable.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - pairwise_scale: Scale for Cauchy prior on baseline pairwise effects.
 *  - difference_scale: Scale for Cauchy prior on group differences.
 *  - difference_selection_alpha, difference_selection_beta: Hyperparameters for difference-selection prior.
 *  - difference_prior: Choice of prior distribution for group differences.
 *  - iter, warmup: Iteration counts.
 *  - na_impute: If true, perform missing data imputation.
 *  - missing_data_indices: Indices of missing observations.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - difference_selection: If true, perform difference selection updates.
 *  - main_effect_indices: Index ranges for main effects.
 *  - pairwise_effect_indices: Index ranges for pairwise effects.
 *  - target_accept: Target acceptance rate for adaptive methods.
 *  - nuts_max_depth: Maximum tree depth for NUTS.
 *  - learn_mass_matrix: If true, adapt mass matrix during warmup.
 *  - projection: Group projection matrix.
 *  - group_membership: Group assignment for each observation.
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - interaction_index_matrix: Lookup table of variable pairs.
 *  - inclusion_probability_master: Prior inclusion probabilities for pairwise effects.
 *  - chain_rngs: Pre-initialized RNG engines (one per chain).
 *  - update_method: Sampler type ("adaptive-metropolis", "hamiltonian-mc", "nuts").
 *  - hmc_num_leapfrogs: Number of leapfrog steps (HMC).
 *
 * Output:
 *  - results: Vector of `bgmCompareChainResult` objects, one per chain, filled in place.
 *
 * Notes:
 *  - Each worker instance is shared across threads but invoked with different
 *    [begin,end) ranges, corresponding to chain indices.
 *  - Per-chain copies of statistics and observations prevent cross-thread mutation.
 *  - Errors are caught locally so one failing chain does not crash the entire run.
 */
struct GibbsCompareChainRunner : public Worker {
  const arma::imat& observations_master;
  const int num_groups;
  const std::vector<arma::imat>& counts_per_category_master;
  const std::vector<arma::imat>& blume_capel_stats_master;
  const std::vector<arma::mat>&  pairwise_stats_master;
  const arma::ivec& num_categories;
  const double main_alpha;
  const double main_beta;
  const double pairwise_scale;
  const double difference_scale;
  const double difference_selection_alpha;
  const double difference_selection_beta;
  const std::string& difference_prior;
  const int iter;
  const int warmup;
  const bool na_impute;
  const arma::imat& missing_data_indices;
  const arma::uvec& is_ordinal_variable;
  const arma::ivec& baseline_category;
  const bool difference_selection;
  const arma::imat& main_effect_indices;
  const arma::imat& pairwise_effect_indices;
  const double target_accept;
  const int nuts_max_depth;
  const bool learn_mass_matrix;
  const arma::mat& projection;
  const arma::ivec& group_membership;
  const arma::imat& group_indices;
  const arma::imat& interaction_index_matrix;
  const arma::mat& inclusion_probability_master;
  // RNG seeds
  const std::vector<SafeRNG>& chain_rngs;
  const UpdateMethod update_method;
  const int hmc_num_leapfrogs;
  ProgressManager& pm;
  // output
  std::vector<bgmCompareChainResult>& results;

  GibbsCompareChainRunner(
    const arma::imat& observations_master,
    int num_groups,
    const std::vector<arma::imat>& counts_per_category_master,
    const std::vector<arma::imat>& blume_capel_stats_master,
    const std::vector<arma::mat>&  pairwise_stats_master,
    const arma::ivec& num_categories,
    double main_alpha,
    double main_beta,
    double pairwise_scale,
    double difference_scale,
    double difference_selection_alpha,
    double difference_selection_beta,
    const std::string& difference_prior,
    int iter,
    int warmup,
    bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    int nuts_max_depth,
    bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    const arma::mat& inclusion_probability_master,
    const std::vector<SafeRNG>& chain_rngs,
    const UpdateMethod update_method,
    const int hmc_num_leapfrogs,
    ProgressManager& pm,
    std::vector<bgmCompareChainResult>& results
  ) :
    observations_master(observations_master),
    num_groups(num_groups),
    counts_per_category_master(counts_per_category_master),
    blume_capel_stats_master(blume_capel_stats_master),
    pairwise_stats_master(pairwise_stats_master),
    num_categories(num_categories),
    main_alpha(main_alpha),
    main_beta(main_beta),
    pairwise_scale(pairwise_scale),
    difference_scale(difference_scale),
    difference_selection_alpha(difference_selection_alpha),
    difference_selection_beta(difference_selection_beta),
    difference_prior(difference_prior),
    iter(iter),
    warmup(warmup),
    na_impute(na_impute),
    missing_data_indices(missing_data_indices),
    is_ordinal_variable(is_ordinal_variable),
    baseline_category(baseline_category),
    difference_selection(difference_selection),
    main_effect_indices(main_effect_indices),
    pairwise_effect_indices(pairwise_effect_indices),
    target_accept(target_accept),
    nuts_max_depth(nuts_max_depth),
    learn_mass_matrix(learn_mass_matrix),
    projection(projection),
    group_membership(group_membership),
    group_indices(group_indices),
    interaction_index_matrix(interaction_index_matrix),
    inclusion_probability_master(inclusion_probability_master),
    chain_rngs(chain_rngs),
    update_method(update_method),
    hmc_num_leapfrogs(hmc_num_leapfrogs),
    pm(pm),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      bgmCompareChainResult out;
      out.chain_id = static_cast<int>(i + 1);
      out.error = false;

      try {
        // per-chain RNG
        SafeRNG rng = chain_rngs[i];

        // make per-chain copies
        std::vector<arma::imat> counts_per_category = counts_per_category_master;
        std::vector<arma::imat> blume_capel_stats = blume_capel_stats_master;
        std::vector<arma::mat>  pairwise_stats = pairwise_stats_master;
        arma::mat inclusion_probability = inclusion_probability_master;
        arma::imat observations = observations_master;

        // run sampler (pure C++)
        bgmCompareOutput result = run_gibbs_sampler_bgmCompare(
          out.chain_id,
          observations,
          num_groups,
          counts_per_category,
          blume_capel_stats,
          pairwise_stats,
          num_categories,
          main_alpha,
          main_beta,
          pairwise_scale,
          difference_scale,
          difference_selection_alpha,
          difference_selection_beta,
          difference_prior,
          iter,
          warmup,
          na_impute,
          missing_data_indices,
          is_ordinal_variable,
          baseline_category,
          difference_selection,
          main_effect_indices,
          pairwise_effect_indices,
          target_accept,
          nuts_max_depth,
          learn_mass_matrix,
          projection,
          group_membership,
          group_indices,
          interaction_index_matrix,
          inclusion_probability,
          rng,
          update_method,
          hmc_num_leapfrogs,
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
 * Runs multiple parallel Gibbs sampling chains for the bgmCompare model.
 *
 * This function is the main entry point from R into the C++ backend for bgmCompare.
 * It launches `num_chains` independent chains in parallel using TBB,
 * each managed by `GibbsCompareChainRunner`.
 *
 * Workflow:
 *  - Initialize a per-chain RNG from the global seed.
 *  - Construct a `GibbsCompareChainRunner` worker with all shared inputs.
 *  - Launch the worker across chains using `parallelFor`.
 *  - Collect results from all chains into an Rcpp::List.
 *
 * Inputs:
 *  - observations: Observation matrix (persons × variables).
 *  - num_groups: Number of groups.
 *  - counts_per_category: Group-level category counts (for ordinal variables).
 *  - blume_capel_stats: Group-level sufficient statistics (for Blume–Capel variables).
 *  - pairwise_stats: Group-level pairwise sufficient statistics.
 *  - num_categories: Number of categories per variable.
 *  - main_alpha, main_beta: Hyperparameters for Beta priors on main effects.
 *  - pairwise_scale: Scale for Cauchy prior on baseline pairwise effects.
 *  - difference_scale: Scale for Cauchy prior on group differences.
 *  - difference_selection_alpha, difference_selection_beta: Hyperparameters for difference-selection prior.
 *  - difference_prior: Choice of prior distribution for group differences.
 *  - iter: Number of post-warmup iterations to draw.
 *  - warmup: Number of warmup iterations.
 *  - na_impute: If true, perform missing data imputation during sampling.
 *  - missing_data_indices: Indices of missing entries in observations.
 *  - is_ordinal_variable: Indicator (1 = ordinal, 0 = Blume–Capel).
 *  - baseline_category: Reference categories for Blume–Capel variables.
 *  - difference_selection: If true, perform difference selection updates.
 *  - main_effect_indices: Index ranges [row_start,row_end] for each variable.
 *  - pairwise_effect_indices: Lookup table mapping (var1,var2) → row in pairwise_effects.
 *  - target_accept: Target acceptance rate for adaptive samplers.
 *  - nuts_max_depth: Maximum tree depth for NUTS.
 *  - learn_mass_matrix: If true, adapt the mass matrix during warmup.
 *  - projection: Group projection matrix (num_groups × (num_groups − 1)).
 *  - group_membership: Group assignment for each observation.
 *  - group_indices: Row ranges [start,end] for each group in observations.
 *  - interaction_index_matrix: Lookup table of variable pairs.
 *  - inclusion_probability: Prior inclusion probabilities for pairwise effects.
 *  - num_chains: Number of chains to run.
 *  - nThreads: Maximum number of threads for parallel execution.
 *  - seed: Base random seed (incremented per chain).
 *  - update_method: Sampler type ("adaptive-metropolis", "hamiltonian-mc", "nuts").
 *  - hmc_num_leapfrogs: Number of leapfrog steps for HMC.
 *
 * Returns:
 *  - Rcpp::List of length `num_chains`, where each element is either:
 *    * An error record (fields: "error", "chain_id"), if the chain failed.
 *    * A result list containing:
 *        - "main_samples": Posterior samples of main effects.
 *        - "pairwise_samples": Posterior samples of pairwise effects.
 *        - "treedepth__", "divergent__", "energy__": NUTS diagnostics.
 *        - "indicator_samples": Inclusion indicators (if selection was enabled).
 *        - "chain_id": Identifier of the chain.
 *
 * Notes:
 *  - Parallel execution is controlled by TBB; `nThreads` limits concurrency.
 *  - Each chain gets its own RNG stream, initialized as `seed + chain_id`.
 *  - This function is called by the exported R function `bgmCompare()`.
 */
// [[Rcpp::export]]
Rcpp::List run_bgmCompare_parallel(
    const arma::imat& observations,
    int num_groups,
    const std::vector<arma::imat>& counts_per_category,
    const std::vector<arma::imat>& blume_capel_stats,
    const std::vector<arma::mat>&  pairwise_stats,
    const arma::ivec& num_categories,
    double main_alpha,
    double main_beta,
    double pairwise_scale,
    double difference_scale,
    double difference_selection_alpha,
    double difference_selection_beta,
    const std::string& difference_prior,
    int iter,
    int warmup,
    bool na_impute,
    const arma::imat& missing_data_indices,
    const arma::uvec& is_ordinal_variable,
    const arma::ivec& baseline_category,
    bool difference_selection,
    const arma::imat& main_effect_indices,
    const arma::imat& pairwise_effect_indices,
    double target_accept,
    int nuts_max_depth,
    bool learn_mass_matrix,
    const arma::mat& projection,
    const arma::ivec& group_membership,
    const arma::imat& group_indices,
    const arma::imat& interaction_index_matrix,
    const arma::mat& inclusion_probability,
    int num_chains,
    int nThreads,
    int seed,
    const std::string& update_method,
    int hmc_num_leapfrogs,
    int progress_type
) {
  std::vector<bgmCompareChainResult> results(num_chains);

  // per-chain seeds
  std::vector<SafeRNG> chain_rngs(num_chains);
  for (int c = 0; c < num_chains; ++c) {
    chain_rngs[c] = SafeRNG(seed + c);
  }

  UpdateMethod update_method_enum = update_method_from_string(update_method);

  // only used to determine the total no. warmup iterations, a bit hacky
  WarmupSchedule warmup_schedule_temp(warmup, difference_selection, (update_method_enum != adaptive_metropolis));
  int total_warmup = warmup_schedule_temp.total_warmup;
  ProgressManager pm(num_chains, iter, total_warmup, 50, progress_type);

  GibbsCompareChainRunner worker(
      observations, num_groups,
      counts_per_category, blume_capel_stats, pairwise_stats,
      num_categories, main_alpha, main_beta, pairwise_scale, difference_scale,
      difference_selection_alpha, difference_selection_beta, difference_prior,
      iter, warmup, na_impute, missing_data_indices, is_ordinal_variable,
      baseline_category, difference_selection, main_effect_indices,
      pairwise_effect_indices, target_accept, nuts_max_depth, learn_mass_matrix,
      projection, group_membership, group_indices, interaction_index_matrix,
      inclusion_probability, chain_rngs, update_method_enum, hmc_num_leapfrogs,
      pm, results
  );

  {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, num_chains, worker);
  }

  // wrap results back into Rcpp::List
  Rcpp::List output(num_chains);
  for (int i = 0; i < num_chains; ++i) {
    if (results[i].error) {
      output[i] = Rcpp::List::create(
        Rcpp::Named("error") = results[i].error_msg,
        Rcpp::Named("chain_id") = results[i].chain_id
      );
    } else {
      const auto& r = results[i].result;
      Rcpp::List chain_out = Rcpp::List::create(
        Rcpp::Named("main_samples") = r.main_samples,
        Rcpp::Named("pairwise_samples") = r.pairwise_samples,
        Rcpp::Named("treedepth__") = r.treedepth_samples,
        Rcpp::Named("divergent__") = r.divergent_samples,
        Rcpp::Named("energy__") = r.energy_samples,
        Rcpp::Named("chain_id") = r.chain_id
      );
      if (r.has_indicator) {
        chain_out["indicator_samples"] = r.indicator_samples;
      }
      chain_out["userInterrupt"] = r.userInterrupt;
      output[i] = chain_out;
    }
  }

  pm.finish();

  return output;
}