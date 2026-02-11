// [[Rcpp::depends(RcppParallel, RcppArmadillo, dqrng, BH)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "math/explog_switch.h"
#include "rng/rng_utils.h"
#include "utils/progress_manager.h"
#include <vector>
#include <string>

using namespace Rcpp;
using namespace RcppParallel;


// ============================================================================
//   MRF Simulation Core Functions (Thread-Safe)
// ============================================================================

/**
 * Simulate observations from an ordinal MRF using Gibbs sampling
 *
 * Thread-safe implementation that can be used both for the standalone
 * mrfSampler() function and for parallel simulation from posterior draws.
 *
 * @param num_states Number of observations to simulate
 * @param num_variables Number of variables
 * @param num_categories Number of categories per variable
 * @param pairwise Symmetric pairwise effect matrix
 * @param main Main effect parameters (variables x max_categories)
 * @param iter Number of Gibbs iterations
 * @param rng Thread-safe RNG instance
 * @return Integer matrix of simulated observations
 */
arma::imat simulate_mrf_ordinal(
    int num_states,
    int num_variables,
    const arma::ivec& num_categories,
    const arma::mat& pairwise,
    const arma::mat& main,
    int iter,
    SafeRNG& rng) {

  arma::imat observations(num_states, num_variables);
  int max_num_categories = arma::max(num_categories);
  arma::vec probabilities(max_num_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double u = 0.0;
  int score = 0;

  // Random (uniform) starting values
  for(int variable = 0; variable < num_variables; variable++) {
    for(int person = 0; person < num_states; person++) {
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < num_categories[variable]; category++) {
        cumsum += 1;
        probabilities[category + 1] = cumsum;
      }

      u = cumsum * runif(rng);

      score = 0;
      while (u > probabilities[score]) {
        score++;
      }
      observations(person, variable) = score;
    }
  }

  // Gibbs sampling iterations
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int variable = 0; variable < num_variables; variable++) {
      for(int person = 0; person < num_states; person++) {
        rest_score = 0.0;
        for(int vertex = 0; vertex < num_variables; vertex++) {
          rest_score += observations(person, vertex) *
            pairwise(vertex, variable);
        }

        cumsum = 1.0;
        probabilities[0] = 1.0;
        for(int category = 0; category < num_categories[variable]; category++) {
          exponent = main(variable, category);
          exponent += (category + 1) * rest_score;
          cumsum += MY_EXP(exponent);
          probabilities[category + 1] = cumsum;
        }

        u = cumsum * runif(rng);

        score = 0;
        while (u > probabilities[score]) {
          score++;
        }
        observations(person, variable) = score;
      }
    }
  }

  return observations;
}

/**
 * Simulate observations from an MRF with Blume-Capel/ordinal mixed variables
 *
 * Thread-safe implementation that can be used both for the standalone
 * mrfSampler() function and for parallel simulation from posterior draws.
 *
 * @param num_states Number of observations to simulate
 * @param num_variables Number of variables
 * @param num_categories Number of categories per variable
 * @param pairwise Symmetric pairwise effect matrix
 * @param main Main effect parameters (variables x max_categories)
 * @param variable_type Type of each variable ("ordinal" or "blume-capel")
 * @param baseline_category Baseline category for Blume-Capel variables
 * @param iter Number of Gibbs iterations
 * @param rng Thread-safe RNG instance
 * @return Integer matrix of simulated observations
 */
arma::imat simulate_mrf_blumecapel(
    int num_states,
    int num_variables,
    const arma::ivec& num_categories,
    const arma::mat& pairwise,
    const arma::mat& main,
    const std::vector<std::string>& variable_type,
    const arma::ivec& baseline_category,
    int iter,
    SafeRNG& rng) {

  arma::imat observations(num_states, num_variables);
  int max_num_categories = arma::max(num_categories);
  arma::vec probabilities(max_num_categories + 1);
  double exponent = 0.0;
  double rest_score = 0.0;
  double cumsum = 0.0;
  double u = 0.0;
  int score = 0;

  // Random (uniform) starting values
  for(int variable = 0; variable < num_variables; variable++) {
    for(int person = 0; person < num_states; person++) {
      cumsum = 1.0;
      probabilities[0] = 1.0;
      for(int category = 0; category < num_categories[variable]; category++) {
        cumsum += 1;
        probabilities[category + 1] = cumsum;
      }

      u = cumsum * runif(rng);

      score = 0;
      while (u > probabilities[score]) {
        score++;
      }
      observations(person, variable) = score;
    }
  }

  // Gibbs sampling iterations
  for(int iteration = 0; iteration < iter; iteration++) {
    for(int variable = 0; variable < num_variables; variable++) {
      for(int person = 0; person < num_states; person++) {
        rest_score = 0.0;
        for(int vertex = 0; vertex < num_variables; vertex++) {
          if(variable_type[vertex] != "blume-capel") {
            rest_score += observations(person, vertex) * pairwise(vertex, variable);
          } else {
            int ref = baseline_category[vertex];
            int obs = observations(person, vertex);
            rest_score += (obs - ref) * pairwise(vertex, variable);
          }
        }

        if(variable_type[variable] == "blume-capel") {
          cumsum = 0.0;
          int ref = baseline_category[variable];
          for(int category = 0; category < num_categories[variable] + 1; category++) {
            const int s = category - ref;
            // The linear term of the Blume-Capel variable
            exponent = main(variable, 0) * s;
            // The quadratic term of the Blume-Capel variable
            exponent += main(variable, 1) * s * s;
            // The pairwise effects
            exponent += rest_score * s;
            cumsum += MY_EXP(exponent);
            probabilities[category] = cumsum;
          }
        } else {
          cumsum = 1.0;
          probabilities[0] = cumsum;
          for(int category = 0; category < num_categories[variable]; category++) {
            exponent = main(variable, category);
            exponent += (category + 1) * rest_score;
            cumsum += MY_EXP(exponent);
            probabilities[category + 1] = cumsum;
          }
        }

        u = cumsum * runif(rng);

        score = 0;
        while (u > probabilities[score]) {
          score++;
        }
        observations(person, variable) = score;
      }
    }
  }

  return observations;
}


// ============================================================================
//   R Interface for mrfSampler()
// ============================================================================

// [[Rcpp::export]]
IntegerMatrix sample_omrf_gibbs(int num_states,
                                int num_variables,
                                IntegerVector num_categories,
                                NumericMatrix pairwise,
                                NumericMatrix main,
                                int iter,
                                int seed) {

  SafeRNG rng(seed);

  // Convert inputs to arma types
  arma::ivec num_categories_arma = Rcpp::as<arma::ivec>(num_categories);
  arma::mat pairwise_arma = Rcpp::as<arma::mat>(pairwise);
  arma::mat main_arma = Rcpp::as<arma::mat>(main);

  // Call the shared implementation
  arma::imat result = simulate_mrf_ordinal(
    num_states,
    num_variables,
    num_categories_arma,
    pairwise_arma,
    main_arma,
    iter,
    rng
  );

  // Check for user interrupt periodically (only in non-parallel context)
  Rcpp::checkUserInterrupt();

  return Rcpp::wrap(result);
}

// [[Rcpp::export]]
IntegerMatrix sample_bcomrf_gibbs(int num_states,
                                  int num_variables,
                                  IntegerVector num_categories,
                                  NumericMatrix pairwise,
                                  NumericMatrix main,
                                  StringVector variable_type_r,
                                  IntegerVector baseline_category,
                                  int iter,
                                  int seed) {

  SafeRNG rng(seed);

  // Convert inputs to arma/std types
  arma::ivec num_categories_arma = Rcpp::as<arma::ivec>(num_categories);
  arma::mat pairwise_arma = Rcpp::as<arma::mat>(pairwise);
  arma::mat main_arma = Rcpp::as<arma::mat>(main);
  arma::ivec baseline_category_arma = Rcpp::as<arma::ivec>(baseline_category);

  std::vector<std::string> variable_type(num_variables);
  for (int i = 0; i < num_variables; i++) {
    variable_type[i] = Rcpp::as<std::string>(variable_type_r[i]);
  }

  // Call the shared implementation
  arma::imat result = simulate_mrf_blumecapel(
    num_states,
    num_variables,
    num_categories_arma,
    pairwise_arma,
    main_arma,
    variable_type,
    baseline_category_arma,
    iter,
    rng
  );

  // Check for user interrupt (only in non-parallel context)
  Rcpp::checkUserInterrupt();

  return Rcpp::wrap(result);
}


// ============================================================================
//   Parallel Simulation for simulate.bgms() with Posterior Draws
// ============================================================================

// Structure to hold individual simulation results
struct SimulationResult {
  arma::imat observations;
  int draw_index;
  bool error;
  std::string error_msg;
};


/**
 * Worker class for parallel simulation across posterior draws
 */
class SimulationWorker : public RcppParallel::Worker {
public:
  // Input data
  const arma::mat& pairwise_samples;
  const arma::mat& main_samples;
  const arma::ivec& draw_indices;
  const int num_states;
  const int num_variables;
  const arma::ivec& num_categories;
  const std::vector<std::string>& variable_type;
  const arma::ivec& baseline_category;
  const int iter;
  const bool has_blume_capel;
  const arma::ivec& main_param_counts;

  // RNGs
  const std::vector<SafeRNG>& draw_rngs;

  // Progress
  ProgressManager& pm;

  // Output
  std::vector<SimulationResult>& results;

  SimulationWorker(
    const arma::mat& pairwise_samples,
    const arma::mat& main_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int num_variables,
    const arma::ivec& num_categories,
    const std::vector<std::string>& variable_type,
    const arma::ivec& baseline_category,
    int iter,
    bool has_blume_capel,
    const arma::ivec& main_param_counts,
    const std::vector<SafeRNG>& draw_rngs,
    ProgressManager& pm,
    std::vector<SimulationResult>& results
  ) :
    pairwise_samples(pairwise_samples),
    main_samples(main_samples),
    draw_indices(draw_indices),
    num_states(num_states),
    num_variables(num_variables),
    num_categories(num_categories),
    variable_type(variable_type),
    baseline_category(baseline_category),
    iter(iter),
    has_blume_capel(has_blume_capel),
    main_param_counts(main_param_counts),
    draw_rngs(draw_rngs),
    pm(pm),
    results(results)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      SimulationResult result;
      result.draw_index = draw_indices[i];
      result.error = false;

      try {
        // Get RNG for this draw
        SafeRNG rng = draw_rngs[i];

        // Reconstruct pairwise matrix from flat vector
        arma::mat pairwise(num_variables, num_variables, arma::fill::zeros);
        int idx = 0;
        for (int col = 0; col < num_variables; col++) {
          for (int row = col + 1; row < num_variables; row++) {
            pairwise(row, col) = pairwise_samples(draw_indices[i] - 1, idx);
            pairwise(col, row) = pairwise_samples(draw_indices[i] - 1, idx);
            idx++;
          }
        }

        // Reconstruct main effect matrix
        int max_main = arma::max(main_param_counts);
        arma::mat main(num_variables, max_main, arma::fill::zeros);
        idx = 0;
        for (int v = 0; v < num_variables; v++) {
          for (int t = 0; t < main_param_counts[v]; t++) {
            main(v, t) = main_samples(draw_indices[i] - 1, idx);
            idx++;
          }
        }

        // Run MRF simulation
        if (has_blume_capel) {
          result.observations = simulate_mrf_blumecapel(
            num_states,
            num_variables,
            num_categories,
            pairwise,
            main,
            variable_type,
            baseline_category,
            iter,
            rng
          );
        } else {
          result.observations = simulate_mrf_ordinal(
            num_states,
            num_variables,
            num_categories,
            pairwise,
            main,
            iter,
            rng
          );
        }

      } catch (const std::exception& e) {
        result.error = true;
        result.error_msg = e.what();
      } catch (...) {
        result.error = true;
        result.error_msg = "Unknown error";
      }

      results[i] = result;

      // Update progress - treating each draw as a "chain" for progress display
      pm.update(0);
    }
  }
};


/**
 * Run parallel simulations across posterior draws
 *
 * @param pairwise_samples Matrix of pairwise samples (ndraws x n_pairwise)
 * @param main_samples Matrix of main/threshold samples (ndraws x n_main)
 * @param draw_indices 1-based indices of which draws to use
 * @param num_states Number of observations to simulate per draw
 * @param num_variables Number of variables
 * @param num_categories Number of categories per variable
 * @param variable_type Type of each variable ("ordinal" or "blume-capel")
 * @param baseline_category Baseline category for each variable
 * @param iter Number of Gibbs iterations per simulation
 * @param nThreads Number of parallel threads
 * @param seed Random seed
 * @param progress_type Progress bar type (0=none, 1=total, 2=per-chain)
 *
 * @return List of simulation results (each is an integer matrix)
 */
// [[Rcpp::export]]
Rcpp::List run_simulation_parallel(
    const arma::mat& pairwise_samples,
    const arma::mat& main_samples,
    const arma::ivec& draw_indices,
    int num_states,
    int num_variables,
    const arma::ivec& num_categories,
    const Rcpp::StringVector& variable_type_r,
    const arma::ivec& baseline_category,
    int iter,
    int nThreads,
    int seed,
    int progress_type) {

  int ndraws = draw_indices.n_elem;

  // Convert variable_type to std::vector<std::string>
  std::vector<std::string> variable_type(num_variables);
  bool has_blume_capel = false;
  for (int i = 0; i < num_variables; i++) {
    variable_type[i] = Rcpp::as<std::string>(variable_type_r[i]);
    if (variable_type[i] == "blume-capel") {
      has_blume_capel = true;
    }
  }

  // Compute number of main parameters per variable
  arma::ivec main_param_counts(num_variables);
  for (int v = 0; v < num_variables; v++) {
    if (variable_type[v] == "blume-capel") {
      main_param_counts[v] = 2;  // linear and quadratic
    } else {
      main_param_counts[v] = num_categories[v];  // K-1 main for K categories
    }
  }

  // Prepare one independent RNG per draw
  std::vector<SafeRNG> draw_rngs(ndraws);
  for (int d = 0; d < ndraws; d++) {
    draw_rngs[d] = SafeRNG(seed + d);
  }

  // Prepare results storage
  std::vector<SimulationResult> results(ndraws);

  // Single-chain progress (we report across all draws as one unit)
  ProgressManager pm(1, ndraws, 0, 50, progress_type);

  // Create worker
  SimulationWorker worker(
    pairwise_samples,
    main_samples,
    draw_indices,
    num_states,
    num_variables,
    num_categories,
    variable_type,
    baseline_category,
    iter,
    has_blume_capel,
    main_param_counts,
    draw_rngs,
    pm,
    results
  );

  // Run in parallel
  {
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, nThreads);
    parallelFor(0, ndraws, worker);
  }

  pm.finish();

  // Convert results to R list
  Rcpp::List output(ndraws);
  for (int i = 0; i < ndraws; i++) {
    if (results[i].error) {
      Rcpp::stop("Error in simulation draw %d: %s",
                 results[i].draw_index, results[i].error_msg.c_str());
    }
    output[i] = Rcpp::wrap(results[i].observations);
  }

  return output;
}
