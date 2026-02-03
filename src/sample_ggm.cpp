#include <vector>
#include <memory>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/global_control.h>

#include "ggm_model.h"
#include "utils/progress_manager.h"
#include "chainResultNew.h"

void run_mcmc_sampler_single_thread(
    ChainResultNew& chain_result,
    BaseModel& model,
    const int no_iter,
    const int no_warmup,
    const int chain_id,
    ProgressManager& pm
) {

    chain_result.chain_id = chain_id + 1;
    size_t i = 0;
    for (size_t iter = 0; iter < no_iter + no_warmup; ++iter) {

        model.do_one_mh_step();

        // update hyperparameters (BetaBinomial, SBM, etc.)

        if (iter >= no_warmup) {

            chain_result.store_sample(i, model.get_vectorized_parameters());
            ++i;
        }

        pm.update(chain_id);
        if (pm.shouldExit()) {
            chain_result.userInterrupt = true;
            break;
        }
    }
}

struct GGMChainRunner : public RcppParallel::Worker {
    std::vector<ChainResultNew>& results_;
    std::vector<std::unique_ptr<BaseModel>>& models_;
    size_t no_iter_;
    size_t no_warmup_;
    int seed_;
    ProgressManager& pm_;

  GGMChainRunner(
      std::vector<ChainResultNew>& results,
      std::vector<std::unique_ptr<BaseModel>>& models,
      const size_t no_iter,
      const size_t no_warmup,
      const int seed,
      ProgressManager& pm
  ) :
    results_(results),
    models_(models),
    no_iter_(no_iter),
    no_warmup_(no_warmup),
    seed_(seed),
    pm_(pm)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {

        ChainResultNew& chain_result = results_[i];
        BaseModel& model = *models_[i];
        model.set_seed(seed_ + i);
      try {

        run_mcmc_sampler_single_thread(chain_result, model, no_iter_, no_warmup_, i, pm_);

      } catch (std::exception& e) {
        chain_result.error = true;
        chain_result.error_msg = e.what();
      } catch (...) {
        chain_result.error = true;
        chain_result.error_msg = "Unknown error";
      }
    }
  }
};

void run_mcmc_sampler_threaded(
    std::vector<ChainResultNew>& results,
    std::vector<std::unique_ptr<BaseModel>>& models,
    const int no_iter,
    const int no_warmup,
    const int seed,
    const int no_threads,
    ProgressManager& pm
) {

    GGMChainRunner runner(results, models, no_iter, no_warmup, seed, pm);
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, no_threads);
    RcppParallel::parallelFor(0, results.size(), runner);
}


std::vector<ChainResultNew> run_mcmc_sampler(
    BaseModel& model,
    const int no_iter,
    const int no_warmup,
    const int no_chains,
    const int seed,
    const int no_threads,
    ProgressManager& pm
) {

    Rcpp::Rcout << "Allocating results objects..." << std::endl;
    std::vector<ChainResultNew> results(no_chains);
    for (size_t c = 0; c < no_chains; ++c) {
        results[c].reserve(model.parameter_dimension(), no_iter);
    }

    if (no_threads > 1) {

        Rcpp::Rcout << "Running multi-threaded MCMC sampling..." << std::endl;
        std::vector<std::unique_ptr<BaseModel>> models;
        models.reserve(no_chains);
        for (size_t c = 0; c < no_chains; ++c) {
            models.push_back(model.clone());  // deep copy via virtual clone
        }
        run_mcmc_sampler_threaded(results, models, no_iter, no_warmup, seed, no_threads, pm);

    } else {

        model.set_seed(seed);
        Rcpp::Rcout << "Running single-threaded MCMC sampling..." << std::endl;
        // TODO: this is actually not correct, each chain should have its own model object
        // now chain 2 continues from chain 1 state
        for (size_t c = 0; c < no_chains; ++c) {
            run_mcmc_sampler_single_thread(results[c], model, no_iter, no_warmup, c, pm);
        }

    }
    return results;
}

Rcpp::List convert_sampler_output_to_ggm_result(const std::vector<ChainResultNew>& results) {

    Rcpp::List output(results.size());
    for (size_t i = 0; i < results.size(); ++i) {

        Rcpp::List chain_i;
        chain_i["chain_id"] = results[i].chain_id;
        if (results[i].error) {
            chain_i["error"] = results[i].error_msg;
        } else {
            chain_i["samples"] = results[i].samples;
            chain_i["userInterrupt"] = results[i].userInterrupt;

        }
        output[i] = chain_i;
    }
    return output;
}

// [[Rcpp::export]]
Rcpp::List sample_ggm(
    const Rcpp::List& inputFromR,
    const arma::mat& prior_inclusion_prob,
    const arma::imat& initial_edge_indicators,
    const int no_iter,
    const int no_warmup,
    const int no_chains,
    const bool edge_selection,
    const int seed,
    const int no_threads,
    const int progress_type
) {

    // should be done dynamically
    // also adaptation method should be specified differently
    // GaussianVariables model(X, prior_inclusion_prob, initial_edge_indicators, edge_selection);
    GaussianVariables model = createGaussianVariablesFromR(inputFromR, prior_inclusion_prob, initial_edge_indicators, edge_selection);

    ProgressManager pm(no_chains, no_iter, no_warmup, 50, progress_type);

    std::vector<ChainResultNew> output = run_mcmc_sampler(model, no_iter, no_warmup, no_chains, seed, no_threads, pm);

    Rcpp::List ggm_result = convert_sampler_output_to_ggm_result(output);

    pm.finish();

    return ggm_result;
}