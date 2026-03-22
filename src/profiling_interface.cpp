#include "mcmc/profiler.h"
#include "Rcpp.h"

// [[Rcpp::export]]
void bgms_profiler_enable() {
    auto& prof = RattleProfiler::instance();
    prof.reset();
    prof.enabled = true;
}

// [[Rcpp::export]]
void bgms_profiler_disable() {
    RattleProfiler::instance().enabled = false;
}

// [[Rcpp::export]]
void bgms_profiler_reset() {
    RattleProfiler::instance().reset();
}

// [[Rcpp::export]]
Rcpp::List bgms_profiler_results() {
    auto& prof = RattleProfiler::instance();

    int n = static_cast<int>(prof.records.size());
    Rcpp::CharacterVector names(n);
    Rcpp::NumericVector total_us(n);
    Rcpp::IntegerVector count(n);
    Rcpp::NumericVector mean_us(n);
    Rcpp::NumericVector min_us(n);
    Rcpp::NumericVector max_us(n);

    int i = 0;
    for(const auto& kv : prof.records) {
        names[i]    = kv.first;
        total_us[i] = kv.second.total_ns / 1000.0;
        count[i]    = static_cast<int>(kv.second.count);
        mean_us[i]  = kv.second.count > 0
            ? (kv.second.total_ns / 1000.0) / kv.second.count
            : 0.0;
        min_us[i] = kv.second.min_ns == INT64_MAX
            ? 0.0
            : kv.second.min_ns / 1000.0;
        max_us[i] = kv.second.max_ns / 1000.0;
        i++;
    }

    Rcpp::DataFrame timings = Rcpp::DataFrame::create(
        Rcpp::Named("name")     = names,
        Rcpp::Named("total_us") = total_us,
        Rcpp::Named("count")    = count,
        Rcpp::Named("mean_us")  = mean_us,
        Rcpp::Named("min_us")   = min_us,
        Rcpp::Named("max_us")   = max_us,
        Rcpp::Named("stringsAsFactors") = false
    );

    return Rcpp::List::create(
        Rcpp::Named("timings")                    = timings,
        Rcpp::Named("nuts_steps")                 = prof.total_nuts_steps,
        Rcpp::Named("total_tree_depth")           = prof.total_tree_depth,
        Rcpp::Named("total_leapfrogs")            = prof.total_leapfrogs,
        Rcpp::Named("total_constrained_leapfrogs") = prof.total_constrained_leapfrogs,
        Rcpp::Named("total_pcg_iterations")       = prof.total_pcg_iterations,
        Rcpp::Named("total_pcg_calls")            = prof.total_pcg_calls,
        Rcpp::Named("total_pcg_constraints")      = prof.total_pcg_constraints
    );
}
