#pragma once

#include <chrono>
#include <cstdint>
#include <map>
#include <string>

// ------------------------------------------------------------------
// TimingRecord
// ------------------------------------------------------------------
// Accumulates wall-clock timing statistics for a named profiling probe.

struct TimingRecord {
    int64_t total_ns = 0;
    int64_t count    = 0;
    int64_t min_ns   = INT64_MAX;
    int64_t max_ns   = 0;
};

// ------------------------------------------------------------------
// RattleProfiler
// ------------------------------------------------------------------
// Lightweight singleton profiler for RATTLE leapfrog and NUTS.
//
// When enabled is false the BGMS_PROF_* macros reduce to a single
// branch on a bool (no clock reads, no map lookups).  When enabled
// is true each probe reads the high-resolution clock and accumulates
// min / max / total / count statistics.
//
// NOT thread-safe.  Use with a single MCMC chain.

class RattleProfiler {
public:
    static RattleProfiler& instance() {
        static RattleProfiler prof;
        return prof;
    }

    bool enabled = false;

    std::map<std::string, TimingRecord> records;

    int64_t total_nuts_steps            = 0;
    int64_t total_tree_depth            = 0;
    int64_t total_leapfrogs             = 0;
    int64_t total_constrained_leapfrogs = 0;
    int64_t total_pcg_iterations        = 0;
    int64_t total_pcg_calls             = 0;
    int64_t total_pcg_constraints       = 0;

    void reset() {
        records.clear();
        total_nuts_steps            = 0;
        total_tree_depth            = 0;
        total_leapfrogs             = 0;
        total_constrained_leapfrogs = 0;
        total_pcg_iterations        = 0;
        total_pcg_calls             = 0;
        total_pcg_constraints       = 0;
    }

    void record(const char* name, int64_t ns) {
        auto& rec = records[name];
        rec.total_ns += ns;
        rec.count++;
        if(ns < rec.min_ns) rec.min_ns = ns;
        if(ns > rec.max_ns) rec.max_ns = ns;
    }

private:
    RattleProfiler() = default;
    RattleProfiler(const RattleProfiler&) = delete;
    RattleProfiler& operator=(const RattleProfiler&) = delete;
};

// ------------------------------------------------------------------
// Profiling macros
// ------------------------------------------------------------------
// BGMS_PROF_START(var)          — capture a time point.
// BGMS_PROF_RECORD(name, var)  — record elapsed ns since var
//                                 (no-op when profiling is disabled).

#define BGMS_PROF_START(var) \
    auto var = std::chrono::high_resolution_clock::now()

#define BGMS_PROF_RECORD(name, var) \
    do { \
        auto& _bgms_prof_ = RattleProfiler::instance(); \
        if(_bgms_prof_.enabled) { \
            _bgms_prof_.record(name, \
                std::chrono::duration_cast<std::chrono::nanoseconds>( \
                    std::chrono::high_resolution_clock::now() - (var)).count()); \
        } \
    } while(0)
