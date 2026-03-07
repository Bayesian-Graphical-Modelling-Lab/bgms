#pragma once

#include <string>

/**
 * SamplerConfig - Configuration for MCMC sampling
 *
 * Holds all settings for the generic MCMC runner, including:
 * - Sampler type selection (Metropolis, NUTS, HMC)
 * - Iteration counts
 * - NUTS/HMC specific parameters
 * - Edge selection settings
 */
struct SamplerConfig {
    /// Sampler type: "adaptive_metropolis", "nuts", or "hmc".
    std::string sampler_type = "adaptive_metropolis";

    /// Number of post-warmup iterations.
    int no_iter = 1000;
    /// Number of warmup iterations.
    int no_warmup = 500;

    /// Maximum NUTS tree depth.
    int max_tree_depth = 10;
    /// Number of leapfrog steps (HMC only).
    int num_leapfrogs = 10;
    /// Initial step size for gradient-based samplers.
    double initial_step_size = 0.1;
    /// Target acceptance rate for dual-averaging adaptation.
    double target_acceptance = 0.8;

    /// Enable spike-and-slab edge selection.
    bool edge_selection = false;

    /// Enable missing-data imputation during sampling.
    bool na_impute = false;

    /// Random seed.
    int seed = 42;

    /// Default constructor.
    SamplerConfig() = default;
};
