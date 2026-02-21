#pragma once

#include <string>

/**
 * SamplerConfig - Configuration for MCMC sampling
 *
 * Holds all settings for the generic MCMC runner, including:
 * - Sampler type selection (MH, NUTS, HMC)
 * - Iteration counts
 * - NUTS/HMC specific parameters
 * - Edge selection settings
 */
struct SamplerConfig {
    // Sampler type: "adaptive_metropolis", "nuts", "hmc"
    std::string sampler_type = "adaptive_metropolis";

    // Iteration counts
    int no_iter = 1000;
    int no_warmup = 500;

    // NUTS/HMC parameters
    int max_tree_depth = 10;
    int num_leapfrogs = 10;  // For HMC only
    double initial_step_size = 0.1;
    double target_acceptance = 0.8;

    // Edge selection settings
    bool edge_selection = false;

    // Missing data imputation
    bool na_impute = false;

    // Random seed
    int seed = 42;

    // Constructor with defaults
    SamplerConfig() = default;
};
