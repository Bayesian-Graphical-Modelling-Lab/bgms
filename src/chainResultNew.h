#pragma once

#include <string>
#include <RcppArmadillo.h>

/**
 * ChainResultNew - Storage for a single MCMC chain's output
 *
 * Holds samples, diagnostics, and error state for one chain.
 * Designed for use with both MH and NUTS/HMC samplers.
 */
class ChainResultNew {

public:
    ChainResultNew() {}

    // Error handling
    bool        error = false;
    bool        userInterrupt = false;
    std::string error_msg;

    // Chain identifier
    int         chain_id = 0;

    // Parameter samples (param_dim × n_iter)
    arma::mat   samples;

    // Edge indicator samples (n_edges × n_iter), only if edge_selection = true
    arma::imat  indicator_samples;
    bool        has_indicators = false;

    // NUTS/HMC diagnostics (n_iter), only if using NUTS/HMC
    arma::ivec  treedepth_samples;
    arma::ivec  divergent_samples;
    arma::vec   energy_samples;
    bool        has_nuts_diagnostics = false;

    /**
     * Reserve storage for samples
     * @param param_dim  Number of parameters per sample
     * @param n_iter     Number of sampling iterations
     */
    void reserve(const size_t param_dim, const size_t n_iter) {
        samples.set_size(param_dim, n_iter);
    }

    /**
     * Reserve storage for edge indicator samples
     * @param n_edges  Number of edges (p * (p - 1) / 2)
     * @param n_iter   Number of sampling iterations
     */
    void reserve_indicators(const size_t n_edges, const size_t n_iter) {
        indicator_samples.set_size(n_edges, n_iter);
        has_indicators = true;
    }

    /**
     * Reserve storage for NUTS diagnostics
     * @param n_iter  Number of sampling iterations
     */
    void reserve_nuts_diagnostics(const size_t n_iter) {
        treedepth_samples.set_size(n_iter);
        divergent_samples.set_size(n_iter);
        energy_samples.set_size(n_iter);
        has_nuts_diagnostics = true;
    }

    /**
     * Store a parameter sample
     * @param iter    Iteration index (0-based)
     * @param sample  Parameter vector
     */
    void store_sample(const size_t iter, const arma::vec& sample) {
        samples.col(iter) = sample;
    }

    /**
     * Store edge indicator sample
     * @param iter        Iteration index (0-based)
     * @param indicators  Edge indicator vector
     */
    void store_indicators(const size_t iter, const arma::ivec& indicators) {
        indicator_samples.col(iter) = indicators;
    }

    /**
     * Store NUTS diagnostics for one iteration
     * @param iter       Iteration index (0-based)
     * @param tree_depth Tree depth from NUTS
     * @param divergent  Whether a divergence occurred
     * @param energy     Final Hamiltonian energy
     */
    void store_nuts_diagnostics(const size_t iter, int tree_depth, bool divergent, double energy) {
        treedepth_samples(iter) = tree_depth;
        divergent_samples(iter) = divergent ? 1 : 0;
        energy_samples(iter) = energy;
    }
};

