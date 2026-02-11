#pragma once

#include <RcppArmadillo.h>
#include <stdexcept>
#include <memory>

// Forward declarations
struct SamplerResult;
struct SafeRNG;

class BaseModel {
public:
    virtual ~BaseModel() = default;

    // Capability queries
    virtual bool has_gradient() const { return false; }
    virtual bool has_adaptive_mh() const { return false; }
    virtual bool has_nuts() const { return has_gradient(); }
    virtual bool has_edge_selection() const { return false; }

    // Core methods (to be overridden by derived classes)
    virtual double logp(const arma::vec& parameters) = 0;

    virtual arma::vec gradient(const arma::vec& parameters) {
        if (!has_gradient()) {
            throw std::runtime_error("Gradient not implemented for this model");
        }
        throw std::runtime_error("Gradient method must be implemented in derived class");
    }

    virtual std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) {
        if (!has_gradient()) {
            throw std::runtime_error("Gradient not implemented for this model");
        }
        return {logp(parameters), gradient(parameters)};
    }

    // For Metropolis-Hastings (model handles parameter groups internally)
    virtual void do_one_mh_step() {
        throw std::runtime_error("do_one_mh_step method must be implemented in derived class");
    }

    // Edge selection (for models with spike-and-slab priors)
    virtual void update_edge_indicators() {
        throw std::runtime_error("update_edge_indicators not implemented for this model");
    }

    virtual arma::vec get_vectorized_parameters() const {
        throw std::runtime_error("get_vectorized_parameters method must be implemented in derived class");
    }

    virtual void set_vectorized_parameters(const arma::vec& parameters) {
        throw std::runtime_error("set_vectorized_parameters method must be implemented in derived class");
    }

    virtual arma::ivec get_vectorized_indicator_parameters() {
        throw std::runtime_error("get_vectorized_indicator_parameters method must be implemented in derived class");
    }

    // Full parameter dimension (for fixed-size output, includes all possible params)
    virtual size_t full_parameter_dimension() const {
        return parameter_dimension();  // Default: same as active dimension
    }

    // Get full vectorized parameters (zeros for inactive, for consistent output)
    virtual arma::vec get_full_vectorized_parameters() const {
        throw std::runtime_error("get_full_vectorized_parameters must be implemented in derived class");
    }

    // Return dimensionality of the active parameter space
    virtual size_t parameter_dimension() const = 0;

    virtual void set_seed(int seed) {
        throw std::runtime_error("set_seed method must be implemented in derived class");
    }

    virtual std::unique_ptr<BaseModel> clone() const {
        throw std::runtime_error("clone method must be implemented in derived class");
    }

    // RNG access for samplers
    virtual SafeRNG& get_rng() {
        throw std::runtime_error("get_rng method must be implemented in derived class");
    }

    // Step size for gradient-based samplers
    virtual void set_step_size(double step_size) { step_size_ = step_size; }
    virtual double get_step_size() const { return step_size_; }

    // Inverse mass matrix for HMC/NUTS
    virtual void set_inv_mass(const arma::vec& inv_mass) { inv_mass_ = inv_mass; }
    virtual const arma::vec& get_inv_mass() const { return inv_mass_; }

    // Get active inverse mass (for models with edge selection, may be subset)
    virtual arma::vec get_active_inv_mass() const { return inv_mass_; }

    // Edge selection activation
    virtual void set_edge_selection_active(bool active) {
        (void)active;  // Default: no-op
    }

    // Initialize graph structure for edge selection
    virtual void initialize_graph() {
        // Default: no-op
    }

    // Missing data imputation
    virtual bool has_missing_data() const { return false; }
    virtual void impute_missing() {
        // Default: no-op
    }

    // Edge prior support: models expose these so the external edge prior
    // class can read current indicators and update inclusion probabilities.
    virtual const arma::imat& get_edge_indicators() const {
        throw std::runtime_error("get_edge_indicators not implemented for this model");
    }
    virtual arma::mat& get_inclusion_probability() {
        throw std::runtime_error("get_inclusion_probability not implemented for this model");
    }
    virtual int get_num_variables() const {
        throw std::runtime_error("get_num_variables not implemented for this model");
    }
    virtual int get_num_pairwise() const {
        throw std::runtime_error("get_num_pairwise not implemented for this model");
    }

protected:
    BaseModel() = default;
    double step_size_ = 0.1;
    arma::vec inv_mass_;
};
