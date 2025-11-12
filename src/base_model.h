#pragma once

#include <RcppArmadillo.h>
#include <stdexcept>
#include <memory>

class BaseModel {
public:
    virtual ~BaseModel() = default;

    // Capability queries
    virtual bool has_gradient() const { return false; }
    virtual bool has_adaptive_mh() const { return false; }

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

    virtual arma::vec get_vectorized_parameters() {
        throw std::runtime_error("get_vectorized_parameters method must be implemented in derived class");
    }

    virtual arma::ivec get_vectorized_indicator_parameters() {
        throw std::runtime_error("get_vectorized_indicator_parameters method must be implemented in derived class");
    }

    // Return dimensionality of the parameter space
    virtual size_t parameter_dimension() const = 0;

    virtual void set_seed(int seed) {
        throw std::runtime_error("set_seed method must be implemented in derived class");
    }

    virtual std::unique_ptr<BaseModel> clone() const {
        throw std::runtime_error("clone method must be implemented in derived class");
    }


protected:
    BaseModel() = default;
};
