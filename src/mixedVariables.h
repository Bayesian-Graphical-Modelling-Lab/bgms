#pragma once

#include <RcppArmadillo.h>
#include <stdexcept>
#include <memory>
#include "base_model.h"


// Forward declaration - this class is work in progress and not yet functional
class MixedVariableTypes : public BaseModel {
public:

    //  Constructor
    MixedVariableTypes(
        Rcpp::List input_from_R,
        const arma::mat& inclusion_probability,
        const arma::imat& initial_edge_indicators,
        const bool edge_selection = true
    )
    {

        instantiate_variable_types(input_from_R);
        // instantiate_variable_interactions();  // TODO: not yet implemented

        dim_ = 0;
        for (const auto& var_type : variable_types_) {
            dim_ += var_type->parameter_dimension();
        }
        // dim_ += interactions_.size();  // TODO: proper dimension calculation
    }

    // Capability queries
    bool has_gradient() const override
    {
        for (const auto& var_type : variable_types_) {
            if (var_type->has_gradient()) {
                return true;
            }
        }
        return false;
    }
    bool has_adaptive_mh() const override
    {
        for (const auto& var_type : variable_types_) {
            if (var_type->has_adaptive_mh()) {
                return true;
            }
        }
        return false;
    }

    // Return dimensionality of the parameter space
    size_t parameter_dimension() const override {
        return dim_;
    }

    arma::vec get_vectorized_parameters() const override {
        arma::vec result(dim_);
        size_t current = 0;
        for (size_t i = 0; i < variable_types_.size(); ++i) {
            arma::vec var_params = variable_types_[i]->get_vectorized_parameters();
            result.subvec(current, current + var_params.n_elem - 1) = var_params;
            current += var_params.n_elem;
        }
        for (size_t i = 0; i < interactions_.size(); ++i) {
            const arma::mat& interactions_mat = interactions_[i];
            for (size_t c = 0; c < interactions_mat.n_cols; ++c) {
                    for (size_t r = 0; r < interactions_mat.n_rows; ++r) {
                    result(current) = interactions_mat(r, c);
                    ++current;
                }
            }
        }
        return result;
    }

    arma::ivec get_vectorized_indicator_parameters() override {
        for (size_t i = 0; i < variable_types_.size(); ++i) {
            auto& [from, to] = indicator_parameters_indices_[i];
            vectorized_indicator_parameters_.subvec(from, to) = variable_types_[i]->get_vectorized_indicator_parameters();
        }
        size_t current = indicator_parameters_indices_.empty() ? 0 : indicator_parameters_indices_.back().second + 1;
        for (size_t i = 0; i < interactions_indicators_.size(); ++i) {
            const arma::imat& indicator_mat = interactions_indicators_[i];
            for (size_t c = 0; c < indicator_mat.n_cols; ++c) {
                    for (size_t r = 0; r < indicator_mat.n_rows; ++r) {
                    vectorized_indicator_parameters_(current) = indicator_mat(r, c);
                    ++current;
                }
            }
        }

        return vectorized_indicator_parameters_;
    }


    double logp(const arma::vec& parameters) override
    {
        double total_logp = 0.0;
        for (size_t i = 0; i < variable_types_.size(); ++i) {
            auto& [from, to] = parameters_indices_[i];
            // need to do some transformation here!
            arma::vec var_params = parameters.subvec(from, to);
            total_logp += variable_types_[i]->logp(var_params);
        }
        // interactions log-probability can be added here if needed
        return total_logp;
    }

    arma::vec gradient(const arma::vec& parameters) override {

        // TODO: only should call the gradient for variable types that have it
        // the rest are assumed to be constant, so have gradient zero
        arma::vec total_gradient = arma::zeros<arma::vec>(parameters.n_elem);
        for (size_t i = 0; i < variable_types_.size(); ++i)
        {
            if (!variable_types_[i]->has_gradient()) {
                continue;
            }
            auto& [from, to] = parameters_indices_[i];
            arma::vec var_params = parameters.subvec(from, to);
            // maybe need to do some transformation here!
            arma::vec var_gradient = variable_types_[i]->gradient(var_params);
            total_gradient.subvec(from, to) = var_gradient;
        }

        return total_gradient;
    }

    std::pair<double, arma::vec> logp_and_gradient(
        const arma::vec& parameters) override {
        if (!has_gradient()) {
            throw std::runtime_error("Gradient not implemented for this model");
        }
        return {logp(parameters), gradient(parameters)};
    }

    void do_one_mh_step() override {
        for (auto& var_type : variable_types_) {
            var_type->do_one_mh_step();
        }
    }

    void set_seed(int seed) override {
        for (auto& var_type : variable_types_) {
            var_type->set_seed(seed);
        }
    }

    std::unique_ptr<BaseModel> clone() const override {
        throw std::runtime_error("clone method not yet implemented for MixedVariableTypes");
    }


private:
    std::vector<std::unique_ptr<BaseModel>> variable_types_;
    std::vector<arma::mat> interactions_;
    std::vector<arma::imat> interactions_indicators_;
    size_t dim_;
    arma::vec vectorized_parameters_;
    arma::ivec vectorized_indicator_parameters_;
    arma::ivec indices_from_;
    arma::ivec indices_to_;
    std::vector<std::pair<size_t, size_t>> parameters_indices_;
    std::vector<std::pair<size_t, size_t>> indicator_parameters_indices_;

    void instantiate_variable_types(const Rcpp::List input_from_R);

};
