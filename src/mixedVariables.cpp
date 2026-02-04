#include "ggm_model.h"
#include "mixedVariables.h"
#include "rng/rng_utils.h"


void MixedVariableTypes::instantiate_variable_types(Rcpp::List input_from_R)
{
    // instantiate variable_types_
    for (Rcpp::List var_type_list : input_from_R) {
        std::string type = Rcpp::as<std::string>(var_type_list["type"]);
        if (type == "Continuous") {
            variable_types_.push_back(std::make_unique<GaussianVariables>(
                Rcpp::as<arma::mat>(var_type_list["observations"]),
                Rcpp::as<arma::mat>(var_type_list["inclusion_probability"]),
                Rcpp::as<arma::imat>(var_type_list["initial_edge_indicators"]),
                Rcpp::as<bool>(var_type_list["edge_selection"])
            ));
        // } else if (type == "Ordinal") {
        //     variable_types_.push_back(std::make_unique<OrdinalVariables>(
        //         var_type_list["observations"],
        //         var_type_list["inclusion_probability"],
        //         var_type_list["initial_edge_indicators"],
        //         var_type_list["edge_selection"]
        //     ));
        // } else if (type == "Blume-Capel") {
        //     variable_types_.push_back(std::make_unique<BlumeCapelVariables>(
        //         var_type_list["observations"],
        //         var_type_list["inclusion_probability"],
        //         var_type_list["initial_edge_indicators"],
        //         var_type_list["edge_selection"]
        //     ));
        // } else if (type == "Count") {
        //     variable_types_.push_back(std::make_unique<CountVariables>(
        //         var_type_list["observations"],
        //         var_type_list["inclusion_probability"],
        //         var_type_list["initial_edge_indicators"],
        //         var_type_list["edge_selection"]
        //     ));
        } else {
            throw std::runtime_error("MixedVariableTypes received an unknown variable type in sublist fro input_from_R: " + type);
        }
    }
}