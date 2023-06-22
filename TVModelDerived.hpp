
#ifndef TIME_VARYING_MODEL_DERIVED
#define TIME_VARYING_MODEL_DERIVED

// Include header files
#include "TVModelBase.hpp"
#include "QuadraturePoints.hpp"
#include "ToolLikelihood.hpp"
#include "Parameters.hpp"

// Include libraries
#include <iostream>
#include <functional>

// Forward declaration of the Base class
//class TVModel::ModelBase;

namespace TVModel{
using T = TypeTraits;

class PowerParameterModel final: public TVModel::ModelBase{
public:
    // Constructor
    PowerParameterModel() = default;
    PowerParameterModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Method for executing the log-likelihood
    T::VariableType compute_loglikelihood() override;



private:
    // Complex data structure
    Params::Parameters parameters;                                          // Class for the parameters

    // Simple data structure
    T::NumberType n_ranges_parameters = 4;                                  // Number of minimum or maximum ranges that have to be provided for this method
    T::NumberType n_parameters;                                             // Number of parameter of the model
    T::VectorNumberType all_n_parameters;
    T::IdNameType name_method = "PowerParameter";                           // Name of this method
    T::VectorXdr & v_parameters = parameters.get_v_parameters();

    // Functions for implementing the likelihood
    // It collects the result
    // std::function<T::VariableType()> ll_pp;
    std::function<T::VariableType(T::VectorXdr)> ll_pp; // T::VariableType, T::IndexType, 
    std::function<T::VariableType(T::VectorXdr, T::SharedPtrType)> ll_group_pp;     // T::VariableType, T::IndexType, 


    // Quadrature nodes and weights
    QuadraturePoints::Points9 points9;
    T::NumberType n_nodes = 9;
    std::array<T::VariableType, 9>& nodes = points9.nodes;
    std::array<T::VariableType, 9>& weights = points9.weights;

    // Other tools
    ToolsLikelihood::Tools tool;
    T::VariableType factor_c = tool.factor_c_pp;

    // Virtual method for computing the number of parameters
    T::NumberType compute_n_parameters() override;

    // Method for building the log-likelihood
    void build_loglikelihood();

    // Virtual method to derive the interval variance of the frailty
    // void compute_sd_frailty() override;

};

}
#endif // TIME_VARYING_MODEL_DERIVED


