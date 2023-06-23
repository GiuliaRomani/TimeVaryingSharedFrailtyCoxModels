
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


namespace TVModel{
using T = TypeTraits;

// Class for implementing PowerParameter Model
class PowerParameterModel final: public TVModel::ModelBase{
public:
    // Constructor
    PowerParameterModel() = default;
    PowerParameterModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Method for executing the log-likelihood
    void optimize_loglikelihood() override;

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

    // Method for extracting the parameters from the vector
    T::TuplePPType extract_parameters(T::VectorXdr v_parameters_) override;

    // Method for building the log-likelihood
    void build_loglikelihood() override;
};

//--------------------------------------------------------------------------------------------------
/*
// Class for implementing PAIK Model
class PaikModel final: public TVModel::ModelBase{
public:
    // Constructor
    PaikModel() = default;
    PaikModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Method for executing the log-likelihood
    void optimize_loglikelihood() override;

private:
    // Complex data structure
    Params::Parameters parameters;                                          // Class for the parameters

    // Simple data structure
    T::NumberType n_ranges_parameters = 5;                                  // Number of minimum or maximum ranges that have to be provided for this method
    T::NumberType n_parameters;                                             // Number of parameter of the model
    T::VectorNumberType all_n_parameters;                        
    T::VectorXdr & v_parameters = parameters.get_v_parameters();

    T::IdNameType name_method = "Paik";                                     // Name of this method

    // Functions for implementing the likelihood
    // std::function<T::VariableType()> ll_pp;
    std::function<T::VariableType(T::VectorXdr)> ll_paik; // T::VariableType, T::IndexType, 
    std::function<T::VariableType(T::VectorXdr, T::SharedPtrType)> ll_group_paik;     // T::VariableType, T::IndexType, 

    // Other tools
    ToolsLikelihood::Tools tool;
    T::VariableType factor_c = tool.factor_c_paik;

    // Virtual method for computing the number of parameters
    T::NumberType compute_n_parameters() override;

    // Method for building the log-likelihood
    void build_loglikelihood() override;

    // Method for building the matrix A_ijk
    T::VariableType compute_A_i();
    T::VectorXdr compute_A_ik();
    T::VectorXdr compute_d_ik();

};
*/

//-----------------------------------------------------------------------------------------

}
#endif // TIME_VARYING_MODEL_DERIVED


