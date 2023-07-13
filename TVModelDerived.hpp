
#ifndef TIME_VARYING_MODEL_DERIVED
#define TIME_VARYING_MODEL_DERIVED

// Include header files
#include "TVModelBase.hpp"
#include "QuadraturePoints.hpp"
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
    void evaluate_loglikelihood() override;

private:
    // Complex data structure
    Params::Parameters parameters;                                          // Class for the parameters

    // Simple data structure
    T::NumberType n_ranges_parameters = 4;                                  // Number of minimum or maximum ranges that have to be provided for this method
    T::NumberType n_parameters;                                             // Number of parameter of the model
    T::VectorNumberType all_n_parameters;
    T::VectorXdr & v_parameters = parameters.get_v_parameters();

    T::IdNameType name_method = "PowerParameter";                           // Name of this method

    // Functions for implementing the likelihood
    std::function<T::VariableType(T::VectorXdr&)> ll_pp;  
    std::function<T::VariableType(T::VectorXdr&, T::SharedPtrType)> ll_group_pp;     
    
    // Function for computing the second derivative wrt one dimension
    std::function<T::VariableType(T::IndexType, T::VectorXdr&)> dd_ll_pp;	


    // Quadrature nodes and weights
    QuadraturePoints::Points9 points9;
    T::NumberType n_nodes = 9;
    std::array<T::VariableType, 9>& nodes = points9.nodes;
    std::array<T::VariableType, 9>& weights = points9.weights;


    // Virtual method for computing the number of parameters
    void compute_n_parameters() override;

    // Method for extracting the parameters from the vector
    T::TuplePPType extract_parameters(const T::VectorXdr&) const ;
    
    // Method for initializing the hessian diagonal and standard error
    void compute_hessian_diagonal(T::VectorXdr&) override;
    void compute_se(T::VectorXdr&) override;

    // Method for building the log-likelihood
    void build_loglikelihood() override;
    void build_dd_loglikelihood() override;
};


//--------------------------------------------------------------------------------------------------
// Class for implementing PAIK Model
class PaikModel final: public TVModel::ModelBase{
public:
    // Constructor
    PaikModel() = default;
    PaikModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Method for executing the log-likelihood
    void evaluate_loglikelihood() override;

private:
    // Complex data structure
    Params::Parameters parameters;                                          // Class for the parameters

    // Simple data structure
    T::NumberType n_ranges_parameters = 5;                                  // Number of minimum or maximum ranges that have to be provided for this method
    T::NumberType n_parameters;                                             // Number of parameter of the model
    T::VectorNumberType all_n_parameters;                        
    T::VectorXdr& v_parameters = parameters.get_v_parameters();

    T::IdNameType name_method = "Paik";                                     // Name of this method

    // Functions for implementing the likelihood
    std::function<T::VariableType(T::VectorXdr&)> ll_paik; 
    std::function<T::VariableType(T::VectorXdr&, T::SharedPtrType)> ll_group_paik;     
    
    // Function for computing the second derivative wrt one dimension
    std::function<T::VariableType(T::IndexType, T::VectorXdr&)> dd_ll_paik;	

    // Virtual method for computing the number of parameters
    void compute_n_parameters() override;

    void compute_hessian_diagonal(T::VectorXdr&);
    void compute_se(T::VectorXdr&);
    
    // Method for extracting the parameters from the vector
    T::TuplePaikType extract_parameters(T::VectorXdr&);
    T::TupleMatrixAType extract_matrixA_variables(T::SharedPtrType, T::VectorXdr&, T::VectorXdr&);
    T::TupleDropoutType extract_dropout_variables(T::SharedPtrType);

    // Method for building the log-likelihood
    void build_loglikelihood() override;
    void build_dd_loglikelihood() override;
};

/*
//-----------------------------------------------------------------------------------------
// Class for implementing LOG FRAILTY Model
class LogFrailtyModel final: public TVModel::ModelBase{
public:
    // Constructor
    LogFrailtyModel() = default;
    LogFrailtyModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Method for executing the log-likelihood
    void evaluate_loglikelihood() override;

private:
    // Complex data structure
    Params::Parameters parameters;                                          // Class for the parameters

    // Simple data structure
    T::NumberType n_ranges_parameters = 5;                                  // Number of minimum or maximum ranges that have to be provided for this method
    T::NumberType n_parameters;                                             // Number of parameter of the model
    T::VectorNumberType all_n_parameters;                        
    T::VectorXdr& v_parameters = parameters.get_v_parameters();

    T::IdNameType name_method = "LogFrailty";                                     // Name of this method

    // Functions for implementing the likelihood
    std::function<T::VariableType(T::VectorXdr&)> ll_lf; 
    std::function<T::VariableType(T::VectorXdr&, T::SharedPtrType)> ll_group_lf;     
    std::function<T::VariableType(T::VariableType, T::SharedPtrType, T::VectorXdr&)> G; 
    std::function<T::VariableType(T::VariableType, T::IndexType, T::VariableType, T::VectorXdr&)> f_ijk;
    
    // Function for computing the second derivative wrt one dimension
    std::function<T::VariableType(T::IndexType, T::VectorXdr&)> dd_ll_lf;

    // Quadrature weights and nodes
    QuadraturePoints::Points10 points10;
    T::NumberType n_nodes = 10;
    std::array<T::VariableType, 10>& nodes = points10.nodes;
    std::array<T::VariableType, 10>& weights = points10.weights;

    // Virtual method for computing the number of parameters
    T::NumberType compute_n_parameters() override;
    
    T::VectorXdr compute_hessian_diagonal(T::VectorXdr&);
    T::VectorXdr compute_se(T::VectorXdr&);
    
    // Method for extracting the parameters from the vector
    T::TupleLFType extract_parameters(T::VectorXdr& v_parameters_) const;
    T::TupleDropoutType extract_dropout_variables(const T::SharedPtrType) const;
    T::VectorXdr extract_time_to_event(const T::SharedPtrType) const;

    // Method for building the log-likelihood
    void build_loglikelihood() override;
    void build_dd_loglikelihood() override;
};

*/

}
#endif // TIME_VARYING_MODEL_DERIVED



