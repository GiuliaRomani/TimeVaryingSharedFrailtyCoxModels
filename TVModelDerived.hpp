
#ifndef TIME_VARYING_MODEL_DERIVED
#define TIME_VARYING_MODEL_DERIVED

#include "TVModelBase.hpp"
#include "TypeTraits.hpp"
#include "QuadraturePoints.hpp"

#include <iostream>

namespace TVModel{
using T = TypeTraits;

class PowerParameterModel: public ModelBase{
public:
    // Constructor
    PowerParameterModel() = default;
    PowerParameterModel(const T::FileNameType& filename1, const T::FileNameType& filename2);



private:
    // Complex data structure
    Params::Parameters parameters;                                          // Class for the parameters
    //T::MatrixXdr 

    // Simple data structure
    T::NumberType n_ranges_parameters = 4;                                  // Number of minimum or maximum ranges that have to be provided for this method
    T::NumberType n_parameters;                                             // Number of parameter of the model
    T::VectorNumberType all_n_parameters;
    T::IdNameType name_method = "PowerParameter";                           // Name of this method
    T::VectorXdr & v_parameters = parameters.get_v_parameters();

    // Functions for implementing the likelihood
    // It collects the result
    std::function<T::VariableType(T::VariableType, T::IndexType, T::VectorXdr) ll_pp; 


    // Quadrature nodes and weights
    QuadraturePoints::Points9 points9;
    T::NumberType n_node = 9;
    std::array<T::VariableType, 9>& nodes = points9.nodes;
    std::array<T::VariableType, 9>& weights = points9.weights;

    // Virtual method for computing the number of parameters
    T::NumberType compute_n_parameters() override;

    // Method for 

    // Virtual method to derive the interval variance of the frailty
    // void compute_sd_frailty() override;

};
}


#endif // TIME_VARYING_MODEL_DERIVED


