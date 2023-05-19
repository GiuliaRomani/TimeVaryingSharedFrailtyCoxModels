#ifndef RESULTS_HPP
#define RESULTS_HPP

#include "TypeTraits.hpp"

using T = TypeTraits;

namespace ResultsMethod{
class Results{
public:
    // Constructor
    Results(const T::IdNameType& name_method_, const T::NumberType n_parameters_,  
            const T::VectorXdr& optimal_parameters_, const T::VariableType optimal_loglikelihood_): 
            name_method(name_method_), 
            n_parameters(n_parameters_), 
            optimal_parameters(optimal_parameters_), 
            optimal_loglikelihood(optimal_loglikelihood_),
            AIC(compute_AIC()) 
            {};

    // Method for printing in a fine way the results
    void print_results() const;

    // Method for getting some general info about the class Results
    void print_help() const;

private:
    T::IdNameType name_method;                          // Name of the applied method

    T::NumberType n_parameters;                         // Number of parameters of the method
    T::VectorXdr optimal_parameters;                    // Parameters obtained by the maximization procedure

    T::VariableType optimal_loglikelihood;              // Value of the maximez log-likelihooh
    T::VariableType AIC;                                // AIKAKE INFORMATION CRITERION


    // Method for compute the AIC 
    T::VariableType compute_AIC() const;
}
} // end namespace



#endif // RESULTS_HPP