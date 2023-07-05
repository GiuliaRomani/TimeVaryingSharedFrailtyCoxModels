#ifndef RESULTS_HPP
#define RESULTS_HPP

#include "TypeTraits.hpp"

namespace ResultsMethod{
using T = TypeTraits;

class Results{
public:
    // Constructor
    Results() = default;
    Results(const T::IdNameType& name_model_, const T::NumberType n_parameters_,  const T::VectorXdr& optimal_parameters_, 
            const T::VariableType optimal_loglikelihood_, const T::VectorXdr& se_);

    // Method for printing in a fine way the results
    void print_results() const;


private:
    T::IdNameType name_model;				// Name of the model
    T::NumberType n_parameters;                         // Number of parameters of the method
    T::VectorXdr optimal_parameters;                    // Parameters obtained by the maximization procedure

    T::VariableType optimal_loglikelihood;              // Value of the maximez log-likelihooh
    T::VariableType AIC;                                // AKAIKE INFORMATION CRITERION
    
    T::VectorXdr se;					// Standard error of the parameters


    // Method for compute the AIC 
    void compute_AIC();
};
} // end namespace



#endif // RESULTS_HPP
