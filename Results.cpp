#include "Results.hpp"

#include <iostream>
#include <iomanip>

namespace ResultsMethod{

Results::Results(const T::IdNameType& name_model_, const T::NumberType n_parameters_,  const T::VectorXdr& optimal_parameters_, 
                 const T::VariableType optimal_loglikelihood_, const T::VectorXdr& se_): 
                 // Initialize variables
                 name_model(name_model_),
                 n_parameters(n_parameters_), 
                 optimal_parameters(optimal_parameters_), 
                 optimal_loglikelihood(optimal_loglikelihood_) {
                    // Initialize the AIC
                    compute_AIC();

                    // Resize thestandard error
                    se.resize(n_parameters);
                       
                    // Initialize them
                    se = se_;

                };


void Results::compute_AIC() {
    AIC = (2 * n_parameters - 2 * optimal_loglikelihood);
};


void Results::print_results() const {
    std::cout << "------------------- Printing results -----------------------" << std::endl;
    std::cout << name_model << " model" << std::endl;
    std::cout << "Number of parameters " << n_parameters << std::endl;
    std::cout << "Log-likelihood = " << std::setprecision(5) << optimal_loglikelihood << std::endl;
    std::cout << "AIC = " << std::setprecision(5) << AIC << std::endl;
    std::cout << "Standard Error of the parameters: \n" << se << std::endl;
};

} // end namespace
