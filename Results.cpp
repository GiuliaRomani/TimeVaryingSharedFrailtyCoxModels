#include "Results.hpp"

#include <iostream>
#include <iomanip>

namespace ResultsMethod{

Results::Results(const T::IdNameType& name_model_, const T::NumberType n_parameters_,  const T::VectorXdr& optimal_parameters_, 
                 const T::VariableType optimal_loglikelihood_, const T::VectorXdr& se_, const T::VectorXdr& sd_frailty_,
                 T::NumberType n_threads_): 
                 // Initialize variables
                 name_model(name_model_),
                 n_parameters(n_parameters_), 
                 optimal_parameters(optimal_parameters_), 
                 optimal_loglikelihood(optimal_loglikelihood_),
                 se(se_),
                 sd_frailty(sd_frailty_),
                 n_threads(n_threads_) {
                    // Initialize the AIC
                    compute_AIC();
                };


void Results::compute_AIC() {
    AIC = (2 * n_parameters - 2 * optimal_loglikelihood);
};


void Results::print_results() const {
    std::cout << "------------------- Printing results -----------------------" << std::endl;
    std::cout << "Name of the model: " << name_model << std::endl;
    if(n_threads == 1)
        std::cout << "Parallel execution: NO" << std::endl;
    else
        std::cout << "Parallel execution with " << n_threads << " threads" << std::endl; 

    std::cout << "Number of parameter: " << n_parameters << std::endl;
    std::cout << "Log-likelihood: " << std::setiosflags(std::ios::scientific) << std::setprecision(3) << optimal_loglikelihood << std::endl;
    std::cout << "AIC:  " << std::setprecision(3) << AIC << std::endl;
    //std::cout << "Standard Error of the parameters: \n" <<std::setprecision(3) << se << std::endl;
    //std::cout << "Standard Deviation of the frailty: \n" << std::setprecision(3) << sd_frailty << std::endl;
};

} // end namespace
