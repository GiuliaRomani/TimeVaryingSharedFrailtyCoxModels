#include "Results.hpp"

namespace ResultsMethod{
T::VariableType Results::compute_AIC() const {
    return (2 * n_parameters - 2 * optimal_loglikelihood);
};

void Results::print_results() const {
    std::cout << "------------------- Printing results -----------------------"
    std::cout << "Method " << name_method << ", with " << n_parameters << "parameters." <<std::endl;
    std::cout << "Log-likelihood = " << std::setprecision(3) << optimal_loglikelihood << std::endl;
    std::cout << "AIC = " << std::setprecision(3) << AIC << std::endl;
};

void Results::print_help() const {
    // ... TODO ...
};

} // end namespace