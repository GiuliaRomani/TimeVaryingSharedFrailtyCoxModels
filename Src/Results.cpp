// Include header files
#include "Results.hpp"

// Include libraries
#include <iostream>
#include <iomanip>

namespace TVSFCM{

Results::Results(const T::IdNameType& name_model_, const T::NumberType n_parameters_,  const T::VectorXdr& optimal_parameters_, 
                 const T::VariableType optimal_loglikelihood_, const T::VectorXdr& se_, const T::VectorXdr& sd_frailty_,
                 T::NumberType n_threads_, T::NumberType chunk_size_, const T::ScheduleType& schedule_type_name_): 
                 // Initialize variables
                 name_model(name_model_),
                 n_parameters(n_parameters_), 
                 optimal_parameters(optimal_parameters_), 
                 optimal_loglikelihood(optimal_loglikelihood_),
                 se(se_),
                 sd_frailty(sd_frailty_),
                 n_threads(n_threads_),
                 chunk_size(chunk_size_),
                 schedule_type_name(schedule_type_name_) {
                    // Initialize the AIC
                    compute_AIC();
                };


void Results::compute_AIC() {
    AIC = (2 * n_parameters - 2 * optimal_loglikelihood);
};

// Method for printing the results of a method call
void Results::print_results() const {
    if(n_threads == 1)
        print_results_noparallel();
    else
        print_results_parallel();
};


void Results::print_results_noparallel() const{
    std::cout << "--------------------- Appplication Summary -----------------------" << std::endl;
    std::cout << std::endl;

    std::cout << "Parallel execution: NO" << std::endl;
    std::cout << std::endl; 
    
    std::cout << "Model executed: " << name_model << std::endl;
    std::cout << std::endl;

    std::cout << "Model results: number of parameters = " << n_parameters << std::endl;
    std::cout << "               log-likelihood       = " << std::setiosflags(std::ios::scientific) << std::setprecision(3) << optimal_loglikelihood << std::endl;
    std::cout << "               AIC                  = " << std::setprecision(3) << AIC << std::endl;
    std::cout << std::endl;

    std::cout << "Standard error of the parameters:" << std::endl;
    for(T::IndexType i = 1; i <= se.size(); ++i){
        std::cout << std::setprecision(3) << se(i-1) << " ";
        if (i % 6 == 0)
            std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Standard deviation of the frailty:" << std::endl;;
    for(T::IndexType i = 1; i <= sd_frailty.size(); ++i){
        std::cout << std::setprecision(3) << sd_frailty(i-1) << " ";
        if (i % 6 == 0)
            std::cout << std::endl;
    }
    std::cout << std::endl;
};

void Results::print_results_parallel() const{
    std::cout << "--------------------- Application Summary -----------------------" << std::endl;
    std::cout << std::endl;

    std::cout << "Parallel execution with: number of threads = " << n_threads << std::endl;
    std::cout << "                         chunk size        = " << chunk_size << std::endl;
    std::cout << "                         schedule type     = " << schedule_type_name << std::endl;
    std::cout << std::endl;

    std::cout << "Model results: number of parameters = " << n_parameters << std::endl;
    std::cout << "               log-likelihood       = " << std::setiosflags(std::ios::scientific) << std::setprecision(3) << optimal_loglikelihood << std::endl;
    std::cout << "               AIC                  = " << std::setprecision(3) << AIC << std::endl;
    std::cout << std::endl;

    std::cout << "Standard error of the parameters:" << std::endl;
    for(T::IndexType i = 1; i <= se.size(); ++i){
        std::cout << std::setprecision(3) << se(i-1) << " ";
        if (i % 6 == 0)
            std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Standard deviation of the frailty:" << std::endl;;
    for(T::IndexType i = 1; i <= sd_frailty.size(); ++i){
        std::cout << std::setprecision(3) << sd_frailty(i-1) << " ";
        if (i % 6 == 0)
            std::cout << std::endl;
    }
    std::cout << std::endl;
};



} // end namespace
