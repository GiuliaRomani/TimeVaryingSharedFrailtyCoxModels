#ifndef RESULTS_HPP
#define RESULTS_HPP

// Include header files
#include "TypeTraits.hpp"

/**
 * Results class contains the results of the application of one time-varying model, 
 * whose name is specified.
*/

namespace TVSFCM{
using T = TypeTraits;

class Results{
public:
    /**
     * Default Constructor
    */
    Results() = default;

    /**
    * Constructor: It initializes the private components of the Results class
    * @param name_model_ The name of the used time-varying model
    * @param n_parameters_ The number of parameters 
    * @param optimal_parameters_ The vector formed by the optimal parameters, with dimension equal to n_parameters_
    * @param optimal_loglikelihood_ The value of the log-likelihood associated to the optimal_parameters and the dataset
    * @param se_ The vector of standard error of the optimal_parameters_
    * @param sd_frailty_ The vector of standard deviation of the frailty
    * @param n_threads_ The number of threads used in the parallel version
    * @param chunk_size_ The number of chunk for the for loop in the parallel version
    * @param schedule_type_name_ The name of type of schedule used in the for loop strategy
    */
    Results(const T::IdNameType& name_model_, const T::NumberType n_parameters_,  const T::VectorXdr& optimal_parameters_, 
            const T::VariableType optimal_loglikelihood_, const T::VectorXdr& se_, const T::VectorXdr& sd_frailty_, 
            T::NumberType n_threads_, T::NumberType chunk_size_, const T::ScheduleType& schedule_type_name_);


    /**
    * Method for printing in a customized way the content of the class. 
    * 
    * According to the application or not of the parallel version, it calls another print method that 
    * prints the results.
    */
    void print_results() const;


private:
    T::IdNameType name_model;				            //! Name of the model

    T::NumberType n_parameters;                         //! Number of parameters of the model
    T::VectorXdr optimal_parameters;                    //! Parameters obtained by the maximization procedure
    T::VariableType optimal_loglikelihood;              //! Value of the maximized log-likelihood
    T::VariableType AIC;                                //! AKAIKE INFORMATION CRITERION
    T::VectorXdr se;					                //! Standard error of the parameters
    T::VectorXdr sd_frailty;                            //! Standard deviation of the frailty

    T::NumberType n_threads;                            //! Number of threads for the parallel version
    T::NumberType chunk_size;                           //! Dimension of the chunk for for loop
    T::ScheduleType schedule_type_name;                 //! Name of the type of for loop schedule

    /**
     * Compute the Akaike Informaion Criterior (AIC)
    */
    void compute_AIC();

    /**
     * Method for printing the results in case of parallel implementation
    */
    void print_results_parallel() const;

    /**
     * Method for printing the results in case of not-parallel implementation
    */
    void print_results_noparallel() const;

};

} // end namespace

#endif // RESULTS_HPP
