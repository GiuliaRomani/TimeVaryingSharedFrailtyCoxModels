#ifndef RESULTS_HPP
#define RESULTS_HPP

#include "TypeTraits.hpp"

namespace ResultsMethod{
using T = TypeTraits;

/**
 * Results class contains the results of the application of one time-varying model, 
 * whose name is specified.
*/

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
    */
    Results(const T::IdNameType& name_model_, const T::NumberType n_parameters_,  const T::VectorXdr& optimal_parameters_, 
            const T::VariableType optimal_loglikelihood_, const T::VectorXdr& se_, const T::VectorXdr& sd_frailty_);


    /**
    * Method for printing in a customized way the content of the class
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


    /**
     * Compute the Akaike Informaion Criterior (AIC)
    */
    void compute_AIC();
};
} // end namespace



#endif // RESULTS_HPP
