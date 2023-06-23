
#ifndef TIME_VARYING_MODEL_BASE
#define TIME_VARYING_MODEL_BASE

// Include all necessary header files
#include "TimeDomain.hpp"
#include "Dataset.hpp"
#include "Results.hpp"

// Include libraries
#include <iostream>

// Class
namespace TVModel{
using T = TypeTraits;

class ModelBase{
public:
    // Base constructor
    ModelBase() = default;
    ModelBase(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Getter 
    T::VectorXdr get_sd_frailty() const {return sd_frailty;};
    T::VectorXdr get_variance_frailty() const {return variance_frailty;};

    // Define and Compute the log-likelihood
    virtual void optimize_loglikelihood() = 0;

    // Print 
    void print_map_groups() const {return database.print_map_groups();};
    void print_n_regressors() const {std::cout << n_regressors << std::endl;};
    void print_n_intervals() const {std::cout << n_intervals << std::endl;};


    // Destructor
    ~ ModelBase() = default;


protected:
    // Complex data structures
    DatasetInfoClass::DatasetInfo database;                                 // Dataset containing the individual covariates
    TimeDomainInfo::TimeDomain time;                                        // Class time 
    ResultsMethod::Results result;                                          // Class for the results

    // Simple variables
    T::NumberType& n_regressors = database.get_n_regressors();                          // Number of regressors
    T::NumberType& n_groups = database.get_n_groups();                                  // Number of groups
    T::NumberType& n_intervals = time.get_n_intervals();                                // Number of intervals

    // Simple data structures
    T::VectorXdr variance_frailty;                                          // Vector for time-interval variance of the frailty
    T::VectorXdr sd_frailty;                                                // Vector for time-interval sd of the frailty


    // Virtual method to compute the number of parameters of each model
    virtual T::NumberType compute_n_parameters() = 0;

    // Virtual method for extracting the parameters from the vector
    virtual T::TuplePPType extract_parameters(T::VectorXdr v_parameters_) = 0;

    // Method for building the log-likelihood
    virtual void build_loglikelihood() = 0;

    // Virtual method to derive the interval variance of the frailty
    //virtual void compute_sd_frailty() = 0;


};

} // end namespace


#endif // TIME_VARYING_MODEL_BASE

