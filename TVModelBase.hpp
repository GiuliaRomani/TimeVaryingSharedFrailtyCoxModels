
#ifndef TIME_VARYING_MODEL_BASE
#define TIME_VARYING_MODEL_BASE

// Include all necessary header files
#include "TimeDomain.hpp"
#include "Dataset.hpp"
#include "Results.hpp"
#include "ToolLikelihood.hpp"

// Include libraries
#include <iostream>

// Class
namespace TVModel{
using T = TypeTraits;
using Dataset = DatasetInfoClass::DatasetInfo;
using Time = TimeDomainInfo::TimeDomain;

class ModelBase: public Time,
		            public Dataset{
public:
    // Base constructor
    ModelBase() = default;
    ModelBase(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Getter 
    // T::VectorXdr get_sd_frailty() const {return sd_frailty;};
    // T::VectorXdr get_variance_frailty() const {return variance_frailty;};

    // Define and Compute the log-likelihood
    virtual void optimize_loglikelihood() = 0;
    virtual void evaluate_loglikelihood(T::VectorXdr&) = 0;

    // virtual void print_extract_parameters() = 0;

    // Print 
    void print_map_groups() const {return Dataset::print_map_groups();};
    void print_n_regressors() const {std::cout << Dataset::n_regressors << std::endl;};
    void print_n_intervals() const {std::cout << Time::n_intervals << std::endl;};


    // Destructor
    ~ ModelBase() = default;


protected:
    // Complex data structures
    ResultsMethod::Results result;                                          // Class for the results
    
    // Base variables for optimization method
    ToolsLikelihood::Tools tool;
    // T::VariableType tol_ll = tool.tol_ll;
    // T::VariableType tol_optim = tool.tol_ll;
    // T::NumberType n_extrarun = tool.n_extrarun;
    T::VariableType h_dd = tool.h_dd;

    // Simple data structures
    T::VectorXdr variance_frailty;                                          // Vector for time-interval variance of the frailty
    T::VectorXdr sd_frailty;                                                // Vector for time-interval sd of the frailty

    // Virtual method to compute the number of parameters of each model
    virtual T::NumberType compute_n_parameters() = 0;

    // Method for building the log-likelihood
    virtual void build_loglikelihood() = 0;
    virtual void build_dd_loglikelihood() = 0;
    

    // Virtual method to derive the interval variance of the frailty
    //virtual void compute_sd_frailty() = 0;

};

} // end namespace
#endif // TIME_VARYING_MODEL_BASE

