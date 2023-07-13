
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

    // Define and Compute the log-likelihood
    virtual void evaluate_loglikelihood() = 0;

    // Destructor
    ~ ModelBase() = default;


protected:
    // Complex data structures
    ResultsMethod::Results result;                                          // Class for the results
    
    // Base variables for optimization method
    ToolsLikelihood::Tools tool;
    T::VariableType h_dd = tool.h_dd;

    // Simple data structures
    T::VectorXdr variance_frailty;                                          // Vector for time-interval variance of the frailty
    T::VectorXdr sd_frailty;                                                // Vector for time-interval sd of the frailty
    
    // Simple data structure for the estimate of the standard error
    T::VectorXdr hessian_diag;
    T::VectorXdr se;

    // Virtual method to compute the number of parameters of each model
    virtual void compute_n_parameters() = 0;

    // Method for building the log-likelihood
    virtual void build_loglikelihood() = 0;
    virtual void build_dd_loglikelihood() = 0;
    
    virtual void compute_hessian_diagonal(T::VectorXdr&) = 0;
    virtual void compute_se(T::VectorXdr&) = 0;
    

    // Virtual method to derive the interval variance of the frailty
    //virtual void compute_sd_frailty() = 0;

};

} // end namespace
#endif // TIME_VARYING_MODEL_BASE

