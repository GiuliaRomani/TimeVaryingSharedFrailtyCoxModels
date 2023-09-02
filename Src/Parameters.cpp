// Include header files
#include "Parameters.hpp"
#include "GetPot"
#include "MyException.hpp"

// Include 
#include <limits>
#include <cmath>
#include <random>
#include <fstream>
#include <string>

namespace TVSFCM{

// Constructor
Parameters::Parameters(const T::FileNameType& filename1_, 
            const T::NumberType& n_parameters_, const T::NumberType& n_intervals_, 
            const T::NumberType& n_regressors_, const T::NumberType& n_ranges_,
            const T::VectorNumberType& all_n_parameters_):
            n_parameters(n_parameters_),
            n_intervals(n_intervals_),
            n_regressors(n_regressors_),
            n_ranges(n_ranges_){
            	// Resize the vector of parameters
            	v_parameters.resize(n_parameters);	
            	
            	// Initialize the vector of the number of parameters
                initialize_all_n_parameters(all_n_parameters_);
                
                // Read data from file
                read_from_file(filename1_);
};

// Method for reading data from file
void Parameters::read_from_file(const T::FileNameType& filename1_){
    GetPot datafile(filename1_.c_str());

    // Resize the dimension of the vectors of ranges according to the number of categories
    range_min_parameters.resize(n_ranges);
    range_max_parameters.resize(n_ranges);

    // Read the min and max range of parameters and then check they are correctly provided
    for(T::IndexType i = 0; i < n_ranges; ++i){
        range_min_parameters[i] = datafile("Parameters/Ranges/range_min", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
        range_max_parameters[i] = datafile("Parameters/Ranges/range_max", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
    }
    check_condition(range_min_parameters, range_max_parameters);

    // Read the parameters from the input file and check they are correctly contained in the ranges
    for(T::IndexType i = 0; i < n_parameters; ++i){
    	v_parameters(i) = datafile("Parameters/Values/params", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
    }
    check_condition(v_parameters);
};

// Method for initializing the vector containing the numerosity of each category
void Parameters::initialize_all_n_parameters(const T::VectorNumberType& all_n_parameters_){
    // Resize the dimension of the vector to the number of categories
    all_n_parameters.resize(n_ranges);

    // Copy the content of the vector and control that the number of parameters is equal to the one provided
    T::NumberType n_parameters_check = 0;
    for(T::IndexType p = 0; p < n_ranges; p++){
        all_n_parameters[p] = all_n_parameters_[p];
        n_parameters_check += all_n_parameters[p];
    }

    // Check the sum of the element in this vector coincides with the n_parameters
    if(n_parameters_check != n_parameters){
        throw MyException("Error in all_n_parameters. Too many or not enough values.");
    }
};

// Method for checking conditions for range bound
void Parameters::check_condition(const T::VectorType& range_min_, const T::VectorType& range_max_) const{
    const T::NumberType& n = range_min_.size();
    
    // Control that the values contained in the ranges vector are not NaN
    for(T::IndexType i = 0; i < n; ++i){
        if(std::isnan(range_min_[i]) || std::isnan(range_max_[i])){
            throw MyException("Either the minimum or maximum range of at least a category is not provided.");
        }

        // Control that the minimum range is smaller than the maximum
        if(range_min_[i] > range_max_[i]){
            T::ExceptionType msg1 = "For category ";
            T::ExceptionType msg2 = msg1.append(std::to_string(i));
            T::ExceptionType msg3 = msg2.append(", min range is greater than max range.");
            throw MyException(msg3);
        }
    }
};

// Method for checking conditions for parameters values
void Parameters::check_condition(const T::VectorXdr& v_parameters_) const{    
    // Check that each single value of the parameter vector is properly defined and contained in its category range
    T::NumberType n = 0;                // Numerosity of a category
    T::IndexType actual_j = 0;          // Index for the parameter vector
    T::VariableType a,b = 0.;           // Min and max range
    
    // Loop over the categories in all_n_parameters
    for(T::IndexType i = 0; i < n_ranges; ++i){
    	n = all_n_parameters[i];
	    a = range_min_parameters[i];
	    b = range_max_parameters[i];
	
        // Loop over the parameters in the category
	    for(T::IndexType j = 0; j < n; ++j){
	    	if(std::isnan(v_parameters(actual_j))){
                throw MyException("At least one parameter is not provided ");
	    	}
	        else if((v_parameters(actual_j) < a) || (v_parameters(actual_j) > b)){
                T::ExceptionType msg1 = "Value of parameter in position ";
                T::ExceptionType msg2 = msg1.append(std::to_string(actual_j));
                T::ExceptionType msg3 = msg2.append(" not in the range.");
                throw MyException(msg3);
	        }
	        actual_j += 1;
	    }
    }
};

} // end namespace

