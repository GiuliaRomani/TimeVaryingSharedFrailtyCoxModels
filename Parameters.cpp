#include "Parameters.hpp"
#include "GetPot"

#include <limits>
#include <cmath>
#include <random>
#include <fstream>

namespace Params{

// Constructor
Parameters::Parameters(const T::FileNameType& filename, 
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
                read_from_file(filename);
};

// Method for reading data from file
void Parameters::read_from_file(const T::FileNameType& filename){
    if(!check_filename(filename))
        std::exit(1);

    GetPot datafile(filename.c_str());

    // Reserve some space for the vectors
    range_min_parameters.reserve(n_ranges);
    range_max_parameters.reserve(n_ranges);

    // Read the min and max range of parameters and then check they are provided
    for(T::NumberType i = 0; i < n_ranges; ++i){
        range_min_parameters.push_back(datafile("Parameters/Ranges/range_min", std::numeric_limits<T::VariableType>::quiet_NaN(), i));
        range_max_parameters.push_back(datafile("Parameters/Ranges/range_max", std::numeric_limits<T::VariableType>::quiet_NaN(), i));
    }
    if(!check_condition(range_min_parameters, range_max_parameters)){
        std::exit(1);
    }
    
    // Read the parameters from the input file
    for(T::NumberType i = 0; i < n_parameters; ++i){
    	v_parameters(i) = datafile("Parameters/Values/params", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
    }
    if(!check_condition(v_parameters)){
        std::exit(1);
    }
};

// Method for initializing the vector containing the number of all parameters
void Parameters::initialize_all_n_parameters(const T::VectorNumberType& all_n_parameters_){
    // Resize the vector to the right dimension
    all_n_parameters.resize(n_ranges);

    // Check that the actual dimension coincides with the one provided
    T::NumberType n_ranges_input = all_n_parameters_.size();
    if((n_ranges != n_ranges_input)){
        std::cerr << "Wrong dimension!" << std::endl;
        std::exit(1);
    }

    T::NumberType n_parameters_check = 0;

    // If the dimension is correct, fill the vector
    for(T::IndexType p = 0; p < n_ranges; p++){
        all_n_parameters[p] = all_n_parameters_[p];
        n_parameters_check += all_n_parameters[p];
    }

    // Check the sum of the element in this vector coincides with the n_parameters
    if(n_parameters_check != n_parameters){
        std::cerr << "Error in all_n_parameters. Too many or not enough values." << std::endl;
        std::exit(1);
    }
};

// Method for checking the filename is correct and exits
T::CheckType Parameters::check_filename(const T::FileNameType& filename_) const{
    std::ifstream check(filename_);
    if(check.fail()){
        std::cerr << "File called " << filename_ << " does not exist." << std::endl;
        return false;
    }

    // The filename is correctly provided
    return true;
};

// Method for checking conditions for range bound
T::CheckType Parameters::check_condition(const T::VectorType& range_min_, const T::VectorType& range_max_) const{
    const T::NumberType& n = range_min_.size();
    
    for(T::NumberType i = 0; i < n; ++i){
        if(std::isnan(range_min_[i]) || std::isnan(range_max_[i])){
            std::cerr << "A parameter range is not provided." << std::endl;
            return false;
        }

        if(range_min_[i] > range_max_[i]){
            std::cerr << "For a parameter, min range is greater than max range" << std::endl;
            return false;
        }
    }

    // All the parameters ranges are provided
    return true;
};

// Method for checking conditions for parameters values
T::CheckType Parameters::check_condition(const T::VectorXdr& v_parameters_) const{    
    // Check that each single value of the parameter vector is properly defined
    T::NumberType n, actual_j = 0;
    T::VariableType a,b = 0.;
    
    for(T::NumberType i = 0; i < n_ranges; ++i){
    	n = all_n_parameters[i];
	    a = range_min_parameters[i];
	    b = range_max_parameters[i];
	
	    for(T::NumberType j = 0; j < n; ++j){
	        if((v_parameters(actual_j) < a) || (v_parameters(actual_j) > b)){
	     	    std::cerr << "Value of parameter in position " << actual_j << " not in the ranges" << std::endl;
		        return false;
	        }
	        actual_j += 1;
	    }
    }

    // All the parameter values are in the correct ranges
    std::cout << "All parameter values have correct assignment." << std::endl;
    return true;
};

} // end namespace

