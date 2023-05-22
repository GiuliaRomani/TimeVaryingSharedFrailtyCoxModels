/*
#include "Parameters.hpp"
#include "GetPot"

#include <limits>
#include <cmath>
#include <random>
#include <fstream>

namespace Params{

// Constructor
Parameters::Parameters(const T::IdNameType& name_method_, const T::NumberType& n_parameters_, 
            const T::NumberType n_intervals_, const T::NumberType n_regressors_, 
            const T::FileNameType& filename):
            name_method(name_method_),
            n_parameters(n_parameters_),
            n_intervals(n_intervals_),
            n_regressors(n_regressors_),
            n_other_parameters(compute_n_other_parameters()){
                // Initialize the number of other parameters
                n_other_parameters = compute_n_other_parameters();

                // Read data from file
                read_from_file(filename);

                // Initialize the vector of parameters
                initialize_v_parameters();   
            };

// Mehtod for computing the number of other parameter
T::NumberType Parameters::compute_n_other_parameters(){
    return (n_parameters - n_regressors - n_intervals);
}

// Method for reading data from file
void Parameters::read_from_file(const T::FileNameType& filename){
    if(!check_filename(filename))
        std::exit(1);

    GetPot datafile(filename.c_str());

    // Read the number of ranges, i.e. the dimension of the vector containing the parameters ranges
    // Then check the value is provided and is correct
    T::NumberType n_ranges = datafile("Parameters/n_ranges", std::numeric_limits<T::IntType>::quiet_NaN());
    if(!check_condition(n_ranges)){
        std::exit(1);
    }

    // Reserve some space for the vectors
    range_min_parameters.reserve(n_ranges);
    range_max_parameters.reserve(n_ranges);

    // Read the min and max range of parameters and then check they are provided
    for(T::NumberType i = 0; i < n_ranges; ++i){
        range_min_parameters.push_back(datafile("Parameters/range_min", std::numeric_limits<T::VariableType>::quiet_NaN(), i));
        range_max_parameters.push_back(datafile("Parameters/range_max", std::numeric_limits<T::VariableType>::quiet_NaN(), i));
    }
    if(!check_condition(range_min_parameters, range_max_parameters)){
        std::exit(1);
    }
};

// Method for initializing the vector of parameters
void Parameters::initialize_v_parameters() {
    // Resize the vector of parameters to the right dimension
    v_parameters.resize(n_parameters);

    // Condense in an array all the different number of sub-parameters
    std::array<T::NumberType, 3> all_n_parameters = {n_intervals, n_regressors, n_other_parameters};

    // Define the index for storing the elements in the dynamic vector
    T::NumberType j_actual = 0;

    // Defie the engine and the distribution for the pseudo-random-number generator
    std::default_random_engine generator;
    for(T::NumberType i = 0; i < 3; ++i){
        T::NumberType & n = all_n_parameters[i];
        T::VariableType & a = range_min_parameters[i];
        T::VariableType & b = range_max_parameters[i];

        std::uniform_real_distribution<> distribution(a, b);
        for(T::NumberType j = 0; j < n; ++j){
            v_parameters(j_actual) = distribution(generator);
            j_actual += 1;
        }   
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

// Method for checking condition for the number of intervals
T::CheckType Parameters::check_condition(const T::NumberType& n_ranges_) const{
    if((n_ranges_ == 0) || std::isnan(n_ranges_)){
        std::cerr << "No ranges provided for the parameters initialization." << std::endl;
        return false; 
    }

    if(n_ranges_ < 3 ){
        std::cerr << "Provided number of ranges not sufficient." << std::endl;
        return false;
    }

    if(n_ranges_ > 3 ){
        std::cerr << "Too many provided number of ranges." << std::endl;
        return false;
    }

    // If the number of ranges is provided
    return true;
};

// Method for checking conditions for time bounds
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



} // end namespace

*/