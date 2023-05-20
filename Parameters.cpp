#include "Parameters.hpp"
#include "GetPot"

#include <limits>
#include <cmath>
#include <random>

namespace Params{

void Parameters::read_from_file(const T::FileNameType& filename){
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
    for(T::NumberType i=0; i < n_ranges; ++i){
        range_min_parameters.push_back(datafile("Parameters/range_min", std::numeric_limits<T::VariableType>::quiet_NaN(), i));
        range_max_parameters.push_back(datafile("Parameters/range_max", std::numeric_limits<T::VariableType>::quiet_NaN(), i));
    }
    if(!check_condition(range_min_parameters, range_max_parameters)){
        std::exit(1);
    }
};

void Parameters::initialize_v_parameters() {
    std::default_random_engine generator;
    fors


    std::uniform_real_distribution<T::VariableType> distribution(a,b);

     distribution(generator);

};

// Method for checking condition for the number of intervals
T::CheckType Parameters::check_condition(const T::NumberType& n_ranges_) const{
    if((n_ranges_ == 0) || std::isnan(n_ranges_)){
        std::cerr << "No ranges provided for the parameters initialization." << std::endl;
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