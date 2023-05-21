#include "TimeDomain.hpp"
#include "GetPot"

#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

namespace TimeDomainInfo{

// Method for reading data from file using GetPot
void TimeDomain::read_from_file(const T::FileNameType& filename){
    // Check filename is correctly provided
    if(!check_filename(filename))
        std::exit(1);

    // If the filename is correctly provided
    GetPot datafile(filename.c_str());

    // Read number of subdivision of time domain and check correctness
    const T::NumberType size_intervals = datafile("TimeDomain/length_vector_intervals", std::numeric_limits<T::IntType>::quiet_NaN());
    if(!check_condition(size_intervals))
        std::exit(1);

    // To check if the vector of time intervals is sorted, to load it into a normal vector
    T::VectorType vector_intervals_copy(size_intervals, 0.);
    for(T::NumberType i = 0; i < size_intervals; ++i){
        vector_intervals_copy[i] = datafile("TimeDomain/vector_intervals", std::numeric_limits<T::TimeType>::quiet_NaN(), i);
    }
    std::sort(vector_intervals_copy.begin(), vector_intervals_copy.end());

    if(!check_condition(vector_intervals_copy)){
        std::exit(1);
    }

    // Initialize the time bounds
    //time_begin = *(vector_intervals_copy.begin());
    //time_end = *(vector_intervals_copy.end() - 1);

    time_begin = vector_intervals_copy[0];
    time_end = vector_intervals_copy[size_intervals - 1];

    // Initialize the vector of the time domain subdivision and check correctness
    vector_intervals.resize(size_intervals);
    for(T::NumberType i = 0; i < size_intervals; ++i)
        vector_intervals(i) = vector_intervals_copy[i];

    // Initialize the number of intervals 
    n_intervals = size_intervals - 1;
};

// Method for checking the filname is correct
T::CheckType TimeDomain::check_filename(const T::FileNameType& filename_) const{
    std::ifstream check(filename_);
    if(check.fail()){
        std::cerr << "File called " << filename_ << " does not exist." << std::endl;
        return false;
    }

    // The filename is correctly provided
    return true;
}

// Method for checking condition for the number of intervals
T::CheckType TimeDomain::check_condition(const T::NumberType& size_int) const{
    if(size_int == 0){
        std::cerr << "Null number of subdivisions of time domain." << std::endl;
        return false; 
    }
    if(std::isnan(size_int)){
        std::cerr << "Number of subdivision of the time domain not provided." << std::endl;
        return false;
    }

    // If number of subdivision of time domain is correctly provided
    return true;
};

// Method for checking conditions for time bounds
T::CheckType TimeDomain::check_condition(const T::VectorType & vector_intervals_) const{
    // If the entire vector is not provided, the first element will be NaN.
    if(std::isnan(*(vector_intervals_.begin()))){
        std::cerr << "List of time intervals is not provided" << std::endl;
        return false;
    }

    // If the vector is provided but with a different dimension from the one 
    // indicated, the programs aborts
    for(const auto & t: vector_intervals_){
        if(std::isnan(t)){
            std::cerr << "Wrong information about the lenght of the time intervals vector" << std::endl;
            return false;
        }
    }

    // If both time bounds are correctly provided
    return true;
};


} // end namespace