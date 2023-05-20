#include "TimeDomain.hpp"
#include "GetPot"

#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>


namespace TimeDomainInfo{

// Method for reading data from file using GetPot
// const T::FileNameType& filename
void TimeDomain::read_from_file(const T::FileNameType& filename){
    GetPot datafile(filename.c_str());

    // Read number of subdivision of time domain and check correctness
    const T::NumberType size_intervals = datafile("time/length_list_intervals", std::numeric_limits<T::IntType>::quiet_NaN());
    if(!check_condition(size_intervals))
        std::exit(1);

    // Initialize the vector of the time domain subdivision and check correctness
    time_intervals.reserve(size_intervals);
    for(T::NumberType i=0; i < size_intervals; ++i){
        time_intervals.push_back(datafile("time/list_intervals", std::numeric_limits<T::TimeType>::quiet_NaN(), i));
    }
    std::sort(time_intervals.begin(), time_intervals.end());
    if(!check_condition(time_intervals)){
        std::exit(1);
    }

    // Initialize the time bounds
    time_begin = time_intervals[0];
    time_end = time_intervals[size_intervals-1];

    // Initialize the number of intervals 
    n_interval = size_intervals - 1;
};

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
T::CheckType TimeDomain::check_condition(const T::VectorType& time_intervals_) const{
    // If the entire vector is not provided, the first element will be NaN.
    if(std::isnan(*(time_intervals.begin()))){
        std::cerr << "List of time intervals is not provided" << std::endl;
        return false;
    }

    // If the vector is provided but with a different dimension from the one 
    // indicated, the programs aborts
    for(const auto& instant: time_intervals){
        if(std::isnan(instant)){
            std::cerr << "Wrong information about the lenght of the time intervals vector" << std::endl;
            return false;
        }
    }

    // If both time bounds are correctly provided
    return true;

};


} // end namespace