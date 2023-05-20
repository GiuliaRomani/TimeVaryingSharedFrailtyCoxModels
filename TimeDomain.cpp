#include "TimeDomain.hpp"

#include <limits>


namespace TimeDomainInfo{

// Define this variable to indicate if a value is not provided (NP)
const inline T::NumberType NP1 = std::numeric_limits<T::NumberType>::quiet_NaN();
const inline T::TimeType NP2 = std::numeric_limits<T::TimeType>::quiet_NaN();

// Method for reading data from file using GetPot
// const T::FileNameType& filename
void TimeDomain::read_from_file(){
    GetPot datafile("DataToolFile.txt");

    // Read number of subdivision of time domain and check correctness
    const T::NumberType size_intervals = datafile("TimeParameters/Intervals/length_list_intervals", NP1);
    if(!check_condition(size_intervals))
        std::exit(1);

    // Initialize the vector of the time domain subdivision
    time_intervals.resize(size_intervals);
    for(T::NumberType i=0; i < size_intervals; ++i){
        time_intervals(i) = datafile("TimeParameters/Intervals/list_intervals", NP2, i);
    }

    // Initialize the time bounds
    time_begin = time_intervals(1);
    time_end = time_intervals(size_intervals,1);
    if(!check_condition(time_begin, time_end)){
        std::exit(1);
    }

    // Initialize the number of intervals 
    n_interval = size_intervals - 1;
};

// Method for checking condition for the number of intervals
T::CheckType TimeDomain::check_condition(const T::NumberType& size_int) const{
    if(size_int == 0){
        std::cerr << "Null number of subdivisions of time domain." << std::endl;
        return false; 
    }
    if(size_int == NP1){
        std::cerr << "Number of subdivision of the time domain not provided." << std::endl;
        return false;
    }

    // If number of subdivision of time domain is correctly provided
    return true;

};

// Method for checking conditions for time bounds
T::CheckType TimeDomain::check_condition(const T::TimeType& begin, const T::TimeType& end) const{

    // usare qualcosa di diverso per controllare il segno tipo is_signed
    // Check for positiviness of time bounds
    if((time_begin < 0) || (time_end < 0)){
        std::cerr << "Negative time bounds" << std::endl;
        std::exit(1);
    }

    // Check for different time bounds
    if((time_begin == time_end)){
        std::cerr << "Equal time bounds" << std::endl;
        std::exit(1);
    }

    if((time_begin == NP2) || (time_end == NP2)){
        std::cerr << "At least one time bound is not provided" << std::endl;
        std::exit(1);
    }

};


} // end namespace