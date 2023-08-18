// Include header files
#include "TimeDomain.hpp"
#include "GetPot"
#include "MyException.hpp"

// Include libraries
#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

namespace TVSFCM{

// Default constructor
TimeDomain::TimeDomain(): n_intervals(0) {
    v_intervals.resize(n_intervals);
};

TimeDomain::TimeDomain(const T::FileNameType& filename1_) { 
    // Read the data from file
    read_from_file(filename1_);
};

// Method for reading data from file using GetPot
void TimeDomain::read_from_file(const T::FileNameType& filename1_){
    GetPot datafile(filename1_.c_str());

    // Read dimension of the time-domain vector
    T::IntType size_v_intervals = datafile("TimeDomain/length_vector_intervals", 0);
    T::NumberType n_intervals = check_condition(size_v_intervals);
    
    // To check if the vector of time intervals is sorted, to load it into a normal vector
    T::VectorType v_intervals(n_intervals, 0.);
    for(T::NumberType i = 0; i < size_v_intervals; ++i){
        v_intervals[i] = datafile("TimeDomain/vector_intervals", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
    }

    // Check correctness of the initialization
    check_condition(v_intervals);

    // Sort the vector in case it is not sorted
    std::sort(v_intervals.begin(), v_intervals.end()); 

    // Initialize the number of intervals as the dimension of the time-domain vector - 1 
    n_intervals -= 1;
};

// Method for checking that the number of elements of the time-domain vector is non negative and null
T::NumberType TimeDomain::check_condition(const T::IntType& size_v_intervals_) const{
    if(size_v_intervals_ < 0)
        throw MyException("Provided negative dimension of the time domain vector.");
    else if(size_v_intervals_ == 0){
        throw MyException("Null or not provided dimension of the time-domain vector.");
    }
    return static_cast<T::NumberType>(size_v_intervals_);
};

// Method for checking conditions for time bounds
void TimeDomain::check_condition(const T::VectorType & v_intervals_) const{
    // If the entire vector is not provided, the first element will be NaN.
    if(std::isnan(*(v_intervals_.begin()))){
        throw MyException("Time-domain vector is not provided.");
    }

    // If the vector is provided but with a different dimension from the one 
    // indicated, the programs throws an exception.
    // If at least one element of the vector is negative, the program throws an exception.
    for(const auto & t: v_intervals_){
        if(std::isnan(t))
            throw MyException("Wrong information about the lenght of the time-domain vector.");
        else if(t < 0)
            throw MyException("At least one element of the time-domain vector is negative.");
    }
};

} // end namespace

