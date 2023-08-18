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
    // Once sure it exists
    read_from_file(filename1_);
};

// Method for reading data from file using GetPot
void TimeDomain::read_from_file(const T::FileNameType& filename1_){
    GetPot datafile(filename1_.c_str());

    // Read number of subdivision of time domain and check correctness
    T::IntType size_intervals_ = datafile("TimeDomain/length_vector_intervals", 0);
    T::NumberType size_intervals = check_condition(size_intervals_);
    
    // To check if the vector of time intervals is sorted, to load it into a normal vector
    T::VectorType v_intervals_copy(size_intervals, 0.);
    for(T::NumberType i = 0; i < size_intervals; ++i){
        v_intervals_copy[i] = datafile("TimeDomain/vector_intervals", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
    }

    // Check correctness of the initialization
    check_condition(v_intervals_copy);

    // Sort the vector in case it is not sorted
    std::sort(v_intervals_copy.begin(), v_intervals_copy.end());

    // Initialize the vector of the time domain subdivision using an Eigen trick
    v_intervals.resize(size_intervals);
    T::MappedVectorType vv_intervals(v_intervals_copy.data(), size_intervals); 
    v_intervals = vv_intervals;  

    // Initialize the number of intervals 
    n_intervals = size_intervals - 1;
};

// Method for checking that the number of elements of the time vector is not null and is really provided
T::NumberType TimeDomain::check_condition(const T::IntType& size_int) const{
    if(size_int < 0)
        throw MyException("Provided negative number of subdivision of time domain.");
    else if(size_int == 0){
        throw MyException("Null or not provided number of subdivisions of time domain.");
    }
    return static_cast<T::NumberType>(size_int);
};

// Method for checking conditions for time bounds
void TimeDomain::check_condition(const T::VectorType & v_intervals_) const{
    // If the entire vector is not provided, the first element will be NaN.
    if(std::isnan(*(v_intervals_.begin()))){
        throw MyException("Vector of time intervals is not provided.");
    }

    // If the vector is provided but with a different dimension from the one 
    // indicated, the programs throws an exception.
    // If at least one element of the vector is negative, the program throws an exception.
    for(const auto & t: v_intervals_){
        if(std::isnan(t))
            throw MyException("Wrong information about the lenght of the time intervals vector");
        else if(t < 0)
            throw MyException("At least one element of the vector of time intervals is negative.");
    }
};


} // end namespace

