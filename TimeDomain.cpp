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

namespace TimeDomainInfo{

// Default constructor
TimeDomain::TimeDomain(): n_intervals(0) {
    v_intervals.resize(n_intervals);
};

TimeDomain::TimeDomain(const T::FileNameType& filename1_) {
    // Check the input file exists
    check_filename(filename1_);
    
    // Once sure it exists
    read_from_file(filename1_);
};

// Method for reading data from file using GetPot
void TimeDomain::read_from_file(const T::FileNameType& filename1_){
    GetPot datafile(filename1_.c_str());

    // Read number of subdivision of time domain and check correctness
    const T::NumberType size_intervals = static_cast<T::NumberType>(datafile("TimeDomain/length_vector_intervals", std::numeric_limits<T::VariableType>::quiet_NaN()));
    check_condition(size_intervals);

    // To check if the vector of time intervals is sorted, to load it into a normal vector
    T::VectorType v_intervals_copy(size_intervals, 0.);
    for(T::NumberType i = 0; i < size_intervals; ++i){
        v_intervals_copy[i] = datafile("TimeDomain/vector_intervals", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
    }
    std::sort(v_intervals_copy.begin(), v_intervals_copy.end());

    // Check correctness of the initialization
    check_condition(v_intervals_copy);

    // Initialize the vector of the time domain subdivision using an Eigen trick
    v_intervals.resize(size_intervals);
    T::MappedVectorType vv_intervals(v_intervals_copy.data(), size_intervals); 
    v_intervals = vv_intervals;  

    // Initialize the number of intervals 
    n_intervals = size_intervals - 1;
};

// Method for checking the filename is correct and exists
void TimeDomain::check_filename(const T::FileNameType& filename1_) const{
    std::ifstream check(filename1_);
    if(check.fail()){
        T::ExceptionType msg1 = "File ";
        T::ExceptionType msg2 = msg1.append((filename1_).c_str());
        T::ExceptionType msg3 = msg2.append(" does not exist.");
        //throw MyException("File provided does not exist.");
        throw MyException(msg3);
    }
};

// Method for checking that the number of elements of the time vector is not null and is really provided
void TimeDomain::check_condition(const T::NumberType& size_int) const{
    if(size_int == 0){
        throw MyException("Null number of subdivisions of time domain.");
    }
    if(std::isnan(size_int)){
        throw MyException("Number of subdivision of the time domain not provided.");
    }
};

// Method for checking conditions for time bounds
void TimeDomain::check_condition(const T::VectorType & v_intervals_) const{
    // If the entire vector is not provided, the first element will be NaN.
    if(std::isnan(*(v_intervals_.begin()))){
        throw MyException("Vector of time intervals is not provided.");
    }

    // If the vector is provided but with a different dimension from the one 
    // indicated, the programs thorws an exception
    for(const auto & t: v_intervals_){
        if(std::isnan(t)){
            throw MyException("Wrong information about the lenght of the time intervals vector");
        }
    }
};


} // end namespace

