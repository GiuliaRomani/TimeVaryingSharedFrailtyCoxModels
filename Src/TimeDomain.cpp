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

//! Constructor
TimeDomain::TimeDomain(const T::FileNameType& filename1_) { 
    //! Read the data from file
    read_from_file(filename1_);
};

//! Method for reading data from file using GetPot
void TimeDomain::read_from_file(const T::FileNameType& filename1_){
    GetPot datafile(filename1_.c_str());

    //! Read dimension of the time-domain vector and check it is properly defined
    T::IntType size_v_intervals = datafile("TimeDomain/length_vector_intervals", 0);
    check_condition(size_v_intervals);
    
    //! Resize the vector of time domain and fill it with its element
    v_intervals.resize(n_intervals);
    for(T::IndexType i = 0; i < n_intervals; ++i){
        v_intervals[i] = datafile("TimeDomain/vector_intervals", std::numeric_limits<T::VariableType>::quiet_NaN(), i);
    }

    //! Check correctness of the initialization of the time-domain vector
    check_condition(v_intervals);

    //! Sort the vector in case it is not sorted
    std::sort(v_intervals.begin(), v_intervals.end()); 

    //! Initialize the number of intervals as the (dimension of the time-domain vector - 1) 
    n_intervals -= 1;
};

//! Method for checking that the number of elements of the time-domain vector is non negative and null.
//! Then convert it from an integer to an unsigned int
void TimeDomain::check_condition(const T::IntType size_v_intervals_){
    if(size_v_intervals_ < 0)
        throw MyException("Provided negative dimension of the time domain vector.");
    else if(size_v_intervals_ == 0){
        throw MyException("Null or not provided dimension of the time-domain vector.");
    }
    n_intervals = static_cast<T::NumberType>(size_v_intervals_);
};

//! Method for checking conditions for time bounds
void TimeDomain::check_condition(const T::VectorTimeType & v_intervals_) const{
    //! If the entire vector is not provided, the first element will be NaN.
    if(std::isnan(*(v_intervals_.begin()))){
        throw MyException("Time-domain vector is not provided.");
    }

    //! If the vector is provided but with a different dimension from the one 
    //! indicated, the programs throws an exception.
    //! If at least one element of the vector is negative, the program throws an exception.
    for(const auto & t: v_intervals_){
        if(std::isnan(t))
            throw MyException("Wrong information about the lenght of the time-domain vector.");
        else if(t < 0)
            throw MyException("At least one element of the time-domain vector is negative.");
    }
};

} // end namespace

