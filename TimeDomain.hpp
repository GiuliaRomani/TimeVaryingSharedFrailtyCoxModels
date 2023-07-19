#ifndef TIMEDOMAIN_HPP
#define TIMEDOMAIN_HPP

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>

/**
 * Class for the temporal domain. 
 * 
 * It reads the variables from the input file using GetPot and controls everything has been defined correctly, 
 * with no mistakes. Otherwise an exception is thrown for each occurred error.
*/

namespace TimeDomainInfo{
using T = TypeTraits;

class TimeDomain{
public:
    /**
     * Default constructor
    */
    TimeDomain();

    /**
     * Constructor for the initialization of the protected components of the class.
     * @param filename1_ Name of the .txt file in which the time variables are kept
    */
    TimeDomain(const T::FileNameType& filename1_);

    /**
     * Default destructor 
    */
    ~ TimeDomain() = default;

protected:
    T::VectorXdr v_intervals;                               //! Vector of time instants, that constitute the intervals of the time domain
    T::NumberType n_intervals;                              //! Number of intervals of the time domain

    /**
     * Method for reading time variables from a file, using GetPot
     * @param filename1_ Name of a .txt file containing time variables
    */
    void read_from_file(const T::FileNameType& filename1_);

    /**
     * Method for checking that the file from which we read the time variables really exists.
     * Otherwise, it throws an exception.
     * @param filename1_ Name of the file .txt containing time variables
    */
    void check_filename(const T::FileNameType& filename1_) const;

    /**
     * Method for checking that the dimension of the vector of time intervals is provided and has not null dimension.
     * Otherwise, it throws an exception.
     * @param size dimension of the vector of time intervals
    */
    void check_condition(const T::NumberType& size) const;

    /**
     * Method for checking that all the elements of the vector are really provided.
     * Otherwise, it throws an exception.
     * @param v_intervals_ Vector of time intervals
    */
    void check_condition(const T::VectorType& v_intervals_) const;

};

} // end namespace

#endif // TIMEDOMAIN_HPP

