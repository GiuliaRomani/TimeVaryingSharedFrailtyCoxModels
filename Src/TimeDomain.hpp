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

namespace TVSFCM{
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
    explicit TimeDomain(const T::FileNameType& filename1_);

    /**
     * Default destructor 
    */
    ~ TimeDomain() = default;

protected:
    T::VectorTimeType v_intervals;                           //! Vector of time instants, constituting the intervals of the time domain. Also called time-domain vector
    T::NumberType n_intervals;                               //! Number of intervals of the time domain

    /**
     * Method for reading time variables from a file, using GetPot
     * @param filename1_ Name of a .txt file containing time variables
    */
    void read_from_file(const T::FileNameType& filename1_);

    /**
     * Method for checking that the dimension of the time-domain vector is provided, it is not null
     * and not negative.
     * Otherwise, it throws an exception.
     * @param size Integer dimension of the time-domain vector
     * @return The same dimension but converted into an unsigned int
    */
    void check_condition(const T::IntType size_v_intervals_);

    /**
     * Method for checking that all the elements of the time-domain vector are really provided.
     * Otherwise, it throws an exception.
     * @param v_intervals_ Vector of time intervals
    */
    void check_condition(const T::VectorTimeType& v_intervals_) const;

};

} // end namespace

#endif // TIMEDOMAIN_HPP

