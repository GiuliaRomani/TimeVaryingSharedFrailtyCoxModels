#ifndef TIMEDOMAIN_HPP
#define TIMEDOMAIN_HPP

#include <iostream>

// Include header files
#include "TypeTraits.hpp"


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
     * @param filename name of the .txt file in which the time variables are kept
    */
    TimeDomain(const T::FileNameType& filename) {read_from_file(filename);};

    /*
    // Print vector of intervals
    void print_v_intervals() const {std::cout << v_intervals << std::endl;};
    void print_n_intervals() const {std::cout << n_intervals << std::endl;};
    T::NumberType get_n_intervals() const {return n_intervals;};
    T::VectorXdr get_v_intervals() const {return v_intervals;};
    */


protected:
    T::VectorXdr v_intervals;                               // Subdivison of the time-domain
    T::NumberType n_intervals;                              // Number of intervals in which the time-domain is splitted

    /**
     * Method for reading time variables from a file, using GetPot
     * @param filename name of a .txt file containing time variables
    */
    void read_from_file(const T::FileNameType& filename);

    /**
     * Method for checking the file from which we read the time variables really exists.
     * @param filename name of the file .txt containing time variables
     * return boolean variable for the status of the checking operation
    */
    T::CheckType check_filename(const T::FileNameType& filename) const;

    /**
     * Method for checking that the vector of time intervals is provided and has not null dimension.
     * @param size dimension of the vector of time intervals
     * return boolean variable for the status of the checking operation
    */
    T::CheckType check_condition(const T::NumberType& size) const;

    // Method for checking conditions for time bounds
    /**
     * Method for checking the right initialization of the vector of time intervals
     * @param v_intervals_ vector of time intervals
     * return boolean variable for the status of the checking operation
    */
    T::CheckType check_condition(const T::VectorType& v_intervals_) const;

};

} // end namespace

#endif // TIMEDOMAIN_HPP

