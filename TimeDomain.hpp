#ifndef TIMEDOMAIN_HPP
#define TIMEDOMAIN_HPP

#include <iostream>

// Include header files
#include "TypeTraits.hpp"


namespace TimeDomainInfo{
using T = TypeTraits;

class TimeDomain{
public:
    // Default Constructor
    TimeDomain();

    // Constructor
    TimeDomain(const T::FileNameType& filename) {read_from_file(filename);};

    // Getter for the Time class variables
    T::TimeType get_time_begin() const {return time_begin;};
    T::TimeType get_time_end() const {return time_end;};
    T::NumberType get_n_intervals() const {return n_intervals;};
    T::VectorXdr get_v_intervals() const {return v_intervals;};

    // Overload 
    T::VectorXdr& get_v_intervals() {return v_intervals;};
    T::NumberType& get_n_intervals() {return n_intervals;};

    // Print vector of intervals
    void print_v_intervals() const {std::cout << v_intervals << std::endl;};


private:
    T::TimeType time_begin;                                 // Left bound of the time-domain
    T::TimeType time_end;                                   // Right bound of the time-domain
    T::VectorXdr v_intervals;                               // Subdivison of the time-domain
    T::NumberType n_intervals;                              // Number of intervals in which the time-domain is splitted

    // Method for reading data from file using GetPot
    void read_from_file(const T::FileNameType& filename);

    // Method for checking the filname is correct
    T::CheckType check_filename(const T::FileNameType& filename) const;

    // Method for checking condition for the number of intervals
    T::CheckType check_condition(const T::NumberType& size) const;

    // Method for checking conditions for time bounds
    T::CheckType check_condition(const T::VectorType& v_intervals_) const;

};

} // end namespace

#endif // TIMEDOMAIN_HPP

