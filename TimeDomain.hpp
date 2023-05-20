#ifndef TIMEDOMAIN_HPP
#define TIMEDOMAIN_HPP

#include "TypeTraits.hpp"

#include <iostream>


namespace TimeDomainInfo{
using T = TypeTraits;

class TimeDomain{
public:
    // Constructor
    TimeDomain(const T::FileNameType& filename) {read_from_file(filename);};

    // Getter for the Time class variables
    T::TimeType get_time_begin() const {return time_begin;};
    T::TimeType get_time_end() const {return time_end;};
    T::NumberType get_n_intervals() const {return n_interval;};
    T::VectorType get_time_intervals() const {return time_intervals;};



private:
    T::TimeType time_begin;                     // Left bound of the time-domain
    T::TimeType time_end;                       // Right bound of the time-domain
    T::VectorType time_intervals;         // Subdivison of the time-domain
    T::NumberType n_interval;                   // Number of intervals in which the time-domain is splitted

    // Method for reading data from file using GetPot
    void read_from_file(const T::FileNameType& filename);
    //void read_from_file();

    // Method for checking condition for the number of intervals
    T::CheckType check_condition(const T::NumberType& size) const;

    // Method for checking conditions for time bounds
    T::CheckType check_condition(const T::VectorType& time_intervals_) const;

};

} // end namespace

#endif // TIMEDOMAIN_HPP