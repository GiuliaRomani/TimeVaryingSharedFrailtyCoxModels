#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "TypeTraits.hpp"

#include <iostream>

namespace Params{
using T = TypeTraits;

class Parameters{
public:
    Parameters(const T::NumberType& n_parameters_, const T::IdNameType& name_method_,
                const T::FileNameType& filename):
                n_parameters(n_parameters_),
                name_method(name_method_) {
                    read_from_file(filename);
                    initialize_v_parameters();
                };
    
    // Getter for the class variables
    T::NumberType get_n_parameters() const {return n_parameters;};
    T::IdNameType get_name_method() const {return name_method;};
    T::VectorXdr get_v_parameters() const {return v_parameters;};
    T::VectorType get_range_min_parameters() const {return range_min_parameters;};
    T::VectorType get_range_max_parameters() const {return range_max_parameters;};

    // Setter for the vector of parameters
    void set_v_parameters(const T::VectorType& v_parameters_);
    void set_v_parameter(const T::IndexType& i, const T::VariableType& value);


private:
    T::NumberType n_parameters;                             // Number of parameters of the method
    T::VectorXdr v_parameters;                              // Vector of parameters
    T::VectorType range_min_parameters;                     // Vector containing the range of the parameters
    T::VectorType range_max_parameters;                     // Vector containing the range of the parameters
    
    T::IdNameType name_method;                              // Name of the method

    // Method for reading the range of the parameters from the file provided
    void read_from_file(const T::FileNameType& filename);   

    // Method for initializing parameter vector with pseudo-random-number generator
    void initialize_v_parameters();

    // Method for checking condition for the number of ranges provide (internal to read_from_file())
    T::CheckType check_condition(const T::NumberType& n_ranges_) const;

    // Method for checking conditions for time bounds
    T::CheckType check_condition(const T::VectorType& range_min_, const T::VectorType& range_max_) const;

};
} // end namespace



#endif // PARAMETERS_HPP
