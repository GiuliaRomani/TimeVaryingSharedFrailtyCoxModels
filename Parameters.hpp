#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "TypeTraits.hpp"

#include <iostream>

namespace Params{
using T = TypeTraits;

class Parameters{
public:
    // Constructor
    Parameters() = default;
    Parameters(const T::FileNameType& filename, 
            const T::NumberType& n_parameters_, const T::NumberType& n_intervals_, 
            const T::NumberType& n_regressors_, const T::NumberType& n_ranges_, 
            const T::VectorNumberType& all_n_parameters_);
    
    // Overload the getter so that it can return a non-const l-value reference
    T::VectorXdr& get_v_parameters() {return v_parameters;};
    T::NumberType& get_n_parameters() {return n_parameters;};

    // Setter for the vector of parameters
    // void set_v_parameters(const T::VectorType& v_parameters_);
    // void set_v_parameter(const T::IndexType& i, const T::VariableType& value);

    // Print hthe vector of parameters
    void print_v_parameters() const {std::cout << v_parameters << std::endl;};


private:
    T::NumberType n_parameters;                             // Number of parameters of the method
    T::NumberType n_intervals;                              // Number of intervals of the time domain
    T::NumberType n_regressors;                             // Number of regressors of the dataset

    T::VectorXdr v_parameters;                              // Vector of parameters
    T::NumberType n_ranges;                                 // Number of ranges to be provided
    T::VectorType range_min_parameters;                     // Vector containing the range of the parameters
    T::VectorType range_max_parameters;                     // Vector containing the range of the parameters
    T::VectorNumberType all_n_parameters;                   // Subdivision of the parameters
    

    // Method for reading the range of the parameters from the file provided
    void read_from_file(const T::FileNameType& filename);   

    // Method for initializing the vector containing the number of all parameters
    void initialize_all_n_parameters(const T::VectorNumberType& all_n_parameters_);

    // Method for initializing parameter vector with pseudo-random-number generator
    void initialize_v_parameters();

    // Method for checking the correctness and existence of the filename provided
    T::CheckType check_filename(const T::FileNameType& filename_) const;

    // Method for checking conditions for ranges bounds
    T::CheckType check_condition(const T::VectorType& range_min_, const T::VectorType& range_max_) const;
    T::CheckType check_condition(const T::VectorXdr& v_parameters_) const;

};
} // end namespace



#endif // PARAMETERS_HPP

