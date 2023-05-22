/*
#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "TypeTraits.hpp"

#include <iostream>

namespace Params{
using T = TypeTraits;

class Parameters{
public:
    // Constructor
    Parameters(const T::IdNameType& name_method_, const T::NumberType& n_parameters_, 
                const T::NumberType n_intervals_, const T::NumberType n_regressors_, 
                const T::FileNameType& filename);
    
    // Getter for the class variables
    T::NumberType get_n_parameters() const {return n_parameters;};
    T::NumberType get_other_parameters() const {return n_other_parameters;};

    T::VectorXdr get_v_parameters() const {return v_parameters;};
    T::VectorType get_range_min_parameters() const {return range_min_parameters;};
    T::VectorType get_range_max_parameters() const {return range_max_parameters;};

    T::IdNameType get_name_method() const {return name_method;};

    // Overload the getter so that it can return a non-const l-value reference
    T::VectorXdr& get_v_parameters() {return v_parameters;};
    T::NumberType& get_n_parameters() {return n_parameters;};

    // Setter for the vector of parameters
    // void set_v_parameters(const T::VectorType& v_parameters_);
    // void set_v_parameter(const T::IndexType& i, const T::VariableType& value);

    // Print hthe vector of parameters
    void print_v_parameters() const {std::cout << v_parameters << std::endl;};


private:
    T::IdNameType name_method;                              // Name of the method

    T::NumberType n_parameters;                             // Number of parameters of the method
    T::NumberType n_intervals;                              // Number of intervals of the time domain
    T::NumberType n_regressors;                             // Number of regressors of the dataset
    T::NumberType n_other_parameters;                       // Number of other parameters of the method

    T::VectorXdr v_parameters;                              // Vector of parameters
    T::VectorType range_min_parameters;                     // Vector containing the range of the parameters
    T::VectorType range_max_parameters;                     // Vector containing the range of the parameters
    

    // Method for computing the number of other parameters
    T::NumberType compute_n_other_parameters();

    // Method for reading the range of the parameters from the file provided
    void read_from_file(const T::FileNameType& filename);   

    // Method for initializing parameter vector with pseudo-random-number generator
    void initialize_v_parameters();

    // Method for checking the correctness and existence of the filename provided
    T::CheckType check_filename(const T::FileNameType& filename_) const;

    // Method for checking condition for the number of ranges provide (internal to read_from_file())
    T::CheckType check_condition(const T::NumberType& n_ranges_) const;

    // Method for checking conditions for time bounds
    T::CheckType check_condition(const T::VectorType& range_min_, const T::VectorType& range_max_) const;

};
} // end namespace



#endif // PARAMETERS_HPP
*/
