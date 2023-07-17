#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include "TypeTraits.hpp"

#include <iostream>

/**
 * Parameters class contains the unknown constrained parameters of the models. 
 * Each model defines its own vector and they differ both in shape and numerosity. 
 * Precisely, its dimension is computed starting from some variables belonging to different classes (Time, Dataset).
 * 
 * Theoretically this vector should be obtained as a result of a maximization procedure, that is not here applied.
 * Thus, the optimized vector is already provided and the correctness of its components is proved,
 * by controlling that each variable is contained in the right range so that no problems arise during the 
 * computation of the log-likelihood.
*/

namespace Params{
using T = TypeTraits;

class Parameters{
public:
    /**
     * Default Constructor
    */
    Parameters() = default;

    /**
     * Constructor
     * @param filename Name of the file where some parameters-related variables are located and extracted
     * @param n_parameters_ Number of model parameters. Once the time-varying model is chosen, the number of parameters is defined
     * @param n_intervals_ Number of intervals of the temporal domain
     * @param n_regressors_ Number of regressors of the dataset
     * @param n_ranges_ Number of different type of parameters constituting the vector
     * @param all_n_parameters Vector containing the numerosity associated to each type of parameter.
    */
    Parameters(const T::FileNameType& filename, 
            const T::NumberType& n_parameters_, const T::NumberType& n_intervals_, 
            const T::NumberType& n_regressors_, const T::NumberType& n_ranges_, 
            const T::VectorNumberType& all_n_parameters_);
    

    /**
     * Default desctructor
    */
   ~ Parameters() = default;


protected:
    T::NumberType n_parameters;                             //! Number of parameters of the method
    T::NumberType n_intervals;                              //! Number of intervals of the time domain
    T::NumberType n_regressors;                             //! Number of regressors of the dataset

    T::VectorXdr v_parameters;                              //! Vector of parameters
    T::NumberType n_ranges;                                 //! Number of ranges to be provided
    T::VectorType range_min_parameters;                     //! Vector containing the range of the parameters
    T::VectorType range_max_parameters;                     //! Vector containing the range of the parameters
    T::VectorNumberType all_n_parameters;                   //! Subdivision of the parameters
 

    /**
     * Method for reading the range of the parameters and the optimized vector of parameters from an input file.
     * Some controls on the meaningfull anc correctness of the initialized values are done.
     * @param filename Name of a .txt file
    */
    void read_from_file(const T::FileNameType& filename);   

    /**
     * Simple method for compying the content of the vector of parameters into another vector.
     * Here we also check that the overall number of parameters coincides with the sum of the element inside this vector
     * @param all_n_parameters Vector whose content must be copied
    */
    void initialize_all_n_parameters(const T::VectorNumberType& all_n_parameters_);

    // Method for checking the correctness and existence of the filename provided
    /**
     * Method for checking the existence of the file passed as input
     * @param filename Name of the file, whose existence in the directory must be checked
     * @return Boolean variable indicating the status of the control (FALSE for not-existing)
    */
    T::CheckType check_filename(const T::FileNameType& filename_) const;

    /**
     * This method controls that the vector of the minimum and maximum range of the parameters has been correctly 
     * filled during its construction (that is, the minimum range is less than or equal to the maximum range)
     * and that no elements are missing. Othwerise, FALSE status is returned.
     * @param range_min_ Vector of the minimum range of parameters
     * @param range_max_ Vector of the maximum range of parameters
     * @return Boolean variable indicating the status of the control 
    */
    T::CheckType check_condition(const T::VectorType& range_min_, const T::VectorType& range_max_) const;
    
    /**
     * This method overloads the previous one based on the input and control that each parameter
     * belong to its correct range, otherwise a negative status value is returned
     * @param v_parameters_ Vector of optimized parameters
     * @return Boolean variable indicating the status of the control
    */
    T::CheckType check_condition(const T::VectorXdr& v_parameters_) const;

};
} // end namespace



#endif // PARAMETERS_HPP

