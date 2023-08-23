#ifndef DATASET_HPP
#define DATASET_HPP

// Include header files
#include "TypeTraits.hpp"
#include "TimeDomain.hpp"

// Include libraries
#include <iostream>

// Class 
namespace TVSFCM{
using T = TypeTraits;

/**
 * Dataset class contains the dataset and other basic variables necessary to the execution of the time-varying models.
*/

class Dataset: public TimeDomain{
public:
    /**
     * Default constructor
    */
    Dataset() = default;

    /**
     * Constructor for the class. It initializes the variables starting from two input files
     * @param filename1 File containing time-variables
     * @param filename2 File containing the dataset
    */
    Dataset(const T::FileNameType& filename1, const T::FileNameType& filename2);
  
    /**
     * Default destructor 
    */
    ~ Dataset() = default;

protected:
    T::NumberType n_individuals;                    //! Number of individuals 
    T::NumberType n_regressors;                     //! Number of regressors
    T::NumberType n_groups;                         //! Number of groups in the cluster variable

    T::MatrixXdr dataset;                           //! Matrix of the dataset (n_individuals, n_regressors)
    T::VectorTimeXdr time_to_event;                 //! Vector of time-to-event (n_individuals)
    T::MatrixTimeXdr e_time;                        //! Matrix of the e_time (n_individuals, n_intervals)
    T::MatrixXdr dropout_intervals;                 //! Matrix of the dropout events (n_individuals, n_intervals)
    T::VectorXdrGroupType dataset_group;            //! Vector of the individual group (n_gindividual)
    T::MapType map_groups;                          //! Map associating to each group the index of individuals belonging to that group (n_groups couples)


    /**
     * Method for reading data from file and initializing the protected data structure 
     * @param filename2_ Name of the .txt file containing the dataset
    */
    void read_from_file(const T::FileNameType& filename2_);

    /**
     * Method for checking that the future dimensions of the dataset are not-negative, otherwise an excpetion is thrown
     * @param n_individuals_ Number of individuals
     * @param n_regressors_ Number of regressors
    */
    void check_condition(const T::IntType n_individuals_, const T::IntType n_regressors_);

    /**
     * Add the row index (of the dataset) of an individual belonging to a precise group, to the vector containing all the individuals of the same group
     * @param name_group Name of the group the individual belongs to
     * @param index_individual Row index (of the dataset) of an individual
    */
    void add_to_map_groups(const T::GroupNameType& name_group, const T::IndexType& index_individual) noexcept;

    /**
     * Initialize the temporal dropout variable d_ijk
    */
    void initialize_dropout_intervals() noexcept;

    /**
     * Initialize the temporal variable e_time using the method e_time_function(...)
    */
    void initialize_e_time() noexcept;

    /**
     * Method for defining the individual and interval e_time variable, through the function defined in Appendix A of the report.
     * @param time_t Individual time-to-event
     * @param v_k Left boundary of the k-th interval
     * @param v_kk Right boundary of the k-th interval
     * @return e_time (e_ijk in Appendix)
    */
    T::TimeType e_time_function(const T::TimeType time_t, const T::TimeType v_k, const T::TimeType v_kk) const noexcept;

    /**
     * Extract the vector of dataset indexes of all individuals belonging to a predefined cluster
     * @param name_group Group/cluster, whose individuals must be extracted
     * @return Shared pointer to the vector of indexes
    */
    T::SharedPtrType extract_individuals_group(const T::GroupNameType& name_group) const noexcept;

    /**
     * Method for printing the number of individuals in each group
    */
    void print_dimension_groups() const noexcept;
};

} // end namespace

#endif // DATASET_HPP
