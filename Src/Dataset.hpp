#ifndef DATASET_HPP
#define DATASET_HPP

// Include header files
#include "TypeTraits.hpp"
#include "TimeDomain.hpp"

// Include libraries
#include <iostream>

/**
 * Dataset class contains the dataset and other basic variables necessary to the execution of the time-varying models.
*/

// Class 
namespace TVSFCM{
using T = TypeTraits;

class Dataset: public TimeDomain{
public:
    /**
     * Default constructor
    */
    Dataset() = default;

    /**
     * Constructor for the class. It initializes the variables starting from two input files
     * @param filename1 file .txt containing time-variables
     * @param filename2 file .txt containing the dataset
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
    T::VectorXdr time_to_event;                     //! Vector of time-to-event (n_individuals)
    T::MatrixXdr e_time;                            //! Matrix of the e_time (n_individuals, n_intervals)
    T::MatrixXdr dropout_intervals;                 //! Matrix of the dropout events (n_individuals, n_intervals)
    T::VectorXdrGroupType dataset_group;            //! Vector of the individual group (n_gindividual)
    T::MapType map_groups;                          //! Map associating to each group the index of individuals belonging to that group (n_groups couples)


    /**
     * Method for reading data from file and initializing the protected data structure 
     * @param filename2_ Name of the .txt file containing the dataset
    */
    void read_from_file(const T::FileNameType& filename2_);

    /**
     * Method for checking the file provided in input really exists
     * @param filename2_ Name of .txt file 
    */
    void check_filename(const T::FileNameType& filename2_) const;

    /**
     * Add the index (of the dataset) of an individual belonging to a precise group, to the vector containing all the individuals of the same group
     * @param name_group Name of the group the individual belong to
     * @param index_individual Index (of the dataset) of the individual
    */
    void add_to_map_groups(const T::GroupNameType& name_group, const T::IndexType& index_individual);

    /**
     * Extract the vector of dataset indexes of all individuals belonging to a predefined cluster
     * @param name_group Group/cluster, whose individuals must be extracted
     * @return Shared pointer to the vector of indexes
    */
    T::SharedPtrType extract_individuals_group(const T::GroupNameType& name_group) const;

    /**
     * Initialize the temporal dropout variable d_ijk
    */
    void initialize_dropout_intervals();

    /**
     * Initialize the temporal variable e_time using the method e_time_function(...)
    */
    void initialize_e_time();

    /**
     * Method for defining the individual and interval e_time variable, according to its definition provided in the reference
     * @param time_t individual time-to-event
     * @param v_k left boundary of the k-th interval
     * @param v_kk right boundary of the k-th interval
     * @return e_time (e_ijk in the reference)
    */
    T::VariableType e_time_function(T::VariableType time_t, T::VariableType v_k, T::VariableType v_kk);

    /**
     * Method for printing the number of individuals in each group
    */
   void print_dimension_groups();
};

} // end namespace

#endif // DATASET_HPP
