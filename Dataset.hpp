#ifndef DATASET_HPP
#define DATASET_HPP

// Include header files
#include "TypeTraits.hpp"
#include "TimeDomain.hpp"

// Include libraries
#include <iostream>

/**
 * DatasetInfo class contains the dataset and other basic variables necessary to the execution of the time-varying model.
*/


// Class 
namespace DatasetInfoClass{
using T = TypeTraits;
using Time = TimeDomainInfo::TimeDomain;

class DatasetInfo: public TimeDomainInfo::TimeDomain{
public:
    /**
     * Default constructor
    */
    DatasetInfo() = default;

    /**
     * Constructor for the class. Initializes all the protected variables
     * @param filename1 file .txt containing time-variables
     * @param filename2 file .txt containing the dataset
    */
    DatasetInfo(const T::FileNameType& filename1, const T::FileNameType& filename2);

    /**
     * Extract the vector of indexes of all individuals of the dataset belonging to a predefined cluster
     * @param name_group group/cluster, whose individuals must be extracted
     * @return shared pointer to the vector of indexes
    */
    std::shared_ptr<T::VectorIndexType> extract_individuals_group(const T::GroupNameType& name_group) const;

    
    // Print methods
    void print_dataset() const {std::cout << dataset << std::endl;};
    void print_dataset_group() const {std::cout << dataset_group << std::endl;};
    void print_time_to_event() const {std::cout << time_to_event << std::endl;};
    void print_dropout_intervals() const {std::cout << dropout_intervals << std::endl;};
    void print_map_groups() const;
    void print_individuals_group(const T::GroupNameType& name_group) const;
    void print_e_time() const {std::cout << e_time << std::endl;};
    void print_n_individuals() const {std::cout << n_individuals << std::endl;};
    void print_n_regressors() const {std::cout << n_regressors << std::endl;};
    void print_n_groups() const {std::cout << n_groups << std::endl;};
    


protected:
    T::NumberType n_individuals;                    // Number of individuals 
    T::NumberType n_regressors;                     // Number of regressors
    T::NumberType n_groups;                         // Number of groups in the cluster variable

    T::MatrixXdr dataset;                           // Matrix of the dataset (n_individuals, n_regressors)
    T::VectorXdr time_to_event;                     // Vector of time-to-event (n_individuals)
    T::MatrixXdr e_time;                            // Matrix of the e_time (n_individuals, n_intervals)
    T::MatrixXdr dropout_intervals;                 // Matrix of the dropout events (n_individuals, n_intervals)
    T::VectorXdrGroupType dataset_group;            // Vector of the individual group (n_gindividual)
    T::MapType map_groups;                          // Map associating to each group the index of individuals belonging to that group (n_groups couples)


    /**
     * Method for reading data from file and initializing the protected data structure 
     * @param filename2 name of the .txt file containing the dataset
    */
    void read_from_file(const T::FileNameType& filename2);

    /**
     * Add the index (of the dataset) of an individual belonging to a precise group, to the vector containing all the individuals of the same group
     * @param name_group name of the group the individual belong to
     * @param index_individual index (of the dataset) of the individual
    */
    void add_to_map_groups(const T::GroupNameType& name_group, const T::IndexType& index_individual);

    /**
     * Initialize the protected variable, dropout interval 
    */
    void initialize_dropout_intervals();

    /**
     * Initialize the protected variable, e_time using the method e_time_function(...)
    */
    void initialize_e_time();

    /**
     * Method for defining the individual and interval e_time variable, according to its definition provided in the references
     * @param time_t individual time-to-event
     * @param v_k left boundary of the k-th interval
     * @param v_kk right boundary of the k-th interval
     * @return e_time (e_ijk in the reference)
    */
    T::VariableType e_time_function(T::VariableType time_t, T::VariableType v_k, T::VariableType v_kk);
};

} // end namespace
#endif // DATASET_HPP
