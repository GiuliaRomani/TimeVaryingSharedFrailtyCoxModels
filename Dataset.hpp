#ifndef DATASET_HPP
#define DATASET_HPP

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>


// Class 
namespace DatasetInfoClass{
using T = TypeTraits;

class DatasetInfo{
public:
    // Constructor
    DatasetInfo() = default;
    DatasetInfo(const T::FileNameType& filename2, T::NumberType n_intervals_, const T::VectorXdr& v_intervals_);

    // Extract the uique pointer to the name in the map of groups
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


protected:
    T::NumberType n_intervals;			            // Number of intervals
    T::VectorXdr v_intervals;		                // Vector of intervals

    T::NumberType n_individuals;                    // Number of individuals 
    T::NumberType n_regressors;                     // Number of regressors
    T::NumberType n_groups;                         // Number of groups in the cluster variable

    T::MatrixXdr dataset;                           // Matrix of the dataset (individual, regressors)
    T::VectorXdr time_to_event;                     // Vector of time-to-event
    T::MatrixXdr e_time;
    T::MatrixXdr dropout_intervals;                 // Matrix of the dropout events
    T::VectorXdrGroupType dataset_group;            // Vector of the individual group
    T::MapType map_groups;                          // Map associating to each group the index of individuals belonging to that group


    // Method for reading data from file
    void read_from_file(const T::FileNameType& filename2);

    // Add name of the group and the unique pointer to the vector containing the indexes of individuals belonging to that group
    void add_to_map_groups(const T::GroupNameType& name_group, const T::IndexType& index_individual);

    // Initialize the dropout_intervals variable
    void initialize_dropout_intervals();

    // Initialize the e_time matrix and define the function to compute the e_time
    void initialize_e_time();
    T::VariableType e_time_function(T::VariableType time_t, T::VariableType v_k, T::VariableType v_kk);
};

} // end namespace
#endif // DATASET_HPP
