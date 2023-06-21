#ifndef DATASET_HPP
#define DATASET_HPP

#include "TypeTraits.hpp"
#include "TimeDomain.hpp"

#include <iostream>

namespace DatasetInfoClass{
using T = TypeTraits;

class DatasetInfo{
public:
    // Constructor
    DatasetInfo() = default;
    DatasetInfo(const T::FileNameType& filename1, const T::FileNameType& filename2);

    // Getter
    T::NumberType get_n_groups() const {return n_groups;};
    T::NumberType get_n_individuals() const {return n_individuals;};
    T::NumberType get_n_regressors() const {return n_regressors;};

    // Overload the same operations
    T::NumberType& get_n_groups() {return n_groups;};
    T::NumberType& get_n_individuals() {return n_individuals;};
    T::NumberType& get_n_regressors() {return n_regressors;};
    T::MapType& get_map_groups() {return map_groups;};

    // Extract the uique pointer to the name in the map of groups
    std::shared_ptr<T::VectorIndexType> extract_individuals_group(const T::GroupNameType& name_group) const;

    // Print methods
    void print_dataset() const {std::cout << dataset << std::endl;};
    void print_dataset_group() const {std::cout << dataset_group << std::endl;};
    void print_time_to_event() const {std::cout << time_to_event << std::endl;};
    void print_dropout_intervals() const {std::cout << dropout_intervals << std::endl;};
    void print_map_groups() const;
    void print_individuals_group(const T::GroupNameType& name_group) const;


private:
    TimeDomainInfo::TimeDomain time;                // Class time 
    T::NumberType& n_intervals = time.get_n_intervals();
    T::VectorXdr& v_intervals = time.get_v_intervals();

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
    T::VariableType e_time_function(T::VariableType time_t, T::IndexType k);
};

} // end namespace



#endif // DATASET_HPP
