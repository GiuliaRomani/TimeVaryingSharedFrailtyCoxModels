#ifndef DATASET_HPP
#define DATASET_HPP

#include "TypeTraits.hpp"
// #include "TimeDomain.hpp"

#include <iostream>

namespace DatasetInfoClass{
using T = TypeTraits;

class DatasetInfo{
public:
    // Constructor
    DatasetInfo(const T::FileNameType& filename1, const T::FileNameType& filename2);

    void print_dataset() const {std::cout << dataset << std::endl;};
    void print_dataset_group() const {std::cout << dataset_group << std::endl;};


private:
    //TimeDomainInfo::TimeDomain time;                // Class time 

    T::NumberType n_individuals;                    // Number of individuals 
    T::NumberType n_regressors;                     // Number of regressors
    //T::NumberType n_groups;                         // Number of groups in the cluster variable

    T::MatrixXdr dataset;                           // Matrix of the dataset (individual, regressors)
    T::VectorXdr time_to_event;                     // Vector of time-to-event
    //T::MatrixXdr dropout_intervals;                 // Matrix of the dropout events
    T::VectorXdrGroupType dataset_group;            // Vector of the individual group

    //T::MapType map_groups;                          // Map associating to each group the index of individuals belonging to that group

    // Method for reading data from file
    void read_from_file(const T::FileNameType& filename2);
};

} // end namespace



#endif // DATASET_HPP