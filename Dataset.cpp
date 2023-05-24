#include "Dataset.hpp"
#include "GetPot"

#include <sstream>
#include <memory>

namespace DatasetInfoClass{
// Constructor
DatasetInfo::DatasetInfo(const T::FileNameType& filename1, const T::FileNameType& filename2){
    // Initialize TimeDomain time object through its constructor
    time = TimeDomainInfo::TimeDomain(filename1);

    // Initialize the number of groups to zero
    n_groups = 0;

    // Initialize the rest of the class
    read_from_file(filename2);

    // Initialize the dropout_intervals variable
    initialize_dropout_intervals();
};

// Method for reading from file
void DatasetInfo::read_from_file(const T::FileNameType& filename2){
    T::CheckType first_line = false;
    T::NumberType i = 0;                                // Number of row
    T::NumberType j = 0;                                // Number of column
    std::ifstream filename(filename2.c_str());

    while(!filename.fail() & !filename.eof()){
        std::string line;
        getline(filename, line);
        std::istringstream value_line(line);

        if(!first_line){
            value_line >> n_individuals >> n_regressors;
            first_line = true;

            dataset_group.resize(n_individuals);
            time_to_event.resize(n_individuals);
            dataset.resize(n_individuals, n_regressors);
        }
        else{
            // Save the group inside the vector
            T::GroupNameType code_group;
            getline(value_line, code_group, ',');
            dataset_group(i) = code_group;
            add_to_map_groups(code_group, i);

            // Read the regressors, one at the time,  and convert them into a double from a string
            std::string value_regressor;
            while(getline(value_line, value_regressor, ',')){
                double value = std::stod(value_regressor.c_str());
                if( j < n_regressors)
                    dataset(i,j) = value;
                else if(j == n_regressors)
                    time_to_event(i) = value;               // The last element is not a regressor but the time-to-event
                j += 1;
            }   
            // Reset j and update i
            j = 0;
            i += 1;
        }
    }
};

// Method for adding individual code group name and index to the map
void DatasetInfo::add_to_map_groups(const T::GroupNameType& name_group, const T::IndexType& index_individual){
    T::MapType::iterator group_position = map_groups.find(name_group);
    if(group_position == map_groups.end()){
        map_groups[name_group] = std::make_unique<T::VectorIndexType>();
        map_groups[name_group]->push_back(index_individual);
        n_groups += 1;
    }
    else{
        group_position->second->push_back(index_individual);
    }
};

// Initialize the dropout_intervals variable
void DatasetInfo::initialize_dropout_intervals(){
    // Get the necessary variables as const reference
    const T::NumberType & n_intervals = time.get_n_intervals();
    const T::VectorXdr & vector_intervals = time.get_vector_intervals();
    std::cout << time.get_n_intervals() << std::endl;

    // Resize the matrix according to the right dimensions and fill it with null elements
    dropout_intervals.resize(n_individuals, n_intervals);
    dropout_intervals.setZero();

    // Fill the matrix according to the condition
    for(T::NumberType i = 0; i < n_individuals; ++i){
        for(T::NumberType k = 0; k < n_intervals; ++k){
            if((time_to_event(i) < vector_intervals(k+1)) & (time_to_event(i) >= vector_intervals(k)))
                dropout_intervals(i,k) = 1;
        }
    } 
}

// Method for print the element of the map
void DatasetInfo::print_map_groups() const{
    for(auto& [name_group, ptr_vector_index]: map_groups){
        std::cout << "In group " << name_group << " individuals with index:"<<std::endl;
        for(auto& index: *ptr_vector_index){
            std::cout << index << std::endl;
        }
    }
};


} // end namespace

