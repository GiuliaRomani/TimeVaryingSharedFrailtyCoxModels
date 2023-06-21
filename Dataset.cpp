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

    // Initialize the e_time matrix
    initialize_e_time();
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
            code_group = code_group.substr(1,4);
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
        map_groups[name_group] = std::make_shared<T::VectorIndexType>();
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
    const T::VectorXdr & vector_intervals = time.get_v_intervals();

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
};

// Initialize the e_time matrix
void DatasetInfo::initialize_e_time(){
    // Resize the matrix
    e_time.resize(n_individuals, n_intervals);

    // Fill the matrix
    for(T::IndexType i = 0; i < n_individuals; ++i){
        T::VariableType time_individual = time_to_event(j);
        for(T::IndexType k = 0; k < n_intervals; ++k){
            e_time(i,j) = e_time_function(time_individual, k);
        }
    }
}

// Define the function to compute the e_time value in the matrix
T::VariableType DatasetInfo::e_time_function(T::VariableType time_t, T::IndexType k){
    if(time_t < v_intervals(k))
        return 0.;
    else if((time_t >= v_intervals(k)) & (time_t < v_intervals(k+1)))
        return (time_t - v_intervals(k));
    else if(time_t >= v_intervals(k+1))
        return (v_intervals(k+1) - v_intervals(k));
};

// Method for print the element of the map
void DatasetInfo::print_map_groups() const{
    for(auto& [name_group, ptr_vector_index]: map_groups){
        std::cout << "In group " << name_group << " individuals with index:"<<std::endl;
        for(auto& index: *ptr_vector_index){
            std::cout << index << std::endl;
        }
    }
};

// Method for ptinting the element of a single group
void DatasetInfo::print_individuals_group(const T::GroupNameType& name_group) const{
    std::shared_ptr<T::VectorIndexType> individuals_group = extract_individuals_group(name_group);
    if(individuals_group == nullptr)
        std::cerr << "No group with this name!" << std::endl;
    else{
        std::cout << "Indexes in group " << name_group << std::endl;
        for(const auto i: *individuals_group)
            std::cout << i << std::endl;
    }
};


// Extract the shared pointer to the name group in the map of groups
std::shared_ptr<T::VectorIndexType> DatasetInfo::extract_individuals_group(const T::GroupNameType& name_group) const{
    T::MapType::const_iterator group_position = map_groups.find(name_group);
    if(group_position == map_groups.cend())
        return nullptr;
    else{
        return group_position->second;
    }
};


} // end namespace

