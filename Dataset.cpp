#include "Dataset.hpp"
#include "GetPot"

#include <sstream>

namespace DatasetInfoClass{
// Constructor
DatasetInfo::DatasetInfo(const T::FileNameType& filename1, const T::FileNameType& filename2){
    // Initialize TimeDomain time object through its constructor
    //time(filename1);

    // Initialize the number of groups to zero
    n_groups = 0;

    // Initialize the rest of the class
    read_from_file(filename2);

};

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


        // Save the first element (group belonging) into the corresponding variable

    }
};

void DatasetInfo::add_to_map_groups(T::GroupNameType& name_group, T::IndexType index_individual){
    T::MapType::iterator group_position = map_groups.find(name_group);
    if(group_position == map_groups.end()){
        map_groups[name_group].push_back(index_individual);
        n_groups += 1;
    }
    else{
        group_position->second.push_back(index_individual);
    }

}


} // end namespace