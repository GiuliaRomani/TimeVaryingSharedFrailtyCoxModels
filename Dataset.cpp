#include "Dataset.hpp"
#include "GetPot"

#include <sstream>

namespace DatasetInfoClass{
// Constructor
DatasetInfo::DatasetInfo(T::FileNameType& filename1, T::FileNameType& filename2){
    // Initialize TimeDomain time object through its constructor
    time(filename1);

    // Initialize the rest of the class
    read_from_file(filename2);

};

void DatasetInfo::read_from_file(T::FileNameType& filename2){
    T::CheckType first_line = false;
    T::NumberType i = 0;
    T::NumberType j = 0;
    std::ifstream filename(filename2.c_str());

    while(!filename.fail() & !filename.eof()){
        std::string line;
        getline(filename, line);
        istringstream value_line(line);

        if(!first_line){
            first_line >> n_individuals >> n_regressors;
            first_line = true;

            dataset_group.resize(n_individuals);
            dataset.resize(n_individuals, n_regressors);
        }
        else{
            getline(value_line, dataset_group(i), " ");
            istringstream value_regressor;
            while(getline(value_line, value_regressor, " ")){
                dataset(i,j) = value_regressor;
                j++;
            }
            i++;
        }


        // Save the first element (group belonging) into the corresponding variable

    }

};


} // end namespace