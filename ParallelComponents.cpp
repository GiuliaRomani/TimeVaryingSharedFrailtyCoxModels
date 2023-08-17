// Include header files
#include "ParallelComponents.hpp"
#include "GetPot"
#include "MyException.hpp"

// Include libraries

namespace TVSFCM{
    
// Constuctor
ParallelComponents::ParallelComponents(const T::FileNameType& filename1_){
    // Check the input file really exists
    check_filename(filename1_);

    // If it exists, extract and store the variables
    read_from_file(filename1_);
};

// Method for reading data from file
void ParallelComponents::read_from_file(const T::FileNameType& filename1_){
    GetPot datafile(filename1_.c_str());

    n_threads = static_cast<T::NumberType>(datafile("ParallelVersion/n_threads", 1));
    chunk_size = static_cast<T::NumberType>(datafile("ParallelVersion/chunk_size", 1));
    schedule_type = static_cast<T::NumberType>(datafile("ParallelVersion/id_schedule_type", 1));
};

// Method for checking the filename is correct and exits
void ParallelComponents::check_filename(const T::FileNameType& filename1_) const{
    std::ifstream check(filename1_);
    if(check.fail()){
        T::ExceptionType msg1 = "File ";
        T::ExceptionType msg2 = msg1.append((filename1_).c_str());
        T::ExceptionType msg3 = msg2.append(" does not exist.");
        throw MyException(msg3);
    }
};

} // end namespace