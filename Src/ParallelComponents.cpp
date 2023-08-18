// Include header files
#include "ParallelComponents.hpp"
#include "GetPot"
#include "MyException.hpp"

// Include libraries

namespace TVSFCM{
    
// Constuctor
ParallelComponents::ParallelComponents(const T::FileNameType& filename1_){
    // If it exists, extract and store the variables
    read_from_file(filename1_);

    // Check that the id of the schedule type exists
    check_schedule_type();

    // Set the name of the schedule_type
    set_schedule_type_name();
};

// Method for reading data from file and check they are positive
void ParallelComponents::read_from_file(const T::FileNameType& filename1_){
    GetPot datafile(filename1_.c_str());

    T::IntType n_threads_, chunk_size_, schedule_type_;

    n_threads_ = datafile("ParallelVersion/n_threads", 1);
    chunk_size_ = datafile("ParallelVersion/chunk_size", 1);
    schedule_type_ = datafile("ParallelVersion/id_schedule_type", 1); 

    // Check the input values are not negative
    check_condition(n_threads_, chunk_size_, schedule_type_);
};

// Method for checking the parallel input variables are positive
void ParallelComponents::check_condition(T::IntType n_threads_, T::IntType chunk_size_, T::IntType schedule_type_){
    if(n_threads_ < 0)
        throw MyException("Provided negative value for number of threads.");
    else if(chunk_size_ < 0)
        throw MyException("Provided negative value for chunk size.");
    else if(schedule_type_ < 0)
        throw MyException("Provided negative value for the type of schedule");

    n_threads = static_cast<T::NumberType>(n_threads_);
    chunk_size = static_cast<T::NumberType>(chunk_size_);
    schedule_type = static_cast<T::NumberType>(schedule_type_);
};

// Method for checking the schedule type exists
void ParallelComponents::check_schedule_type() const{
    if(schedule_type < 1 or schedule_type > 4){
        T::ExceptionType msg1 = "Schedule id ";
        T::ExceptionType msg2 = msg1.append(std::to_string(schedule_type));
        T::ExceptionType msg3 = msg2.append(" does not exist.");
        throw MyException(msg3);
    }
};

// Method for defining the name associated to the id_schedule
void ParallelComponents::set_schedule_type_name(){
    if(schedule_type == 1)
        schedule_type_name = "static";
    else if(schedule_type == 2)
        schedule_type_name = "dynamic";
    else if(schedule_type == 3)
        schedule_type_name = "guided";
    else   
        schedule_type_name = "auto";
};

} // end namespace