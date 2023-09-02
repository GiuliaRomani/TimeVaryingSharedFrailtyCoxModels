# ifndef PARALLELCOMPONENTS_HPP
# define PARALLELCOMPONENTS_HPP

// Include header files
#include "TypeTraits.hpp"

// Include 
#include <iostream>

namespace TVSFCM{
using T = TypeTraits;

/**
 * ParallelComponents Class contains the necessary variables for the execution of the 
 * parallel version of the log-likelihood function. 
 */

class ParallelComponents{
public: 
    /**
     * Default constructor
    */
    ParallelComponents() = default;

    /**
     * Constructor
     * @param filename1_ Name of the file where the variables related to parallel methods are stored
    */
    explicit ParallelComponents(const T::FileNameType& filename1_);

    /**
     * Default destructor
    */
    ~ ParallelComponents() = default;

protected:
    T::NumberType n_threads;                                //! Number of threads for the omp parallel version
    T::NumberType chunk_size;                               //! Number of chunk size for the for loop
    T::IdType schedule_type;                                //! Id number for the for loop scheduling
    T::ScheduleType schedule_type_name;                     //! Name of the for loop scheduling

    /**
     * Method for reading data from file and initializing the protected data structures 
     * @param filename1_ Name of the .txt file containing the data
    */
    void read_from_file(const T::FileNameType& filename1_);

    /**
     * Method for checking that the parallel variables read from file are positive. Otherwise an exception is thrown.
     * @param n_threads_
     * @param chunk_size_
     * @param schedule_type_
    */
    void check_condition(T::IntType n_threads_, T::IntType chunk_size_, T::IntType schedule_type_);

    /**
     * Method for checking that the id of the schedule type exists, otherwise an exception is thrown.
    */
    void check_schedule_type() const;

    /**
     * Method for initializing the name of the adopted scheduling strategy
    */
    void set_schedule_type_name() noexcept;
};

} // end namespace

#endif // PARALLELCOMPONENTS_HPP
