# ifndef PARALLELCOMPONENTS_HPP
# define PARALLELcOMPONENTS_HPP

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>

namespace TVSFCM{
using T = TypeTraits;

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


    // n_threads(n_threads_), chunk_size(chunk_size_), schedule_type(schedule_type_){};

protected:
    T::NumberType n_threads;                                //! Number of threads for the omp parallel version
    T::NumberType chunk_size;                               //! Number of chunk size for the for loop
    T::NumberType schedule_type;                            //! Id number for the for loop scheduling

    /**
     * Method for checking the existence of the file where the variables should be contained
     * @param filename1_ Name of the .txt file, whose existence in the directory needs to be checked
    */
    void check_filename(const T::FileNameType& filename1_) const;

    /**
     * Method for reading data from file and initializing the protected data structure 
     * @param filename1_ Name of the .txt file containing the data
    */
    void read_from_file(const T::FileNameType& filename1_);
};

} // end namespace

#endif // PARALLELCOMPONENTS_HPP