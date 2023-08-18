#ifndef TYPETRAITS_HPP
#define TYPETRAITS_HPP

// Include libraries
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <tuple>
#include <Eigen/Dense>

/**
 * Struct containing the type aliases used in the project.
*/

namespace TVSFCM{

struct TypeTraits{
    using VariableType = double;                    //! Type used to indicate the basic variable
    using NumberType = unsigned int;                //! Type used to indicate a positive integer number of elements
    using IndexType = unsigned int;                 //! Type used to indicate any index of a matrix or of a vector
    using CheckType = bool;                         //! Type used for returning the status of a check/control

    using FileNameType = std::string;               //! Type used to indicate the name of the file passed to constructors

    using ExceptionType = std::string;              //! Type of the content of an exception

    using ScheduleType = std::string;               //! Type indicating the schedule type used in the OpenMP for loop

    using IdType = unsigned int;                    //! Type of the numeric Id of a method
    using IdNameType = std::string;                 //! Type of the name of a method
    
    using GroupNameType = std::string;              //! Type used to indicate the name of the clusters/groups in which the individuals are divided

    using VectorXdr = Eigen::Matrix<VariableType, Eigen::Dynamic, 1>;                      //! Type used for any dynamic vector
    using MatrixXdr = Eigen::Matrix<VariableType, Eigen::Dynamic, Eigen::Dynamic>;         //! Type used for any dynamic matrix Eigen::RowMajor
    using VectorXdrGroupType = Eigen::Matrix<GroupNameType, Eigen::Dynamic, 1>;            //! Type used for collecting all the groups
    using MappedVectorType = Eigen::Map<VectorXdr>; 

    using VectorType = std::vector<VariableType>;                                           //! Type used for any vector of VariableType type
    using VectorIndexType = std::vector<IndexType>;                                         //! Type used to for any vector of IndexType type
    using VectorNumberType = std::vector<NumberType>;                                       //! Type used for storing integers


    /**
     * Tye of the factory. It is a map, where the first element indicates the numeric id of one model and the second element its name. 
    */
    using FactoryType = std::map<IdType, IdNameType>;     

    /**
     * Type used to indicate a shared pointed to a vector of indexes, where this vector containes
     * the indexes of the dataset, of the individuals belonging to a cluster.
    */
    using SharedPtrType = std::shared_ptr<VectorIndexType>;

    /**
     * Type used to define a map, where the first element is the name of the group and the second 
     * element is a shared pointer to the individuals belonging to that cluster
    */
    using MapType = std::map<GroupNameType, SharedPtrType>;

    /**
     * Type used to indicate a tuple associated to the "CSFM with Power Parameter". It contains the type of the different 
     * categories that constitutes the vector of parameters.
    */
    using TuplePPType = std::tuple<VectorXdr, VectorXdr, VectorXdr, VariableType>;

    /**
     * Type used to indicate a tuple associated to the "Adapted Paik et al.'s model". It contains the type of the different 
     * categories that constitutes the vector of parameters.
    */
    using TuplePaikType = std::tuple<VectorXdr, VectorXdr, VariableType, VariableType, VariableType, VectorXdr>;

    /**
     * Type used to indicate a tuple associated to the "Stochastic Time-Dependent CSFM ". It contains the type of the different 
     * categories that constitutes the vector of parameters.
    */
    using TupleLFType = std::tuple<VectorXdr, VectorXdr, VariableType, VariableType, VariableType, VariableType, VariableType>;

    /**
     * Type used to indicate the tuple for the dropout variables
    */
    using TupleDropoutType = std::tuple<MatrixXdr, VectorXdr, VariableType>;

    /**
     * Type used to indicate the tuple for the variables of the dataset
    */
    using TupleMatrixAType = std::tuple<MatrixXdr, VectorXdr, VariableType>;

};

} // end namespace


#endif //TYPETRAITS_HPP
