#ifndef TYPETRAITS_HPP
#define TYPETRAITS_HPP

#include <string>
#include <map>
#include <vector>
#include <memory>

#include <Eigen/Dense>

struct TypeTraits{
    // Basic types
    using VariableType = double;                    // Type used to indicate the basic type variables
    using NumberType = unsigned int;                // Type used to indicate a positive integer number of elements
    using IntType = int;                            // Type used for an integer elements
    using IndexType = unsigned int;                 // Type used to indicate any index of a matrix or of a vector
    using SizeType = std::size_t;                    // Type used for a size
    using CheckType = bool;                        // Type used to indicate any check condition

    // Type for file name
    using FileNameType = std::string;               // Type used to indicate the name of the file passed to constructors

    // Types for methods
    using IdType = unsigned int;                    // Type for the numeric Id of a method
    using IdNameType = std::string;                 // Type for the name of a method
    
    // Type for groups
    using GroupNameType = std::string;              // Type used to indicate the name of the group an individual belongs to

    // Basic containers types
    using VectorXdr = Eigen::Matrix<VariableType, Eigen::Dynamic, 1>;                                        // Type used for any dynamic vector
    using MatrixXdr = Eigen::Matrix<VariableType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;         // Type used for any dynamic matrix
    using VectorType = std::vector<VariableType>;                               // Type used for any vector of VariableType type
    using VectorIndexType = std::vector<IndexType>;                             // Type used to for any vector of IndexType type
    using VectorXdrGroupType = Eigen::Matrix<GroupNameType, Eigen::Dynamic, 1>;    // Type used for collecting all the group belonging


    // Type for time variables
    using TimeType = double;                        // Type used for the time variables
    using TimeIntervalType = VectorXdr;            // Type used for the time interval variable

    // Other containers types
    // Type of the factory used. It is a map, where the first element indicates the numeric id of the method and the second element
    // is a string indicating the name of the method
    using FactoryType = std::map<IdType, IdNameType>;     

    // Type used for mapping the group name with the indexes of the individuals belonging to that group
    using MapType = std::map<GroupNameType, std::unique_ptr<VectorIndexType>>;


};


#endif //TYPETRAITS_HPP