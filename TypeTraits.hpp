#ifndef TYPETRAITS_HPP
#define TYPETRAITS_HPP

#include <string>
#include <map>
#include <vector>

#include <Eigen/dense>

struct TypeTraits{
    using NumberType = unsigned int;                // Type used to indicate a positive integer number of elements
    using IdType = unsigned int;                    // Type for the numeric Id of a method
    using IdNameType = std::string;                 // Type for the name of a method
    using GroupNameType = std::string;                  // Type used to indicate the name of the group an individual belongs to

    using VariableType = double;                    // Type used to indicate the basic type variables
    using TimeType = double;                        // Type used for the time variables
    using IntervalType = VectorXd;                  // Type used for the time interval variable

    using VectorXd = Eigen::Vector<VariableType, Dynamic, 1>;                               // Type used for any dynamic vector
    using MatrixXd = Eigen::Matrix<VariableType, Dynamic, Dynamic, RowMajor>;               // Type used for any dynamic matrix

    using VectorType = std::vector<VariableType>;                                           // Type used for any vector containing variableType object
    using IndexType = unsigned int;                     // Type used to indicate any index of a matrix or of a vector
    using VectorIndexType = std::vector<IndexType>;     // Type used to indicate a vector of indexes


    // Type of the factory used. It is a map, where the first element indicates the numeric id of the method and the second element
    // is a string indicating the name of the method
    using FactoryType = std::map<IdType, IdNameType>;     

    // Type used for mapping the group name with the indexes of the individuals belonging to that group
    using MapType = std::map<GroupNameType, VectorIndexType>;


};


#endif //TYPETRAITS_HPP