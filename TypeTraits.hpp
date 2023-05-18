#ifndef TYPETRAITS_HPP
#define TYPETRAITS_HPP

#include <string>
#include <map>

#include <Eigen/dense>

struct TypeTraits{
    using NumberType = unsigned int;                // Type used to indicate a positive integer number of elements
    using IdType = unsigned int;                    // Type for the numeric Id of a method
    using IdNameType = std::string;                 // Type for the name of a method

    using VariableType = double;                    // Type used to indicate the basic type variables
    using TimeType = double;                        // Type used for the time variables
    using IntervalType = Eigen::VectorXd;           // Type used for the time interval variable

    using VectorXd = Eigen::VectorXd;               // Type used for any dynamic vector
    using MatrixXd = Eigen::MatrixXd;               // Type used for any dynamic matrix


    // Type of the factory used. It is a map, where the first element indicates the numeric id of the method and the second element
    // is a string indicating the name of the method
    using FactoryType = std::map<IdType, IdNameType>       


};


#endif //TYPETRAITS_HPP