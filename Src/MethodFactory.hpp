#ifndef METHODFACTORY_HPP
#define METHODFACTORY_HPP

// Include header files
#include "ModelBase.hpp"
#include "ModelDerived.hpp"
#include "MyException.hpp"

// Include libraries
#include <memory>
#include <iostream>

namespace TVSFCM{
using T = TypeTraits;

/**
 * This method inserts in a FactoryType object the available Time-Varying Shared Frailty Cox Models, 
 * each represented by a couple composed of (numerical_id, name) of the model itself. 
 * @return Initialized FactoryType object
*/
T::FactoryType RegisterMethods() noexcept{
    T::FactoryType MethodFactory;

    MethodFactory.insert(std::pair(1, "AdaptedPaikeaM"));
    MethodFactory.insert(std::pair(2, "CSFM with Power Parameter"));
    MethodFactory.insert(std::pair(3, "Stochastic Time-Dependent CSFM"));

    return MethodFactory;
};

/**
 * This method controls that the id of the method is positive. 
 * If so, convert it from an integer to an unsigned integer. 
 * Otherwise, it throws an exception.
*/
T::IdType check_id_type(T::IntType id_) {
    if(id_ < 0)
        throw MyException("Provided negative id method.");

    return static_cast<T::IdType>(id_);
};

/**
 * This method prints the numerical id and the name of the time-vaying model stored in the FactoryType object.
 * 
 * If no models are contained, it throws an exception.
 * @param MethodFactory FactoryMethods object, already intialized
*/
void PrintMethods(const T::FactoryType & MethodFactory){
    if(MethodFactory.empty()){
        throw MyException("MethodFactory object is empty!");
    }

    std::cout << "Registered methods in the factory:" << std::endl;
    for (const auto & [id, name]: MethodFactory){
        std::cout << id << " " << name << std::endl;
    }
    std::cout << std::endl;
};

/**
 * This method returns a unique pointer to the base class, initialized with one of the time-varying models
 * present in the FactoryType object.
 * 
 * If the selected model does not exist, an exception is thrown.
 * 
 * @param id Numeric id of the model the user want to apply. 
 * @param filename1_ Name of the first file from which the parameters have to be extracted
 * @param filename2_ Name of the second file from which the parameters have to be extracted
 * @return Unique pointer to base class ModelBase, initialized with an object of the derived class, that could be AdaptedPaikeaM, 
 * CSFMwithPowerParameter and StochasticTimeDependentCSFM. If a wrong numeric id is provided, it throws an exception.
*/
std::unique_ptr<ModelBase> MakeLikelihoodModel(const T::IdType id, const T::FileNameType& filename1_, const T::FileNameType& filename2_) {
    switch(id){
        case 1:  return std::make_unique<AdaptedPaikeaM>(filename1_, filename2_);
        case 2:  return std::make_unique<CSFMwithPowerParameter>(filename1_, filename2_);
        case 3:  return std::make_unique<StochasticTimeDependentCSFM>(filename1_, filename2_);            
        default:  throw MyException("Not existent or not provided id method!");
    };
};

} // end namespace

#endif // METHODFACTORY_HPP



