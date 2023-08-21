#ifndef METHODFACTORY_HPP
#define METHODFACTORY_HPP

// Include header files
#include "TVModelBase.hpp"
#include "TVModelDerived.hpp"
#include "MyException.hpp"

// Include libraries
#include <memory>
#include <iostream>

namespace TVSFCM{
using T = TypeTraits;

/**
 * This method insert in a FactoryType object the available time-varying models, each represented by a couple composed of (numerical_id, name) of the model itself. 
 * @return Initialized FactoryType object
*/
T::FactoryType RegisteredMethods() noexcept{
    T::FactoryType mapMethod;

    mapMethod.insert(std::pair(1, "AdaptedPaikeaM"));
    mapMethod.insert(std::pair(2, "CSFM with Power Parameter"));
    mapMethod.insert(std::pair(3, "Stochastic Time-Dependent CSFM"));

    return mapMethod;
};

// Check that the id of the method is positive and then convert it
T::IdType check_id_type(T::IntType id_) {
    if(id_ < 0)
        throw MyException("Provided negative id method.");

    return static_cast<T::IdType>(id_);
};

/**
 * This method print the numerical id and the name of the time-vaying model stored in the FactoryType object.
 * 
 * If no models are contained, it throws an exception.
 * @param FactoryMethods FactoryMethods object, already intialized
*/
void PrintMethods(const T::FactoryType & FactoryMethods){
    if(FactoryMethods.empty()){
        throw MyException("FactoryMethods object is empty!");
    }

    std::cout << "Registered methods in the factory:" << std::endl;
    for (const auto & [id, name]: FactoryMethods){
        std::cout << id << " " << name << std::endl;
    }
    std::cout << std::endl;
};

/**
 * This method return a unique pointer to one of the time-varying model present in the FactoryType object.
 * 
 * If the selected model does not exist, an exception is thrown.
 * 
 * @param id Numeric id of the model the user want to apply. 
 * @param args Variadic template for the arguments of the model call.
 * @return Unique pointer to base class TVModel::ModelBase, initialized with an object of the derived class, that could be TVModel::PaikModel, 
 * TVModel::PowerParameterModel, TVModel::LogFrailtyModel. If a wrong numeric id is provided, it throws an exception.
*/
template<class ... Args>
std::unique_ptr<ModelBase> MakeLikelihoodMethod(const T::IdType id, Args && ... args) {
    switch(id){
        case 1:  return std::make_unique<AdaptedPaikeaM>(std::forward<Args>(args) ...);
        case 2:  return std::make_unique<CSFMwithPowerParameter>(std::forward<Args>(args) ...);
        case 3:  return std::make_unique<StochasticTimeDependentCSFM>(std::forward<Args>(args) ...);            
        default:  throw MyException("Not existent or not provided id method!");
    };
};

} // end namespace

#endif // METHODFACTORY_HPP



