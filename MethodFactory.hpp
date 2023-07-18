
#ifndef METHODFACTORY_HPP
#define METHODFACTORY_HPP

// Include header files
#include "TVModelBase.hpp"
#include "TVModelDerived.hpp"
#include "MyException.hpp"

// Include libraries
#include <memory>
#include <iostream>


using T = TypeTraits;

/**
 * This method insert in a FactoryType object the available time-varying models, each represented by a couple composed of (numerical_id, name) of the model itself. 
 * @return Initialized FactoryType object
*/
T::FactoryType RegisteredMethods() {
    T::FactoryType mapMethod;

    mapMethod.insert(std::pair(1, "Paik"));
    mapMethod.insert(std::pair(2, "Power Parameter"));
    mapMethod.insert(std::pair(3, "Stochastic Frailty"));

    return mapMethod;
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
std::unique_ptr<TVModel::ModelBase> MakeLikelihoodMethod(const T::IdType id, Args && ... args){
    switch(id){
        case 1:  return std::make_unique<TVModel::PaikModel>(std::forward<Args>(args) ...);
        case 2:  return std::make_unique<TVModel::PowerParameterModel>(std::forward<Args>(args) ...);
        case 3:  return std::make_unique<TVModel::LogFrailtyModel>(std::forward<Args>(args) ...);            
        default:  throw MyException("Not existent id method!");
    };
};


#endif // METHODFACTORY_HPP



