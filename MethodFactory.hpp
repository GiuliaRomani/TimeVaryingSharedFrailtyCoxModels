
#ifndef METHODFACTORY_HPP
#define METHODFACTORY_HPP

// Include header files
#include "TVModelBase.hpp"
#include "TVModelDerived.hpp"

// Include libraries
#include <memory>
#include <iostream>
#include <exception>


using T = TypeTraits;

/**
 * This method insert in a FactoryType object the available time-varying models, each represented by a couple composed of (numerical_id, name) of the model itself. 
 * 
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
 * This method print the numerical id and the name of the time-vaying model stored in the FactoryType object
 * 
 * @param FactoryMethods reference to const FactoryType object
*/
void PrintMethods(const T::FactoryType & FactoryMethods){
    std::cout << "Registered methods in the factory:" << std::endl;
    for (const auto & [id, name]: FactoryMethods){
        std::cout << id << " " << name << std::endl;
    }
    std::cout << std::endl;
};


// This function receives the id of the method the user wants to use and some parameter (through a variadic template)
// and returns a unique pointer to the base class, but initialised with a derived class object with the arguments passed.
// A valid id is controlled in the main() function. 
/**
 * This method return a unique pointer to one of the time-varying model present in the FactoryType object.
 * 
 * @param id Numeric id of the model the user want to apply. If the correct id is not provided, an exception is thrown.
 * @param args Variadic template for the arguments of the model call.
 * 
 * @return Unique pointer to base class TVModel::ModelBase, initialized with an object of the derived class, that could be TVModel::PaikModel, 
 * TVModel::PowerParameterModel, TVModel::LogFrailtyModel. If a wrong numeric id is provided, it throws an exception.
*/
template<class ... Args>
std::unique_ptr<TVModel::ModelBase> MakeLikelihoodMethod(const T::IdType id, Args && ... args){
    switch(id){
        case 1:  return std::make_unique<TVModel::PaikModel>(std::forward<Args>(args) ...);
        case 2:  return std::make_unique<TVModel::PowerParameterModel>(std::forward<Args>(args) ...);
        case 3:  return std::make_unique<TVModel::LogFrailtyModel>(std::forward<Args>(args) ...);            
        default:  throw std::runtime_error("Not existent id method!");
    };
};


#endif // METHODFACTORY_HPP



