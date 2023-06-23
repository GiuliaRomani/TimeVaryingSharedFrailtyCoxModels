
#ifndef METHODFACTORY_HPP
#define METHODFACTORY_HPP

// Include header files
#include "TVModelBase.hpp"
#include "TVModelDerived.hpp"

// Include libraries
#include <memory>
#include <iostream>


using T = TypeTraits;

// This function simply inserts the available methods in a FactoryType objet and return the object created
T::FactoryType RegisteredMethods() {
    T::FactoryType mapMethod;

    mapMethod.insert(std::pair(1, "Paik"));
    mapMethod.insert(std::pair(2, "Power Parameter"));
    mapMethod.insert(std::pair(3, "Stochastic Frailty"));

    return mapMethod;
};

// For each method available in the FactoryType object, its id and name are printed.
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
template<class ... Args>
std::unique_ptr<TVModel::ModelBase> MakeLikelihoodMethod(const T::IdType id, Args && ... args){
    switch(id){
        case 1:  return std::make_unique<TVModel::PaikModel>(std::forward<Args>(args) ...);
        case 2:  return std::make_unique<TVModel::PowerParameterModel>(std::forward<Args>(args) ...);
        //case 3:  return std::make_unique<StochasticFrailty>(std::forward<Args>(args) ...);            
        default:  throw "Not existent id method!";
    };
};


#endif // METHODFACTORY_HPP



