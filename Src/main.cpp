// Include header files
#include "TypeTraits.hpp"
#include "MethodFactory.hpp"
#include "GetPot"
#include "MyException.hpp"

// Include libraries
#include <iostream>
#include <utility>
#include <chrono>
#include <iomanip>

// Type alias
using T = TVSFCM::TypeTraits;

int main(int argc, char *argv[]){

    // Check the number of input files
    try{
        if(argc != 3){
            throw TVSFCM::MyException("Not enough files provided.");  
        }

        // If there are two input .txt files, then save them
        T::FileNameType filename1 = argv[1];
        T::FileNameType filename2 = argv[2];

        // Time-Varying Shared Frailty Cox Model
        GetPot datafile(filename1.c_str());
        T::IntType id_ = datafile("Model/id_model", 0);
        T::IdType id = TVSFCM::check_id_type(id_);

        // Create the object factory map
        static T::FactoryType methods(TVSFCM::RegisteredMethods());

        // Call the desired method
        std::unique_ptr<TVSFCM::ModelBase> ptrMethod = TVSFCM::MakeLikelihoodMethod(id, filename1.c_str(), filename2.c_str());

        // Measure time elapsed
        const auto start = std::chrono::steady_clock::now();
        ptrMethod -> evaluate_loglikelihood();
        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<T::VariableType> elapsed_seconds = end - start;
        std::cout << std::endl;
        std::cout << "Elapsed time: " << resetiosflags(std::ios::scientific) << elapsed_seconds.count() << "s" << std::endl; 

    }
    catch(const TVSFCM::MyException& e){
        std::cout << e.what() << std::endl;
    };
    
    return 0;
}
