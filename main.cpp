#include "QuadraturePoints.hpp"
//#include "TimeDomain.hpp"
//#include "Parameters.hpp"
#include "TypeTraits.hpp"
//#include "Results.hpp"
#include "Dataset.hpp"

#include <iostream>
#include <utility>

int main(){
    using T = TypeTraits;
    /*
    // Prove Time class for errors
    TimeDomainInfo::TimeDomain time(".txt");
    std::cout << time.get_time_begin() << std::endl;
    std::cout << time.get_n_intervals() << std::endl;
    time.print_v_intervals();
    */
    
    /*
    // Prova Parameters class and Results class for errors 
    Params::Parameters params("Paik", 20, 8, 8, "DataToolFile.txt");
    T::VectorXdr & v_params = params.get_v_parameters();
    T::NumberType & n_params = params.get_n_parameters();
    T::VariableType opt_likelihoood = -10.0;

    ResultsMethod::Results results("Paik", n_params, v_params, opt_likelihoood);
    results.print_results();
    */

    // Prova Dataset class for errors
    DatasetInfoClass::DatasetInfo database("DataToolFile.txt", "DataIndividualsFile.txt");
    database.print_dataset();
    database.print_dataset_group();
    

    return 0;
}

