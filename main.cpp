#include "TypeTraits.hpp"
#include "MethodFactory.hpp"
#include "GetPot"
//#include "QuadraturePoints.hpp"
//#include "TimeDomain.hpp"
//#include "Parameters.hpp"
//#include "Results.hpp"
//#include "Dataset.hpp"
//#include "TVModelBase.hpp"
//#include "TVModelDerived.hpp"

#include <iostream>
#include <utility>
#include <chrono>

int main(){
    using T = TypeTraits;

    /*
    // Prove QuadraturePoints formula
    QuadraturePoints::Points9 points9;
    std::array<T::VariableType, 9>& nodes = points9.nodes;
    for(const auto& n: nodes)
        std::cout << n << std::endl;
    */
    
    
    // Prove Time class for errors
    // std::cout << "Direct initialization" << std::endl;
    // TimeDomainInfo::TimeDomain time("DataToolFile.txt");
    // time.print_v_intervals();
    // time.print_n_intervals();
    
    
    
    /*
    std::cout << "Using a different initialization" << std::endl;
    TimeDomainInfo::TimeDomain time_alternative(TimeDomainInfo::TimeDomain("DataToolFile.txt"));
    time.print_v_intervals();
    */
    
    
    /*
    // Prova Parameters class and Results class for errors 
    // T::VectorNumberType all_n_parameters{8,5,7,1};
    T::VectorNumberType all_n_parameters{8, 5, 1, 1, 8};
    Params::Parameters params("DataToolFile.txt", 23, 8, 5, 5, all_n_parameters);
    T::VectorXdr & v_params = params.get_v_parameters();
    T::NumberType & n_params = params.get_n_parameters();
    std::cout << "Number of parameters: "<< n_params << std::endl;
    std::cout << v_params << std::endl;
    */
    
    

    /*
    ResultsMethod::Results results("Paik", n_params, v_params, opt_likelihoood);
    results.print_results();
    */

    /*
    // Prova Dataset class for errors
    // DatasetInfoClass::DatasetInfo database("DataIndividualsFile.txt", time.get_n_intervals(), time.get_v_intervals());
    DatasetInfoClass::DatasetInfo database("DatasetYear2018.txt", time.get_n_intervals(), time.get_v_intervals());
    //database.print_dataset();
    database.print_n_regressors();
    database.print_n_individuals();
    database.print_n_groups();
    //database.print_dataset_group();
    //database.print_map_groups();
    //database.print_dropout_intervals();
    //database.print_individuals_group("EngC");
    //database.print_e_time();
    */

    /*
    // Prova TVModelBase fo errors
    TVModel::ModelBase modelbase("DataToolFile.txt", "DataIndividualsFile.txt");
    modelbase.print_map_groups();
    modelbase.print_n_regressors();
    modelbase.print_n_intervals();
    */
    
    // Time-Varying Shared Frailty Cox Model
    GetPot datafile("DataToolFile.txt");
    T::IdType id = datafile("Model/id_model",2);

    static T::FactoryType methods(RegisteredMethods());
    PrintMethods(methods);

    try{
    //std::unique_ptr<TVModel::ModelBase> ptrMethod = MakeLikelihoodMethod(id, "DataToolFile.txt", "DataIndividualsFile.txt");
    std::unique_ptr<TVModel::ModelBase> ptrMethod = MakeLikelihoodMethod(id, "DataToolFile.txt", "DatasetYear2018.txt");

    // Measure time elapsed
    const auto start = std::chrono::steady_clock::now();

    ptrMethod -> evaluate_loglikelihood();

    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<T::VariableType> elapsed_seconds = end - start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s" << std::endl; 
    }
    catch(const std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    };
    
    
    
    
    
    return 0;
}

