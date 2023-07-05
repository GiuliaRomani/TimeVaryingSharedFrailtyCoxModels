#include "TypeTraits.hpp"
#include "MethodFactory.hpp"
//#include "QuadraturePoints.hpp"
//#include "TimeDomain.hpp"
//#include "Parameters.hpp"
//#include "Results.hpp"
//#include "Dataset.hpp"
//#include "TVModelBase.hpp"
//#include "TVModelDerived.hpp"

#include <iostream>
#include <utility>

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
    T::VectorNumberType all_n_parameters{3,4,2,1};
    Params::Parameters params("DataToolFile.txt", 4, 10, 3, 4, all_n_parameters);
    T::VectorXdr & v_params = params.get_v_parameters();
    T::NumberType & n_params = params.get_n_parameters();
    T::VariableType opt_likelihoood = -10.0;
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
    

    
    // Prove the PowerParameterModel
    static T::FactoryType methods(RegisteredMethods());
    PrintMethods(methods);
    T::IdType id = 2;

    //std::unique_ptr<TVModel::ModelBase> ptrMethod = MakeLikelihoodMethod(id, "DataToolFile.txt", "DataIndividualsFile.txt");
    std::unique_ptr<TVModel::ModelBase> ptrMethod = MakeLikelihoodMethod(id, "DataToolFile.txt", "DatasetYear2018.txt");
    ptrMethod -> optimize_loglikelihood();
    // ptrMethod -> print_extract_parameters();
    
    
    
    return 0;
}

