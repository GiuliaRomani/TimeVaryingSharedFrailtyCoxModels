
#include "TVModelDerived.hpp"

namespace TVModel{
using T = TypeTraits;

// Constructor
PowerParameterModel::PowerParameterModel(const T::FileNameType& filename1, const T::FileNameType& filename2):
        // Base constructor for base class
        ModelBase(filename1, filename2) {
            // Initialize the number of parameters
            n_parameters = compute_n_parameters();

            // Initialize the vector of number of parameters
            all_n_parameters = {n_intervals, n_regressors, n_intervals - 1, 1};

            // Construct the class parameters
            parameters = Params::Parameters(filename1, n_ranges_parameters, n_parameters, 
                                        n_intervals, n_regressors, all_n_parameters);

};
        
// Virtual method for computing the number of parameters
T::NumberType PowerParameterModel::compute_n_parameters() {
    return (2*n_intervals + n_regressors);
};

/*
// Virtual method to derive the interval variance of the frailty
void PowerParameterModel::compute_sd_frailty() {
    // Extract the useful parameter from the vector of parameters (i.e. gamma vector and sigma)
    T::VectorXdr& v_parameters = parameters.get_v_parameters();
    T::VariableType sigma = v_parameters.tail(1);

    T::IndexType begin_gamma = n_intervals + n_regressors;
    T::IndexType end_gamma = 2 * n_intervals + n_regressors - 2;
    T::VectorXdr v_gamma(n_intervals);
    v_gamma(0) = 1;
    v_gamma.block(1,0,n_intervals-1,1) = v_parameters.block(begin_gamma,0, n_intervals-1,1);

    for(const auto& p: v_gamma){
        variance_frailty(p) = (sigma * v_gamma(p)) * (sigma * v_gamma(p));
        sd_frailty(p) = sigma * v_gamma(p);
    };
};
*/

// Implement the function ll_pp
auto ll_pp = [this] (T::VariableType x, T::IndexType index, T::VectorXdr v_parameters_){
    T::VariableType log_likelihood = 0;
    T::MapType::iterator it_map = database.get_map_groups().begin();
    T::MapType::iterator it_map_end = database.get_map_groups().end();

    // For each group, compute the likelihood and then sum them
    for(auto it_map; it_map != it_map_end; ++it_map){
        // All the indexes in a group
        T::SharedPtrType indexes_group = it_map->second;

        //log_likelihood_group = ll_group_pp(x, index, v_parameters_, indexes_group);
        log_likelihood += log_likelihood_group;
    }

    // Subtract the constant term
    log_likelihood -= (n_groups/2)*log(M_PI);
    return log_likelihood;
}



} // end namespace


