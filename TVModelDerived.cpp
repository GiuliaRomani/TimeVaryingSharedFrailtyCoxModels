
#include "TVModelDerived.hpp"

#include <cmath>

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
            parameters = Params::Parameters(filename1, n_parameters, n_intervals, n_regressors, 
                                            n_ranges_parameters, all_n_parameters);

            // Build the log-likelihood
            build_loglikelihood();

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


// Method for building the log-likelihood
void PowerParameterModel::build_loglikelihood(){

    // Implement the function ll_pp
    // T::VariableType x, T::IndexType index, 
    ll_pp = [this] (T::VectorXdr v_parameters_){
    T::VariableType log_likelihood = 0;
    T::MapType::iterator it_map = database.get_map_groups().begin();
    T::MapType::iterator it_map_end = database.get_map_groups().end();

    // For each group, compute the likelihood and then sum them
    for(it_map; it_map != it_map_end; ++it_map){
        // All the indexes in a group
        T::SharedPtrType indexes_group = it_map->second;

        T::VariableType log_likelihood_group = ll_group_pp(v_parameters_, indexes_group); //x, index, 
        log_likelihood += log_likelihood_group;
        //log_likelihood += (*indexes_group)[0];
    }

    // Subtract the constant term
    log_likelihood -= (n_groups/2)*log(M_PI);
    return log_likelihood;
    };

    
    // Implement the function ll_group_pp
    // T::VariableType x, T::IndexType index, 
    ll_group_pp = [this] (T::VectorXdr v_parameters_, T::SharedPtrType indexes_group_){

        // Extract the whole dataset as const reference
        const auto dataset = database.get_dataset();
        const auto dropout_intervals = database.get_dropout_intervals();
        const auto e_time = database.get_e_time();

        // Extract single parameters from the vector
        // v_parameters_(index) = x;

        T::VectorXdr phi = v_parameters_.head(n_intervals);                 // block(0,0,n_intervals,1);  
        T::VectorXdr betar = v_parameters_.block(n_intervals, 0, n_regressors,1);
        T::VectorXdr gammak(n_intervals);
        gammak(0) = 1.;
        gammak.block(1,0,n_intervals-1,1) = v_parameters_.block(n_intervals+n_regressors,0,n_intervals-1,1);

        T::VariableType sigma = v_parameters_(n_parameters - 1);
        sigma = sqrt(sigma);

        // Compute the necessary
        T::VariableType loglik1, partial1 = 0;
        for(const auto &i: *(indexes_group_)){
            for(T::IndexType k = 0; k < n_intervals; ++k){
                loglik1 += (dataset.row(i) * betar + phi(k)) * dropout_intervals(i,k);
                partial1 += dropout_intervals(i,k) * gammak(k);
            }
        }

        T::VariableType loglik2 = 0;
        for(T::IndexType q = 0; q < n_nodes; ++q){
            T::VariableType node = nodes[q];
            T::VariableType weight = weights[q];
            T::VariableType partial2 = 0;

            for(const auto &i: *(indexes_group_)){
                T::VariableType dataset_betar = dataset.row(i) * betar;
                for(T::IndexType k = 0; k < n_intervals; ++k){
                    partial2 += exp(dataset_betar + phi(k) + sqrt(2) * sigma * gammak(k) * node) * e_time(i,k);
                }
            }
            loglik2 += factor_c * weight * exp(sqrt(2) * sigma * node * partial1 - partial2);
        }
        loglik2 = log(loglik2);
        return (loglik1 + loglik2);
    };
};


// Method for executing the log-likelihood
T::VariableType PowerParameterModel::compute_loglikelihood(){
    T::VectorXdr& v_parameters = parameters.get_v_parameters();
    return ll_pp(v_parameters);
};



} // end namespace


