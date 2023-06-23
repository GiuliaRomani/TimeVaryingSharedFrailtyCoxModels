// Include header files
#include "TVModelDerived.hpp"

// Include libraries
#include <cmath>
#include <tuple>

namespace TVModel{
using T = TypeTraits;

//-------------------------------------------------------------------------------------------------------------
// Implementations for the PoweParameter model

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


T::TuplePPType PowerParameterModel::extract_parameters(T::VectorXdr v_parameters_){
    // Extract parameters from the vector
    T::VectorXdr phi = v_parameters_.head(n_intervals);                 // block(0,0,n_intervals,1);  
    T::VectorXdr betar = v_parameters_.block(n_intervals, 0, n_regressors,1);
    T::VectorXdr gammak(n_intervals);
    gammak(0) = 1.;
    gammak.block(1,0,n_intervals-1,1) = v_parameters_.block(n_intervals+n_regressors,0,n_intervals-1,1);

    T::VariableType sigma = v_parameters_(n_parameters - 1);
    sigma = sqrt(sigma);

    return std::make_tuple(phi, betar, gammak, sigma);
};


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
        auto [phi, betar, gammak, sigma] = extract_parameters(v_parameters_);

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
void PowerParameterModel::optimize_loglikelihood(){
    T::VectorXdr& v_parameters = parameters.get_v_parameters();
    T::VariableType optimal_ll_pp = ll_pp(v_parameters);

    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_parameters, optimal_ll_pp);
    result.print_results();
};

/*
//-------------------------------------------------------------------------------------------------------
// Implementations for the Paik Model
// Constructor
PaikModel::PaikModel(const T::FileNameType& filename1, const T::FileNameType& filename2):
        // Base constructor for base class
        ModelBase(filename1, filename2) {
            // Initialize the number of parameters
            n_parameters = compute_n_parameters();

            // Initialize the vector of number of parameters
            all_n_parameters = {n_intervals, n_regressors, 1, 1, n_intervals};

            // Construct the class parameters
            parameters = Params::Parameters(filename1, n_parameters, n_intervals, n_regressors, 
                                            n_ranges_parameters, all_n_parameters);

            // Build the log-likelihood
            build_loglikelihood();

};
        
// Virtual method for computing the number of parameters
T::NumberType PaikModel::compute_n_parameters() {
    return (2*n_intervals + n_regressors + 2);
};

// Method for building the log-likelihood
void PaikModel::build_loglikelihood(){
    // Implement the function ll_paik
    ll_paik = [this] (T::VectorXdr v_parameters_){              // T::VariableType x, T::IndexType index, 
        T::VariableType log_likelihood = 0;
        T::MapType::iterator it_map = database.get_map_groups().begin();
        T::MapType::iterator it_map_end = database.get_map_groups().end();

        // For each group, compute the likelihood and then sum them
        for(it_map; it_map != it_map_end; ++it_map){
            // All the indexes in a group
            T::SharedPtrType indexes_group = it_map->second;

            T::VariableType log_likelihood_group = ll_group_paik(v_parameters_, indexes_group); //x, index, 
            log_likelihood += log_likelihood_group;
            //log_likelihood += (*indexes_group)[0];
        }

        // Subtract the constant term
        log_likelihood -= (n_groups/2)*log(M_PI);
        return log_likelihood;
    };

    
    // Implement the function ll_group_pp
    ll_group_paik = [this] (T::VectorXdr v_parameters_, T::SharedPtrType indexes_group_){     // T::VariableType x, T::IndexType index, 
        // Extract the whole dataset as const reference
        const auto dataset = database.get_dataset();
        const auto dropout_intervals = database.get_dropout_intervals();
        const auto e_time = database.get_e_time();

        // Extract single parameters from the vector
        // v_parameters_(index) = x;
        T::VectorXdr phi = v_parameters_.head(n_intervals);                 // block(0,0,n_intervals,1);  
        T::VectorXdr betar = v_parameters_.block(n_intervals, 0, n_regressors,1);
        T::VariableType mu1 = v_parameters_(n_intervals + n_regressors);
        T::VariableType mu2 = 1 - mu1;
        T::VariableType nu = v_parameters_(n_intervals + n_regressors + 1);
        T::VectorXdr gammak = v_parameters_.tail(n_intervals);

        // Compute some necessary components
        T::NumberType n_individuals = (*indexes_group_).size();
        T::MatrixXdr A_ijk(n_individuals, n_intervals);
        T::VariableType A_i;
        T::MatrixXdr A_ik(1,n_intervals);
        T::VectorXdr d_ik(n_intervals);

        // Compute the first component of the likelihood
        T::VariableType loglik1 = 0;
        for(const auto &i: *(indexes_group_)){
            T::VariableType dataset_betar = dataset.row(i) * betar;
            for(T::IndexType k = 0; k < n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * dropout_intervals(i,k);
                A_ijk(i,k) = e_time(i,k) * exp(dataset_betar + phi(k));
            }
        }
        A_ik = A_ijk.colwise().sum();
        //A_i = A_ik.colwise().sum();
        loglik1 -= (mu1/nu) * log(1 + nu * A_i);

        // Compute the second line of the formula
        T::VariableType loglik2 = 0;
        for(T::IndexType k = 0; k < n_intervals; ++k){
            loglik2 -= (mu2/gammak(k)) * log(1 + gammak(k) * A_ik(k));
        }

        // Compute the third line of the formula
        T::VariableType loglik3 = 0;
        for(T::IndexType k = 0; k < n_intervals; ++k){
            T::VariableType loglik4 = 0;
            T::NumberType d_ik_size = d_ik(k);
            for(T::IndexType l = 0; l <= d_ik_size; ++l){
                T::VariableType coeff_bin = 2; // TODO
                T::VariableType tgamma1 = tgamma(d_ik(k) + mu2/gammak(k) - l);
                T::VariableType tgamma2 = tgamma(mu1/nu + l);
                T::VariableType tgamma3 = tgamma(mu2/gammak(k));
                T::VariableType tgamma4 = tgamma(mu1/nu);
                T::VariableType pow1 = pow(A_ik(k) + 1/gammak(k), l - d_ik(k));
                T::VariableType pow2 = pow(A_i + 1/nu, l);
                loglik4 += coeff_bin * tgamma1 * tgamma2 * pow1 * factor_c/ (tgamma3 * tgamma4 * pow2);                 
            }
            loglik3 += log(loglik4);
        }
    return (loglik1 + loglik2 + loglik3);
    };
};


// Method for executing the log-likelihood
void PowerParameterModel::optimize_loglikelihood(){
    T::VectorXdr& v_parameters = parameters.get_v_parameters();
    T::VariableType optimal_ll_pp = ll_pp(v_parameters);

    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_parameters, optimal_ll_pp);
    result.print_results();
};

*/



} // end namespace


