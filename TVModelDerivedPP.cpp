
// Include header files
#include "TVModelDerived.hpp"

// Include libraries
#include <cmath>
#include <tuple>
#include <algorithm>
#include <random>

namespace TVModel{
using T = TypeTraits;

//-------------------------------------------------------------------------------------------------------------
// Implementations for the PoweParameter model

// Constructor
PowerParameterModel::PowerParameterModel(const T::FileNameType& filename1, const T::FileNameType& filename2):
        // Constructor for base class
        ModelBase(filename1, filename2) {
            // Constructor for current class
            
            // Initialize the number of parameters
            compute_n_parameters();
            
            // Resize standard error vector and hessian diagonal according to the number of parameters
            hessian_diag.resize(n_parameters);
            se.resize(n_parameters);

            // Initialize the vector of number of parameters
            all_n_parameters = {Dataset::n_intervals, Dataset::n_regressors, Dataset::n_intervals - 1, 1};

            // Construct the class parameters
            parameters = Params::Parameters(filename1, n_parameters, Dataset::n_intervals, Dataset::n_regressors, 
                                            n_ranges_parameters, all_n_parameters);

            // Build the log-likelihood and the one-directional second derivative
            build_loglikelihood();
            build_dd_loglikelihood();
};
        
// Virtual method for computing the number of parameters
void PowerParameterModel::compute_n_parameters() {
    n_parameters = (2 * Dataset::n_intervals + Dataset::n_regressors);
};

// Method for extracting the parameters from the vector of parameters
T::TuplePPType PowerParameterModel::extract_parameters(const T::VectorXdr& v_parameters_) const{
    // Extract parameters from the vector
    T::VectorXdr phi = v_parameters_.head(Dataset::n_intervals);                 // block(0,0,n_intervals,1);  
    T::VectorXdr betar = v_parameters_.block(Dataset::n_intervals, 0, Dataset::n_regressors, 1);
    T::VectorXdr gammak(Dataset::n_intervals);
    gammak(0) = 1.;
    gammak.block(1,0,Dataset::n_intervals-1,1) = v_parameters_.block(Dataset::n_intervals + Dataset::n_regressors, 0, Dataset::n_intervals - 1, 1);

    T::VariableType sigma = v_parameters_(n_parameters - 1);
    sigma = sqrt(sigma);

    return std::make_tuple(phi, betar, gammak, sigma);
};


// Method for building the log-likelihood
void PowerParameterModel::build_loglikelihood(){
    // Implement the function ll_pp

    ll_pp = [this] (T::VectorXdr& v_parameters_){
        T::VariableType log_likelihood_group, log_likelihood = 0;
        T::SharedPtrType indexes_group = nullptr;

        // For each group, compute the likelihood and then sum them
        T::MapType::iterator it_map = Dataset::map_groups.begin();
        T::MapType::iterator it_map_end = Dataset::map_groups.end();
        for(; it_map != it_map_end; ++it_map){
            // All the indexes of a group
            indexes_group = it_map->second;

            log_likelihood_group = ll_group_pp(v_parameters_, indexes_group); 
            log_likelihood += log_likelihood_group;

            indexes_group = nullptr;
        }

        // Subtract the constant term
        log_likelihood -= ((Dataset::n_groups)/2)*log(M_PI);
        return log_likelihood;
    };

    
    // Implement the function ll_group_pp
    ll_group_pp = [this] (T::VectorXdr& v_parameters_, T::SharedPtrType indexes_group_){

        // Extract single parameters from the vector
        auto [phi, betar, gammak, sigma ] = extract_parameters(v_parameters_);

        // Compute the necessary
        T::VariableType dataset_betar, loglik1 = 0;
        T::VariableType partial1 = 0;
        for(const auto &i: *(indexes_group_)){
            dataset_betar = Dataset::dataset.row(i) * betar;
            for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * Dataset::dropout_intervals(i,k);
                partial1 += Dataset::dropout_intervals(i,k) * gammak(k);
            }
        }

        T::VariableType node, weight, loglik2 = 0;
        T::VariableType partial2 = 0;
        for(T::IndexType q = 0; q < n_nodes; ++q){
            node = nodes[q];
            weight = weights[q];
            partial2 = 0;

            for(const auto &i: *indexes_group_){
                dataset_betar = Dataset::dataset.row(i) * betar;
                for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
                    partial2 += (exp(dataset_betar + phi(k) + sqrt(2) * sigma * gammak(k) * node)) * Dataset::e_time(i,k);
                }
            }
            loglik2 += weight * exp(sqrt(2) * sigma * node * partial1 - partial2);
        }
        loglik2 = log(loglik2);
        return (loglik1 + loglik2);
    };
};

// Method for building the second derivative of the function wrt one direction
void PowerParameterModel::build_dd_loglikelihood(){
    // Implement the function dd_ll_pp
    dd_ll_pp = [this] (T::IndexType index_, T::VectorXdr& v_parameters_){
        T::VariableType value = v_parameters_(index_);
        T::VariableType valueplush = value + h_dd;
        T::VariableType valueminush = value - h_dd;
        
        T::VectorXdr v_parameters_plus = v_parameters_; 
        T::VectorXdr v_parameters_minus = v_parameters_;
        v_parameters_plus(index_) = valueplush;
        v_parameters_minus(index_) = valueminush;
        
        T::VariableType result = (ll_pp(v_parameters_plus) + ll_pp(v_parameters_minus) - 2*ll_pp(v_parameters_))/(h_dd * h_dd);
        return result;
    };
};

// Method for computing the second derivtaive of the log-likelihood
void PowerParameterModel::compute_hessian_diagonal(T::VectorXdr& v_parameters_){  
    // Initialize the hessian diagonal matrix 
    for(T::IndexType i = 0; i < n_parameters; ++i){
        hessian_diag(i) = dd_ll_pp(i, v_parameters_);
    }
};

// compute the standard error of the parameters
void PowerParameterModel::compute_se(T::VectorXdr& v_parameters_){
     // Initialize the diagonal of the hessian matrix
     compute_hessian_diagonal(v_parameters_);
     
     // Define an element that stores the information element
     T::VariableType information_element;
     
     // Initialize the standard error vector
     for(T::IndexType i = 0; i < n_parameters; ++i){
         information_element = -hessian_diag(i);
         se(i) = 1/(sqrt(information_element));
     }
};


// Method for builfing the result, provided the optimal vector
void PowerParameterModel::evaluate_loglikelihood(){
    T::VariableType optimal_ll_pp = ll_pp(v_parameters);

    // Initialize the standard error of the parameters
    compute_se(v_parameters);

    // Store the results in the class
    result = ResultsMethod::Results(name_method, n_parameters, v_parameters, optimal_ll_pp, se);
    result.print_results();
};

} // end namespace



