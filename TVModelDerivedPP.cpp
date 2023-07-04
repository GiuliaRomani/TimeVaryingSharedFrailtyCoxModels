
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
        // Base constructor for base class
        ModelBase(filename1, filename2) {
            // Initialize the number of parameters
            n_parameters = compute_n_parameters();
            hessian_diag.resize(n_parameters);
            se.resize(n_parameters);

            // Initialize the vector of number of parameters
            all_n_parameters = {Time::n_intervals, Dataset::n_regressors, Time::n_intervals - 1, 1};

            // Construct the class parameters
            parameters = Params::Parameters(filename1, n_parameters, Time::n_intervals, Dataset::n_regressors, 
                                            n_ranges_parameters, all_n_parameters);

            // Build the log-likelihood and the one-directional second derivative
            build_loglikelihood();
            build_dd_loglikelihood();
};
        
// Virtual method for computing the number of parameters
T::NumberType PowerParameterModel::compute_n_parameters() {
    return (2 * Time::n_intervals + Dataset::n_regressors);
};

T::TuplePPType PowerParameterModel::extract_parameters(T::VectorXdr& v_parameters_){
    // Extract parameters from the vector
    T::VectorXdr phi = v_parameters_.head(Time::n_intervals);                 // block(0,0,n_intervals,1);  
    T::VectorXdr betar = v_parameters_.block(Time::n_intervals, 0, Dataset::n_regressors, 1);
    T::VectorXdr gammak(Time::n_intervals);
    gammak(0) = 1.;
    gammak.block(1,0,Time::n_intervals-1,1) = v_parameters_.block(Time::n_intervals + Dataset::n_regressors, 0, Time::n_intervals - 1, 1);

    T::VariableType sigma = v_parameters_(n_parameters - 1);
    sigma = sqrt(sigma);

    return std::make_tuple(phi, betar, gammak, sigma);
};


// Method for building the log-likelihood
void PowerParameterModel::build_loglikelihood(){

    // Implement the function ll_pp
    // T::VariableType x, T::IndexType index, 
    ll_pp = [this] (T::VectorXdr& v_parameters_){
        T::VariableType log_likelihood_group, log_likelihood = 0;
        T::SharedPtrType indexes_group = nullptr;

        // For each group, compute the likelihood and then sum them
        T::MapType::iterator it_map = Dataset::map_groups.begin();
        T::MapType::iterator it_map_end = Dataset::map_groups.end();
        for(; it_map != it_map_end; ++it_map){
            // All the indexes in a group
            indexes_group = it_map->second;

            log_likelihood_group = ll_group_pp(v_parameters_, indexes_group); 
            log_likelihood += log_likelihood_group;

            //indexes_group = nullptr;
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
        T::VariableType loglik1 = 0;
        T::VariableType partial1 = 0;
        T::VariableType dataset_betar = 0;
        for(const auto &i: *(indexes_group_)){
            dataset_betar = Dataset::dataset.row(i) * betar;
            for(T::IndexType k = 0; k < Time::n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * Dataset::dropout_intervals(i,k);
                partial1 += Dataset::dropout_intervals(i,k) * gammak(k);
            }
        }

        T::VariableType loglik2 = 0;
        T::VariableType partial2 = 0;
        T::VariableType node, weight = 0;
        for(T::IndexType q = 0; q < n_nodes; ++q){
            node = nodes[q];
            weight = weights[q];
            partial2 = 0;

            for(const auto &i: *indexes_group_){
                dataset_betar = Dataset::dataset.row(i) * betar;
                for(T::IndexType k = 0; k < Time::n_intervals; ++k){
                    partial2 += (exp(dataset_betar + phi(k) + sqrt(2) * sigma * gammak(k) * node)) * Dataset::e_time(i,k);
                }
            }
            loglik2 += factor_c * weight * exp(sqrt(2) * sigma * node * partial1 - partial2);
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
T::VectorXdr PowerParameterModel::compute_hessian_diagonal(T::VectorXdr& v_parameters_){
    T::VectorXdr hessian_diag(n_parameters);
    
    for(T::IndexType i = 0; i < n_parameters; ++i){
        hessian_diag(i) = dd_ll_pp(i, v_parameters_);
    }
    
    return hessian_diag;
};

// compute the standard error of the parameters
T::VectorXdr PowerParameterModel::compute_se(T::VectorXdr& v_parameters_){
     T::VectorXdr hessian_diag = compute_hessian_diagonal(v_parameters_);
     T::VectorXdr se(n_parameters);
     T::VariableType information_element;
     
     for(T::IndexType i = 0; i < n_parameters; ++i){
         information_element = -hessian_diag(i);
         se(i) = 1/(sqrt(information_element));
     }
     
     std::cout << se << std::endl;
     return se;
};


// Method for executing the log-likelihood
void PowerParameterModel::evaluate_loglikelihood(T::VectorXdr& v_parameters_){
    T::VariableType optimal_ll_pp = ll_pp(v_parameters);
       
    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_parameters, optimal_ll_pp);
    result.print_results();
};


// Method for executing the log-likelihood
void PowerParameterModel::optimize_loglikelihood(){
    // T::VectorXdr& v_parameters = parameters.get_v_parameters();
    T::VectorType optimal_parameters{-4.767, -1.904, -2.113, -6.430, -1.891, -2.384, -3.2470,
                                        -3.8260, -0.1499, -0.1106, 0.1307,
                                        -0.0490, -1.3909, 4.3153, 8.2398, 1.0167, 6.3769,
                                         8.7448, 1.9151, 3.320, 0.104};
    using MappedVectorType = Eigen::Map<T::VectorXdr>; 
    MappedVectorType v_parameters(optimal_parameters.data(), n_parameters); 
    T::VectorXdr v_opt_parameters = v_parameters;                              

    compute_se(v_opt_parameters);

    T::VariableType optimal_ll_pp = ll_pp(v_opt_parameters);
    
    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_opt_parameters, optimal_ll_pp);
    result.print_results();
};


} // end namespace



