
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
// Implementations for the LogFrailty model
// Constructor
LogFrailtyModel::LogFrailtyModel(const T::FileNameType& filename1, const T::FileNameType& filename2):
        // Constructor for base classes
        ModelBase(filename1, filename2),
        Parameters(filename1, Dataset::n_intervals + Dataset::n_regressors + 3, 
                    Dataset::n_intervals, Dataset::n_regressors, 5, 
                    {Dataset::n_intervals, Dataset::n_regressors, 1, 1, 1}) {

            // Initialize the number of parameters
            compute_n_parameters();

            // Resize the vectors according to the number of parameters
            hessian_diag.resize(n_parameters);
            se.resize(n_parameters);

            // Build the log-likelihood
            build_loglikelihood();
            build_dd_loglikelihood();

};
        
// Virtual method for computing the number of parameters
void LogFrailtyModel::compute_n_parameters() {
    n_parameters = (Dataset::n_intervals + Dataset::n_regressors + 3);
};


T::TupleLFType LogFrailtyModel::extract_parameters(T::VectorXdr& v_parameters_) const{
    // Extract parameters from the vector
    T::VectorXdr phi = v_parameters_.head(Dataset::n_intervals);                 // block(0,0,n_intervals,1);  
    T::VectorXdr betar = v_parameters_.block(Dataset::n_intervals, 0, Dataset::n_regressors,1);
    T::VariableType lambda1 = v_parameters_(Dataset::n_intervals + Dataset::n_regressors);
    T::VariableType lambda2 = v_parameters_(Dataset::n_intervals + Dataset::n_regressors+1);
    T::VariableType angle_alpha = v_parameters_(Dataset::n_intervals + Dataset::n_regressors+2);

    // Compute the original variable 
    T::VariableType cos_angle = cos(angle_alpha);
    T::VariableType sin_angle = sin(angle_alpha);
    T::VariableType sigma2c = lambda1 * (cos_angle * cos_angle) + lambda2 * (sin_angle * sin_angle);
    T::VariableType sigmacb = (lambda1 - lambda2) * sin_angle * cos_angle;
    T::VariableType sigma2b = lambda1 * (sin_angle * sin_angle) + lambda2 * (cos_angle * cos_angle);

    // Compute other variables
    T::VariableType gammas = sigmacb / sigma2b;
    T::VariableType sigma2r = sigma2c - sigma2b * gammas * gammas;

    return std::make_tuple(phi, betar, sigma2c, sigmacb, sigma2b, gammas, sigma2r);
};

T::TupleDropoutType LogFrailtyModel::extract_dropout_variables(const T::SharedPtrType indexes_group_) const{
    // Define the variables 
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr d_ijk(n_individuals_group, Dataset::n_intervals);
    T::VectorXdr d_ij(n_individuals_group);
    T::VariableType d_i;

    T::IndexType index = 0;
    // Initialize them
    for(const auto &i: *indexes_group_){
    	for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
    		d_ijk(index,k) = Dataset::dropout_intervals(i,k);
    	}
    	index += 1;
    }    
    d_ij = d_ijk.rowwise().sum();
    d_i = d_ijk.sum();

    return std::tuple(d_ijk, d_ij, d_i);
};

T::VectorXdr LogFrailtyModel::extract_time_to_event(const T::SharedPtrType indexes_group_) const{
    // Define the variables 
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::VectorXdr time_to_event_group(n_individuals_group);
    T::IndexType index = 0;
    for(const auto& i: *indexes_group_){
    	time_to_event_group(index) = Dataset::time_to_event(i);
        index += 1;
    }
    return time_to_event_group;
};

// Method for building the log-likelihood
void LogFrailtyModel::build_loglikelihood(){

    // Implement the function ll_lf
    ll_lf = [this] (T::VectorXdr& v_parameters_){
        T::VariableType log_likelihood_group, log_likelihood = 0;
        T::SharedPtrType indexes_group = nullptr;

        // For each group, compute the likelihood and then sum them
        T::MapType::iterator it_map = Dataset::map_groups.begin();
        T::MapType::iterator it_map_end = Dataset::map_groups.end();
        for(; it_map != it_map_end; ++it_map){
            // All the indexes in a group
            indexes_group = it_map->second;

            // Compute the log-likelihood related to a group
            log_likelihood_group = ll_group_lf(v_parameters_, indexes_group);
            log_likelihood += log_likelihood_group;
            
            indexes_group = nullptr;
        }

        // Subtract the constant term
        log_likelihood -= (Dataset::n_groups)*log(M_PI);
        return log_likelihood;
    };

    // Implement the function ll_group_lf
    ll_group_lf = [this] (T::VectorXdr& v_parameters_, T::SharedPtrType indexes_group_){

        // Extract single parameters from the vector
        // v_parameters_(index) = x;
        auto [phi, betar, sigma2c, sigmacb, sigma2b, gammas, sigma2r] = extract_parameters(v_parameters_);
        auto [d_ijk, d_ij, d_i] = extract_dropout_variables(indexes_group_);

        // Compute the first component of the likelihood
        T::VariableType dataset_betar, loglik1 = 0;
        for(const auto &i: *indexes_group_){
            dataset_betar = Dataset::dataset.row(i) * betar;
            for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * Dataset::dropout_intervals(i,k);
            }
        }

        // Compute the second line of the likelihood
        T::VariableType loglik2 = 0;
        T::VariableType weight, node, exp1, G1;
        T::VariableType sigmar2 = sqrt(2*sigma2r);
        for(T::IndexType q = 0; q < n_nodes; ++q){
            weight = weights[q];
            node = nodes[q];

            exp1 = exp(node * d_i * sigmar2);
            G1 = G(node, indexes_group_, v_parameters_);
            loglik2 += weight * exp1 * G1;
        }
        loglik2 = log(loglik2);

        T::VariableType result = loglik1 + loglik2;
        return result;
    };

    // Implement the function G
    G = [this] (T::VariableType z, const T::SharedPtrType& indexes_group_, T::VectorXdr& v_parameters_){

        // Extract parameters and variables from the vectors
        auto [phi, betar, sigma2c, sigmacb, sigma2b, gammas, sigma2r] = extract_parameters(v_parameters_);
        auto d_ij = std::get<1>(extract_dropout_variables(indexes_group_));
        auto d_i = std::get<2>(extract_dropout_variables(indexes_group_));
        auto time_to_event_group(extract_time_to_event(indexes_group_));

        // Define some useful variables
        T::VariableType partial1, partial = 0.;
        T::VariableType weight, node;
        T::VariableType dataset_betar, time_to_event_i;
        T::VariableType arg1, arg2, arg3, arg4, arg5, res_f;

        arg1 = gammas * d_i;
        arg2 = d_ij.dot(time_to_event_group);
        for(T::NumberType u = 0; u < n_nodes; ++u ){
            partial1 = 0.;
            node = nodes[u];
            weight = weights[u];
            arg3 = sqrt(2 * sigma2b) * node;
            for(const auto &i: *indexes_group_){
                dataset_betar = Dataset::dataset.row(i) * betar;
                time_to_event_i = Dataset::time_to_event(i);
                for (T::NumberType k = 0; k < Dataset::n_intervals; ++k){
                    res_f = f_ijk(arg3, k, time_to_event_i, v_parameters_);
                    partial1 += exp(dataset_betar) * res_f;
                }
            }
            arg4 = sqrt(2 * sigma2r) * z + arg3 * gammas;
            arg5 = arg3 * (arg1 + arg2) - partial1 * (exp(arg4)) / arg3;
            partial += weight * exp(arg5); 
        }
        return partial;
    };

    f_ijk = [this] (T::VariableType b, T::IndexType kkk, T::VariableType time_to_i, T::VectorXdr& v_parameters_){
        // Extract the baseline components from the vector of parameters
        T::VectorXdr phi = std::get<0>(extract_parameters(v_parameters_));
        const auto& v_intervals = Dataset::v_intervals;

        // Define some useful variables
        T::VariableType exp1, exp2, exp3;
        T::VariableType result;

        // Check conditions
        if(time_to_i < v_intervals[kkk])
            result = 0.;
        else if((time_to_i >= v_intervals[kkk]) & (time_to_i < v_intervals[kkk+1])){
            exp1 = exp(phi(kkk));
            exp2 = exp(b * time_to_i);
            exp3 = exp(b * v_intervals[kkk]);
            result = exp1 * (exp2 - exp3);
        }
        else if(time_to_i >= v_intervals[kkk+1]){
            exp1 = exp(phi(kkk));
            exp2 = exp(b * v_intervals[kkk+1]);
            exp3 = exp(b * v_intervals[kkk]);
            result = exp1 * (exp2 - exp3);
        }
        return result;
    };
};

// Method for building the second derivative of the function wrt one direction
void LogFrailtyModel::build_dd_loglikelihood(){
    // Implement the function dd_ll_pp
    dd_ll_lf = [this] (T::IndexType index_, T::VectorXdr& v_parameters_){
        T::VariableType value = v_parameters_(index_);
        T::VariableType valueplush = value + h_dd;
        T::VariableType valueminush = value - h_dd;
        
        T::VectorXdr v_parameters_plus = v_parameters_; 
        T::VectorXdr v_parameters_minus = v_parameters_;
        v_parameters_plus(index_) = valueplush;
        v_parameters_minus(index_) = valueminush;
        
        T::VariableType result = (ll_lf(v_parameters_plus) + ll_lf(v_parameters_minus) - 2*ll_lf(v_parameters_))/(h_dd * h_dd);
        return result;
    };
};

// Method for computing the second derivtaive of the log-likelihood
void LogFrailtyModel::compute_hessian_diagonal(T::VectorXdr& v_parameters_){
    for(T::IndexType i = 0; i < n_parameters; ++i){
        hessian_diag(i) = dd_ll_lf(i, v_parameters_);
    }
};

// compute the standard error of the parameters
void LogFrailtyModel::compute_se(T::VectorXdr& v_parameters_){
     compute_hessian_diagonal(v_parameters_);
     T::VariableType information_element;
     
     for(T::IndexType i = 0; i < n_parameters; ++i){
         information_element = -hessian_diag(i);
         se(i) = 1/(sqrt(information_element));
     }
};

// Method for executing the log-likelihood
void LogFrailtyModel::evaluate_loglikelihood(){
    T::VariableType optimal_ll_lf = ll_lf(v_parameters);

    compute_se(v_parameters);
       
    // Store the results in the class
    result = ResultsMethod::Results(name_method, n_parameters, v_parameters, optimal_ll_lf, se);
    result.print_results();
};


} // end namespace


