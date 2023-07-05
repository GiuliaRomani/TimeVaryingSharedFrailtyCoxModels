/*
// Include header files
#include "TVModelDerived.hpp"
#include "SupportFunctions.cpp"

// Include libraries
#include <cmath>
#include <tuple>
#include <algorithm>
#include <random>


namespace TVModel{
using T = TypeTraits;

// Implementations for the Paik Model
// Constructor
PaikModel::PaikModel(const T::FileNameType& filename1, const T::FileNameType& filename2):
        // Base constructor for base class
        ModelBase(filename1, filename2) {
            // Initialize the number of parameters
            n_parameters = compute_n_parameters();
            hessian_diag.resize(n_parameters);
            se.resize(n_parameters);

            // Initialize the vector of number of parameters
            all_n_parameters = {Time::n_intervals, Dataset::n_regressors, 1, 1, Time::n_intervals};

            // Construct the class parameters
            parameters = Params::Parameters(filename1, n_parameters, Time::n_intervals, Dataset::n_regressors, 
                                            n_ranges_parameters, all_n_parameters);

            // Build the log-likelihood
            build_loglikelihood();
            build_dd_loglikelihood();

};
        
// Virtual method for computing the number of parameters
T::NumberType PaikModel::compute_n_parameters() {
    return (2 * Time::n_intervals + Dataset::n_regressors + 2);
};

// Virtual method for extracting the parameters fromt the vector
T::TuplePaikType PaikModel::extract_parameters(T::VectorXdr& v_parameters_){
    // Extract parameters from the vector
    T::VectorXdr phi = v_parameters_.head(Time::n_intervals);                 // block(0,0,n_intervals,1);  
    T::VectorXdr betar = v_parameters_.block(Time::n_intervals, 0, Dataset::n_regressors,1);
    T::VariableType mu1 = v_parameters_(Time::n_intervals + Dataset::n_regressors);
    T::VariableType mu2 = 1 - mu1;
    T::VariableType nu = v_parameters_(Time::n_intervals + Dataset::n_regressors + 1);
    T::VectorXdr gammak = v_parameters_.tail(Time::n_intervals);

    return std::make_tuple(phi, betar, mu1, mu2, nu, gammak);
};

T::TupleMatrixAType PaikModel::extract_matrixA_variables(T::SharedPtrType indexes_group_, T::VectorXdr& phi_, T::VectorXdr& betar_){
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr A_ijk(n_individuals_group, Time::n_intervals);
    T::VectorXdr A_ik(Time::n_intervals);
    T::VariableType A_i;

    T::IndexType index = 0;
    T::VariableType dataset_betar = 0;
    for(const auto &i: *(indexes_group_)){
        dataset_betar = Dataset::dataset.row(i) * betar_;
        for(T::IndexType k = 0; k < Time::n_intervals; ++k){
            A_ijk(index,k) = Dataset::e_time(i,k) * exp(dataset_betar + phi_(k));
        }
        index += 1;
    }

    A_ik = A_ijk.colwise().sum();
    A_i = A_ik.sum();

    return std::make_tuple(A_ijk, A_ik, A_i);
};

T::TupleDropoutType PaikModel::extract_dropout_variables(T::SharedPtrType indexes_group_){
    // Define the variables 
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr d_ijk(n_individuals_group, Time::n_intervals);
    T::VectorXdr d_ik(Time::n_intervals);
    T::VariableType d_i;

    // Initialize them
    T::IndexType index = 0;
    for(const auto &i: (*indexes_group_)){
    	for(T::IndexType k = 0; k < Time::n_intervals; ++k){
    		d_ijk(index,k) = Dataset::dropout_intervals(i,k);
    	}
    	index += 1;
    }   

    d_ik = d_ijk.colwise().sum();
    d_i = d_ijk.sum();

    return std::tuple(d_ijk, d_ik, d_i);
};

// Method for building the log-likelihood
void PaikModel::build_loglikelihood(){
    // Implement the function ll_paik
    ll_paik = [this] (T::VectorXdr& v_parameters_){              // T::VariableType x, T::IndexType index, 
        T::VariableType log_likelihood_group, log_likelihood = 0;
        T::SharedPtrType indexes_group = nullptr;

        // For each group, compute the likelihood and then sum them
        T::MapType::iterator it_map = Dataset::map_groups.begin();
        T::MapType::iterator it_map_end = Dataset::map_groups.end();
        for(; it_map != it_map_end; ++it_map){
            // All the indexes in a group
            indexes_group = it_map->second;

            log_likelihood_group = ll_group_paik(v_parameters_, indexes_group); 
            log_likelihood += log_likelihood_group;

            indexes_group = nullptr;
        }
        return log_likelihood;
    };

    
    // Implement the function ll_group_pp
    ll_group_paik = [this] (T::VectorXdr& v_parameters_, T::SharedPtrType indexes_group_){     

        // Extract single parameters from the vector
        auto [phi, betar, mu1, mu2, nu, gammak] = extract_parameters(v_parameters_);
        auto [A_ijk, A_ik, A_i] = extract_matrixA_variables(indexes_group_, phi, betar);
        auto [d_ijk, d_ik, d_i] = extract_dropout_variables(indexes_group_);

        // Compute the first component of the likelihood
        T::VariableType loglik1 = 0;
        T::VariableType dataset_betar = 0;
        for(const auto &i: *(indexes_group_)){
            dataset_betar = Dataset::dataset.row(i) * betar;
            for(T::IndexType k = 0; k < Time::n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * Dataset::dropout_intervals(i,k);
            }
        }
        loglik1 -= (mu1/nu) * log(1 + nu * A_i);

        // Compute the second line of the formula
        T::VariableType loglik2 = 0;
        for(T::IndexType k = 0; k < Time::n_intervals; ++k){
            loglik2 -= (mu2/gammak(k)) * log(1 + gammak(k) * A_ik(k));
        }

        // Compute the third line of the formula
        T::VariableType loglik4, loglik3 = 0.;
        T::VariableType coeff_bin, tgamma1, tgamma2, tgamma3, tgamma4, tgamma5, tgamma6, tgamma7, pow1, pow2 = 0.;
        T::NumberType d_ik_size = 0;

        tgamma4 = tgamma(mu1/nu);
        tgamma5 = (A_i + 1/nu);
        for(T::IndexType k = 0; k < Time::n_intervals; ++k){
            loglik4 = 0.;
            d_ik_size = d_ik(k);
            tgamma3 = tgamma(mu2/gammak(k));
            tgamma6 = (d_ik_size + mu2/gammak(k));
            tgamma7 = (A_ik(k) + 1/gammak(k));

            for(T::IndexType l = 0; l <= d_ik_size; ++l){
                coeff_bin = binom(d_ik_size, l);
                tgamma1 = tgamma(tgamma6 - l);
                tgamma2 = tgamma(mu1/nu + l);
                pow1 = pow(tgamma7, l - d_ik_size);
                pow2 = pow(tgamma5, l);
                loglik4 += coeff_bin * tgamma1 * tgamma2 * pow1 * factor_c/ (tgamma3 * tgamma4 * pow2);                 
            }
            loglik3 += log(loglik4);
        }

    return (loglik1 + loglik2 + loglik3);
    };
};

// Method for building the second derivative of the function wrt one direction
void PaikModel::build_dd_loglikelihood(){
    // Implement the function dd_ll_pp
    dd_ll_paik = [this] (T::IndexType index_, T::VectorXdr& v_parameters_){
        T::VariableType value = v_parameters_(index_);
        T::VariableType valueplush = value + h_dd;
        T::VariableType valueminush = value - h_dd;
        
        T::VectorXdr v_parameters_plus = v_parameters_; 
        T::VectorXdr v_parameters_minus = v_parameters_;
        v_parameters_plus(index_) = valueplush;
        v_parameters_minus(index_) = valueminush;
        
        T::VariableType result = (ll_paik(v_parameters_plus) + ll_paik(v_parameters_minus) - 2*ll_paik(v_parameters_))/(h_dd * h_dd);
        return result;
    };
};

// Method for computing the second derivtaive of the log-likelihood
T::VectorXdr PaikModel::compute_hessian_diagonal(T::VectorXdr& v_parameters_){    
    for(T::IndexType i = 0; i < n_parameters; ++i){
        hessian_diag(i) = dd_ll_paik(i, v_parameters_);
    }
    
    return hessian_diag;
};

// compute the standard error of the parameters
T::VectorXdr PaikModel::compute_se(T::VectorXdr& v_parameters_){
     hessian_diag = compute_hessian_diagonal(v_parameters_);
     T::VariableType information_element;
     
     for(T::IndexType i = 0; i < n_parameters; ++i){
         information_element = -hessian_diag(i);
         se(i) = 1/(sqrt(information_element));
     }
     
     std::cout << se << std::endl;
     return se;
};

// Method for executing the log-likelihood
void PaikModel::evaluate_loglikelihood(T::VectorXdr& v_parameters_){
    T::VariableType optimal_ll_paik = ll_paik(v_parameters);
       
    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_parameters, optimal_ll_paik);
    result.print_results();
};


// Method for executing the log-likelihood
void PaikModel::optimize_loglikelihood(){
    // T::VectorXdr& v_parameters = parameters.get_v_parameters();
    T::VectorType optimal_parameters{-5.099, -3.33, -4.8521, -6.962, -4.017, -5.315, -3.916, -4.913, 
                                        -1.386e-1, -1.134e-1, 1.37e-1, -4.661e-2, -1.385, 1.974e-6,
                                        2.00, 6.413e-3, 7.6376e-3, 3.507e-2, 1.145, 1.299e-2, 2.767e-1, 
                                        1.65761e-1, 1.5344e-1 };
    using MappedVectorType = Eigen::Map<T::VectorXdr>; 
    MappedVectorType v_parameters(optimal_parameters.data(), n_parameters); 
    T::VectorXdr v_opt_parameters = v_parameters;                              

    // compute_se(v_opt_parameters);

    T::VariableType optimal_ll_paik = ll_paik(v_opt_parameters);
    // Ideal value = -2153.992
    
    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_opt_parameters, optimal_ll_paik);
    result.print_results();
};


} // end namespace
*/