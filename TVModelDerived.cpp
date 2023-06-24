// Include header files
#include "TVModelDerived.hpp"
#include "SupportFunctions.cpp"

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


T::TuplePPType PowerParameterModel::extract_parameters(T::VectorXdr& v_parameters_){
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
    ll_pp = [this] (T::VectorXdr& v_parameters_){
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
    ll_group_pp = [this] (T::VectorXdr& v_parameters_, T::SharedPtrType& indexes_group_){

        // Extract single parameters from the vector
        // v_parameters_(index) = x;
        auto [phi, betar, gammak, sigma] = extract_parameters(v_parameters_);

        // Compute the necessary
        T::VariableType loglik1, partial1 = 0;
        for(const auto &i: *(indexes_group_)){
            T::VariableType dataset_betar = dataset.row(i) * betar;
            for(T::IndexType k = 0; k < n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * dropout_intervals(i,k);
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

// Virtual method for extracting the parameters fromt the vector
T::TuplePaikType PaikModel::extract_parameters(T::VectorXdr& v_parameters_){
    // Extract parameters from the vector
    T::VectorXdr phi = v_parameters_.head(n_intervals);                 // block(0,0,n_intervals,1);  
    T::VectorXdr betar = v_parameters_.block(n_intervals, 0, n_regressors,1);
    T::VariableType mu1 = v_parameters_(n_intervals + n_regressors);
    T::VariableType mu2 = 1 - mu1;
    T::VariableType nu = v_parameters_(n_intervals + n_regressors + 1);
    T::VectorXdr gammak = v_parameters_.tail(n_intervals);

    return std::make_tuple(phi, betar, mu1, mu2, nu, gammak);
};

T::TupleMatrixAType PaikModel::extract_matrixA_variables(T::SharedPtrType& indexes_group_, T::VectorXdr& phi_, T::VectorXdr& betar_){
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr A_ijk(n_individuals_group, n_intervals);
    T::VectorXdr A_ik(n_intervals);
    T::VariableType A_i;

    T::IndexType index = 0;
    for(const auto &i: *indexes_group_){
        T::VariableType dataset_betar = dataset.row(i) * betar_;
        for(T::IndexType k = 0; k < n_intervals; ++k){
            A_ijk(index,k) = e_time(i,k) * exp(dataset_betar + phi_(k));
        }
        index += 1;
    }

    A_ik = A_ijk.colwise().sum();
    A_i = A_ik.sum();

    return std::make_tuple(A_ijk, A_ik, A_i);
};

T::TupleDropoutType PaikModel::extract_dropout_variables(T::SharedPtrType& indexes_group_){
    // Define the variables 
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr d_ijk(n_individuals_group, n_intervals);
    T::VectorXdr d_ij(n_individuals_group);
    T::VariableType d_i;

    // Initialize them
    T::IndexType index = 0;
    for(const auto &i: (*indexes_group_)){
    	for(T::IndexType k = 0; k < n_intervals; ++k){
    		d_ijk(index,k) = dropout_intervals(i,k);
    	}
    	index += 1;
    }    
    d_ij = d_ijk.rowwise().sum();
    d_i = d_ijk.sum();

    return std::tuple(d_ijk, d_ij, d_i);
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
    ll_group_paik = [this] (T::VectorXdr& v_parameters_, T::SharedPtrType& indexes_group_){     // T::VariableType x, T::IndexType index, 

        // Extract single parameters from the vector
        // v_parameters_(index) = x;
        auto [phi, betar, mu1, mu2, nu, gammak] = extract_parameters(v_parameters_);
        auto [A_ijk, A_ik, A_i] = extract_matrixA_variables(indexes_group_, phi, betar);
        auto [d_ijk, d_ik, d_i] = extract_dropout_variables(indexes_group_);

        // Compute the first component of the likelihood
        T::VariableType loglik1 = 0;
        for(const auto &i: *(indexes_group_)){
            T::VariableType dataset_betar = dataset.row(i) * betar;
            for(T::IndexType k = 0; k < n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * dropout_intervals(i,k);
            }
        }
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
                T::VariableType coeff_bin = binom(d_ik_size, l);
                T::VariableType tgamma1 = tgamma(d_ik_size + mu2/gammak(k) - l);
                T::VariableType tgamma2 = tgamma(mu1/nu + l);
                T::VariableType tgamma3 = tgamma(mu2/gammak(k));
                T::VariableType tgamma4 = tgamma(mu1/nu);
                T::VariableType pow1 = pow(A_ik(k) + 1/gammak(k), l - d_ik_size);
                T::VariableType pow2 = pow(A_i + 1/nu, l);
                loglik4 += coeff_bin * tgamma1 * tgamma2 * pow1 * factor_c/ (tgamma3 * tgamma4 * pow2);                 
            }
            loglik3 += log(loglik4);
        }
    return (loglik1 + loglik2 + loglik3);
    };
};



// Method for executing the log-likelihood
void PaikModel::optimize_loglikelihood(){
    T::VectorXdr& v_parameters = parameters.get_v_parameters();
    T::VariableType optimal_ll_pp = ll_paik(v_parameters);

    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_parameters, optimal_ll_pp);
    result.print_results();
};
*/

/*
//-------------------------------------------------------------------------------------------------------------
// Implementations for the LogFrailty model
// Constructor
LogFrailtyModel::LogFrailtyModel(const T::FileNameType& filename1, const T::FileNameType& filename2):
        // Base constructor for base class
        ModelBase(filename1, filename2) {
            // Initialize the number of parameters
            n_parameters = compute_n_parameters();

            // Initialize the vector of number of parameters
            all_n_parameters = {n_intervals, n_regressors, 1, 1, 1};

            // Construct the class parameters
            parameters = Params::Parameters(filename1, n_parameters, n_intervals, n_regressors, 
                                            n_ranges_parameters, all_n_parameters);

            // Build the log-likelihood
            build_loglikelihood();

};
        
// Virtual method for computing the number of parameters
T::NumberType LogFrailtyModel::compute_n_parameters() {
    return (n_intervals + n_regressors + 3);
};


T::TupleLFType LogFrailtyModel::extract_parameters(T::VectorXdr& v_parameters_) const{
    // Extract parameters from the vector
    T::VectorXdr phi = v_parameters_.head(n_intervals);                 // block(0,0,n_intervals,1);  
    T::VectorXdr betar = v_parameters_.block(n_intervals, 0, n_regressors,1);
    T::VariableType lambda1 = v_parameters_(n_intervals + n_regressors);
    T::VariableType lambda2 = v_parameters_(n_intervals + n_regressors+1);
    T::VariableType angle_alpha = v_parameters_(n_intervals + n_regressors+2);

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

T::TupleDropoutType LogFrailtyModel::extract_dropout_variables(const T::SharedPtrType& indexes_group_) const{
    // Define the variables 
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr d_ijk(n_individuals_group, n_intervals);
    T::VectorXdr d_ij(n_individuals_group);

    T::VariableType d_i;
    T::IndexType index = 0;

    // Initialize them
    for(const auto &i: (*indexes_group_)){
    	for(T::IndexType k = 0; k < n_intervals; ++k){
    		d_ijk(index,k) = dropout_intervals(i,k);
    	}
    	index += 1;
    }    
    d_ij = d_ijk.rowwise().sum();
    d_i = d_ijk.sum();

    return std::tuple(d_ijk, d_ij, d_i);
};

T::VectorXdr LogFrailtyModel::extract_time_to_event(const T::SharedPtrType& indexes_group_) const{
    // Define the variables 
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::VectorXdr time_to_event_group(n_individuals_group);
    T::IndexType index = 0;
    
    for(const auto& i: (*indexes_group_)){
    	time_to_event_group(index) = time_to_event(i);
        index += 1;
    }
    
    return time_to_event_group;
};

// Method for building the log-likelihood
void LogFrailtyModel::build_loglikelihood(){

    // Implement the function ll_pp
    // T::VariableType x, T::IndexType index, 
    ll_lf = [this] (T::VectorXdr& v_parameters_){
        T::VariableType log_likelihood = 0;
        T::MapType::iterator it_map = database.get_map_groups().begin();
        T::MapType::iterator it_map_end = database.get_map_groups().end();

        // For each group, compute the likelihood and then sum them
        for(it_map; it_map != it_map_end; ++it_map){
            // All the indexes in a group
            T::SharedPtrType indexes_group = it_map->second;

            // Compute the log-likelihood related to a group
            T::VariableType log_likelihood_group = ll_group_lf(v_parameters_, indexes_group); //x, index, 
            log_likelihood += log_likelihood_group;
        }

        // Subtract the constant term
        log_likelihood -= n_groups*log(M_PI);
        return log_likelihood;
    };

    
    // Implement the function ll_group_pp
    // T::VariableType x, T::IndexType index, 
    ll_group_lf = [this] (T::VectorXdr& v_parameters_, T::SharedPtrType& indexes_group_){

        // Extract single parameters from the vector
        // v_parameters_(index) = x;
        auto [phi, betar, sigma2c, sigmacb, sigma2b, gammas, sigma2r] = extract_parameters(v_parameters_);
        auto [d_ijk, d_ij, d_i] = extract_dropout_variables(indexes_group_);
        auto time_to_event_group(extract_time_to_event(indexes_group_));

        // Compute the first component of the likelihood
        T::VariableType loglik1 = 0;
        for(const auto &i: (*indexes_group_)){
            T::VariableType dataset_betar = dataset.row(i) * betar;
            for(T::IndexType k = 0; k < n_intervals; ++k){
                loglik1 += (dataset_betar + phi(k)) * dropout_intervals(i,k);
            }
        }

        // Compute the second line of the likelihood
        T::VariableType loglik2 = 0;
        T::VariableType weight, node, exp1, G1;
        for(T::IndexType q = 0; q < n_nodes; ++q){
            weight = weights[q];
            node = nodes[q];
            exp1 = exp(sqrt(2 * sigma2r) * node * d_i);
            G1 = G(node, indexes_group_, v_parameters_);
            loglik2 += weight * exp1 * G1 * factor_c;
        }
        loglik2 = log(loglik2);

        return (loglik1 + loglik2);
    };

    // Implement the function G
    G = [this] (T::VariableType z, const T::SharedPtrType& indexes_group_, T::VectorXdr& v_parameters_){

        // Extract parameters and variables from the vectors
        auto [phi, betar, sigma2c, sigmacb, sigma2b, gammas, sigma2r] = extract_parameters(v_parameters_);
        auto d_ij = std::get<1>(extract_dropout_variables(indexes_group_));
        auto d_i = std::get<2>(extract_dropout_variables(indexes_group_));
        auto time_to_event_group(extract_time_to_event(indexes_group_));

        // Define some useful variables
        T::VariableType partial1, partial2, partial, partial3;
        T::VariableType node, weight;
        T::VariableType arg_f, res_f, arg_exp1, arg_exp2, arg_exp3, exp1;
        T::VariableType time_to_event_i, dataset_betar;

        partial1 = gammas * d_i;
        partial2 = d_ij.dot(time_to_event_group);
        partial = 0;
        for(T::IndexType u = 0; u < n_nodes; ++u){
            node = nodes[u];
            weight = weights[u];
            partial3 = 0;
            arg_f = sqrt(2 * sigma2b) * node;
            for(const auto &i: *indexes_group_){
                time_to_event_i = time_to_event(i);
                dataset_betar = dataset.row(i) * betar;
                for(T::IndexType kk = 0; kk < n_intervals; ++kk){
                    res_f = f_ijk(arg_f, kk, time_to_event_i, v_parameters_);
                    partial3 += res_f * exp(dataset_betar);
                }
            }
            arg_exp1 = arg_f * (partial1 - partial2);
            arg_exp2 = sqrt(2 * sigma2r) * z;
            arg_exp3 = arg_f * gammas;
            exp1 = exp(arg_exp2  + arg_exp3);
            partial += weight * exp(arg_exp1 - partial3 * exp1 / arg_f);
        }
        return partial;
    };

    f_ijk = [this] (T::VariableType b, T::IndexType kkk, T::VariableType time_to_i, T::VectorXdr& v_parameters_){
        // Extract the baseline components from the vector of parameters
        T::VectorXdr phi = std::get<0>(extract_parameters(v_parameters_));
        const auto& v_intervals = time.get_v_intervals();

        // Define some useful variables
        T::VariableType exp1, exp2, exp3;

        // Check conditions
        if(time_to_i < v_intervals[kkk])
            return 0.;
        else if((time_to_i >= v_intervals[kkk]) & (time_to_i < v_intervals[kkk+1])){
            exp1 = exp(phi(kkk));
            exp2 = exp(b * time_to_i);
            exp3 = exp(b * v_intervals[kkk]);
            return (exp1 * (exp2 - exp3));
        }
        else if(time_to_i >= v_intervals[kkk+1]){
            exp1 = exp(phi(kkk));
            exp2 = exp(b * v_intervals[kkk+1]);
            exp3 = exp(b * v_intervals[kkk]);
            return (exp1 * (exp2 - exp3));
        }
    };
};

// Method for executing the log-likelihood
void LogFrailtyModel::optimize_loglikelihood(){
    T::VectorXdr& v_parameters = parameters.get_v_parameters();
    T::VariableType optimal_ll_lf = ll_lf(v_parameters);

    // Store the results in the class
    result = ResultsMethod::Results(n_parameters, v_parameters, optimal_ll_lf);
    result.print_results();
};
*/






} // end namespace


