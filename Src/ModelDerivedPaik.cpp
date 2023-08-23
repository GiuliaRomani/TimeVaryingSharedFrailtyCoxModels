
// Include header files
#include "ModelDerived.hpp"
#include "MathFunctions.cpp"

// Include libraries
#include <cmath>
#include <tuple>
#include <algorithm>
#include <random>

/**
 * Implementation of the methods declared in the class "AdaptedPaikeaM"
*/

namespace TVSFCM{
using T = TypeTraits;

//! Constructor
AdaptedPaikeaM::AdaptedPaikeaM(const T::FileNameType& filename1_, const T::FileNameType& filename2_):
        //! Construction of the base classes
        ModelBase(filename1,_ filename2_),
        Parameters(filename1_, 2 * Dataset::n_intervals + Dataset::n_regressors + 2, 
                    Dataset::n_intervals, Dataset::n_regressors, 5, 
                    {Dataset::n_intervals, Dataset::n_regressors, 1, 1, Dataset::n_intervals}){
            
            //! Initialize the number of parameters
            compute_n_parameters();

            //! Resize standard error vector and hessian diagonal according to the number of parameters
            hessian_diag.resize(n_parameters);
            se.resize(n_parameters);

            //! Build the log-likelihood functions
            build_loglikelihood();
            build_dd_loglikelihood();

            //! If the number of declared threads is > 1, then build also the log-likelihood function for the parallel version
            if(n_threads > 1)
                build_loglikelihood_parallel();

};
        
//! Virtual method for computing the number of parameters
void AdaptedPaikeaM::compute_n_parameters() noexcept{
    n_parameters = (2 * Dataset::n_intervals + Dataset::n_regressors + 2);
};

//! Virtual method for extracting the parameters from the vector of parameters
T::TuplePaikType AdaptedPaikeaM::extract_parameters(T::VectorXdr& v_parameters_) noexcept{
    T::VectorXdr phi = v_parameters_.head(Dataset::n_intervals);                  
    T::VectorXdr betar = v_parameters_.block(Dataset::n_intervals, 0, Dataset::n_regressors,1);
    T::VariableType mu1 = v_parameters_(Dataset::n_intervals + Dataset::n_regressors);
    T::VariableType mu2 = 1 - mu1;
    T::VariableType nu = v_parameters_(Dataset::n_intervals + Dataset::n_regressors + 1);
    T::VectorXdr gammak = v_parameters_.tail(Dataset::n_intervals);

    return std::make_tuple(phi, betar, mu1, mu2, nu, gammak);
};

//! Method for extracting the variables related to dataset and to a precise group
T::TupleMatrixAType AdaptedPaikeaM::extract_matrixA_variables(T::SharedPtrType indexes_group_, T::VectorXdr& phi_, T::VectorXdr& betar_) noexcept{
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr A_ijk(n_individuals_group, Dataset::n_intervals);
    T::VectorXdr A_ik(Dataset::n_intervals);
    T::VariableType A_i;

    //! Compute A_ijk
    T::IndexType index = 0;
    T::VariableType dataset_betar = 0;
    for(const auto &i: *(indexes_group_)){
        dataset_betar = Dataset::dataset.row(i) * betar_;
        for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
            A_ijk(index,k) = Dataset::e_time(i,k) * exp(dataset_betar + phi_(k));
        }
        index += 1;
    }

    //! Compute A_ik and A_i
    A_ik = A_ijk.colwise().sum();
    A_i = A_ik.sum();

    return std::make_tuple(A_ijk, A_ik, A_i);
};

//! Method for extracting the dropout variables related to a group
T::TupleDropoutType AdaptedPaikeaM::extract_dropout_variables(T::SharedPtrType indexes_group_) noexcept{
    //! Define the variables 
    T::NumberType n_individuals_group = (*indexes_group_).size();
    T::MatrixXdr d_ijk(n_individuals_group, Dataset::n_intervals);
    T::VectorXdr d_ik(Dataset::n_intervals);
    T::VariableType d_i;

    //! Compute d_ijk
    T::IndexType index = 0;
    for(const auto &i: (*indexes_group_)){
    	for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
    		d_ijk(index,k) = Dataset::dropout_intervals(i,k);
    	}
    	index += 1;
    }   

    //! Compute d_ik and d_i
    d_ik = d_ijk.colwise().sum();
    d_i = d_ijk.sum();

    return std::tuple(d_ijk, d_ik, d_i);
};

//! Method for building the overall and group log-likelihood functions
void AdaptedPaikeaM::build_loglikelihood() noexcept{
    //! Implement the function ll_paik
    ll_paik = [this] (T::VectorXdr& v_parameters_){              
        T::VariableType log_likelihood_group, log_likelihood = 0;

        //! For each group, compute the log-likelihood and then sum them
        T::MapType::iterator it_map = Dataset::map_groups.begin();
        T::MapType::iterator it_map_end = Dataset::map_groups.end();
        for(; it_map != it_map_end; ++it_map){
            //! All the indexes in a group
            const auto& indexes_group = it_map->second;

            log_likelihood_group = ll_group_paik(v_parameters_, indexes_group); 
            log_likelihood += log_likelihood_group;
        }
        return log_likelihood;
    };

    
    //! Implement the function ll_group_paik
    ll_group_paik = [this] (T::VectorXdr& v_parameters_, T::SharedPtrType indexes_group_){     
        //! Extract single parameters from the vector using the methods for extracting the variables and the parameters
        auto [phi, betar, mu1, mu2, nu, gammak] = extract_parameters(v_parameters_);
        auto [A_ijk, A_ik, A_i] = extract_matrixA_variables(indexes_group_, phi, betar);
        auto [d_ijk, d_ik, d_i] = extract_dropout_variables(indexes_group_);

        //! Compute the first term of the log-likelihood and then subtract the second term
	    T::VariableType  dataset_betar, loglik1 = 0.;
	    for(const auto &i: *indexes_group_){
	        dataset_betar = Dataset::dataset.row(i) * betar;
	        for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
	            loglik1 += (dataset_betar + phi(k)) * Dataset::dropout_intervals(i,k);
	        }
	    }
	    loglik1 -= (mu1/nu) * log(1 + nu * A_i);

        //! Compute the third term of the formula
	    T::VariableType loglik2 = 0.;
	    for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
	        loglik2 -= (mu2/gammak(k)) * log(1 + gammak(k) * A_ik(k));
	    }
	
        //! Compute the rest of the formula (second line)
        T::VariableType result, loglik4, loglik3 = 0.;
        T::VariableType gamma_res1, gamma_res2, gamma_res3, gamma_res4 = 0.;
        T::VariableType arg1, arg2, arg3 = 0.;
        T::VariableType exp1 = 0.;
        T::VariableType pow1, pow2 = 0.;
        T::NumberType d_ik_size, coeff_binom = 0;
        T::VariableType actual_gammak = 0.;
        T::VariableType ll, d_ik_sizel = 0.;
        
        gamma_res1 = tgamma(mu1/nu);
        arg1 = (A_i + 1/nu);
        for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
            loglik4 = 0.;
            d_ik_size = d_ik(k);
            actual_gammak = gammak(k);
            gamma_res2 = tgamma(mu2/actual_gammak);
            arg2 = (d_ik_size + mu2/actual_gammak);
            arg3 = (A_ik(k) + 1/actual_gammak);
            for(T::IndexType l = 0; l <= d_ik_size; ++l){
                coeff_binom = binom(static_cast<T::NumberType>(d_ik_size), l);
                gamma_res3 = tgamma(arg2 - l);
                gamma_res4 = tgamma(mu1/nu + l);

                //! Cast for negative power
                ll = static_cast<T::VariableType>(l);
                d_ik_sizel = static_cast<T::VariableType>(d_ik_size);
                exp1 = ll - d_ik_sizel;
                pow1 = pow(arg3, exp1);
                pow2 = pow(arg1, ll);
                loglik4 += coeff_binom * gamma_res3 / (gamma_res2 * gamma_res1 * pow2) * gamma_res4 * pow1;
            }
            loglik3 += log(loglik4);           
        }
    	result = loglik1 + loglik2 + loglik3;
    	return (result);
    };
};

//! Build the overall log-likelihood function in the parallel version
void AdaptedPaikeaM::build_loglikelihood_parallel() noexcept {
    ll_paik_parallel = [this] (T::VectorXdr& v_parameters_){
        T::VariableType log_likelihood = 0;             //! Overall log-likelihood value
        T::IdType id = 0;                               //! Id of the thread executing an iteration

        //! Loop over the map through an iterator
        T::MapType::iterator it_map_begin = Dataset::map_groups.begin();
        T::MapType::iterator it_map_end = Dataset::map_groups.end();
        T::MapType::iterator it_map;

    //! Parallel region
    //! If you want to print the iterations execution order, uncomment the related lines and add (id) to (firstprivate)
    omp_set_schedule(omp_sched_t(ParallelComponents::schedule_type), ParallelComponents::chunk_size);
    #pragma omp parallel for num_threads(ParallelComponents::n_threads) firstprivate(it_map, id) schedule(runtime) reduction(+:log_likelihood)
        for(T::IndexType j = 0; j < n_groups; ++j){
            if(it_map != it_map_end){
                it_map = it_map_begin;
                std::advance(it_map, j);
                const auto& indexes_group = it_map->second;

                log_likelihood += ll_group_paik(v_parameters_, indexes_group);

                id = omp_get_thread_num();
                std::cout << "Iteration " << j << " executed by thread " << id << " out of " << n_threads << std::endl;
            }           
        }
        return log_likelihood;
    };
};

//! Method for building the second derivative of the log-likelihood function wrt one direction.
//! Even if we use more than one thread, the derivative is not computed with the parallel version of the log-likelihood but with the serial version.
void AdaptedPaikeaM::build_dd_loglikelihood() noexcept{
    //! Implement the function dd_ll_pp
    dd_ll_paik = [this] (T::IndexType index_, T::VectorXdr& v_parameters_){
        T::VariableType value = v_parameters_(index_);
        T::VariableType valueplush = value + h_dd;
        T::VariableType valueminush = value - h_dd;
        
        T::VectorXdr v_parameters_plus = v_parameters_; 
        T::VectorXdr v_parameters_minus = v_parameters_;
        v_parameters_plus(index_) = valueplush;
        v_parameters_minus(index_) = valueminush;

        T::VariableType result = 0.;
        result = (ll_paik(v_parameters_plus) + ll_paik(v_parameters_minus) - 2*ll_paik(v_parameters_))/(h_dd * h_dd);

        return result;
    };
};

//! Method for computing the diagonal of the hessian matrix
void AdaptedPaikeaM::compute_hessian_diagonal(T::VectorXdr& v_parameters_) noexcept{    
    for(T::IndexType i = 0; i < n_parameters; ++i){
        hessian_diag(i) = dd_ll_paik(i, v_parameters_);
    }
};

//! Compute the standard error of the parameters
void AdaptedPaikeaM::compute_se(T::VectorXdr& v_parameters_) noexcept{
    compute_hessian_diagonal(v_parameters_);
    T::VariableType information_element;
     
     for(T::IndexType i = 0; i < n_parameters; ++i){
         information_element = -hessian_diag(i);
         se(i) = 1/(sqrt(information_element));
     }
};

//! Compute the standrd deviation of the frailty
void AdaptedPaikeaM::compute_sd_frailty(T::VectorXdr& v_parameters_) noexcept{
    T::TuplePaikType extracted_parameters = extract_parameters(v_parameters_);
    auto mu1 = std::get<2>(extracted_parameters);
    auto mu2 = std::get<3>(extracted_parameters);
    auto nu = std::get<4>(extracted_parameters);
    auto gammak = std::get<5>(extracted_parameters);

    for(T::IndexType k = 0; k < Dataset::n_intervals; ++k){
        variance_frailty(k) = mu1 * nu + mu2 * gammak(k);
        sd_frailty(k) = std::sqrt(variance_frailty(k));
    }
}

//! Method for executing the log-likelihood
void AdaptedPaikeaM::evaluate_loglikelihood() noexcept{
    //! To print the numerosity of each group
    //Dataset::print_dimension_groups();

    //! According to the number of threads, compute the log-likelihood in the parallel or in the serial way
    T::VariableType optimal_ll_paik = 0.;
    if(n_threads == 1)
        optimal_ll_paik = ll_paik(v_parameters);
    else
        optimal_ll_paik = ll_paik_parallel(v_parameters);

    //! Compute the standard error of the parameters
    //compute_se(v_parameters);

    //! Compute the standard deviation fo the frailty
    //compute_sd_frailty(v_parameters);
       
    //! Store the results in the class
    result = Results(name_method, n_parameters, v_parameters, optimal_ll_paik, se, sd_frailty, 
                    ParallelComponents::n_threads, ParallelComponents::chunk_size, ParallelComponents::schedule_type_name);
};

} // end namespace
