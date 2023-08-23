#ifndef MODEL_DERIVED_HPP
#define MODEL_DERIVED_HPP

// Include header files
#include "ModelBase.hpp"
#include "QuadraturePoints.hpp"
#include "Parameters.hpp"

// Include libraries
#include <iostream>
#include <functional>

/**
 * For each model, a class publicly derived from the ModelBase class and from the Parameters class is created. This implements one of the three
 *  Time-Varying Shared Frailty Cox model.
 * 
 * This class is thus initialized only after having initialized the other two base classes.
*/


namespace TVSFCM{
using T = TypeTraits;

//! Class for implementing the Adapted Paik et al.'s Model
class AdaptedPaikeaM final: public ModelBase, 
                            public Parameters{
public:
    /**
     * Default Constructor
    */
    AdaptedPaikeaM() = default;

    /**
     * Constructor that initializes firstly the base classes and then the model variables.
     * @param filename1_ Name of the first file
     * @param filename2_ Name of the second file
    */
    AdaptedPaikeaM(const T::FileNameType& filename1_, const T::FileNameType& filename2_);

    /**
     * Method for computing the value of the log-likelihood, given the optimal vector of parameters and the dataset.
     * 
     * Eventually, it initializes the result class.
    */
    void evaluate_loglikelihood() noexcept override;

    /**
     * Virtual destructor
    */
    virtual ~ AdaptedPaikeaM() = default;

private:
    T::NumberType n_parameters;                                             //! Number of parameter of the model       
    T::VectorXdr& v_parameters = Parameters::v_parameters;                  //! Vector of parameters

    T::IdNameType name_method = "Adapted Paik eaM";                         //! Name of this method

    /**
     * This function is later built as a lambda function. 
     * It implement the model overall log-likelihood as has been defined in the report. 
     * @param v_parameters_ Vector of parameters, where the log-likelihood has to be evaluated
     * @return Evaluation of the log-likelihood at the input, for all the clusters
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_paik; 

    /**
     * This function is later built as a lambda function.
     * It implelemnts the model overall log-likelihood as has been defined in the report, 
     * but using the parallel computing provided by OpenMP, with a number of threads provided by the user.
     * If no threads or just 1 is provided, no parallel implementation is applied.
     * @param v_parameters_ Vector of parameters, where the log-likelihood has to be evaluated
     * @return The evaluation of the log-likelihood
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_paik_parallel; 

    /**
     * This function compute the group log-likelihood for all the individuals belonging to a cluster/group/faculty.
     * @param v_parameters Vector of parameters
     * @param indexes_group_ Shared pointer to the vector of indexes of the dataset, of the students belonging to the cluster
     * @return Evaluation of the log-likelihood in this group
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_, T::SharedPtrType indexes_group_)> ll_group_paik;     
    
    /**
     * This function, implemented as a lambda function, compute the approximation of the second derivative of the 
     * log-likelihood function, along a single dimension.
     * @param index_ Index of the direction with respect to compute the approximation of the second derivative
     * @param v_parameters_ Vector of parameters
     * @return Evaluation of the directional second directive at the vector of parameters
    */
    std::function<T::VariableType(T::IndexType index_, T::VectorXdr& v_parameters_)> dd_ll_paik;	

    /**
     * Method for computing thr number of parameters of this model. It initalizes the private variable n_parameters.
    */
    void compute_n_parameters() noexcept override;

    /**
     * Method for computing the diagonal of the hessian matrix of the log-likelihood function.
     * @param v_parameters_ Vector of parameters, where the hessian must be evaluated
    */
    void compute_hessian_diagonal(T::VectorXdr& v_parameters_) noexcept override;

    /**
     * Method for computing the standard error of the parameters.
     * @param v_parameters_ Vector of parameers, whose standard error must be produced
    */
    void compute_se(T::VectorXdr& v_parameters_) noexcept override;

    /**
     * Method for computing the standard deviation of the frailty.
     * @param v_parameters_ Vector of parameters
    */
    void compute_sd_frailty(T::VectorXdr& v_parameters_) noexcept override;
    
    /**
     * Method for extracting single parameters of the model or vectors of similar parameters (category).
     * @param v_parameters_ Vector of parameters
     * @return A customized model tuple containing all the extracted parameters
    */
    T::TuplePaikType extract_parameters(T::VectorXdr& v_parameters_) noexcept;

    /**
     * Method for extracting the variables related to the dataset, for each group/cluster of individuals.
     * @param indexes_group_ Vector of indexes of the dataset, of the indivividuals belonging to a cluster/group/faculty
     * @param phi_ Vector of baseline parameters
     * @param betar_ Vector of regressor coefficients
     * @return A customized tuple containing the variables to be extracted
    */
    T::TupleMatrixAType extract_matrixA_variables(T::SharedPtrType indexes_group_, T::VectorXdr& phi_, T::VectorXdr& betar_) noexcept;

    /**
     * Method for extracting the variables related to the dropout, for each group/cluster of individuals.
     * @param indexes_group_ Vector of indexes of the dataset, of the individuals belonging to a cluster/group/faculty
     * @return A customized tuple containing the dropout variables to be extracted
    */
    T::TupleDropoutType extract_dropout_variables(T::SharedPtrType indexes_group_) noexcept;

    /**
     * Method for building the model log-likelihood function. It simply contains the definition
     * of the functions ll_paik, ll_group_paik.
    */
    void build_loglikelihood() noexcept override;

    /**
     * Method for building the model log-likelihood function in parallel. It simply contains the 
     * definition of the parallel verison of the function ll_paik_parallel.
    */
    void build_loglikelihood_parallel() noexcept override;

    /**
     * Method for building the second derivative of the log-likelihood function. It simply contains the 
     * definition of the funciton dd_ll_paik.
     * It uses a centered finite difference scheme, with accuracy of the secon order.
    */
    void build_dd_loglikelihood() noexcept override;
};
//--------------------------------------------------------------------------------------------------

//! Class for implementing the Centre-Specific Frailty Model with Power Parameter
class CSFMwithPowerParameter final: public ModelBase,
                                    public Parameters{

public:
    /**
     * Deafult constructor
    */
    CSFMwithPowerParameter() = default;

    /**
     * Constructor that initializes firstly the base classes and then the model variables.
     * @param filename1_ Name of the first file
     * @param filename2_ Name of the second file
    */
    CSFMwithPowerParameter(const T::FileNameType& filename1, const T::FileNameType& filename2);

    /**
     * Method for computing the value of the log-likelihood, given the optimal vector of parameters. 
     * 
     * Eventually, it initializes the result class.
    */
    void evaluate_loglikelihood() noexcept override;

    /**
     * Virtual destructor 
    */
    virtual ~ CSFMwithPowerParameter() = default;

private:
    T::NumberType n_parameters;                                             //! Number of model parameters         
    T::VectorXdr & v_parameters = Parameters::v_parameters;                 //! Vector of parameters

    T::IdNameType name_method = "CSFM with Power Parameter";                           //! Name of this method

    /**
     * This function is later built as a lambda function. 
     * It implememts the model log-likelihood as has been defined in the report. 
     * @param v_parameters Vector of parameters, where the log-likelihood has to be evaluated
     * @return Evaluation of the log-likelihood at the input, for all the clusters
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_pp;  

    /**
     * This function is later built as a lambda function.
     * It implements the model log-likelihood as has been defined in the report, but in a parallel version
     * using OpenMP and a number of threads provided in input by the user. 
     * 
     * The default provided value is 1, that implies no paralellization. A greater value corresponds to parallelization.
     * 
     * @param v_parameters Vector of parameters, where the log-likelihood is evaluated
     * @return Evaluation of the log-likelihood at the input, for all clusters
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_pp_parallel;  

    /**
     * This function compute the log-likelihood for all the individuals belonging to a cluster/group/faculty.
     * @param v_parameters Vector of parameters
     * @param indexes_group_ Shared pointer to the vector of indexes of the dataset, of the students belonging to the cluster
     * @return Evaluation of the log-likelihood in this group
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_, T::SharedPtrType indexes_group_)> ll_group_pp;     
    
    /**
     * This function, implemented as a lambda function, compute the approximation of the second derivative of the 
     * log-likelihood function, along a single dimension.
     * @param index_ Index of the direction with respect to compute the approximation of the second derivative
     * @param v_parameters_ Vector of parameters
     * @return Evaluation of the directional second directive at the vector of parameters
    */
    std::function<T::VariableType(T::IndexType index_, T::VectorXdr& v_parameters_)> dd_ll_pp;	

    //! Quadrature points
    QuadraturePoints9 points9;                                      //! Struct for the nine points quadrature formula
    T::NumberType n_nodes = 9;                                              //! Number of points 
    std::array<T::VariableType, 9>& nodes = points9.nodes;                  //! Array of nodes
    std::array<T::VariableType, 9>& weights = points9.weights;              //! Array of weights


    /**
     * Method for computing thr number of parameters of this model. It initalizes the private variable n_parameters.
    */
    void compute_n_parameters() noexcept override;

    /**
     * Method for extracting each single parameter of the model or a single vector of similar parameters.
     * @param v_parameters_ Vector of parameters
     * @return A customized model tuple containing all the extracted parameters
    */
    T::TuplePPType extract_parameters(const T::VectorXdr& v_parameters_) const noexcept;
    
    /**
     * Method for computing the diagonal of the hessian matrix of the log-likelihood function.
     * @param v_parameters_ Vector of parameters, where the hessian must be evaluated
    */
    void compute_hessian_diagonal(T::VectorXdr& v_parameters_) noexcept override;

    /**
     * Method for computing the standard error of the parameters.
     * @param v_parameters_ Vector of parameters, whose standard error must be produced
    */
    void compute_se(T::VectorXdr& v_parameters) noexcept override;

    /**
     * Method for computing the standard deviation of the frailty.
     * @param v_parameters_ Vector of parameters
    */
    void compute_sd_frailty(T::VectorXdr& v_parameters_) noexcept override;

    /**
     * Method for building the model log-likelihood function. It contains the definition of the 
     * functions ll_pp and ll_group_pp.
    */
    void build_loglikelihood() noexcept override;

    /**
     * Virtual pure method for constructing the log-likelihood in the parallel version. It contains 
     * the definition of the function ll_pp_parallel.
    */
    void build_loglikelihood_parallel() noexcept override;

    /**
     * Method for building the second derivative of the log-likelihood function. It contains
     * the definition of the function dd_ll_pp.
     * It uses a centered finite difference scheme, with accuracy of the secon order.
    */
    void build_dd_loglikelihood() noexcept override;
};


//-----------------------------------------------------------------------------------------
//! Class for implementing the Stochastic Time-Dependent Centre-Specific Frailty Model
class StochasticTimeDependentCSFM final: public ModelBase,
                                        public Parameters{
public:
    /**
     * Default constructor
    */
    StochasticTimeDependentCSFM() = default;

    /**
     * Constructor that initializes firstly the base classes and then the model variables.
     * @param filename1_ Name of the first file
     * @param filename2_ Name of the second file
    */
    StochasticTimeDependentCSFM(const T::FileNameType& filename1, const T::FileNameType& filename2);

    /**
     * Method for computing the value of the log-likelihood, given the optimal vector of parameters. 
     * 
     * Eventually, it initializes the result class.
    */
    void evaluate_loglikelihood() noexcept override;

    /**
     * Virtual destructor
    */
    virtual ~ StochasticTimeDependentCSFM() = default;

private:
    T::NumberType n_parameters;                                             //! Number of parameter of the model 
    T::VectorXdr& v_parameters = Parameters::v_parameters;                  //! Vector of parameters

    T::IdNameType name_method = "Stochastic Time Dependent CSFM";           //! Name of this method

    /**
     * This function is later built as a lambda function. 
     * It implememts the model log-likelihood as has been defined in the report.
     * @param v_parameters Vector of parameters, where the log-likelihood has to be evaluated
     * @return Evaluation of the log-likelihood at the input, for all the clusters
    */    
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_lf; 

    /**
     * This function is later built as a lambda function. 
     * It implements the model log-likelihood as has been defined in the report, in a parallel version
     * using OpenMP and a number of threads specified by the user. If no threads or just 1 are provided, no parallelization is used.
     * @param v_parameters Vector of parameters, where the log-likelihood has to be evaluated
     * @return Evaluation of the log-likelihood at the input, for all the clusters
    */    
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_lf_parallel; 

    /**
     * This function compute the log-likelihood for all the individuals belonging to a cluster/group/faculty.
     * @param v_parameters Vector of parameters
     * @param indexes_group_ Shared pointer to the vector of indexes of the dataset, of the students belonging to the cluster
     * @return Evaluation of the log-likelihood in this group
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_, T::SharedPtrType indexes_)> ll_group_lf;     

    /**
     * This function, implemented as a lambda function, builds the function G  defined in the report.
     * @param z Argument of the function
     * @param indexes_group_ Vector of the indexes of the dataset, of the individuals in a cluster
     * @param v_parameters_ Vector of parameters
     * @return Evaluation of the function
    */
    std::function<T::VariableType(T::VariableType z, T::SharedPtrType indexes_group_, T::VectorXdr& v_parameters_)> G; 

    /**
     * This function, implemented as a lambda function, builds the f_ijk function defined in the report
     * @param b Argument of the function
     * @param kkk Index referring to an element of the time-domain
     * @param time_to_i Time-to-event of an individual
     * @param v_parameters_ Vector of parameters
     * @return Evaluation of the function
    */
    std::function<T::VariableType(T::VariableType b, T::IndexType kkk, T::VariableType time_to_i, T::VectorXdr& v_parameters_)> f_ijk;
    
    /**
     * This function, implemented as a lambda function, compute the approximation of the second derivative of the 
     * log-likelihood function, along a single dimension.
     * @param index_ Index of the direction with respect to compute the approximation of the second derivative
     * @param v_parameters_ Vector of parameters
     * @return Evaluation of the directional second directive at the vector of parameters
    */
    std::function<T::VariableType(T::IndexType, T::VectorXdr&)> dd_ll_lf;

    //! Quadrature Points
    QuadraturePoints10 points10;                                    //! Struct for the ten points of the quadrature rule
    T::NumberType n_nodes = 10;                                             //! Number of points
    std::array<T::VariableType, 10>& nodes = points10.nodes;                //! Array of nodes
    std::array<T::VariableType, 10>& weights = points10.weights;            //! Array of weights

    /**
     * Method for computing thr number of parameters of this model. It initalizes the private variable n_parameters.
    */
    void compute_n_parameters() noexcept override;
    
    /**
     * Method for computing the diagonal of the hessian matrix of the log-likelihood function.
     * @param v_parameters_ Vector of parameters, where the hessian must be evaluated
    */
    void compute_hessian_diagonal(T::VectorXdr& v_parameters_) noexcept override;

    /**
     * Method for computing the standard error of the parameters.
     * @param v_parameters_ Vector of parameers, whose standard error must be produced
    */
    void compute_se(T::VectorXdr& v_parameters_) noexcept override;

    /**
     * Method for computing the standard deviation of the frailty.
     * @param v_parameters_ Vector of parameters
    */
    void compute_sd_frailty(T::VectorXdr& v_parameters_) noexcept override;
    
    /**
     * Method for extracting each single parameter of the model or a single vector of similar parameters.
     * @param v_parameters_ Vector of parameters
     * @return A customized model tuple containing all the extracted parameters
    */
    T::TupleLFType extract_parameters(T::VectorXdr& v_parameters_) const noexcept;

    /**
     * Method for extracting the variables related to the dropout, for each group/cluster of individuals.
     * @param indexes_group_ Vector of indexes of the dataset, of the individuals belonging to a cluster/group/faculty
     * @return A customized tuple containing the dropout variables to be extracted
    */
    T::TupleDropoutType extract_dropout_variables(const T::SharedPtrType indexes_group_) const noexcept;

    /**
     * Method for extracting the time-to-event variable of individuals belonging to a cluster.
     * @param indexes_group_ Vector of indexes fo the dataset, of individuals belonging to the cluster
     * @return Vector of time-to-event
    */
    T::VectorXdr extract_time_to_event(const T::SharedPtrType indexes_group) const noexcept;

    /**
     * Method for building the model log-likelihood function. It contains the definition
     * of the functions ll_lf, ll_group_lf, G, f_ijk.
    */
    void build_loglikelihood() noexcept override;

    /**
     * Method for building the model log-likelihood function in parallel. It contains the
     * parallel definition of the function ll_lf.
    */
    void build_loglikelihood_parallel() noexcept override;

    /**
     * Method for building the second derivative of the log-likelihood function. 
     * It uses a centered finite difference scheme, with accuracy of the secon order.
    */
    void build_dd_loglikelihood() noexcept override;
};

}
#endif // MODEL_DERIVED_HPP



