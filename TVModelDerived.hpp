
#ifndef TIME_VARYING_MODEL_DERIVED
#define TIME_VARYING_MODEL_DERIVED

// Include header files
#include "TVModelBase.hpp"
#include "QuadraturePoints.hpp"
#include "Parameters.hpp"

// Include libraries
#include <iostream>
#include <functional>

/**
 * For each model, a class publicly derived from the ModelBase class is created. This implements one of the three
 *  Time-Varying Shared Frailty Cox model.
 * 
 * It is composed of the Parameters class, that depends on the model we choose and on both the dataset and the time domain. 
 * This class is thus initialized only after having initialized the other two base classes.
*/


namespace TVModel{
using T = TypeTraits;
using Parameters = Params::Parameters;

// Class for implementing PowerParameter Model
class PowerParameterModel final: public TVModel::ModelBase,
                                 public Parameters{
public:
    /**
     * Deafult constructor
    */
    PowerParameterModel() = default;

    /**
     * Constructor that initializes firstly the base classes and then compute the number of the parameters
     * of this model, to be able to resize all the other variables and instantiate an object of the class Parameters.
     * Finally, the model log-likelihood and its multi-directional second derivative are built.
    */
    PowerParameterModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    /**
     * Method for computing the value of the log-likelihood, given the optimal vector of parameters. 
     * 
     * Eventually, it initializes the result class.
    */
    void evaluate_loglikelihood() override;

    /**
     * Virtual destructor 
    */
    virtual ~ PowerParameterModel() = default;

private:
    T::NumberType n_parameters;                                             //! Number of model parameters         
    T::VectorXdr & v_parameters = Parameters::v_parameters;                 //! Vector of parameters

    T::IdNameType name_method = "PowerParameter";                           //! Name of this method

    /**
     * This function is later built as a lambda function. 
     * It implememts the model log-likelihood as has been defined in the reference and in the report.
     * 
     * Inside this function, other methods or functions are called: 
     * 
     * 1) extract_parameters for extracting the vector of parameters
     * 2) ll_group_pp to compute the log-likelihood related to a cluster/group
     * 
     * @param v_parameters Vector of parameters, where the log-likelihood has to be evaluated
     * @return Evaluation of the log-likelihood at the input, for all the clusters
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_pp;  

    /**
     * This function is later built as a lambda function.
     * It implements the model log-likelihood as has been defined in the reference, but in a parallel version
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
    
    // Function for computing the second derivative wrt one dimension
    /**
     * This function, implemented as a lambda function, compute the approximation of the second derivative of the 
     * log-likelihood function, along a single dimension.
     * @param index_ Index of the direction with respect to compute the approximation of the second derivative
     * @param v_parameters_ Vector of parameters
     * @return Evaluation of the directional second directive at the vector of parameters
    */
    std::function<T::VariableType(T::IndexType index_, T::VectorXdr& v_parameters_)> dd_ll_pp;	


    QuadraturePoints::Points9 points9;                                      //! Struct for the nine points quadrature formula
    T::NumberType n_nodes = 9;                                              //! Number of points 
    std::array<T::VariableType, 9>& nodes = points9.nodes;                  //! Array of nodes
    std::array<T::VariableType, 9>& weights = points9.weights;              //! Array of weights


    /**
     * Method for computing thr number of parameters of this model. It initalizes the private variable n_parameters
    */
    void compute_n_parameters() override;

    /**
     * Method for extracting each single parameter of the model or a single vector of similar parameters
     * @param v_parameters_ Vector of parameters
     * @return A customized model tuple containing all the extracted parameters
    */
    T::TuplePPType extract_parameters(const T::VectorXdr& v_parameters_) const ;
    
    /**
     * Method for computing the diagonal of the hessian matrix of the log-likelihood function.
     * @param v_parameters_ Vector of parameters, where the hessian must be evaluated
    */
    void compute_hessian_diagonal(T::VectorXdr& v_parameters_) override;

    /**
     * Method for computing the standard error of the parameters
     * @param v_parameters_ Vector of parameters, whose standard error must be produced
    */
    void compute_se(T::VectorXdr& v_parameters) override;

    /**
     * Method for computing the standard deviation of the frailty 
     * @param v_parameters_ Vector of parameters
    */
    void compute_sd_frailty(T::VectorXdr& v_parameters_) override;

    /**
     * Method for building the model log-likelihood function.
    */
    void build_loglikelihood() override;

    /**
     * Virtual pure method for constructing the log-likelihood in the parallel version
    */
    void build_loglikelihood_parallel() override;

    /**
     * Method for building the second derivative of the log-likelihood function. 
     * It uses a centered finite difference scheme, with accuracy of the secon order.
    */
    void build_dd_loglikelihood() override;
};


//--------------------------------------------------------------------------------------------------
// Class for implementing PAIK Model
class PaikModel final: public TVModel::ModelBase, 
                       public Parameters{
public:
    /**
     * Default Constructor
    */
    PaikModel() = default;

    /**
     * Constructor that initializes firstly the base classes and then compute the number of the parameters
     * of this model, to be able to resize all the other variables and instantiate an object of the class Parameters.
     * Finally, the model log-likelihood and its multi-directional second derivative are built.
    */
    PaikModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    /**
     * Method for computing the value of the log-likelihood, given the optimal vector of parameters. 
     * 
     * Eventually, it initializes the result class.
    */
    void evaluate_loglikelihood() override;

    virtual ~ PaikModel() = default;

private:
    T::NumberType n_parameters;                                             //! Number of parameter of the model       
    T::VectorXdr& v_parameters = Parameters::v_parameters;                  //! Vector of parameters

    T::IdNameType name_method = "Paik";                                     //! Name of this method

    /**
     * This function is later built as a lambda function. 
     * It implememts the model log-likelihood as has been defined in the reference and in the report.
     * 
     * Inside this function, other methods or functions are called: 
     * 
     * 1) extract_parameters for extracting the vector of parameters
     * 2) ll_group_pp to compute the log-likelihood related to a cluster/group
     * 
     * @param v_parameters_ Vector of parameters, where the log-likelihood has to be evaluated
     * @return Evaluation of the log-likelihood at the input, for all the clusters
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_paik; 

    /**
     * This function is later built as a lmabda function.
     * It implelemnts the model log-likelihood as has been defined in the reference and in the report, 
     * but using parallel implementation provided by OpenMP with a number of threads specified by the user. 
     * If no threads or just 1 is provided, no parallel implementation is chosen.
     * 
     * @param v_parameters_ Vector of parameters, where the log-likelihood has to be evaluated
     * @return The evaluation of the log-likelihood
    */
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_paik_parallel; 

    /**
     * This function compute the log-likelihood for all the individuals belonging to a cluster/group/faculty.
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
     * Method for computing thr number of parameters of this model. It initalizes the private variable n_parameters
    */
    void compute_n_parameters() override;

    /**
     * Method for computing the diagonal of the hessian matrix of the log-likelihood function.
     * @param v_parameters_ Vector of parameters, where the hessian must be evaluated
    */
    void compute_hessian_diagonal(T::VectorXdr& v_parameters_);

    /**
     * Method for computing the standard error of the parameters
     * @param v_parameters_ Vector of parameers, whose standard error must be produced
    */
    void compute_se(T::VectorXdr& v_parameters_);

    /**
     * Method for computing the standard deviation of the frailty 
     * @param v_parameters_ Vector of parameters
    */
    void compute_sd_frailty(T::VectorXdr& v_parameters_) override;
    
    /**
     * Method for extracting each single parameter of the model or a single vector of similar parameters
     * @param v_parameters_ Vector of parameters
     * @return A customized model tuple containing all the extracted parameters
    */
    T::TuplePaikType extract_parameters(T::VectorXdr& v_parameters_);

    /**
     * Method for extracting the variables related to the dataset, for each group/cluster of individuals
     * @param indexes_group_ Vector of indexes of the dataset, of the indivividuals belonging to a cluster/group/faculty
     * @param phi_ Vector of baseline parameters
     * @param betar_ Vector of regressor coefficients
     * @return A customized tuple containing the variables to be extracted
    */
    T::TupleMatrixAType extract_matrixA_variables(T::SharedPtrType indexes_group_, T::VectorXdr& phi_, T::VectorXdr& betar_);

    /**
     * Method for extracting the variables related to the dropout, for each group/cluster of individuals
     * @param indexes_group_ Vector of indexes of the dataset, of the individuals belonging to a cluster/group/faculty
     * @return A customized tuple containing the dropout variables to be extracted
    */
    T::TupleDropoutType extract_dropout_variables(T::SharedPtrType indexes_group_);

    /**
     * Method for building the model log-likelihood function.
    */
    void build_loglikelihood() override;

    /**
     * Method for building the model log-likelihood function in parallel
    */
    void build_loglikelihood_parallel() override;

    /**
     * Method for building the second derivative of the log-likelihood function. 
     * It uses a centered finite difference scheme, with accuracy of the secon order.
    */
    void build_dd_loglikelihood() override;
};


//-----------------------------------------------------------------------------------------
// Class for implementing LOG FRAILTY Model
class LogFrailtyModel final: public TVModel::ModelBase,
                             public Parameters{
public:
    /**
     * Default constructor
    */
    LogFrailtyModel() = default;

    /**
     * Constructor that initializes firstly the base classes and then compute the number of the parameters
     * of this model, to be able to resize all the other variables and instantiate an object of the class Parameters.
     * Finally, the model log-likelihood and its multi-directional second derivative are built.
    */
    LogFrailtyModel(const T::FileNameType& filename1, const T::FileNameType& filename2);

    /**
     * Method for computing the value of the log-likelihood, given the optimal vector of parameters. 
     * 
     * Eventually, it initializes the result class.
    */
    void evaluate_loglikelihood() override;

    /**
     * Virtual destructor
    */
    virtual ~ LogFrailtyModel() = default;

private:
    T::NumberType n_parameters;                                             //! Number of parameter of the model 
    T::VectorXdr& v_parameters = Parameters::v_parameters;                  //! Vector of parameters

    T::IdNameType name_method = "LogFrailty";                               //! Name of this method

    /**
     * This function is later built as a lambda function. 
     * It implememts the model log-likelihood as has been defined in the reference and in the report.
     * 
     * Inside this function, other methods or functions are called: 
     * 
     * 1) extract_parameters for extracting the vector of parameters
     * 2) ll_group_pp to compute the log-likelihood related to a cluster/group
     * 
     * @param v_parameters Vector of parameters, where the log-likelihood has to be evaluated
     * @return Evaluation of the log-likelihood at the input, for all the clusters
    */    
    std::function<T::VariableType(T::VectorXdr& v_parameters_)> ll_lf; 

    /**
     * This function is later built as a lambda function. 
     * It implements the model log-likelihood as has been defined in the reference and in the report, in a parallel version
     * using OpenMP and a number of threads specified by the user. If no threads or just 1 are provided, no parallelization is used.
     * 
     * Inside this function, other methods or functions are called: 
     * 
     * 1) extract_parameters for extracting the vector of parameters
     * 2) ll_group_pp to compute the log-likelihood related to a cluster/group
     * 
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
     * This function, implemented as a lambda function, builds the function G  defined in the reference and report.
     * @param z Argument of the function
     * @param indexes_group_ Vector of the indexes of the dataset, of the individuals in a cluster
     * @param v_parameters_ Vector of parameters
     * @return Evaluation of the function
    */
    std::function<T::VariableType(T::VariableType z, T::SharedPtrType indexes_group_, T::VectorXdr& v_parameters_)> G; 

    /**
     * This function, implemented as a lambda function, builds the f_ijk function defined in the reference and in the report
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


    QuadraturePoints::Points10 points10;                                    //! Struct for the ten points of the quadrature rule
    T::NumberType n_nodes = 10;                                             //! Number of points
    std::array<T::VariableType, 10>& nodes = points10.nodes;                //! Array of nodes
    std::array<T::VariableType, 10>& weights = points10.weights;            //! Array of weights

    /**
     * Method for computing thr number of parameters of this model. It initalizes the private variable n_parameters
    */
    void compute_n_parameters() override;
    
    /**
     * Method for computing the diagonal of the hessian matrix of the log-likelihood function.
     * @param v_parameters_ Vector of parameters, where the hessian must be evaluated
    */
    void compute_hessian_diagonal(T::VectorXdr& v_parameters_);

    /**
     * Method for computing the standard error of the parameters
     * @param v_parameters_ Vector of parameers, whose standard error must be produced
    */
    void compute_se(T::VectorXdr& v_parameters_);

    /**
     * Method for computing the standard deviation of the frailty 
     * @param v_parameters_ Vector of parameters
    */
    void compute_sd_frailty(T::VectorXdr& v_parameters_) override;
    
    /**
     * Method for extracting each single parameter of the model or a single vector of similar parameters
     * @param v_parameters_ Vector of parameters
     * @return A customized model tuple containing all the extracted parameters
    */
    T::TupleLFType extract_parameters(T::VectorXdr& v_parameters_) const;

    /**
     * Method for extracting the variables related to the dropout, for each group/cluster of individuals
     * @param indexes_group_ Vector of indexes of the dataset, of the individuals belonging to a cluster/group/faculty
     * @return A customized tuple containing the dropout variables to be extracted
    */
    T::TupleDropoutType extract_dropout_variables(const T::SharedPtrType indexes_group_) const;

    /**
     * Method for extracting the time-to-event variable of individuals belonging to a cluster
     * @param indexes_group_ Vector of indexes fo the dataset, of individuals belonging to the cluster
     * @return Vector of time-to-event
    */
    T::VectorXdr extract_time_to_event(const T::SharedPtrType indexes_group) const;

    /**
     * Method for building the model log-likelihood function.
    */
    void build_loglikelihood() override;

    /**
     * Method for building the model log-likelihood function in parallel
    */
    void build_loglikelihood_parallel() override;

    /**
     * Method for building the second derivative of the log-likelihood function. 
     * It uses a centered finite difference scheme, with accuracy of the secon order.
    */
    void build_dd_loglikelihood() override;
};



}
#endif // TIME_VARYING_MODEL_DERIVED



