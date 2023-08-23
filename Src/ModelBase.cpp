// Include header files
#include "ModelBase.hpp"
#include "GetPot"
#include "MyException.hpp"

namespace TVSFCM{

ModelBase::ModelBase(const T::FileNameType& filename1, const T::FileNameType& filename2): 
	    Dataset(filename1, filename2),
        ParallelComponents(filename1){

    //! Initialize GetPot and read the discretization step. Then check its positiveness
    GetPot datafile(filename1.c_str());
    h_dd = datafile("DiscretizationStep/h_dd", 1e-3);
    check_condition(h_dd);

    //! Initialize the frailty variance and the standard deviation to the number of intervals of the time-domain
    variance_frailty.resize(Dataset::n_intervals);
    sd_frailty.resize(Dataset::n_intervals);
    
    //! Initialize to 0 the dimension of the hessian matrix and the vector of standard error, because they depend
    //! on the parameter vector
    hessian_diag.resize(0);
    se.resize(0);
};

//! Method for checking that the discretization step is positive
void ModelBase::check_condition(T::VariableType h_dd_) const{
    if(h_dd_ < 0)
        throw MyException("Provided negative discretization step.");
};

//! Method for printing the results of the model application
void ModelBase::print_results() const noexcept{
    result.print_results();
};

} // end namespace

