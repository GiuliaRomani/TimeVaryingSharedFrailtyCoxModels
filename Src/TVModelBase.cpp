// Include header files
#include "TVModelBase.hpp"
#include "GetPot"
#include "MyException.hpp"

namespace TVSFCM{

ModelBase::ModelBase(const T::FileNameType& filename1, const T::FileNameType& filename2): 
		Dataset(filename1, filename2),
        ParallelComponents(filename1){

    GetPot datafile(filename1.c_str());
    h_dd = datafile("DiscretizationStep/h_dd", 1e-3);
    check_condition(h_dd);

    // Initialize the other simple data structire
    variance_frailty.resize(Dataset::n_intervals);
    sd_frailty.resize(Dataset::n_intervals);
    
    // Correct dimension will be set up for each method
    hessian_diag.resize(0);
    se.resize(0);
};


void ModelBase::check_condition(T::VariableType h_dd_){
    if(h_dd_ < 0)
        throw MyException("Provided negative discretization step.");
};

} // end namespace

