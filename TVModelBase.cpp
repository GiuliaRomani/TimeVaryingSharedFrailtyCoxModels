// Include header files
#include "TVModelBase.hpp"

namespace TVModel{

ModelBase::ModelBase(const T::FileNameType& filename1, const T::FileNameType& filename2): 
		Dataset(filename1, filename2){

    // Initialize the other simple data structire
    variance_frailty.resize(Dataset::n_intervals);
    sd_frailty.resize(Dataset::n_intervals);
    
    // Correct dimension will be set up for each method
    hessian_diag.resize(0);
    se.resize(0);
};


} // end namespace

