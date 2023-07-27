// Include header files
#include "TVModelBase.hpp"
#include "GetPot"
#include "MyException.hpp"

namespace TVSFCM{

ModelBase::ModelBase(const T::FileNameType& filename1, const T::FileNameType& filename2): 
		Dataset(filename1, filename2){

    // Before extracting the discretization step from the input file, check that it exists.
    check_filename(filename1);
    GetPot datafile(filename1.c_str());
    h_dd = datafile("DiscretizationStep/h_dd",1e-3);
    n_threads = static_cast<T::NumberType>(datafile("NumberThreads/n_threads", 1));

    // Initialize the other simple data structire
    variance_frailty.resize(Dataset::n_intervals);
    sd_frailty.resize(Dataset::n_intervals);
    
    // Correct dimension will be set up for each method
    hessian_diag.resize(0);
    se.resize(0);
};

// Method for checking the filename is correct and exists
void ModelBase::check_filename(const T::FileNameType& filename_) const{
    std::ifstream check(filename_);
    if(check.fail()){
        T::ExceptionType msg1 = "File ";
        T::ExceptionType msg2 = msg1.append((filename_).c_str());
        T::ExceptionType msg3 = msg2.append(" does not exist.");
        //throw MyException("File provided does not exist.");
        throw MyException(msg3);
    }
};


} // end namespace

