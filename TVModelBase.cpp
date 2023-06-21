
#include "TVModelBase.hpp"

namespace TVModel{
using T = TypeTraits;

ModelBase::ModelBase(const T::FileNameType& filename1, const T::FileNameType& filename2){
    // Construct the complex data structures using their constructor
    database = DatasetInfoClass::DatasetInfo(filename1, filename2);
    time = TimeDomainInfo::TimeDomain(filename1);
    //result = ResultsMethod::Results();

    // Initialize the other simple data structire
    variance_frailty.resize(n_intervals);
    sd_frailty.resize(n_intervals);
};


} // end namespace

