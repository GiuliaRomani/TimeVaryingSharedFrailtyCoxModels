#ifndef TOOL_LIKELIHOOD
#define TOOL_LIKELIHOOD

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>

namespace ToolsLikelihood{
using T = TypeTraits;

struct Tools{
    T::VariableType factor_c_pp = 1;
    T::VariableType factor_c_paik = 1;
    T::VariableType factor_c_lf = 1;
    
    T::VariableType tol_ll = 1e-4;
    T::VariableType tol_optim = 1e-4;
    
    T::NumberType n_extrarun = 10;
};

}
#endif // TOOL_LIKELIHOOD
