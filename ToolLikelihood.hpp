#ifndef TOOL_LIKELIHOOD
#define TOOL_LIKELIHOOD

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>

namespace ToolsLikelihood{
using T = TypeTraits;

struct Tools{
    T::VariableType factor_c_pp = 1e+100;
    T::VariableType factor_c_paik = 1e+100;
    T::VariableType factor_c_lf = 1e+100;
};

}
#endif // TOOL_LIKELIHOOD