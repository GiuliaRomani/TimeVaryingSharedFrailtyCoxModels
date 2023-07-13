#ifndef TOOL_LIKELIHOOD
#define TOOL_LIKELIHOOD

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>

namespace ToolsLikelihood{
using T = TypeTraits;

struct Tools{
    T::VariableType h_dd = 1e-3;
};

}
#endif // TOOL_LIKELIHOOD
