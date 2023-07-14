#ifndef TOOL_LIKELIHOOD
#define TOOL_LIKELIHOOD

// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>

/**
 * Tools is a struct containing some variables, that could be used insise the time-varying models.
 * Actually, it only contains the discretization step for the computation of the second derivative of the log-likelihood functions
*/
namespace ToolsLikelihood{
using T = TypeTraits;

struct Tools{
    T::VariableType h_dd = 1e-3;            // Discretization step for the second derivative
};

}
#endif // TOOL_LIKELIHOOD
