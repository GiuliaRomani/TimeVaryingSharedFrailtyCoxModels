// Include header files
#include "TypeTraits.hpp"

// Include libraries
#include <iostream>
#include <cmath>

namespace TVSFCM{
using T = TypeTraits;

/**
 * The function computes the natual logarithm of the binomial coefficient between two integer numbers,
 * using the fact that the gamma function of an (integer + 1) corresponds to the factorial of that integer.
 * @param n First term of the binomial coefficient
 * @param k Second term of the binomial coefficient
 * @return Logarithm of the binomial coefficient
*/
inline T::VariableType logbinom(T::NumberType n, T::NumberType k) noexcept{
    return std::lgamma(n+1) - std::lgamma(n-k+1) - std::lgamma(k+1);
};

/**
 * Once the logarithmic binomial coefficient has been computed, we evaluate its exponential to get the the pure binomial coefficient.
 * Inside the function, the logbinom(n,k) is called
 * @param n First term of the binomial coefficient
 * @param k Second term of the binomial coefficient
 * @return Binomial coefficient
*/
inline T::VariableType binom(T::NumberType n, T::NumberType k) noexcept{
    return std::exp(logbinom(n,k));
};

} // end namespace

