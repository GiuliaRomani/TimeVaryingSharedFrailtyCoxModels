#include "TypeTraits.hpp"

#include <iostream>
#include <cmath>

using T = TypeTraits;


// These two functions are used to compute the binomial coefficient 
inline T::VariableType logbinom(T::NumberType n, T::NumberType k){
    return std::lgamma(n+1) - std::lgamma(n-k+1) - std::lgamma(k+1);
};

inline T::VariableType binom(T::NumberType n, T::NumberType k){
    return std::exp(logbinom(n,k));
};


