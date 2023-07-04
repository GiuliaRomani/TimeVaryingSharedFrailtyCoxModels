#include "TypeTraits.hpp"

#include <iostream>
#include <cmath>

using T = TypeTraits;

/*
T::VariableType factorial(T::VariableType n){
    if(n <= 1)
        return 1;
    else{
        return (n * factorial(n-1));
    }
};

T::VariableType binom(T::NumberType a, T::NumberType b){
    return factorial(a)/(factorial(b) * factorial(a - b));
};
*/



// These two functions are used to compute the binomial coefficient 
inline double logbinom(double n, double k){
    return std::lgamma(n+1) - std::lgamma(n-k+1) - std::lgamma(k+1);
};

inline double binom(double n, double k){
    return std::exp(logbinom(n,k));
};

