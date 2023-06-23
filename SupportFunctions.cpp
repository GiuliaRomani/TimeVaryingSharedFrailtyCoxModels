#include <iostream>
#include <cmath>

// These two functions are used to cpmute the binomial coefficient 
inline double logbinom(double n, double k){
    return std::lgamma(n+1) - std::lgamma(n-k+1) - std::lgamma(k+1);
};

inline double binom(double n, double k){
    return std::exp(logbinom(n,k));
};
