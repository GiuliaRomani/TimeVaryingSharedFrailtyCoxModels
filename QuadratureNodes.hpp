#ifndef QUADRATURENODES_HPP
#define QUADRATURENODES_HPP

#include <Eigen/Dense>

#include "TypeTraits.hpp"

using T = TypeTraits;

namespace QuadratureNodes{
    // QuadratureFormula with 9 points
    struct Points9{
        T::NumberType n = 9;

    };

    // QuadratureFormula with 10 points
    struct Points10{
        T::NumberType n = 10;

    };
};


#endif //QUADRATURENODES_HPP