#ifndef QUADRATURENODES_HPP
#define QUADRATURENODES_HPP

#include "TypeTraits.hpp"


using T = TypeTraits;

namespace QuadratureNodes{
    // QuadratureFormula with 9 points
    struct Points9{
        T::NumberType n = 9;
        // T::VectorXd nodes << ;
        // T::VectorXd weights << ;

    };

    // QuadratureFormula with 10 points
    struct Points10{
        T::NumberType n = 10;
        // T::VectorXd nodes{};
        // T::VectorXd weights{};

    };
};


#endif //QUADRATURENODES_HPP