#ifndef QUADRATURENODES_HPP
#define QUADRATURENODES_HPP

#include "TypeTraits.hpp"

#include <array>

namespace TVSFCM{
using T = TypeTraits;

/**
 * Two structs containing each one a predefined number of nodes and weights, for the application
 * of the Gauss-Hermite quadrature formula.
 * Each one is composed of two arrays, one for the nodes and the other for the weights, and the number of points.
*/

    //! QuadratureFormula with 9 points
    struct QuadraturePoints9{
        T::NumberType n = 9;
        std::array<T::VariableType, 9> nodes{0., -0.723551018752838, 0.723551018752838, 
                                            -1.46855328921667, 1.46855328921667,-2.26658058453184, 
                                            2.26658058453184, -3.19099320178153, 3.19099320178153};
        std::array<T::VariableType, 9> weights{0.7202352156061, 0.43265155900261, 0.4326515590026, 
                                                0.08847452739438, 0.08847452739438, 0.004943624275537, 
                                                0.004943624275537,3.960697726326e-05,3.960697726326e-05};
    };

    //! QuadratureFormula with 10 points
    struct QuadraturePoints10{
        T::NumberType n = 10;
        std::array<T::VariableType, 10> nodes{0.342901327223705,-0.342901327223705,1.03661082978951,
                                            -1.03661082978951,1.75668364929988,-1.75668364929988,
                                            2.53273167423279,-2.53273167423279,3.43615911883774,
                                            -3.43615911883774};
        std::array<T::VariableType, 10> weights{0.6108626337353,0.6108626337353,0.2401386110823,
                                            0.2401386110823,0.03387439445548,0.03387439445548,
                                            0.001343645746781,0.001343645746781,7.640432855233e-06,
                                            7.640432855233e-06};
    };
};


#endif //QUADRATURENODES_HPP