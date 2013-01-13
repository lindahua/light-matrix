/**
 * @file test_simd_math.cpp
 *
 * @brief Unit testing for elementary math on SIMD vectors
 *
 * @author Dahua Lin
 */


#include "test_simd_math_base.h"
#include <light_mat/math/simd_math.h>

using namespace lmat;
using namespace lmat::test;

const int TTimes = 10;

LTEST_INIT_AUTOSUITE

// power

DEFINE_MATH_TPACK2( pow,   4,  0.5,  3.0, 0.5, 2.0 )
DEFINE_MATH_TPACK2( hypot, 4, -5.0,  5.0, -5.0, 5.0 )
DEFINE_MATH_TPACK1( cbrt,  4,  1.0, 20.0 )

// exp & log

DEFINE_MATH_TPACK1( exp,   4, -1.0,  2.0 )
DEFINE_MATH_TPACK1( log,   4,  1.0, 10.0 )
DEFINE_MATH_TPACK1( log10, 4,  1.0, 50.0 )

DEFINE_MATH_TPACK1( exp2,  4, -1.0,  2.0 )
DEFINE_MATH_TPACK1( log2,  4,  1.0, 10.0 )
DEFINE_MATH_TPACK1( expm1, 4, -0.5,  0.5 )
DEFINE_MATH_TPACK1( log1p, 4, -0.5,  0.5 )

DEFINE_MATH_TPACK2( xlogy, 4, -1.0, 1.0, 0.0, 1.0 )

// trigonometry

DEFINE_MATH_TPACK1( sin, 4, -3.0, 3.0 )
DEFINE_MATH_TPACK1( cos, 4, -3.0, 3.0 )
DEFINE_MATH_TPACK1( tan, 4,  0.0, 1.2 )

DEFINE_MATH_TPACK1( asin, 4, -1.0, 1.0 )
DEFINE_MATH_TPACK1( acos, 4, -1.0, 1.0 )
DEFINE_MATH_TPACK1( atan, 4, -2.0, 2.0 )
DEFINE_MATH_TPACK2( atan2, 4, -3.0, 3.0, -3.0, 3.0 )

// hyperbolic

DEFINE_MATH_TPACK1( sinh, 4, -2.0, 2.0 )
DEFINE_MATH_TPACK1( cosh, 4, -2.0, 2.0 )
DEFINE_MATH_TPACK1( tanh, 4, -2.0, 2.0 )

DEFINE_MATH_TPACK1( asinh, 4, -5.0, 5.0 )
DEFINE_MATH_TPACK1( acosh, 4,  1.0, 4.0 )
DEFINE_MATH_TPACK1( atanh, 4, -0.9, 0.9 )

