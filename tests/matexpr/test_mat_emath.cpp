/**
 * @file test_mat_emath.cpp
 *
 * Unit test of elementary functions on matrices
 *
 * @author Dahua Lin
 */

#include "matfun_test_base.h"

#include <light_mat/matexpr/mat_emath.h>

using namespace lmat;
using namespace lmat::test;


// power functions

DEFINE_MATFUN_TESTS2(pow, 0.5, 2.0, 0.5, 2.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(cbrt, -10.0, 10.0, 1.0e-14)
DEFINE_MATFUN_TESTS2(hypot, -3.0, 3.0, -3.0, 3.0, 1.0e-14)

// exp & log

DEFINE_MATFUN_TESTS1(exp, -1.0, 3.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(log, 0.5, 10.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(log10, 0.5, 20.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(xlogx, -1.0, 3.0, 1.0e-14)
DEFINE_MATFUN_TESTS2(xlogy, -1.0, 3.0, 1.0, 4.0, 1.0e-14)

DEFINE_MATFUN_TESTS1(exp2, -1.0, 3.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(log2, 0.5, 10.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(expm1, -1.0, 1.0, 1.0e-15)
DEFINE_MATFUN_TESTS1(log1p, 0.1, 1.0, 1.0e-15)

// trigonometry

DEFINE_MATFUN_TESTS1(sin, -3.0, 3.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(cos, -3.0, 3.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(tan, -1.0, 1.0, 1.0e-14)

DEFINE_MATFUN_TESTS1(asin, -1.0, 1.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(acos, -1.0, 1.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(atan, -1.0, 1.0, 1.0e-14)
DEFINE_MATFUN_TESTS2(atan2, -3.0, 3.0, -3.0, 3.0, 1.0e-14)

// hyperbolic

DEFINE_MATFUN_TESTS1(sinh, -2.0, 2.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(cosh, -2.0, 2.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(tanh, -2.0, 2.0, 1.0e-14)

DEFINE_MATFUN_TESTS1(asinh, -2.0, 2.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(acosh, 1.0, 5.0, 1.0e-14)
DEFINE_MATFUN_TESTS1(atanh, -1.0, 1.0, 1.0e-14)


