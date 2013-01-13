/**
 * @file test_mat_special.cpp
 *
 * @brief Unit testing of special functions on matrices
 *
 * @author Dahua Lin
 */

#include "matfun_test_base.h"

#include <light_mat/matexpr/mat_special.h>

using namespace lmat;
using namespace lmat::test;


// gauss related functions

DEFINE_MATFUN_TESTS1(erf, -2.0, 2.0, 1.0e-13)
DEFINE_MATFUN_TESTS1(erfc, -2.0, 2.0, 1.0e-13)
DEFINE_MATFUN_TESTS1(norminv, 0.0, 1.0, 1.0e-13)

// gamma related functions

DEFINE_MATFUN_TESTS1(lgamma, 1.0, 5.0, 1.0e-13)
DEFINE_MATFUN_TESTS1(tgamma, 1.0, 5.0, 1.0e-12)
