/**
 * @file test_simd_special.cpp
 *
 * @brief Unit testing of SIMD special functions
 *
 * @author Dahua Lin
 */

#include "test_simd_math_base.h"
#include <light_mat/math/simd_math.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

const int TTimes = 10;


// error functions

DEFINE_MATH_TPACK1( erf, 4, -2.0, 2.0 )
DEFINE_MATH_TPACK1( erfc, 4, -2.0, 2.0 )

// norminv

DEFINE_MATH_TPACK1( norminv, 4, 0.0, 1.0 )

BEGIN_MAIN_SUITE
	ADD_TPACK( simd_erf )
	ADD_TPACK( simd_erfc )
	ADD_TPACK( simd_norminv )
END_MAIN_SUITE


