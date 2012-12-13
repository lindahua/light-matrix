/**
 * @file simd_test_base.h
 *
 * @brief
 *
 * @author Dahua Lin
 */

#ifndef SIMD_TEST_BASE_H_
#define SIMD_TEST_BASE_H_

#include "test_base.h"
#include <light_mat/math/simd_debug.h>

#define ASSERT_SIMD_EQ( a, b ) \
	if (!::lmat::math::test_equal(a, b)) throw ::ltest::assertion_failure(__FILE__, __LINE__, #a " == " #b)


#endif /* SIMD_TEST_BASE_H_ */
