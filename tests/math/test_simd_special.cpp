/**
 * @file test_simd_special.cpp
 *
 * @brief Unit testing of SIMD special functions
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/math/simd_special.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

using lmat::math::simd_pack;
using lmat::math::sse_t;
using lmat::math::avx_t;

const int TTimes = 10;

inline double randunif(double LB, double UB)
{
	double u = (double)std::rand() / RAND_MAX;
	return LB + (UB - LB) * u;
}


/************************************************
 *
 *  Macros to define test cases
 *
 ************************************************/

#define DEFINE_MATH_TEST_1( Name, SIMD, ulp, LB, UB ) \
	T_CASE( SIMD##_special, Name##_##SIMD ) { \
		typedef simd_pack<T, SIMD##_t> pack_t; \
		const unsigned int width = pack_t::pack_width; \
		T a_src[width]; \
		T r0[width]; \
		for (int t = 0; t < TTimes; ++t) { \
			for (unsigned i = 0; i < width; ++i) { \
				a_src[i] = T(randunif(LB, UB)); \
				r0[i] = math::Name(a_src[i]); \
			} \
			pack_t a; a.load_u(a_src); \
			pack_t r = math::Name(a); \
			pack_t r_emu = math::internal::Name##_emulate(a); \
			ASSERT_SIMD_ULP(r, r0, ulp); \
			ASSERT_SIMD_ULP(r_emu, r0, 1); \
		} \
	}


#ifdef LMAT_HAS_AVX

#define DEFINE_MATH_TPACK( Name ) \
	BEGIN_TPACK( simd_##Name ) \
		ADD_T_CASE( sse_special, Name##_sse, float ) \
		ADD_T_CASE( sse_special, Name##_sse, double ) \
		ADD_T_CASE( avx_special, Name##_avx, float ) \
		ADD_T_CASE( avx_special, Name##_avx, double ) \
	END_TPACK
#else
#define DEFINE_MATH_TPACK( Name ) \
	BEGIN_TPACK( simd_##Name ) \
		ADD_T_CASE( sse_special, Name##_sse, float ) \
		ADD_T_CASE( sse_special, Name##_sse, double ) \
	END_TPACK
#endif

/************************************************
 *
 *  specific cases
 *
 ************************************************/

#ifdef LMAT_HAS_CXX11_MATH

// error functions

DEFINE_MATH_TEST_1( erf,  sse, 4, -2.0, 2.0 )
DEFINE_MATH_TEST_1( erf,  avx, 4, -2.0, 2.0 )
DEFINE_MATH_TPACK( erf )

DEFINE_MATH_TEST_1( erfc, sse, 4, -2.0, 2.0 )
DEFINE_MATH_TEST_1( erfc, avx, 4, -2.0, 2.0 )
DEFINE_MATH_TPACK( erfc )

// gamma functions

DEFINE_MATH_TEST_1( lgamma, sse, 4, 1.0, 5.0 )
DEFINE_MATH_TEST_1( lgamma, avx, 4, 1.0, 5.0 )
DEFINE_MATH_TPACK( lgamma )

DEFINE_MATH_TEST_1( tgamma, sse, 4, 1.0, 5.0 )
DEFINE_MATH_TEST_1( tgamma, avx, 4, 1.0, 5.0 )
DEFINE_MATH_TPACK( tgamma )

#endif


BEGIN_MAIN_SUITE
	ADD_TPACK( simd_erf )
	ADD_TPACK( simd_erfc )
	ADD_TPACK( simd_lgamma )
	ADD_TPACK( simd_tgamma )
END_MAIN_SUITE


