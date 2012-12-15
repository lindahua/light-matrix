/**
 * @file test_avx_math.cpp
 *
 * @brief Unit testing of AVX math functions
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/math/avx_math.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

using lmat::math::simd_pack;
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

#define DEFINE_AVX_MATH_TEST_1( Name, ulp, LB, UB ) \
	T_CASE( avx_math, Name ) { \
		typedef simd_pack<T, avx_t> pack_t; \
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
	} \
	BEGIN_TPACK( avx_math_##Name ) \
		ADD_T_CASE( avx_math, Name, float ) \
		ADD_T_CASE( avx_math, Name, double ) \
	END_TPACK

#define DEFINE_AVX_MATH_TEST_2( Name, ulp, LBa, UBa, LBb, UBb ) \
	T_CASE( avx_math, Name ) { \
		typedef simd_pack<T, avx_t> pack_t; \
		const unsigned int width = pack_t::pack_width; \
		T a_src[width]; \
		T b_src[width]; \
		T r0[width]; \
		for (int t = 0; t < TTimes; ++t) { \
			for (unsigned i = 0; i < width; ++i) { \
				a_src[i] = T(randunif(LBa, UBa)); \
				b_src[i] = T(randunif(LBb, UBb)); \
				r0[i] = math::Name(a_src[i], b_src[i]); \
			} \
			pack_t a; a.load_u(a_src); \
			pack_t b; b.load_u(b_src); \
			pack_t r = math::Name(a, b); \
			pack_t r_emu = math::internal::Name##_emulate(a, b); \
			ASSERT_SIMD_ULP(r, r0, ulp); \
			ASSERT_SIMD_ULP(r_emu, r0, 1); \
		} \
	} \
	BEGIN_TPACK( avx_math_##Name ) \
		ADD_T_CASE( avx_math, Name, float ) \
		ADD_T_CASE( avx_math, Name, double ) \
	END_TPACK


#define ADD_MATH_TPACK( Name ) ADD_TPACK( avx_math_##Name )


/************************************************
 *
 *  Specific cases
 *
 ************************************************/

// C++ 03

// power, exp, and log

DEFINE_AVX_MATH_TEST_2( pow,   4,  0.5,  3.0, 0.5, 2.0 )
DEFINE_AVX_MATH_TEST_1( exp,   4, -1.0,  2.0 )
DEFINE_AVX_MATH_TEST_1( log,   4,  1.0, 10.0 )
DEFINE_AVX_MATH_TEST_1( log10, 4,  1.0, 50.0 )

// trigonometry

DEFINE_AVX_MATH_TEST_1( sin, 4, -3.0, 3.0 )
DEFINE_AVX_MATH_TEST_1( cos, 4, -3.0, 3.0 )
DEFINE_AVX_MATH_TEST_1( tan, 4,  0.0, 1.2 )

// inverse trigonometry

DEFINE_AVX_MATH_TEST_1( asin, 4, -1.0, 1.0 )
DEFINE_AVX_MATH_TEST_1( acos, 4, -1.0, 1.0 )
DEFINE_AVX_MATH_TEST_1( atan, 4, -2.0, 2.0 )
DEFINE_AVX_MATH_TEST_2( atan2, 4, -3.0, 3.0, -3.0, 3.0 )

// hyperbolic

DEFINE_AVX_MATH_TEST_1( sinh, 4, -2.0, 2.0 )
DEFINE_AVX_MATH_TEST_1( cosh, 4, -2.0, 2.0 )
DEFINE_AVX_MATH_TEST_1( tanh, 4, -2.0, 2.0 )

// C++ 11

#ifdef LMAT_HAS_CXX11_MATH

// hypot & cbrt

DEFINE_AVX_MATH_TEST_2( hypot, 4, -5.0,  5.0, -5.0, 5.0 )
DEFINE_AVX_MATH_TEST_1( cbrt,  4,  1.0, 20.0 )

// extended exp & log

DEFINE_AVX_MATH_TEST_1( exp2,  4, -1.0,  2.0 )
DEFINE_AVX_MATH_TEST_1( log2,  4,  1.0, 10.0 )
DEFINE_AVX_MATH_TEST_1( expm1, 4, -0.5,  0.5 )
DEFINE_AVX_MATH_TEST_1( log1p, 4, -0.5,  0.5 )

// inverse hyperbolic

DEFINE_AVX_MATH_TEST_1( asinh, 4, -5.0, 5.0 )
DEFINE_AVX_MATH_TEST_1( acosh, 4,  1.0, 4.0 )
DEFINE_AVX_MATH_TEST_1( atanh, 4, -0.9, 0.9 )

// error functions

DEFINE_AVX_MATH_TEST_1( erf,  4, -2.0, 2.0 )
DEFINE_AVX_MATH_TEST_1( erfc, 4, -2.0, 2.0 )

// gamma functions

DEFINE_AVX_MATH_TEST_1( lgamma, 4, 1.0, 5.0 )
DEFINE_AVX_MATH_TEST_1( tgamma, 4, 1.0, 5.0 )

#endif


BEGIN_MAIN_SUITE

	ADD_MATH_TPACK( pow )
	ADD_MATH_TPACK( exp )
	ADD_MATH_TPACK( log )
	ADD_MATH_TPACK( log10 )

	ADD_MATH_TPACK( sin )
	ADD_MATH_TPACK( cos )
	ADD_MATH_TPACK( tan )

	ADD_MATH_TPACK( asin )
	ADD_MATH_TPACK( acos )
	ADD_MATH_TPACK( atan )
	ADD_MATH_TPACK( atan2 )

	ADD_MATH_TPACK( sinh )
	ADD_MATH_TPACK( cosh )
	ADD_MATH_TPACK( tanh )

#ifdef LMAT_HAS_CXX11_MATH

	ADD_MATH_TPACK( hypot )
	ADD_MATH_TPACK( cbrt )

	ADD_MATH_TPACK( exp2 )
	ADD_MATH_TPACK( log2 )
	ADD_MATH_TPACK( expm1 )
	ADD_MATH_TPACK( log1p )

	ADD_MATH_TPACK( asinh )
	ADD_MATH_TPACK( acosh )
	ADD_MATH_TPACK( atanh )

	ADD_MATH_TPACK( erf )
	ADD_MATH_TPACK( erfc )

	ADD_MATH_TPACK( lgamma )
	ADD_MATH_TPACK( tgamma )

#endif

END_MAIN_SUITE

