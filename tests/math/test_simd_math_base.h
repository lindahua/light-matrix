/**
 * @file test_simd_math_base.h
 *
 * @brief The basis for testing SIMD math
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef TEST_SIMD_MATH_BASE_H_
#define TEST_SIMD_MATH_BASE_H_

#include "../simd/simd_test_base.h"


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

#define DEFINE_MATH_CASE1( Name, SIMD, ulp, LB, UB ) \
	T_CASE( simd_math, Name##_##SIMD ) { \
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
			ASSERT_SIMD_ULP(r, r0, ulp); \
		} \
	}

#define DEFINE_MATH_CASE2( Name, SIMD, ulp, LBa, UBa, LBb, UBb ) \
	T_CASE( simd_math, Name##_##SIMD ) { \
		typedef simd_pack<T, SIMD##_t> pack_t; \
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
			ASSERT_SIMD_ULP(r, r0, ulp); \
		} \
	}


#ifdef LMAT_HAS_AVX

#define DEFINE_MATH_TPACK1( Name, ulp, LB, UB ) \
	DEFINE_MATH_CASE1( Name,  sse, ulp, LB, UB ) \
	DEFINE_MATH_CASE1( Name,  avx, ulp, LB, UB ) \
	BEGIN_TPACK( simd_##Name ) \
		ADD_T_CASE( simd_math, Name##_sse, float ) \
		ADD_T_CASE( simd_math, Name##_sse, double ) \
		ADD_T_CASE( simd_math, Name##_avx, float ) \
		ADD_T_CASE( simd_math, Name##_avx, double ) \
	END_TPACK

#define DEFINE_MATH_TPACK2( Name, ulp, LBa, UBa, LBb, UBb ) \
	DEFINE_MATH_CASE2( Name,  sse, ulp, LBa, UBa, LBb, UBb ) \
	DEFINE_MATH_CASE2( Name,  avx, ulp, LBa, UBa, LBb, UBb ) \
	BEGIN_TPACK( simd_##Name ) \
		ADD_T_CASE( simd_math, Name##_sse, float ) \
		ADD_T_CASE( simd_math, Name##_sse, double ) \
		ADD_T_CASE( simd_math, Name##_avx, float ) \
		ADD_T_CASE( simd_math, Name##_avx, double ) \
	END_TPACK

#else

#define DEFINE_MATH_TPACK1( Name, ulp, LB, UB ) \
	DEFINE_MATH_CASE1( Name,  sse, ulp, LB, UB ) \
	BEGIN_TPACK( simd_##Name ) \
		ADD_T_CASE( simd_math, Name##_sse, float ) \
		ADD_T_CASE( simd_math, Name##_sse, double ) \
	END_TPACK

#define DEFINE_MATH_TPACK2( Name, ulp, LBa, UBa, LBb, UBb ) \
	DEFINE_MATH_CASE2( Name,  sse, ulp, LBa, UBa, LBb, UBb ) \
	BEGIN_TPACK( simd_##Name ) \
		ADD_T_CASE( simd_math, Name##_sse, float ) \
		ADD_T_CASE( simd_math, Name##_sse, double ) \
	END_TPACK

#endif


#define ADD_MATH_TPACK( Name ) ADD_TPACK( simd_##Name )


#endif
