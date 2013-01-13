/**
 * @file test_sse_reduce.cpp
 *
 * @brief Unit testing of SSE reduction
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/simd/sse_reduce.h>

using namespace lmat;
using namespace lmat::test;


T_CASE( sse_sum )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T s0 = T(0);

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		s0 += a_src[i];
	}

	pack_t a; a.load_u(a_src);

	ASSERT_EQ( sum(a), s0 );
}

T_CASE( sse_max )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T s0 = T(-1000);

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		if (a_src[i] > s0) s0 = a_src[i];
	}

	pack_t a; a.load_u(a_src);

	ASSERT_EQ( maximum(a), s0 );
}


T_CASE( sse_min )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T s0 = T(1000);

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(-3 - (int)i);
		if (a_src[i] < s0) s0 = a_src[i];
	}

	pack_t a; a.load_u(a_src);

	ASSERT_EQ( minimum(a), s0 );
}


SIMPLE_CASE( sse_booltest_f32 )
{
	const unsigned int M = 15;
	for (unsigned i = 0; i <= M; ++i)
	{
		bool b0 = i & 1;
		bool b1 = i & 2;
		bool b2 = i & 4;
		bool b3 = i & 8;

		sse_f32bpk pk(b0, b1, b2, b3);

		bool all_t = (i == M);
		bool all_f = (i == 0);
		bool any_t = !all_f;
		bool any_f = !all_t;

		ASSERT_EQ( all_true (pk), all_t );
		ASSERT_EQ( all_false(pk), all_f );
		ASSERT_EQ( any_true (pk), any_t );
		ASSERT_EQ( any_false(pk), any_f );

		__m128i vi = _mm_castps_si128(pk);
		ASSERT_EQ( internal::testz_sse2(vi), all_f);
		ASSERT_EQ( internal::testc_sse2(vi), all_t);
	}
}


SIMPLE_CASE( sse_booltest_f64 )
{
	const unsigned int M = 3;
	for (unsigned i = 0; i <= M; ++i)
	{
		bool b0 = i & 1;
		bool b1 = i & 2;

		sse_f64bpk pk(b0, b1);

		bool all_t = (i == M);
		bool all_f = (i == 0);
		bool any_t = !all_f;
		bool any_f = !all_t;

		ASSERT_EQ( all_true (pk), all_t );
		ASSERT_EQ( all_false(pk), all_f );
		ASSERT_EQ( any_true (pk), any_t );
		ASSERT_EQ( any_false(pk), any_f );

		__m128i vi = _mm_castpd_si128(pk);
		ASSERT_EQ( internal::testz_sse2(vi), all_f);
		ASSERT_EQ( internal::testc_sse2(vi), all_t);
	}
}

LTEST_INIT_AUTOSUITE

AUTO_TPACK( sse_stats )
{
	ADD_T_CASE_FP( sse_sum )
	ADD_T_CASE_FP( sse_max )
	ADD_T_CASE_FP( sse_min )
}

AUTO_TPACK( sse_booltest )
{
	ADD_SIMPLE_CASE( sse_booltest_f32 )
	ADD_SIMPLE_CASE( sse_booltest_f64 )
}


