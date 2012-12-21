/**
 * @file test_sse_reduce.cpp
 *
 * @brief Unit testing of SSE reduction
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/math/sse_reduce.h>

using namespace lmat;
using namespace lmat::test;

using lmat::math::simd_pack;
using lmat::math::simd_bpack;
using lmat::math::sse_t;

#define DEF_TPACK( pname, tname ) \
	BEGIN_TPACK( pname##_##tname ) \
		ADD_T_CASE( pname, tname, float ) \
		ADD_T_CASE( pname, tname, double ) \
	END_TPACK


T_CASE( sse_reduce, sum )
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

	ASSERT_EQ( math::sum(a), s0 );
}

T_CASE( sse_reduce, max )
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

	ASSERT_EQ( math::maximum(a), s0 );
}


T_CASE( sse_reduce, min )
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

	ASSERT_EQ( math::minimum(a), s0 );
}


SIMPLE_CASE( sse_reduce, booltest_f32 )
{
	const unsigned int M = 15;
	for (unsigned i = 0; i <= M; ++i)
	{
		bool b0 = i & 1;
		bool b1 = i & 2;
		bool b2 = i & 4;
		bool b3 = i & 8;

		math::sse_f32bpk pk(b0, b1, b2, b3);

		bool all_t = (i == M);
		bool all_f = (i == 0);
		bool any_t = !all_f;
		bool any_f = !all_t;

		ASSERT_EQ( math::all_true (pk), all_t );
		ASSERT_EQ( math::all_false(pk), all_f );
		ASSERT_EQ( math::any_true (pk), any_t );
		ASSERT_EQ( math::any_false(pk), any_f );

		__m128i vi = _mm_castps_si128(pk);
		ASSERT_EQ( math::internal::testz_sse2(vi), all_f);
		ASSERT_EQ( math::internal::testc_sse2(vi), all_t);
	}
}


SIMPLE_CASE( sse_reduce, booltest_f64 )
{
	const unsigned int M = 3;
	for (unsigned i = 0; i <= M; ++i)
	{
		bool b0 = i & 1;
		bool b1 = i & 2;

		math::sse_f64bpk pk(b0, b1);

		bool all_t = (i == M);
		bool all_f = (i == 0);
		bool any_t = !all_f;
		bool any_f = !all_t;

		ASSERT_EQ( math::all_true (pk), all_t );
		ASSERT_EQ( math::all_false(pk), all_f );
		ASSERT_EQ( math::any_true (pk), any_t );
		ASSERT_EQ( math::any_false(pk), any_f );

		__m128i vi = _mm_castpd_si128(pk);
		ASSERT_EQ( math::internal::testz_sse2(vi), all_f);
		ASSERT_EQ( math::internal::testc_sse2(vi), all_t);
	}
}



DEF_TPACK( sse_reduce, sum )
DEF_TPACK( sse_reduce, max )
DEF_TPACK( sse_reduce, min )

BEGIN_TPACK( sse_reduce_booltest )
	ADD_SIMPLE_CASE( sse_reduce, booltest_f32 )
	ADD_SIMPLE_CASE( sse_reduce, booltest_f64 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( sse_reduce_sum )
	ADD_TPACK( sse_reduce_max )
	ADD_TPACK( sse_reduce_min )
	ADD_TPACK( sse_reduce_booltest )
END_MAIN_SUITE

