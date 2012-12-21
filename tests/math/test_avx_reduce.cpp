/**
 * @file test_avx_reduce.cpp
 *
 * Unit testing of AVX reduction
 * 
 * @author Dahua Lin 
 */


#include "simd_test_base.h"
#include <light_mat/math/avx_reduce.h>

using namespace lmat;
using namespace lmat::test;

using lmat::math::simd_pack;
using lmat::math::simd_bpack;
using lmat::math::avx_t;

#define DEF_TPACK( pname, tname ) \
	BEGIN_TPACK( pname##_##tname ) \
		ADD_T_CASE( pname, tname, float ) \
		ADD_T_CASE( pname, tname, double ) \
	END_TPACK


T_CASE( avx_reduce, sum )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_reduce, max )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_reduce, min )
{
	typedef simd_pack<T, avx_t> pack_t;
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


SIMPLE_CASE( avx_reduce, booltest_f32 )
{
	const unsigned int M = 255;
	for (unsigned i = 0; i <= M; ++i)
	{
		bool b0 = i & 1;
		bool b1 = i & 2;
		bool b2 = i & 4;
		bool b3 = i & 8;
		bool b4 = i & 16;
		bool b5 = i & 32;
		bool b6 = i & 64;
		bool b7 = i & 128;

		math::avx_f32bpk pk(b0, b1, b2, b3, b4, b5, b6, b7);

		bool all_t = (i == M);
		bool all_f = (i == 0);
		bool any_t = !all_f;
		bool any_f = !all_t;

		ASSERT_EQ( math::all_true (pk), all_t );
		ASSERT_EQ( math::all_false(pk), all_f );
		ASSERT_EQ( math::any_true (pk), any_t );
		ASSERT_EQ( math::any_false(pk), any_f );
	}
}

SIMPLE_CASE( avx_reduce, booltest_f64 )
{
	const unsigned int M = 15;
	for (unsigned i = 0; i <= M; ++i)
	{
		bool b0 = i & 1;
		bool b1 = i & 2;
		bool b2 = i & 4;
		bool b3 = i & 8;

		math::avx_f64bpk pk(b0, b1, b2, b3);

		bool all_t = (i == M);
		bool all_f = (i == 0);
		bool any_t = !all_f;
		bool any_f = !all_t;

		ASSERT_EQ( math::all_true (pk), all_t );
		ASSERT_EQ( math::all_false(pk), all_f );
		ASSERT_EQ( math::any_true (pk), any_t );
		ASSERT_EQ( math::any_false(pk), any_f );
	}
}



DEF_TPACK( avx_reduce, sum )
DEF_TPACK( avx_reduce, max )
DEF_TPACK( avx_reduce, min )

BEGIN_TPACK( avx_reduce_booltest )
	ADD_SIMPLE_CASE( avx_reduce, booltest_f32 )
	ADD_SIMPLE_CASE( avx_reduce, booltest_f64 )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( avx_reduce_sum )
	ADD_TPACK( avx_reduce_max )
	ADD_TPACK( avx_reduce_min )
	ADD_TPACK( avx_reduce_booltest )
END_MAIN_SUITE


