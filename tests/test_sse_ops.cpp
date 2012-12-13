/**
 * @file test_sse_ops.cpp
 *
 * @brief Unit testing of sse_ops.h
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/math/sse_ops.h>
#include <light_mat/math/math_base.h>

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



T_CASE( sse_arith, add )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T b_src[width];

	T r1[width];
	T r2[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = T(2 * i + 3);

		r1[i] = a_src[i] + b_src[i];
		r2[i] = r1[i] + b_src[i];
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	pack_t r = a + b;
	ASSERT_SIMD_EQ(r, r1);

	r += b;
	ASSERT_SIMD_EQ(r, r2);
}


T_CASE( sse_arith, sub )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T b_src[width];

	T r1[width];
	T r2[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = T(2 * i + 3);

		r1[i] = a_src[i] - b_src[i];
		r2[i] = r1[i] - b_src[i];
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	pack_t r = a - b;
	ASSERT_SIMD_EQ(r, r1);

	r -= b;
	ASSERT_SIMD_EQ(r, r2);
}


T_CASE( sse_arith, mul )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T b_src[width];

	T r1[width];
	T r2[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = T(2 * i + 3);

		r1[i] = a_src[i] * b_src[i];
		r2[i] = r1[i] * b_src[i];
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	pack_t r = a * b;
	ASSERT_SIMD_EQ(r, r1);

	r *= b;
	ASSERT_SIMD_EQ(r, r2);
}

T_CASE( sse_arith, div )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T b_src[width];

	T r1[width];
	T r2[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = T(2 * i + 3);

		r1[i] = a_src[i] / b_src[i];
		r2[i] = r1[i] / b_src[i];
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	pack_t r = a / b;
	ASSERT_SIMD_ULP(r, r1, 1);

	r /= b;
	ASSERT_SIMD_ULP(r, r2, 3);
}


T_CASE( sse_arith, neg )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		if (i % 2 == 0) a_src[i] = - a_src[i];

		r1[i] = - a_src[i];
	}

	pack_t a; a.load_u(a_src);

	pack_t r = -a;
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_arith, abs )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		if (i % 2 == 0) a_src[i] = - a_src[i];

		r1[i] = math::abs(a_src[i]);
	}

	pack_t a; a.load_u(a_src);

	pack_t r = math::abs(a);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_arith, max )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T b_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 2);
		b_src[i] = T(3 * i + 1);

		r1[i] = math::max(a_src[i], b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	pack_t r = math::max(a, b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_arith, min )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T b_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 2);
		b_src[i] = T(3 * i + 1);

		r1[i] = math::min(a_src[i], b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	pack_t r = math::min(a, b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_arith, sqr )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 2);
		if (i % 2 == 0) a_src[i] = - a_src[i];

		r1[i] = math::sqr(a_src[i]);
	}

	pack_t a; a.load_u(a_src);

	pack_t r = math::sqr(a);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_arith, cube )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 2);
		if (i % 2 == 0) a_src[i] = - a_src[i];

		r1[i] = math::cube(a_src[i]);
	}

	pack_t a; a.load_u(a_src);

	pack_t r = math::cube(a);
	ASSERT_SIMD_EQ(r, r1);
}

T_CASE( sse_arith, sqrt )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 2);
		r1[i] = math::sqrt(a_src[i]);
	}

	pack_t a; a.load_u(a_src);

	pack_t r = math::sqrt(a);
	ASSERT_SIMD_ULP(r, r1, 1);
}


T_CASE( sse_arith, rcp )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 2);
		if (i % 2 == 0) a_src[i] = - a_src[i];

		r1[i] = math::rcp(a_src[i]);
	}

	pack_t a; a.load_u(a_src);

	pack_t r = math::rcp(a);
	ASSERT_SIMD_ULP(r, r1, 1);
}

T_CASE( sse_arith, rsqrt )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 2);
		r1[i] = math::rsqrt(a_src[i]);
	}

	pack_t a; a.load_u(a_src);

	pack_t r = math::rsqrt(a);
	ASSERT_SIMD_ULP(r, r1, 1);
}


T_CASE( sse_comp, eq )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, sse_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = i % 2 == 0 ? a_src[i] : a_src[i] + T(1);

		r1[i] = -bint(a_src[i] == b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a == b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_comp, ne )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, sse_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = i % 2 == 0 ? a_src[i] : a_src[i] + T(1);

		r1[i] = -bint(a_src[i] != b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a != b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_comp, gt )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, sse_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = i % 2 == 0 ? a_src[i] - T(1) : a_src[i] + T(1);

		r1[i] = -bint(a_src[i] > b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a > b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_comp, ge )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, sse_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = i % 2 == 0 ? a_src[i] : a_src[i] + T(1);

		r1[i] = -bint(a_src[i] >= b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a >= b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_comp, lt )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, sse_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = i % 2 == 0 ? a_src[i] - T(1) : a_src[i] + T(1);

		r1[i] = -bint(a_src[i] < b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a < b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( sse_comp, le )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, sse_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = i % 2 == 0 ? a_src[i] : a_src[i] + T(1);

		r1[i] = -bint(a_src[i] <= b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a <= b);
	ASSERT_SIMD_EQ(r, r1);
}


DEF_TPACK( sse_arith, add )
DEF_TPACK( sse_arith, sub )
DEF_TPACK( sse_arith, mul )
DEF_TPACK( sse_arith, div )
DEF_TPACK( sse_arith, neg )
DEF_TPACK( sse_arith, abs )
DEF_TPACK( sse_arith, max )
DEF_TPACK( sse_arith, min )
DEF_TPACK( sse_arith, sqr )
DEF_TPACK( sse_arith, cube )
DEF_TPACK( sse_arith, sqrt )
DEF_TPACK( sse_arith, rcp )
DEF_TPACK( sse_arith, rsqrt )

DEF_TPACK( sse_comp, eq )
DEF_TPACK( sse_comp, ne )
DEF_TPACK( sse_comp, gt )
DEF_TPACK( sse_comp, ge )
DEF_TPACK( sse_comp, lt )
DEF_TPACK( sse_comp, le )

BEGIN_MAIN_SUITE
	ADD_TPACK( sse_arith_add )
	ADD_TPACK( sse_arith_sub )
	ADD_TPACK( sse_arith_mul )
	ADD_TPACK( sse_arith_div )
	ADD_TPACK( sse_arith_neg )
	ADD_TPACK( sse_arith_abs )
	ADD_TPACK( sse_arith_max )
	ADD_TPACK( sse_arith_min )
	ADD_TPACK( sse_arith_sqr )
	ADD_TPACK( sse_arith_cube )
	ADD_TPACK( sse_arith_sqrt )
	ADD_TPACK( sse_arith_rcp )
	ADD_TPACK( sse_arith_rsqrt )

	ADD_TPACK( sse_comp_eq )
	ADD_TPACK( sse_comp_ne )
	ADD_TPACK( sse_comp_gt )
	ADD_TPACK( sse_comp_ge )
	ADD_TPACK( sse_comp_lt )
	ADD_TPACK( sse_comp_le )
END_MAIN_SUITE



