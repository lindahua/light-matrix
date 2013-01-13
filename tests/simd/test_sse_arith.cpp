/**
 * @file test_sse_arith.cpp
 *
 * Unit testing of arithmetics on SSE packs
 * 
 * @author Dahua Lin 
 */

#include "simd_test_base.h"
#include <light_mat/simd/sse_arith.h>
#include <light_mat/math/math_base.h>

using namespace lmat;
using namespace lmat::test;


T_CASE( sse_add )
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


T_CASE( sse_sub )
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


T_CASE( sse_mul )
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

T_CASE( sse_div )
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


T_CASE( sse_neg )
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


T_CASE( sse_fma )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	T a_src[width];
	T b_src[width];
	T c_src[width];

	T r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		a_src[i] = T(i + 1);
		b_src[i] = T(2 * i + 3);
		c_src[i] = T(5) - T(i);

		r1[i] = math::fma(a_src[i], b_src[i], c_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);
	pack_t c; c.load_u(c_src);

	pack_t r = math::fma(a, b, c);
	ASSERT_SIMD_EQ(r, r1);
}



T_CASE( sse_abs )
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


T_CASE( sse_sqr )
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


T_CASE( sse_cube )
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

T_CASE( sse_sqrt )
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


T_CASE( sse_rcp )
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

T_CASE( sse_rsqrt )
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


T_CASE( sse_max )
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


T_CASE( sse_min )
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

template<typename T> struct sse_cond_tbody;

template<>
struct sse_cond_tbody<float>
{
	static void run()
	{
		typedef simd_pack<float, sse_t> pack_t;
		typedef simd_bpack<float, sse_t> bpack_t;

		bpack_t b(true, false, false, true);
		pack_t x(1.0f, 2.0f, 3.0f, 4.0f);
		pack_t y(5.0f, 6.0f, 7.0f, 8.0f);

		float r0[4] = {1.0f, 6.0f, 7.0f, 4.0f};
		pack_t r1 = math::cond(b, x, y);
		pack_t r2 = internal::cond_sse2(b, x, y);

		ASSERT_SIMD_EQ( r1, r0 );
		ASSERT_SIMD_EQ( r2, r0 );
	}
};

template<>
struct sse_cond_tbody<double>
{
	static void run()
	{
		typedef simd_pack<double, sse_t> pack_t;
		typedef simd_bpack<double, sse_t> bpack_t;

		bpack_t b(true, false);
		pack_t x(1.0, 2.0);
		pack_t y(5.0, 6.0);

		double r0[2] = {1.0, 6.0};
		pack_t r1 = math::cond(b, x, y);
		pack_t r2 = internal::cond_sse2(b, x, y);

		ASSERT_SIMD_EQ( r1, r0 );
		ASSERT_SIMD_EQ( r2, r0 );
	}
};


T_CASE( sse_cond )
{
	sse_cond_tbody<T>::run();
}

LTEST_INIT_AUTOSUITE

AUTO_TPACK( sse_arith )
{
	ADD_T_CASE_FP( sse_add )
	ADD_T_CASE_FP( sse_sub )
	ADD_T_CASE_FP( sse_mul )
	ADD_T_CASE_FP( sse_div )
	ADD_T_CASE_FP( sse_neg )
	ADD_T_CASE_FP( sse_fma )
}

AUTO_TPACK( sse_spower )
{
	ADD_T_CASE_FP( sse_abs )
	ADD_T_CASE_FP( sse_sqr )
	ADD_T_CASE_FP( sse_cube )

	ADD_T_CASE_FP( sse_rcp )
	ADD_T_CASE_FP( sse_sqrt )
	ADD_T_CASE_FP( sse_rsqrt )
}

AUTO_TPACK( sse_minmax )
{
	ADD_T_CASE_FP( sse_max )
	ADD_T_CASE_FP( sse_min )
}

AUTO_TPACK( sse_cond )
{
	ADD_T_CASE_FP( sse_cond )
}



