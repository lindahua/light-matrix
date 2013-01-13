/**
 * @file test_avx_arith.cpp
 *
 * Test of arithmetics on AVX packs
 * 
 * @author Dahua Lin 
 */

#include "simd_test_base.h"
#include <light_mat/simd/avx_arith.h>
#include <light_mat/math/math_base.h>

using namespace lmat;
using namespace lmat::test;


T_CASE( avx_add )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_sub )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_mul )
{
	typedef simd_pack<T, avx_t> pack_t;
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

T_CASE( avx_div )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_neg )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_fma )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_abs )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_sqr )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_cube )
{
	typedef simd_pack<T, avx_t> pack_t;
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

T_CASE( avx_sqrt )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_rcp )
{
	typedef simd_pack<T, avx_t> pack_t;
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

T_CASE( avx_rsqrt )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_max )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_min )
{
	typedef simd_pack<T, avx_t> pack_t;
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


template<typename T> struct avx_cond_tbody;

template<> struct avx_cond_tbody<float>
{
	static void run()
	{
		typedef simd_pack<float, avx_t> pack_t;
		typedef simd_bpack<float, avx_t> bpack_t;

		bpack_t b(true, false, true, false, true, true, false, false);
		pack_t x( 1.f,  2.f,  3.f,  4.f,  5.f,  6.f,  7.f,  8.f);
		pack_t y(-1.f, -2.f, -3.f, -4.f, -5.f, -6.f, -7.f, -8.f);

		float r0[8] = {1.f, -2.f, 3.f, -4.f, 5.f, 6.f, -7.f, -8.f};

		ASSERT_SIMD_EQ( math::cond(b, x, y), r0 );
	}
};

template<> struct avx_cond_tbody<double>
{
	static void run()
	{
		typedef simd_pack<double, avx_t> pack_t;
		typedef simd_bpack<double, avx_t> bpack_t;

		bpack_t b(true, false, true, false);
		pack_t x( 1.0,  2.0,  3.0,  4.0);
		pack_t y(-1.0, -2.0, -3.0, -4.0);

		double r0[8] = {1.0, -2.0, 3.0, -4.0};
		ASSERT_SIMD_EQ( math::cond(b, x, y), r0 );
	}
};


T_CASE( avx_cond )
{
	avx_cond_tbody<T>::run();
}

LTEST_INIT_AUTOSUITE

AUTO_TPACK( avx_arith )
{
	ADD_T_CASE_FP( avx_add )
	ADD_T_CASE_FP( avx_sub )
	ADD_T_CASE_FP( avx_mul )
	ADD_T_CASE_FP( avx_div )
	ADD_T_CASE_FP( avx_neg )
	ADD_T_CASE_FP( avx_fma )
}

AUTO_TPACK( avx_spower )
{
	ADD_T_CASE_FP( avx_abs )
	ADD_T_CASE_FP( avx_sqr )
	ADD_T_CASE_FP( avx_cube )

	ADD_T_CASE_FP( avx_rcp )
	ADD_T_CASE_FP( avx_sqrt )
	ADD_T_CASE_FP( avx_rsqrt )
}

AUTO_TPACK( avx_minmax )
{
	ADD_T_CASE_FP( avx_max )
	ADD_T_CASE_FP( avx_min )
}

AUTO_TPACK( avx_cond )
{
	ADD_T_CASE_FP( avx_cond )
}





