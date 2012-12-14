/**
 * @file test_avx_ops.cpp
 *
 * Unit testing of basic operations on AVX packs
 * 
 * @author Dahua Lin 
 */

#include "simd_test_base.h"
#include <light_mat/math/avx_ops.h>
#include <light_mat/math/math_base.h>

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


T_CASE( avx_arith, add )
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


T_CASE( avx_arith, sub )
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


T_CASE( avx_arith, mul )
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

T_CASE( avx_arith, div )
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


T_CASE( avx_arith, neg )
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


T_CASE( avx_arith, abs )
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


T_CASE( avx_arith, max )
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


T_CASE( avx_arith, min )
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


T_CASE( avx_arith, sqr )
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


T_CASE( avx_arith, cube )
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

T_CASE( avx_arith, sqrt )
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


T_CASE( avx_arith, rcp )
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

T_CASE( avx_arith, rsqrt )
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



T_CASE( avx_comp, eq )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		int dv = (int)(i % 3) - 1;

		a_src[i] = T(i + 1);
		b_src[i] = a_src[i] + T(dv);

		r1[i] = -bint(a_src[i] == b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a == b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( avx_comp, ne )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		int dv = (int)(i % 3) - 1;

		a_src[i] = T(i + 1);
		b_src[i] = a_src[i] + T(dv);

		r1[i] = -bint(a_src[i] != b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a != b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( avx_comp, gt )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		int dv = (int)(i % 3) - 1;

		a_src[i] = T(i + 1);
		b_src[i] = a_src[i] + T(dv);

		r1[i] = -bint(a_src[i] > b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a > b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( avx_comp, ge )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		int dv = (int)(i % 3) - 1;

		a_src[i] = T(i + 1);
		b_src[i] = a_src[i] + T(dv);

		r1[i] = -bint(a_src[i] >= b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a >= b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( avx_comp, lt )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		int dv = (int)(i % 3) - 1;

		a_src[i] = T(i + 1);
		b_src[i] = a_src[i] + T(dv);

		r1[i] = -bint(a_src[i] < b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a < b);
	ASSERT_SIMD_EQ(r, r1);
}


T_CASE( avx_comp, le )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;

	T a_src[width];
	T b_src[width];

	bint r1[width];

	for (unsigned i = 0; i < width; ++i)
	{
		int dv = (int)(i % 3) - 1;

		a_src[i] = T(i + 1);
		b_src[i] = a_src[i] + T(dv);

		r1[i] = -bint(a_src[i] <= b_src[i]);
	}

	pack_t a; a.load_u(a_src);
	pack_t b; b.load_u(b_src);

	bpack_t r = (a <= b);
	ASSERT_SIMD_EQ(r, r1);
}


template<typename T> struct avx_not_tbody;

template<>
struct avx_not_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(true, false, false, true, false, true, true, false);
		bint r1[8] = {0, -1, -1, 0, -1, 0, 0, -1};

		ASSERT_SIMD_EQ(~a, r1);
	}
};

template<>
struct avx_not_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(true, false, false, true);
		bint r1[4] = {0, -1, -1, 0};

		ASSERT_SIMD_EQ(~a, r1);
	}
};

T_CASE( avx_logical, not )
{
	avx_not_tbody<T>::run();
}


template<typename T> struct avx_and_tbody;

template<>
struct avx_and_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true, true, true, false, false);
		bpack_t b(false, true, false, true, true, false, true, false);
		bint r[8] = {0, 0, 0, -1, -1, 0, 0, 0};

		ASSERT_SIMD_EQ(a & b, r);
	}
};

template<>
struct avx_and_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r[4] = {0, 0, 0, -1};

		ASSERT_SIMD_EQ(a & b, r);
	}
};

T_CASE( avx_logical, and )
{
	avx_and_tbody<T>::run();
}


template<typename T> struct avx_or_tbody;

template<>
struct avx_or_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true, true, true, false, false);
		bpack_t b(false, true, false, true, true, false, true, false);
		bint r[8] = {0, -1, -1, -1, -1, -1, -1, 0};

		ASSERT_SIMD_EQ(a | b, r);
	}
};

template<>
struct avx_or_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r[4] = {0, -1, -1, -1};

		ASSERT_SIMD_EQ(a | b, r);
	}
};

T_CASE( avx_logical, or )
{
	avx_or_tbody<T>::run();
}


template<typename T> struct avx_eq_tbody;

template<>
struct avx_eq_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true, true, true, false, false);
		bpack_t b(false, true, false, true, true, false, true, false);
		bint r[8] = {-1, 0, 0, -1, -1, 0, 0, -1};

		ASSERT_SIMD_EQ(a == b, r);
	}
};

template<>
struct avx_eq_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r[4] = {-1, 0, 0, -1};

		ASSERT_SIMD_EQ(a == b, r);
	}
};

T_CASE( avx_logical, eq )
{
	avx_eq_tbody<T>::run();
}


template<typename T> struct avx_ne_tbody;

template<>
struct avx_ne_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true, true, true, false, false);
		bpack_t b(false, true, false, true, true, false, true, false);
		bint r[8] = {0, -1, -1, 0, 0, -1, -1, 0};

		ASSERT_SIMD_EQ(a != b, r);
	}
};

template<>
struct avx_ne_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r[4] = {0, -1, -1, 0};

		ASSERT_SIMD_EQ(a != b, r);
	}
};

T_CASE( avx_logical, ne )
{
	avx_ne_tbody<T>::run();
}



template<typename T> struct avx_fpclassify_tbody;

template<>
struct avx_fpclassify_tbody<float>
{
	static void run()
	{
		typedef std::numeric_limits<float> lim_t;

		typedef simd_pack<float, avx_t> pack_t;
		typedef simd_bpack<float, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		float v_zero = 0.f;
		float v_neg_zero = -0.f;
		float v_one = 1.f;
		float v_neg_one = -1.f;
		float v_inf = lim_t::infinity();
		float v_neg_inf = -lim_t::infinity();
		float v_nan = lim_t::quiet_NaN();
		float v_neg_nan = -lim_t::quiet_NaN();

		pack_t a1(
				v_zero, v_neg_zero, v_one, v_neg_one,
				v_inf, v_neg_inf, v_nan, v_neg_nan);

		bint is_neg_r1   [8] = { 0, -1, 0, -1, 0, -1, 0, -1 };
		bint is_finite_r1[8] = { -1, -1, -1, -1, 0, 0, 0, 0 };
		bint is_inf_r1   [8] = { 0, 0, 0, 0, -1, -1, 0, 0 };
		bint is_nan_r1   [8] = { 0, 0, 0, 0, 0, 0, -1, -1 };

		ASSERT_SIMD_EQ( math::is_neg(a1), is_neg_r1  );
		ASSERT_SIMD_EQ( math::is_finite(a1), is_finite_r1  );
		ASSERT_SIMD_EQ( math::is_inf(a1), is_inf_r1  );
		ASSERT_SIMD_EQ( math::is_nan(a1), is_nan_r1  );
	}
};


template<>
struct avx_fpclassify_tbody<double>
{
	static void run()
	{
		typedef std::numeric_limits<double> lim_t;

		typedef simd_pack<double, avx_t> pack_t;
		typedef simd_bpack<double, avx_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		double v_zero = 0.f;
		double v_neg_zero = -0.f;
		double v_one = 1.f;
		double v_neg_one = -1.f;
		double v_inf = lim_t::infinity();
		double v_neg_inf = -lim_t::infinity();
		double v_nan = lim_t::quiet_NaN();
		double v_neg_nan = -lim_t::quiet_NaN();

		pack_t a1(v_zero, v_neg_zero, v_one, v_neg_one);
		pack_t a2(v_inf, v_neg_inf, v_nan, v_neg_nan);

		bint is_neg_r1[4] = { 0, -1, 0, -1 };
		bint is_neg_r2[4] = { 0, -1, 0, -1 };

		bint is_finite_r1[4] = { -1, -1, -1, -1 };
		bint is_finite_r2[4] = { 0, 0, 0, 0 };

		bint is_inf_r1[4] = { 0, 0, 0, 0 };
		bint is_inf_r2[4] = { -1, -1, 0, 0 };

		bint is_nan_r1[4] = { 0, 0, 0, 0 };
		bint is_nan_r2[4] = { 0, 0, -1, -1 };

		ASSERT_SIMD_EQ( math::is_neg(a1), is_neg_r1  );
		ASSERT_SIMD_EQ( math::is_neg(a2), is_neg_r2  );

		ASSERT_SIMD_EQ( math::is_finite(a1), is_finite_r1  );
		ASSERT_SIMD_EQ( math::is_finite(a2), is_finite_r2  );

		ASSERT_SIMD_EQ( math::is_inf(a1), is_inf_r1  );
		ASSERT_SIMD_EQ( math::is_inf(a2), is_inf_r2  );

		ASSERT_SIMD_EQ( math::is_nan(a1), is_nan_r1  );
		ASSERT_SIMD_EQ( math::is_nan(a2), is_nan_r2  );
	}
};

T_CASE( avx_fpclassify, all )
{
	avx_fpclassify_tbody<T>::run();
}


DEF_TPACK( avx_arith, add )
DEF_TPACK( avx_arith, sub )
DEF_TPACK( avx_arith, mul )
DEF_TPACK( avx_arith, div )
DEF_TPACK( avx_arith, neg )
DEF_TPACK( avx_arith, abs )
DEF_TPACK( avx_arith, max )
DEF_TPACK( avx_arith, min )
DEF_TPACK( avx_arith, sqr )
DEF_TPACK( avx_arith, cube )
DEF_TPACK( avx_arith, sqrt )
DEF_TPACK( avx_arith, rcp )
DEF_TPACK( avx_arith, rsqrt )

DEF_TPACK( avx_comp, eq )
DEF_TPACK( avx_comp, ne )
DEF_TPACK( avx_comp, gt )
DEF_TPACK( avx_comp, ge )
DEF_TPACK( avx_comp, lt )
DEF_TPACK( avx_comp, le )

DEF_TPACK( avx_logical, not )
DEF_TPACK( avx_logical, and )
DEF_TPACK( avx_logical, or )
DEF_TPACK( avx_logical, eq )
DEF_TPACK( avx_logical, ne )

DEF_TPACK( avx_fpclassify, all )


BEGIN_MAIN_SUITE
	ADD_TPACK( avx_arith_add )
	ADD_TPACK( avx_arith_sub )
	ADD_TPACK( avx_arith_mul )
	ADD_TPACK( avx_arith_div )
	ADD_TPACK( avx_arith_neg )
	ADD_TPACK( avx_arith_abs )
	ADD_TPACK( avx_arith_max )
	ADD_TPACK( avx_arith_min )
	ADD_TPACK( avx_arith_sqr )
	ADD_TPACK( avx_arith_cube )
	ADD_TPACK( avx_arith_sqrt )
	ADD_TPACK( avx_arith_rcp )
	ADD_TPACK( avx_arith_rsqrt )

	ADD_TPACK( avx_comp_eq )
	ADD_TPACK( avx_comp_ne )
	ADD_TPACK( avx_comp_gt )
	ADD_TPACK( avx_comp_ge )
	ADD_TPACK( avx_comp_lt )
	ADD_TPACK( avx_comp_le )

	ADD_TPACK( avx_logical_not )
	ADD_TPACK( avx_logical_and )
	ADD_TPACK( avx_logical_or )
	ADD_TPACK( avx_logical_eq )
	ADD_TPACK( avx_logical_ne )

	ADD_TPACK( avx_fpclassify_all )
END_MAIN_SUITE

