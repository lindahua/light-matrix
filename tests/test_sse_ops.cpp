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



template<typename T> struct sse_not_tbody;

template<>
struct sse_not_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(true, false, false, true);
		bint r1[4] = {0, -1, -1, 0};

		ASSERT_SIMD_EQ(~a, r1);
	}
};

template<>
struct sse_not_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(true, false);
		bint r1[4] = {0, -1};

		ASSERT_SIMD_EQ(~a, r1);
	}
};

T_CASE( sse_logical, not )
{
	sse_not_tbody<T>::run();
}


template<typename T> struct sse_and_tbody;

template<>
struct sse_and_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r1[4] = {0, 0, 0, -1};

		ASSERT_SIMD_EQ(a & b, r1);
	}
};

template<>
struct sse_and_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a1(false, false);
		bpack_t b1(false, true);
		bint r1[4] = {0, 0};

		bpack_t a2(true, true);
		bpack_t b2(false, true);
		bint r2[4] = {0, -1};

		ASSERT_SIMD_EQ(a1 & b1, r1);
		ASSERT_SIMD_EQ(a2 & b2, r2);
	}
};

T_CASE( sse_logical, and )
{
	sse_and_tbody<T>::run();
}


template<typename T> struct sse_or_tbody;

template<>
struct sse_or_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r1[4] = {0, -1, -1, -1};

		ASSERT_SIMD_EQ(a | b, r1);
	}
};

template<>
struct sse_or_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a1(false, false);
		bpack_t b1(false, true);
		bint r1[4] = {0, -1};

		bpack_t a2(true, true);
		bpack_t b2(false, true);
		bint r2[4] = {-1, -1};

		ASSERT_SIMD_EQ(a1 | b1, r1);
		ASSERT_SIMD_EQ(a2 | b2, r2);
	}
};

T_CASE( sse_logical, or )
{
	sse_or_tbody<T>::run();
}


template<typename T> struct sse_eq_tbody;

template<>
struct sse_eq_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r1[4] = {-1, 0, 0, -1};

		ASSERT_SIMD_EQ(a == b, r1);
	}
};

template<>
struct sse_eq_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a1(false, false);
		bpack_t b1(false, true);
		bint r1[4] = {-1, 0};

		bpack_t a2(true, true);
		bpack_t b2(false, true);
		bint r2[4] = {0, -1};

		ASSERT_SIMD_EQ(a1 == b1, r1);
		ASSERT_SIMD_EQ(a2 == b2, r2);
	}
};

T_CASE( sse_logical, eq )
{
	sse_eq_tbody<T>::run();
}


template<typename T> struct sse_ne_tbody;

template<>
struct sse_ne_tbody<float>
{
	static void run()
	{
		typedef simd_bpack<float, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a(false, false, true, true);
		bpack_t b(false, true, false, true);
		bint r1[4] = {0, -1, -1, 0};

		ASSERT_SIMD_EQ(a != b, r1);
	}
};

template<>
struct sse_ne_tbody<double>
{
	static void run()
	{
		typedef simd_bpack<double, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		bpack_t a1(false, false);
		bpack_t b1(false, true);
		bint r1[4] = {0, -1};

		bpack_t a2(true, true);
		bpack_t b2(false, true);
		bint r2[4] = {-1, 0};

		ASSERT_SIMD_EQ(a1 != b1, r1);
		ASSERT_SIMD_EQ(a2 != b2, r2);
	}
};

T_CASE( sse_logical, ne )
{
	sse_ne_tbody<T>::run();
}



template<typename T> struct sse_fpclassify_tbody;

template<>
struct sse_fpclassify_tbody<float>
{
	static void run()
	{
		typedef std::numeric_limits<float> lim_t;

		typedef simd_pack<float, sse_t> pack_t;
		typedef simd_bpack<float, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		float v_zero = 0.f;
		float v_neg_zero = -0.f;
		float v_one = 1.f;
		float v_neg_one = -1.f;
		float v_inf = lim_t::infinity();
		float v_neg_inf = -lim_t::infinity();
		float v_nan = lim_t::quiet_NaN();
		float v_neg_nan = -lim_t::quiet_NaN();

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


template<>
struct sse_fpclassify_tbody<double>
{
	static void run()
	{
		typedef std::numeric_limits<double> lim_t;

		typedef simd_pack<double, sse_t> pack_t;
		typedef simd_bpack<double, sse_t> bpack_t;
		typedef typename bpack_t::bint_type bint;

		double v_zero = 0.0;
		double v_neg_zero = -0.0;
		double v_one = 1.0;
		double v_neg_one = -1.0;
		double v_inf = lim_t::infinity();
		double v_neg_inf = -lim_t::infinity();
		double v_nan = lim_t::quiet_NaN();
		double v_neg_nan = -lim_t::quiet_NaN();

		pack_t a1(v_zero, v_neg_zero);
		pack_t a2(v_one, v_neg_one);
		pack_t a3(v_inf, v_neg_inf);
		pack_t a4(v_nan, v_neg_nan);

		bint is_neg_r1[2] = { 0, -1 };
		bint is_neg_r2[2] = { 0, -1 };
		bint is_neg_r3[2] = { 0, -1 };
		bint is_neg_r4[2] = { 0, -1 };

		bint is_finite_r1[2] = { -1, -1 };
		bint is_finite_r2[2] = { -1, -1 };
		bint is_finite_r3[2] = { 0, 0 };
		bint is_finite_r4[2] = { 0, 0 };

		bint is_inf_r1[2] = { 0, 0 };
		bint is_inf_r2[2] = { 0, 0 };
		bint is_inf_r3[2] = { -1, -1 };
		bint is_inf_r4[2] = { 0, 0 };

		bint is_nan_r1[2] = { 0, 0 };
		bint is_nan_r2[2] = { 0, 0 };
		bint is_nan_r3[2] = { 0, 0 };
		bint is_nan_r4[2] = { -1, -1 };

		ASSERT_SIMD_EQ( math::is_neg(a1), is_neg_r1  );
		ASSERT_SIMD_EQ( math::is_neg(a2), is_neg_r2  );
		ASSERT_SIMD_EQ( math::is_neg(a3), is_neg_r3  );
		ASSERT_SIMD_EQ( math::is_neg(a4), is_neg_r4  );

		ASSERT_SIMD_EQ( math::is_finite(a1), is_finite_r1  );
		ASSERT_SIMD_EQ( math::is_finite(a2), is_finite_r2  );

		ASSERT_SIMD_EQ( math::is_finite(a3), is_finite_r3  );
		ASSERT_SIMD_EQ( math::is_finite(a4), is_finite_r4  );

		ASSERT_SIMD_EQ( math::is_inf(a1), is_inf_r1  );
		ASSERT_SIMD_EQ( math::is_inf(a2), is_inf_r2  );
		ASSERT_SIMD_EQ( math::is_inf(a3), is_inf_r3 );
		ASSERT_SIMD_EQ( math::is_inf(a4), is_inf_r4  );

		ASSERT_SIMD_EQ( math::is_nan(a1), is_nan_r1  );
		ASSERT_SIMD_EQ( math::is_nan(a2), is_nan_r2  );
		ASSERT_SIMD_EQ( math::is_nan(a3), is_nan_r3  );
		ASSERT_SIMD_EQ( math::is_nan(a4), is_nan_r4  );
	}
};

T_CASE( sse_fpclassify, all )
{
	sse_fpclassify_tbody<T>::run();
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

DEF_TPACK( sse_logical, not )
DEF_TPACK( sse_logical, and )
DEF_TPACK( sse_logical, or )
DEF_TPACK( sse_logical, eq )
DEF_TPACK( sse_logical, ne )

DEF_TPACK( sse_fpclassify, all )

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

	ADD_TPACK( sse_logical_not )
	ADD_TPACK( sse_logical_and )
	ADD_TPACK( sse_logical_or )
	ADD_TPACK( sse_logical_eq )
	ADD_TPACK( sse_logical_ne )

	ADD_TPACK( sse_fpclassify_all )
END_MAIN_SUITE

