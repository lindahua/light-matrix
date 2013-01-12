/**
 * @file test_sse_pred.cpp
 *
 * Unit testing of predicates on SSE
 * 
 * @author Dahua Lin 
 */


#include "simd_test_base.h"
#include <light_mat/simd/sse_pred.h>
#include <light_mat/math/math_base.h>

using namespace lmat;
using namespace lmat::test;


T_CASE( sse_eq )
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


T_CASE( sse_ne )
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


T_CASE( sse_gt )
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


T_CASE( sse_ge )
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


T_CASE( sse_lt )
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


T_CASE( sse_le )
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
		bint r1[2] = {0, -1};

		ASSERT_SIMD_EQ(~a, r1);
	}
};

T_CASE( sse_logical_not )
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
		bint r1[2] = {0, 0};

		bpack_t a2(true, true);
		bpack_t b2(false, true);
		bint r2[2] = {0, -1};

		ASSERT_SIMD_EQ(a1 & b1, r1);
		ASSERT_SIMD_EQ(a2 & b2, r2);
	}
};

T_CASE( sse_logical_and )
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
		bint r1[2] = {0, -1};

		bpack_t a2(true, true);
		bpack_t b2(false, true);
		bint r2[2] = {-1, -1};

		ASSERT_SIMD_EQ(a1 | b1, r1);
		ASSERT_SIMD_EQ(a2 | b2, r2);
	}
};

T_CASE( sse_logical_or )
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

T_CASE( sse_logical_eq )
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

T_CASE( sse_logical_ne )
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

		ASSERT_SIMD_EQ( math::signbit(a1), is_neg_r1  );
		ASSERT_SIMD_EQ( math::signbit(a2), is_neg_r2  );

		ASSERT_SIMD_EQ( math::isfinite(a1), is_finite_r1  );
		ASSERT_SIMD_EQ( math::isfinite(a2), is_finite_r2  );

		ASSERT_SIMD_EQ( math::isinf(a1), is_inf_r1  );
		ASSERT_SIMD_EQ( math::isinf(a2), is_inf_r2  );

		ASSERT_SIMD_EQ( math::isnan(a1), is_nan_r1  );
		ASSERT_SIMD_EQ( math::isnan(a2), is_nan_r2  );
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

		ASSERT_SIMD_EQ( math::signbit(a1), is_neg_r1  );
		ASSERT_SIMD_EQ( math::signbit(a2), is_neg_r2  );
		ASSERT_SIMD_EQ( math::signbit(a3), is_neg_r3  );
		ASSERT_SIMD_EQ( math::signbit(a4), is_neg_r4  );

		ASSERT_SIMD_EQ( math::isfinite(a1), is_finite_r1  );
		ASSERT_SIMD_EQ( math::isfinite(a2), is_finite_r2  );

		ASSERT_SIMD_EQ( math::isfinite(a3), is_finite_r3  );
		ASSERT_SIMD_EQ( math::isfinite(a4), is_finite_r4  );

		ASSERT_SIMD_EQ( math::isinf(a1), is_inf_r1  );
		ASSERT_SIMD_EQ( math::isinf(a2), is_inf_r2  );
		ASSERT_SIMD_EQ( math::isinf(a3), is_inf_r3 );
		ASSERT_SIMD_EQ( math::isinf(a4), is_inf_r4  );

		ASSERT_SIMD_EQ( math::isnan(a1), is_nan_r1  );
		ASSERT_SIMD_EQ( math::isnan(a2), is_nan_r2  );
		ASSERT_SIMD_EQ( math::isnan(a3), is_nan_r3  );
		ASSERT_SIMD_EQ( math::isnan(a4), is_nan_r4  );
	}
};

T_CASE( sse_fpclassify )
{
	sse_fpclassify_tbody<T>::run();
}

AUTO_TPACK( sse_comp )
{
	ADD_T_CASE_FP( sse_eq )
	ADD_T_CASE_FP( sse_ne )
	ADD_T_CASE_FP( sse_gt )
	ADD_T_CASE_FP( sse_ge )
	ADD_T_CASE_FP( sse_lt )
	ADD_T_CASE_FP( sse_le )
}

AUTO_TPACK( sse_logical )
{
	ADD_T_CASE_FP( sse_logical_not )
	ADD_T_CASE_FP( sse_logical_and )
	ADD_T_CASE_FP( sse_logical_or )
	ADD_T_CASE_FP( sse_logical_eq )
	ADD_T_CASE_FP( sse_logical_ne )
}

AUTO_TPACK( sse_fpclassify )
{
	ADD_T_CASE_FP( sse_fpclassify )
}


