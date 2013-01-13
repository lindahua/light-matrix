/**
 * @file test_avx_pred.cpp
 *
 * Unit testing of predicates on AVX packs
 * 
 * @author Dahua Lin 
 */

#include "simd_test_base.h"
#include <light_mat/simd/avx_pred.h>
#include <light_mat/math/math_base.h>

using namespace lmat;
using namespace lmat::test;


T_CASE( avx_eq )
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


T_CASE( avx_ne )
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


T_CASE( avx_gt )
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


T_CASE( avx_ge )
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


T_CASE( avx_lt )
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


T_CASE( avx_le )
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

T_CASE( avx_logical_not )
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

T_CASE( avx_logical_and )
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

T_CASE( avx_logical_or )
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

T_CASE( avx_logical_eq )
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

T_CASE( avx_logical_ne )
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

		ASSERT_SIMD_EQ( math::signbit(a1), is_neg_r1  );
		ASSERT_SIMD_EQ( math::isfinite(a1), is_finite_r1  );
		ASSERT_SIMD_EQ( math::isinf(a1), is_inf_r1  );
		ASSERT_SIMD_EQ( math::isnan(a1), is_nan_r1  );
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

T_CASE( avx_fpclassify )
{
	avx_fpclassify_tbody<T>::run();
}


LTEST_INIT_AUTOSUITE

AUTO_TPACK( avx_comp )
{
	ADD_T_CASE_FP( avx_eq )
	ADD_T_CASE_FP( avx_ne )
	ADD_T_CASE_FP( avx_gt )
	ADD_T_CASE_FP( avx_ge )
	ADD_T_CASE_FP( avx_lt )
	ADD_T_CASE_FP( avx_le )
}

AUTO_TPACK( avx_logical )
{
	ADD_T_CASE_FP( avx_logical_not )
	ADD_T_CASE_FP( avx_logical_and )
	ADD_T_CASE_FP( avx_logical_or )
	ADD_T_CASE_FP( avx_logical_eq )
	ADD_T_CASE_FP( avx_logical_ne )
}

AUTO_TPACK( avx_fpclassify )
{
	ADD_T_CASE_FP( avx_fpclassify )
}



