/**
 * @file test_ewise_accum.cpp
 *
 * @brief Unit testing of element-wise accumulation
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#define DEFAULT_M_VALUE 13;
#define DEFAULT_N_VALUE 9;

#include "multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/ewise_eval.h>

#include <light_mat/math/sse_ops.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_ops.h>
#endif

using namespace lmat;
using namespace lmat::test;

struct my_sum_kernel
{
	template<typename T>
	LMAT_ENSURE_INLINE
	void operator() (const T& x, T& s) const
	{
		s += x;
	}
};

struct my_max_kernel
{
	template<typename T>
	LMAT_ENSURE_INLINE
	void operator() (const T& x, T& s) const
	{
		s = math::max(s, x);
	}
};

struct my_min_kernel
{
	template<typename T>
	LMAT_ENSURE_INLINE
	void operator() (const T& x, T& s) const
	{
		s = math::min(s, x);
	}
};


template<typename U, int M, int N>
void test_linear_accum()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<cont, double, M, N>::cmat_t smat_t;

	mat_host<cont, double, M, N> src(m, n);
	src.fill_rand();

	smat_t smat = src.get_cmat();
	matrix_shape<M, N> shape(m, n);

	// sum

	double v_sum = 10.0;
	for (index_t i = 0; i < m * n; ++i) v_sum += smat[i];

	double r_sum = 10.0;
	linear_ewise(U(), shape).apply(my_sum_kernel(),
			in_(smat), in_out_(r_sum, atags::sum()));

	ASSERT_APPROX(r_sum, v_sum, 1.0e-12);

	// max

	double v_max = 0.5;
	for (index_t i = 0; i < m * n; ++i) if (v_max < smat[i]) v_max = smat[i];

	double r_max = 0.5;
	linear_ewise(U(), shape).apply(my_max_kernel(),
			in_(smat), in_out_(r_max, atags::max()));

	ASSERT_EQ( v_max, r_max );

	// min

	double v_min = 0.5;
	for (index_t i = 0; i < m * n; ++i) if (v_min > smat[i]) v_min = smat[i];

	double r_min = 0.5;
	linear_ewise(U(), shape).apply(my_min_kernel(),
			in_(smat), in_out_(r_min, atags::min()));

	ASSERT_EQ( v_min, r_min );

}


MN_CASE( ewise_accum, linear_scalar )
{
	test_linear_accum<atags::scalar, M, N>();
}

MN_CASE( ewise_accum, linear_sse )
{
	test_linear_accum<atags::simd<double, lmat::math::sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, linear_avx )
{
	test_linear_accum<atags::simd<double, lmat::math::avx_t>, M, N>();
}

#endif

BEGIN_TPACK( accum_linear_scalar )
	ADD_MN_CASE_3X3( ewise_accum, linear_scalar, DM, DN );
END_TPACK

BEGIN_TPACK( accum_linear_sse )
	ADD_MN_CASE_3X3( ewise_accum, linear_sse, DM, DN );
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( accum_linear_avx )
	ADD_MN_CASE_3X3( ewise_accum, linear_avx, DM, DN );
END_TPACK

#endif

BEGIN_MAIN_SUITE
	ADD_TPACK( accum_linear_scalar )
	ADD_TPACK( accum_linear_sse )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_linear_avx )
#endif
END_MAIN_SUITE




