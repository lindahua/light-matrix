/**
 * @file test_ewise_accum.cpp
 *
 * @brief Unit testing of element-wise accumulation
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#define DEFAULT_M_VALUE 13
#define DEFAULT_N_VALUE 9

#include "../multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/ewise_eval.h>

#include <light_mat/math/sse_ops.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_ops.h>
#endif

using namespace lmat;
using namespace lmat::test;

template<typename T>
struct my_sum_kernel
{
	typedef T value_type;

	LMAT_ENSURE_INLINE
	void operator() (const T& x, T& s) const
	{
		s += x;
	}

	template<typename Kind>
	LMAT_ENSURE_INLINE
	void operator() (const math::simd_pack<T, Kind>& x, math::simd_pack<T, Kind>& s) const
	{
		s += x;
	}
};


template<typename T>
struct my_max_kernel
{
	typedef T value_type;

	LMAT_ENSURE_INLINE
	void operator() (const T& x, T& s) const
	{
		s = math::max(s, x);
	}

	template<typename Kind>
	LMAT_ENSURE_INLINE
	void operator() (const math::simd_pack<T, Kind>& x, math::simd_pack<T, Kind>& s) const
	{
		s = math::max(s, x);
	}
};


template<typename T>
struct my_min_kernel
{
	typedef T value_type;

	LMAT_ENSURE_INLINE
	void operator() (const T& x, T& s) const
	{
		s = math::min(s, x);
	}

	template<typename Kind>
	LMAT_ENSURE_INLINE
	void operator() (const math::simd_pack<T, Kind>& x, math::simd_pack<T, Kind>& s) const
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
	ewise(my_sum_kernel<double>(), U())(shape, in_(smat), in_out_(r_sum, atags::sum()));

	ASSERT_APPROX(r_sum, v_sum, 1.0e-12);

	// max

	double v_max = 0.5;
	for (index_t i = 0; i < m * n; ++i) if (v_max < smat[i]) v_max = smat[i];

	double r_max = 0.5;
	ewise(my_max_kernel<double>(), U())(shape, in_(smat), in_out_(r_max, atags::max()));

	ASSERT_EQ( v_max, r_max );

	// min

	double v_min = 0.5;
	for (index_t i = 0; i < m * n; ++i) if (v_min > smat[i]) v_min = smat[i];

	double r_min = 0.5;
	ewise(my_min_kernel<double>(), U())(shape, in_(smat), in_out_(r_min, atags::min()));

	ASSERT_EQ( v_min, r_min );
}


template<typename U, int M, int N>
void test_percol_accum()
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
	percol(ewise(my_sum_kernel<double>(), U()), shape,
			in_(smat), in_out_(r_sum, atags::sum()));

	ASSERT_APPROX(r_sum, v_sum, 1.0e-12);

	// max

	double v_max = 0.5;
	for (index_t i = 0; i < m * n; ++i) if (v_max < smat[i]) v_max = smat[i];

	double r_max = 0.5;
	percol(ewise(my_max_kernel<double>(), U()), shape,
			in_(smat), in_out_(r_max, atags::max()));

	ASSERT_EQ( v_max, r_max );

	// min

	double v_min = 0.5;
	for (index_t i = 0; i < m * n; ++i) if (v_min > smat[i]) v_min = smat[i];

	double r_min = 0.5;
	percol(ewise(my_min_kernel<double>(), U()), shape,
			in_(smat), in_out_(r_min, atags::min()));

	ASSERT_EQ( v_min, r_min );
}


template<typename U, typename DTag, int M, int N>
void test_accum_colwise()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<cont, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<DTag, double, 1, N>::mat_t dmat_t;
	typedef dense_row<double, N> rrow_t;

	mat_host<cont, double, M, N> src(m, n);
	src.fill_rand();
	mat_host<DTag, double, 1, N> dst(1, n);

	smat_t smat = src.get_cmat();
	dmat_t drow = dst.get_mat();
	rrow_t rrow(n);

	matrix_shape<M, N> shape(m, n);

	// sum

	for (index_t j = 0; j < n; ++j)
	{
		rrow[j] = double(2 * j + 1);
		drow[j] = rrow[j];
	}

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rrow[j] += smat(i, j);
		}
	}

	percol(ewise(my_sum_kernel<double>(), U()), shape,
			in_(smat), in_out_(drow, atags::colwise_sum()));

	ASSERT_MAT_APPROX( 1, n, drow, rrow, 1.0e-12 );

	// max

	for (index_t j = 0; j < n; ++j)
	{
		rrow[j] = 0.5;
		drow[j] = rrow[j];
	}

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (rrow[j] < smat(i, j))
				rrow[j] = smat(i, j);
		}
	}

	percol(ewise(my_max_kernel<double>(), U()), shape,
			in_(smat), in_out_(drow, atags::colwise_max()));

	ASSERT_MAT_EQ( 1, n, drow, rrow );


	// min

	for (index_t j = 0; j < n; ++j)
	{
		rrow[j] = 0.5;
		drow[j] = rrow[j];
	}

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (rrow[j] > smat(i, j))
				rrow[j] = smat(i, j);
		}
	}

	percol(ewise(my_min_kernel<double>(), U()), shape,
			in_(smat), in_out_(drow, atags::colwise_min()));

	ASSERT_MAT_EQ( 1, n, drow, rrow );
}


template<typename U, typename DTag, int M, int N>
void test_accum_rowwise()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<cont, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<DTag, double, M, 1>::mat_t dmat_t;
	typedef dense_col<double, M> rcol_t;

	mat_host<cont, double, M, N> src(m, n);
	src.fill_rand();
	mat_host<DTag, double, M, 1> dst(m, 1);

	smat_t smat = src.get_cmat();
	dmat_t dcol = dst.get_mat();
	rcol_t rcol(m);

	matrix_shape<M, N> shape(m, n);

	// sum

	for (index_t i = 0; i < m; ++i)
	{
		rcol[i] = double(2 * i + 1);
		dcol[i] = rcol[i];
	}

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rcol[i] += smat(i, j);
		}
	}

	percol(ewise(my_sum_kernel<double>(), U()), shape,
			in_(smat), in_out_(dcol, atags::rowwise_sum()));

	ASSERT_MAT_APPROX( m, 1, dcol, rcol, 1.0e-12 );

	// max

	for (index_t i = 0; i < m; ++i)
	{
		rcol[i] = 0.5;
		dcol[i] = rcol[i];
	}

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (rcol[i] < smat(i, j))
				rcol[i] = smat(i, j);
		}
	}

	percol(ewise(my_max_kernel<double>(), U()), shape,
			in_(smat), in_out_(dcol, atags::rowwise_max()));

	ASSERT_MAT_EQ( m, 1, dcol, rcol );

	// min

	for (index_t i = 0; i < m; ++i)
	{
		rcol[i] = 0.5;
		dcol[i] = rcol[i];
	}

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (rcol[i] > smat(i, j))
				rcol[i] = smat(i, j);
		}
	}

	percol(ewise(my_min_kernel<double>(), U()), shape,
			in_(smat), in_out_(dcol, atags::rowwise_min()));

	ASSERT_MAT_EQ( m, 1, dcol, rcol );
}



// accum linear

MN_CASE( ewise_accum, linear_scalar )
{
	test_linear_accum<atags::scalar, M, N>();
}

MN_CASE( ewise_accum, linear_sse )
{
	test_linear_accum<atags::simd<sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, linear_avx )
{
	test_linear_accum<atags::simd<avx_t>, M, N>();
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


// accum percol

MN_CASE( ewise_accum, percol_scalar )
{
	test_percol_accum<atags::scalar, M, N>();
}

MN_CASE( ewise_accum, percol_sse )
{
	test_percol_accum<atags::simd<sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, percol_avx )
{
	test_percol_accum<atags::simd<avx_t>, M, N>();
}

#endif

BEGIN_TPACK( accum_percol_scalar )
	ADD_MN_CASE_3X3( ewise_accum, percol_scalar, DM, DN );
END_TPACK

BEGIN_TPACK( accum_percol_sse )
	ADD_MN_CASE_3X3( ewise_accum, percol_sse, DM, DN );
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( accum_percol_avx )
	ADD_MN_CASE_3X3( ewise_accum, percol_avx, DM, DN );
END_TPACK

#endif


// accum colwise

MN_CASE( ewise_accum, colwise_scalar_cont )
{
	test_accum_colwise<atags::scalar, cont, M, N>();
}

MN_CASE( ewise_accum, colwise_sse_cont )
{
	test_accum_colwise<atags::simd<sse_t>, cont, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, colwise_avx_cont )
{
	test_accum_colwise<atags::simd<avx_t>, cont, M, N>();
}

#endif

MN_CASE( ewise_accum, colwise_scalar_bloc )
{
	test_accum_colwise<atags::scalar, bloc, M, N>();
}

MN_CASE( ewise_accum, colwise_sse_bloc )
{
	test_accum_colwise<atags::simd<sse_t>, bloc, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, colwise_avx_bloc )
{
	test_accum_colwise<atags::simd<avx_t>, bloc, M, N>();
}

#endif

MN_CASE( ewise_accum, colwise_scalar_grid )
{
	test_accum_colwise<atags::scalar, grid, M, N>();
}

MN_CASE( ewise_accum, colwise_sse_grid )
{
	test_accum_colwise<atags::simd<sse_t>, grid, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, colwise_avx_grid )
{
	test_accum_colwise<atags::simd<avx_t>, grid, M, N>();
}

#endif


BEGIN_TPACK( accum_colwise_scalar_cont )
	ADD_MN_CASE_3X3( ewise_accum, colwise_scalar_cont, DM, DN )
END_TPACK

BEGIN_TPACK( accum_colwise_sse_cont )
	ADD_MN_CASE_3X3( ewise_accum, colwise_sse_cont, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( accum_colwise_avx_cont )
	ADD_MN_CASE_3X3( ewise_accum, colwise_avx_cont, DM, DN )
END_TPACK

#endif

BEGIN_TPACK( accum_colwise_scalar_bloc )
	ADD_MN_CASE_3X3( ewise_accum, colwise_scalar_bloc, DM, DN )
END_TPACK

BEGIN_TPACK( accum_colwise_sse_bloc )
	ADD_MN_CASE_3X3( ewise_accum, colwise_sse_bloc, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( accum_colwise_avx_bloc )
	ADD_MN_CASE_3X3( ewise_accum, colwise_avx_bloc, DM, DN )
END_TPACK

#endif

BEGIN_TPACK( accum_colwise_scalar_grid )
	ADD_MN_CASE_3X3( ewise_accum, colwise_scalar_grid, DM, DN )
END_TPACK

BEGIN_TPACK( accum_colwise_sse_grid )
	ADD_MN_CASE_3X3( ewise_accum, colwise_sse_grid, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( accum_colwise_avx_grid )
	ADD_MN_CASE_3X3( ewise_accum, colwise_avx_grid, DM, DN )
END_TPACK

#endif



// accum rowwise

MN_CASE( ewise_accum, rowwise_scalar_cont )
{
	test_accum_rowwise<atags::scalar, cont, M, N>();
}

MN_CASE( ewise_accum, rowwise_sse_cont )
{
	test_accum_rowwise<atags::simd<sse_t>, cont, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, rowwise_avx_cont )
{
	test_accum_rowwise<atags::simd<avx_t>, cont, M, N>();
}

#endif

MN_CASE( ewise_accum, rowwise_scalar_bloc )
{
	test_accum_rowwise<atags::scalar, bloc, M, N>();
}

MN_CASE( ewise_accum, rowwise_sse_bloc )
{
	test_accum_rowwise<atags::simd<sse_t>, bloc, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( ewise_accum, rowwise_avx_bloc )
{
	test_accum_rowwise<atags::simd<avx_t>, bloc, M, N>();
}

#endif

MN_CASE( ewise_accum, rowwise_scalar_grid )
{
	test_accum_rowwise<atags::scalar, grid, M, N>();
}


BEGIN_TPACK( accum_rowwise_scalar_cont )
	ADD_MN_CASE_3X3( ewise_accum, rowwise_scalar_cont, DM, DN )
END_TPACK

BEGIN_TPACK( accum_rowwise_sse_cont )
	ADD_MN_CASE_3X3( ewise_accum, rowwise_sse_cont, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( accum_rowwise_avx_cont )
	ADD_MN_CASE_3X3( ewise_accum, rowwise_avx_cont, DM, DN )
END_TPACK

#endif

BEGIN_TPACK( accum_rowwise_scalar_bloc )
	ADD_MN_CASE_3X3( ewise_accum, rowwise_scalar_bloc, DM, DN )
END_TPACK

BEGIN_TPACK( accum_rowwise_sse_bloc )
	ADD_MN_CASE_3X3( ewise_accum, rowwise_sse_bloc, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( accum_rowwise_avx_bloc )
	ADD_MN_CASE_3X3( ewise_accum, rowwise_avx_bloc, DM, DN )
END_TPACK

#endif

BEGIN_TPACK( accum_rowwise_scalar_grid )
	ADD_MN_CASE_3X3( ewise_accum, rowwise_scalar_grid, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE

	// accum linear

	ADD_TPACK( accum_linear_scalar )
	ADD_TPACK( accum_linear_sse )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_linear_avx )
#endif

	// accum percol

	ADD_TPACK( accum_percol_scalar )
	ADD_TPACK( accum_percol_sse )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_percol_avx )
#endif

	// accum colwise

	ADD_TPACK( accum_colwise_scalar_cont )
	ADD_TPACK( accum_colwise_sse_cont )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_colwise_avx_cont )
#endif

	ADD_TPACK( accum_colwise_scalar_bloc )
	ADD_TPACK( accum_colwise_sse_bloc )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_colwise_avx_bloc )
#endif

	ADD_TPACK( accum_colwise_scalar_grid )
	ADD_TPACK( accum_colwise_sse_grid )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_colwise_avx_grid )
#endif

	// accum rowwise

	ADD_TPACK( accum_rowwise_scalar_cont )
	ADD_TPACK( accum_rowwise_sse_cont )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_rowwise_avx_cont )
#endif

	ADD_TPACK( accum_rowwise_scalar_bloc )
	ADD_TPACK( accum_rowwise_sse_bloc )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( accum_rowwise_avx_bloc )
#endif

	ADD_TPACK( accum_rowwise_scalar_grid )

END_MAIN_SUITE




