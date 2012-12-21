/**
 * @file test_percol_ewise.cpp
 *
 * @brief Unit testing of percol ewise evaluation
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#define DEFAULT_M_VALUE 13
#define DEFAULT_N_VALUE 9

#include "multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/math/math_functors.h>
#include <light_mat/mateval/ewise_eval.h>

#include <light_mat/math/sse_ops.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_ops.h>
#endif


using namespace lmat;
using namespace lmat::test;


// test cases

template<typename STag, typename DTag, typename U, int M, int N>
void test_percol_ewise()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<STag, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<DTag, double, M, N>::mat_t dmat_t;

	mat_host<STag, double, M, N> src(m, n);
	src.fill_lin();
	mat_host<DTag, double, M, N> dst(m, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<M, N> shape(m, n);

	copy_kernel<double> cpy_kernel;
	accum_kernel<double> upd_kernel;

	percol(ewise(cpy_kernel, U()), shape, in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, n, smat, dmat);

	dense_matrix<double, M, N> rmat(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) rmat(i, j) = smat(i, j) + dmat(i, j);
	}

	percol(ewise(upd_kernel, U()), shape, in_out_(dmat), in_(smat));

	ASSERT_MAT_EQ(m, n, dmat, rmat);
}


template<typename DTag, typename U, int M, int N>
void test_percol_ewise_single()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<DTag, double, M, N>::mat_t dmat_t;

	const double v = 2.56;
	mat_host<DTag, double, M, N> dst(m, n);

	dmat_t dmat = dst.get_mat();

	dense_matrix<double, M, N> rmat(m, n);
	fill(rmat, v);

	matrix_shape<M, N> shape(m, n);

	percol(ewise(copy_kernel<double>(), U()), shape, in_(v, atags::single()), out_(dmat));

	ASSERT_MAT_EQ(m, n, dmat, rmat);
}


template<typename DTag, typename U, int M, int N>
void test_percol_ewise_repcol()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<cont, double, M, 1>::cmat_t col_t;
	typedef typename mat_host<DTag, double, M, N>::mat_t dmat_t;

	mat_host<cont, double, M, 1> src(m, 1);
	src.fill_lin();
	mat_host<DTag, double, M, N> dst(m, n);

	col_t col = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	dense_matrix<double, M, N> rmat(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rmat(i, j) = col[i];
		}
	}

	matrix_shape<M, N> shape(m, n);

	percol(ewise(copy_kernel<double>(), U()), shape, in_(col, atags::repcol()), out_(dmat));

	ASSERT_MAT_EQ(m, n, dmat, rmat);
}


template<typename DTag, typename U, int M, int N>
void test_percol_ewise_reprow()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<cont, double, 1, N>::cmat_t row_t;
	typedef typename mat_host<DTag, double, M, N>::mat_t dmat_t;

	mat_host<cont, double, 1, N> src(1, n);
	src.fill_lin();
	mat_host<DTag, double, M, N> dst(m, n);

	row_t row = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	dense_matrix<double, M, N> rmat(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rmat(i, j) = row[j];
		}
	}

	matrix_shape<M, N> shape(m, n);

	percol(ewise(copy_kernel<double>(), U()), shape, in_(row, atags::reprow()), out_(dmat));

	ASSERT_MAT_EQ(m, n, dmat, rmat);
}


// Specific test cases

#define DEFINE_PERCOL_EWISE_SCALAR_TEST( STag, DTag ) \
		MN_CASE( percol_ewise, scalar_##STag##_##DTag  ) { \
			test_percol_ewise<STag, DTag, atags::scalar, M, N>(); } \
		BEGIN_TPACK( percol_ewise_scalar_##STag##_##DTag ) \
			ADD_MN_CASE_3X3( percol_ewise, scalar_##STag##_##DTag, DM, DN ) \
		END_TPACK

#define DEFINE_PERCOL_EWISE_SIMD_TEST( SKindName, STag, DTag ) \
		MN_CASE( percol_ewise, SKindName##_##STag##_##DTag  ) { \
			test_percol_ewise<STag, DTag, atags::simd<SKindName##_t>, M, N>(); } \
		BEGIN_TPACK( percol_ewise_##SKindName##_##STag##_##DTag ) \
			ADD_MN_CASE_3X3( percol_ewise, scalar_##STag##_##DTag, DM, DN ) \
		END_TPACK

DEFINE_PERCOL_EWISE_SCALAR_TEST( cont, cont )
DEFINE_PERCOL_EWISE_SCALAR_TEST( cont, bloc )
DEFINE_PERCOL_EWISE_SCALAR_TEST( cont, grid )

DEFINE_PERCOL_EWISE_SCALAR_TEST( bloc, cont )
DEFINE_PERCOL_EWISE_SCALAR_TEST( bloc, bloc )
DEFINE_PERCOL_EWISE_SCALAR_TEST( bloc, grid )

DEFINE_PERCOL_EWISE_SCALAR_TEST( grid, cont )
DEFINE_PERCOL_EWISE_SCALAR_TEST( grid, bloc )
DEFINE_PERCOL_EWISE_SCALAR_TEST( grid, grid )

DEFINE_PERCOL_EWISE_SIMD_TEST( sse, cont, cont )
DEFINE_PERCOL_EWISE_SIMD_TEST( sse, cont, bloc )
DEFINE_PERCOL_EWISE_SIMD_TEST( sse, bloc, cont )
DEFINE_PERCOL_EWISE_SIMD_TEST( sse, bloc, bloc )

#ifdef LMAT_HAS_AVX

DEFINE_PERCOL_EWISE_SIMD_TEST( avx, cont, cont )
DEFINE_PERCOL_EWISE_SIMD_TEST( avx, cont, bloc )
DEFINE_PERCOL_EWISE_SIMD_TEST( avx, bloc, cont )
DEFINE_PERCOL_EWISE_SIMD_TEST( avx, bloc, bloc )

#endif


// Testing single reader

MN_CASE( percol_ewise, scalar_single_cont )
{
	test_percol_ewise_single<cont, atags::scalar, M, N>();
}

MN_CASE( percol_ewise, sse_single_cont )
{
	test_percol_ewise_single<cont, atags::simd<sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( percol_ewise, avx_single_cont )
{
	test_percol_ewise_single<cont, atags::simd<avx_t>, M, N>();
}

#endif

BEGIN_TPACK( percol_ewise_scalar_single_cont )
	ADD_MN_CASE_3X3( percol_ewise, scalar_single_cont, DM, DN )
END_TPACK

BEGIN_TPACK( percol_ewise_sse_single_cont )
	ADD_MN_CASE_3X3( percol_ewise, sse_single_cont, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( percol_ewise_avx_single_cont )
	ADD_MN_CASE_3X3( percol_ewise, avx_single_cont, DM, DN )
END_TPACK

#endif



// Testing repcol reader

MN_CASE( percol_ewise, scalar_repcol_cont )
{
	test_percol_ewise_repcol<cont, atags::scalar, M, N>();
}

MN_CASE( percol_ewise, sse_repcol_cont )
{
	test_percol_ewise_repcol<cont, atags::simd<sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( percol_ewise, avx_repcol_cont )
{
	test_percol_ewise_repcol<cont, atags::simd<avx_t>, M, N>();
}

#endif


BEGIN_TPACK( percol_ewise_scalar_repcol_cont )
	ADD_MN_CASE_3X3( percol_ewise, scalar_repcol_cont, DM, DN )
END_TPACK

BEGIN_TPACK( percol_ewise_sse_repcol_cont )
	ADD_MN_CASE_3X3( percol_ewise, sse_repcol_cont, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( percol_ewise_avx_repcol_cont )
	ADD_MN_CASE_3X3( percol_ewise, avx_repcol_cont, DM, DN )
END_TPACK

#endif



// Testing reprow reader

MN_CASE( percol_ewise, scalar_reprow_cont )
{
	test_percol_ewise_reprow<cont, atags::scalar, M, N>();
}

MN_CASE( percol_ewise, sse_reprow_cont )
{
	test_percol_ewise_reprow<cont, atags::simd<sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( percol_ewise, avx_reprow_cont )
{
	test_percol_ewise_reprow<cont, atags::simd<avx_t>, M, N>();
}

#endif


BEGIN_TPACK( percol_ewise_scalar_reprow_cont )
	ADD_MN_CASE_3X3( percol_ewise, scalar_reprow_cont, DM, DN )
END_TPACK

BEGIN_TPACK( percol_ewise_sse_reprow_cont )
	ADD_MN_CASE_3X3( percol_ewise, sse_reprow_cont, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( percol_ewise_avx_reprow_cont )
	ADD_MN_CASE_3X3( percol_ewise, avx_reprow_cont, DM, DN )
END_TPACK

#endif



BEGIN_MAIN_SUITE
	ADD_TPACK( percol_ewise_scalar_cont_cont )
	ADD_TPACK( percol_ewise_scalar_cont_bloc )
	ADD_TPACK( percol_ewise_scalar_cont_grid )
	ADD_TPACK( percol_ewise_scalar_bloc_cont )
	ADD_TPACK( percol_ewise_scalar_bloc_bloc )
	ADD_TPACK( percol_ewise_scalar_bloc_grid )
	ADD_TPACK( percol_ewise_scalar_grid_cont )
	ADD_TPACK( percol_ewise_scalar_grid_bloc )
	ADD_TPACK( percol_ewise_scalar_grid_grid )

	ADD_TPACK( percol_ewise_sse_cont_cont )
	ADD_TPACK( percol_ewise_sse_cont_bloc )
	ADD_TPACK( percol_ewise_sse_bloc_cont )
	ADD_TPACK( percol_ewise_sse_bloc_bloc )

#ifdef LMAT_HAS_AVX
	ADD_TPACK( percol_ewise_avx_cont_cont )
	ADD_TPACK( percol_ewise_avx_cont_bloc )
	ADD_TPACK( percol_ewise_avx_bloc_cont )
	ADD_TPACK( percol_ewise_avx_bloc_bloc )
#endif

	ADD_TPACK( percol_ewise_scalar_single_cont )
	ADD_TPACK( percol_ewise_sse_single_cont )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( percol_ewise_avx_single_cont )
#endif

	ADD_TPACK( percol_ewise_scalar_repcol_cont )
	ADD_TPACK( percol_ewise_sse_repcol_cont )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( percol_ewise_avx_repcol_cont )
#endif

	ADD_TPACK( percol_ewise_scalar_reprow_cont )
	ADD_TPACK( percol_ewise_sse_reprow_cont )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( percol_ewise_avx_reprow_cont )
#endif
END_MAIN_SUITE




