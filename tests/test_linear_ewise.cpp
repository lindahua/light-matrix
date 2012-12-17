/**
 * @file test_linear_ewise.cpp
 *
 * @brief Test Linear element-wise accesses
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

// core functions

struct my_update_kernel
{
	template<typename T>
	LMAT_ENSURE_INLINE
	void operator() (const T& s, T& d) const
	{
		d += s;
	}
};


template<typename U, int M, int N>
void test_linear_ewise_cont_cont()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<cont, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<cont, double, M, N>::mat_t dmat_t;

	mat_host<cont, double, M, N> src(m, n);
	src.fill_lin();
	mat_host<cont, double, M, N> dst(m, n);

	dense_matrix<double, M, N> rmat(m, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<M, N> shape(m, n);

	linear_ewise(U(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, n, smat, dmat);

	for (index_t i = 0; i < m * n; ++i) rmat[i] = smat[i] + dmat[i];

	linear_ewise(U(), shape).apply(
			my_update_kernel(), in_(smat), in_out_(dmat));

	ASSERT_MAT_EQ(m, n, dmat, rmat);
}

template<typename U, typename STag, typename DTag, int M>
void test_linear_ewise_col()
{
	const index_t m = M == 0 ? DM : M;

	typedef typename mat_host<STag, double, M, 1>::cmat_t smat_t;
	typedef typename mat_host<DTag, double, M, 1>::mat_t dmat_t;

	mat_host<STag, double, M, 1> src(m, 1);
	src.fill_lin();
	mat_host<DTag, double, M, 1> dst(m, 1);

	dense_matrix<double, M, 1> rmat(m, 1);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<M, 1> shape(m, 1);

	linear_ewise(U(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, 1, smat, dmat);

	for (index_t i = 0; i < m; ++i) rmat[i] = smat[i] + dmat[i];

	linear_ewise(U(), shape).apply(
			my_update_kernel(), in_(smat), in_out_(dmat));

	ASSERT_MAT_EQ(m, 1, dmat, rmat);
}

template<typename U, typename STag, typename DTag, int N>
void test_linear_ewise_row()
{
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<STag, double, 1, N>::cmat_t smat_t;
	typedef typename mat_host<DTag, double, 1, N>::mat_t dmat_t;

	mat_host<STag, double, 1, N> src(1, n);
	src.fill_lin();
	mat_host<DTag, double, 1, N> dst(1, n);

	dense_matrix<double, 1, N> rmat(1, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<1, N> shape(1, n);

	linear_ewise(U(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(1, n, smat, dmat);

	for (index_t i = 0; i < n; ++i) rmat[i] = smat[i] + dmat[i];

	linear_ewise(U(), shape).apply(
			my_update_kernel(), in_(smat), in_out_(dmat));

	ASSERT_MAT_EQ(1, n, dmat, rmat);
}


template<typename U, int M, int N>
void test_linear_ewise_single_cont()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<cont, double, M, N>::mat_t dmat_t;

	double v = 2.56;
	mat_host<cont, double, M, N> dst(m, n);

	dmat_t dmat = dst.get_mat();

	dense_matrix<double, M, N> rmat(m, n);
	fill(rmat, v);

	matrix_shape<M, N> shape(m, n);

	linear_ewise(U(), shape).apply(
			copy_kernel(), in_(v, atags::single()), out_(dmat));

	ASSERT_MAT_EQ(m, n, dmat, rmat);
}



// Specific test cases


MN_CASE( linear_ewise, scalar_cont_cont  )
{
	test_linear_ewise_cont_cont<atags::scalar, M, N>();
}

N_CASE( linear_ewise, scalar_cont_stepcol  )
{
	test_linear_ewise_col<atags::scalar, cont, grid, N>();
}

N_CASE( linear_ewise, scalar_stepcol_cont  )
{
	test_linear_ewise_col<atags::scalar, grid, cont, N>();
}

N_CASE( linear_ewise, scalar_stepcol_stepcol  )
{
	test_linear_ewise_col<atags::scalar, grid, grid, N>();
}

N_CASE( linear_ewise, scalar_cont_steprow  )
{
	test_linear_ewise_row<atags::scalar, cont, bloc, N>();
}

N_CASE( linear_ewise, scalar_steprow_cont  )
{
	test_linear_ewise_row<atags::scalar, bloc, cont, N>();
}

N_CASE( linear_ewise, scalar_steprow_steprow  )
{
	test_linear_ewise_row<atags::scalar, bloc, bloc, N>();
}


MN_CASE( linear_ewise, sse_cont_cont  )
{
	test_linear_ewise_cont_cont<atags::simd<double, sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( linear_ewise, avx_cont_cont  )
{
	test_linear_ewise_cont_cont<atags::simd<double, avx_t>, M, N>();
}

#endif


MN_CASE( linear_ewise, scalar_single_cont )
{
	test_linear_ewise_single_cont<atags::scalar, M, N>();
}

MN_CASE( linear_ewise, sse_single_cont )
{
	test_linear_ewise_single_cont<atags::simd<double, sse_t>, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( linear_ewise, avx_single_cont )
{
	test_linear_ewise_single_cont<atags::simd<double, avx_t>, M, N>();
}

#endif


// Test packs

BEGIN_TPACK( linear_ewise_scalar_cont_cont )
	ADD_MN_CASE_3X3( linear_ewise, scalar_cont_cont, DM, DN )
END_TPACK

BEGIN_TPACK( linear_ewise_scalar_cont_stepcol )
	ADD_N_CASE_3( linear_ewise, scalar_cont_stepcol, DM )
END_TPACK

BEGIN_TPACK( linear_ewise_scalar_stepcol_cont )
	ADD_N_CASE_3( linear_ewise, scalar_stepcol_cont, DM )
END_TPACK

BEGIN_TPACK( linear_ewise_scalar_stepcol_stepcol )
	ADD_N_CASE_3( linear_ewise, scalar_stepcol_stepcol, DM )
END_TPACK

BEGIN_TPACK( linear_ewise_scalar_cont_steprow )
	ADD_N_CASE_3( linear_ewise, scalar_cont_steprow, DN )
END_TPACK

BEGIN_TPACK( linear_ewise_scalar_steprow_cont )
	ADD_N_CASE_3( linear_ewise, scalar_steprow_cont, DN )
END_TPACK

BEGIN_TPACK( linear_ewise_scalar_steprow_steprow )
	ADD_N_CASE_3( linear_ewise, scalar_steprow_steprow, DN )
END_TPACK


BEGIN_TPACK( linear_ewise_sse_cont_cont )
	ADD_MN_CASE_3X3( linear_ewise, sse_cont_cont, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( linear_ewise_avx_cont_cont )
	ADD_MN_CASE_3X3( linear_ewise, avx_cont_cont, DM, DN )
END_TPACK

#endif

BEGIN_TPACK( linear_ewise_scalar_single_cont )
	ADD_MN_CASE_3X3( linear_ewise, scalar_single_cont, DM, DN )
END_TPACK

BEGIN_TPACK( linear_ewise_sse_single_cont )
	ADD_MN_CASE_3X3( linear_ewise, sse_single_cont, DM, DN )
END_TPACK

#ifdef LMAT_HAS_AVX

BEGIN_TPACK( linear_ewise_avx_single_cont )
	ADD_MN_CASE_3X3( linear_ewise, avx_single_cont, DM, DN )
END_TPACK

#endif


BEGIN_MAIN_SUITE
	ADD_TPACK( linear_ewise_scalar_cont_cont )
	ADD_TPACK( linear_ewise_scalar_cont_stepcol )
	ADD_TPACK( linear_ewise_scalar_stepcol_cont )
	ADD_TPACK( linear_ewise_scalar_stepcol_stepcol )
	ADD_TPACK( linear_ewise_scalar_cont_steprow )
	ADD_TPACK( linear_ewise_scalar_steprow_cont )
	ADD_TPACK( linear_ewise_scalar_steprow_steprow )

	ADD_TPACK( linear_ewise_sse_cont_cont )
#ifdef LMAT_HAS_AVX
	ADD_TPACK( linear_ewise_avx_cont_cont )
#endif

	ADD_TPACK( linear_ewise_scalar_single_cont )
	ADD_TPACK( linear_ewise_sse_single_cont )
	ADD_TPACK( linear_ewise_avx_single_cont )
END_MAIN_SUITE


