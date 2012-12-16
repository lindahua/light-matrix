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

using namespace lmat;
using namespace lmat::test;


// core functions

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

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<M, N> shape(m, n);

	linear_ewise(U(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, n, smat, dmat);
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

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<M, 1> shape(m, 1);

	linear_ewise(U(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, 1, smat, dmat);
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

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<1, N> shape(1, n);

	linear_ewise(U(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(1, n, smat, dmat);
}


MN_CASE( linear_ewise, cont_cont  )
{
	test_linear_ewise_cont_cont<atags::scalar, M, N>();
}


N_CASE( linear_ewise, cont_stepcol  )
{
	test_linear_ewise_col<atags::scalar, cont, grid, N>();
}


N_CASE( linear_ewise, stepcol_cont  )
{
	test_linear_ewise_col<atags::scalar, grid, cont, N>();
}

N_CASE( linear_ewise, stepcol_stepcol  )
{
	test_linear_ewise_col<atags::scalar, grid, grid, N>();
}


N_CASE( linear_ewise, cont_steprow  )
{
	test_linear_ewise_row<atags::scalar, cont, bloc, N>();
}


N_CASE( linear_ewise, steprow_cont  )
{
	test_linear_ewise_row<atags::scalar, bloc, cont, N>();
}


N_CASE( linear_ewise, steprow_steprow  )
{
	test_linear_ewise_row<atags::scalar, bloc, bloc, N>();
}


BEGIN_TPACK( linear_ewise_cont_cont )
	ADD_MN_CASE_3X3( linear_ewise, cont_cont, DM, DN )
END_TPACK

BEGIN_TPACK( linear_ewise_cont_stepcol )
	ADD_N_CASE_3( linear_ewise, cont_stepcol, DM )
END_TPACK

BEGIN_TPACK( linear_ewise_stepcol_cont )
	ADD_N_CASE_3( linear_ewise, stepcol_cont, DM )
END_TPACK

BEGIN_TPACK( linear_ewise_stepcol_stepcol )
	ADD_N_CASE_3( linear_ewise, stepcol_stepcol, DM )
END_TPACK

BEGIN_TPACK( linear_ewise_cont_steprow )
	ADD_N_CASE_3( linear_ewise, cont_steprow, DN )
END_TPACK

BEGIN_TPACK( linear_ewise_steprow_cont )
	ADD_N_CASE_3( linear_ewise, steprow_cont, DN )
END_TPACK

BEGIN_TPACK( linear_ewise_steprow_steprow )
	ADD_N_CASE_3( linear_ewise, steprow_steprow, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( linear_ewise_cont_cont )
	ADD_TPACK( linear_ewise_cont_stepcol )
	ADD_TPACK( linear_ewise_stepcol_cont )
	ADD_TPACK( linear_ewise_stepcol_stepcol )
	ADD_TPACK( linear_ewise_cont_steprow )
	ADD_TPACK( linear_ewise_steprow_cont )
	ADD_TPACK( linear_ewise_steprow_steprow )
END_MAIN_SUITE


