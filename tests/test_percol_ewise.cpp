/**
 * @file test_percol_ewise.cpp
 *
 * @brief Unit testing of percol ewise evaluation
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

	percol_ewise(U(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, n, smat, dmat);
}


#define DEFINE_PERCOL_EWISE_SCALAR_TEST( STag, DTag ) \
		MN_CASE( percol_ewise, scalar_##STag##_##DTag  ) { \
			test_percol_ewise<STag, DTag, atags::scalar, M, N>(); } \
		BEGIN_TPACK( percol_ewise_scalar_##STag##_##DTag ) \
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
END_MAIN_SUITE




