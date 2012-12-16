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

MN_CASE( percol_ewise, cont_cont  )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	typedef typename mat_host<cont, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<cont, double, M, N>::mat_t dmat_t;

	mat_host<cont, double, M, N> src(m, n);
	src.fill_lin();
	mat_host<cont, double, M, N> dst(m, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<M, N> shape(m, n);

	percol_ewise(atags::scalar(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, n, smat, dmat);
}


BEGIN_TPACK( percol_ewise_cont_cont )
	ADD_MN_CASE_3X3( percol_ewise, cont_cont, DM, DN )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( percol_ewise_cont_cont )
END_MAIN_SUITE
