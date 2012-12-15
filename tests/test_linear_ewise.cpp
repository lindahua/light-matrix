/**
 * @file test_linear_ewise.cpp
 *
 * @brief Test Linear element-wise accesses
 *
 * @author Dahua Lin
 */


#include "test_base.h"
#include "multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/ewise_eval.h>

using namespace lmat;
using namespace lmat::test;


// Auxiliary functions

MN_CASE( linear_ewise, cont_cont  )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	typedef typename mat_host<cont, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<cont, double, M, N>::mat_t dmat_t;

	mat_host<STag, double, M, N> src(m, n);
	src.fill_lin();
	mat_host<DTag, double, M, N> dst(m, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	matrix_shape<M, N> shape(m, n);

	linear_ewise(atags::scalar(), shape).apply(
			copy_kernel(), in_(smat), out_(dmat));

	ASSERT_MAT_EQ(m, n, smat, dmat);
}



BEGIN_TPACK( linear_ewise_cont_cont )
	ADD_MN_CASE( linear_ewise, cont_cont, 0, 0 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( linear_ewise_cont_cont )
END_MAIN_SUITE


