/**
 * @file test_mat_fill.cpp
 *
 * @brief Unit testing for matrix filling
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include "../multimat_supp.h"

#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/matrix/ref_block.h>
#include <light_mat/matrix/ref_grid.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;


// Auxiliary classes

template<class Mat, typename T>
bool verify_all_equal(const Mat& mat, const T& v)
{
	const index_t m = mat.nrows();
	const index_t n = mat.ncolumns();

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (mat(i, j) != v) return false;
		}
	}

	return true;
}


template<typename DTag, int M, int N>
void test_matrix_zero()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	typedef typename mat_host<DTag, double, M, N>::mat_t mat_t;
	mat_host<DTag, double, M, N> dst(m, n);
	dst.fill_lin();

	mat_t dmat = dst.get_mat();
	zero(dmat);

	ASSERT_TRUE( verify_all_equal(dmat, 0.0) );
}


template<typename DTag, int M, int N>
void test_matrix_fill()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	typedef typename mat_host<DTag, double, M, N>::mat_t mat_t;
	mat_host<DTag, double, M, N> dst(m, n);

	mat_t dmat = dst.get_mat();
	double v = 12;
	fill(dmat, v);

	ASSERT_TRUE( verify_all_equal(dmat, v) );
}


MN_CASE( mat_zero, cont )
{
	test_matrix_zero<cont, M, N>();
}

MN_CASE( mat_zero, bloc )
{
	test_matrix_zero<bloc, M, N>();
}

MN_CASE( mat_zero, grid )
{
	test_matrix_zero<grid, M, N>();
}


MN_CASE( mat_fill, cont )
{
	test_matrix_fill<cont, M, N>();
}

MN_CASE( mat_fill, bloc )
{
	test_matrix_fill<bloc, M, N>();
}

MN_CASE( mat_fill, grid )
{
	test_matrix_fill<grid, M, N>();
}


BEGIN_TPACK( mat_zero_cont )
	ADD_MN_CASE_3X3( mat_zero, cont, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_zero_bloc )
	ADD_MN_CASE_3X3( mat_zero, bloc, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_zero_grid )
	ADD_MN_CASE_3X3( mat_zero, grid, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_fill_cont )
	ADD_MN_CASE_3X3( mat_fill, cont, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_fill_bloc )
	ADD_MN_CASE_3X3( mat_fill, bloc, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_fill_grid )
	ADD_MN_CASE_3X3( mat_fill, grid, 3, 4 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_zero_cont )
	ADD_TPACK( mat_zero_bloc )
	ADD_TPACK( mat_zero_grid )

	ADD_TPACK( mat_fill_cont )
	ADD_TPACK( mat_fill_bloc )
	ADD_TPACK( mat_fill_grid )
END_MAIN_SUITE





