/**
 * @file test_direct_trans.cpp
 *
 * @brief Unit testing of direct transposition
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include "../multimat_supp.h"

#include <light_mat/matrix/matrix_transpose.h>


template<class Tag1, class Tag2, int M, int N>
void test_direct_trans()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	// M x N ==> N x M

	typedef typename mat_host<Tag1, double, M, N>::cmat_t cmat_t;
	typedef typename mat_host<Tag2, double, N, M>::mat_t mat_t;

	mat_host<Tag1, double, M, N> src(m, n);
	mat_host<Tag2, double, N, M> dst(n, m);

	src.fill_lin();

	cmat_t smat = src.get_cmat();
	mat_t dmat = dst.get_mat();

	transpose(smat, dmat);

	dense_matrix<double, N, M> rmat(n, m);
	for (index_t i = 0; i < m; ++i)
	{
		for (index_t j = 0; j < n; ++j) rmat(j, i) = smat(i, j);
	}

	ASSERT_MAT_EQ(n, m, dmat, rmat);

	// N x M ==> M x N

	typedef typename mat_host<Tag1, double, N, M>::cmat_t cmat2_t;
	typedef typename mat_host<Tag2, double, M, N>::mat_t mat2_t;

	mat_host<Tag1, double, N, M> src2(n, m);
	mat_host<Tag2, double, M, N> dst2(m, n);

	src2.fill_lin();

	cmat2_t smat2 = src2.get_cmat();
	mat2_t dmat2 = dst2.get_mat();

	transpose(smat2, dmat2);

	dense_matrix<double, M, N> rmat2(m, n);
	for (index_t i = 0; i < m; ++i)
	{
		for (index_t j = 0; j < n; ++j) rmat2(i, j) = smat2(j, i);
	}

	ASSERT_MAT_EQ(m, n, dmat2, rmat2);

}


LTEST_INIT_AUTOSUITE

#define TEST_DIRECT_TRANS( sform, dform, name ) \
	MN_CASE( direct_trans_##name ) \
	{ test_direct_trans<sform, dform, M, N>(); } \
	AUTO_TPACK( direct_trans_##name ) { \
	ADD_MN_CASE_3X3( direct_trans_##name, DM, DN ) \
	}

TEST_DIRECT_TRANS( cont, cont, mat_to_mat )
TEST_DIRECT_TRANS( cont, bloc, mat_to_blk )
TEST_DIRECT_TRANS( cont, grid, mat_to_grid )

TEST_DIRECT_TRANS( bloc, cont, blk_to_mat )
TEST_DIRECT_TRANS( bloc, bloc, blk_to_blk )
TEST_DIRECT_TRANS( bloc, grid, blk_to_grid )

TEST_DIRECT_TRANS( grid, cont, grid_to_mat )
TEST_DIRECT_TRANS( grid, bloc, grid_to_blk )
TEST_DIRECT_TRANS( grid, grid, grid_to_grid )



