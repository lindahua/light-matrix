/**
 * @file test_direct_trans.cpp
 *
 * @brief Unit testing of direct transposition
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include "multimat_supp.h"

#include <light_mat/matexpr/matrix_transpose.h>


template<
	template<typename T1, int M1, int N1> class ClassT1,
	template<typename T2, int M2, int N2> class ClassT2,
	int M, int N
>
void test_direct_trans()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	const index_t max_dim = m > n ? m : n;

	dblock<double> sblk(LDim * max_dim, zero());
	dblock<double> dblk(LDim * max_dim, zero());

	for (index_t i = 0; i < sblk.nelems(); ++i) sblk[i] = double(i + 1);

	// M x N ==> N x M

	typedef typename mat_maker<ClassT1, M, N>::cmat_t cmat_t;
	typedef typename mat_maker<ClassT2, N, M>::mat_t mat_t;

	cmat_t smat = mat_maker<ClassT1, M, N>::get_cmat(sblk.ptr_data(), m, n);
	mat_t dmat = mat_maker<ClassT2, N, M>::get_mat(dblk.ptr_data(), n, m);

	direct_transpose(smat, dmat);

	dense_matrix<double, N, M> rmat(n, m);
	for (index_t i = 0; i < m; ++i)
	{
		for (index_t j = 0; j < n; ++j) rmat(j, i) = smat(i, j);
	}

	ASSERT_MAT_EQ(n, m, dmat, rmat);

	// N x M ==> M x N

	typedef typename mat_maker<ClassT1, N, M>::cmat_t cmat2_t;
	typedef typename mat_maker<ClassT2, M, N>::mat_t mat2_t;

	dblk = zero();

	cmat2_t smat2 = mat_maker<ClassT1, N, M>::get_cmat(sblk.ptr_data(), n, m);
	mat2_t dmat2 = mat_maker<ClassT2, M, N>::get_mat(dblk.ptr_data(), m, n);

	direct_transpose(smat2, dmat2);

	dense_matrix<double, M, N> rmat2(m, n);
	for (index_t i = 0; i < m; ++i)
	{
		for (index_t j = 0; j < n; ++j) rmat2(i, j) = smat2(j, i);
	}

	ASSERT_MAT_EQ(m, n, dmat2, rmat2);

}


#define TEST_DIRECT_TRANS( sform, dform, name ) \
	MN_CASE( direct_trans, name ) \
	{ test_direct_trans<sform, dform, M, N>(); } \
	BEGIN_TPACK( direct_trans_##name ) \
	ADD_MN_CASE_3X3( direct_trans, name, DM, DN ) \
	END_TPACK

TEST_DIRECT_TRANS( ref_matrix, ref_matrix, mat_to_mat )
TEST_DIRECT_TRANS( ref_matrix, ref_block, mat_to_blk )
TEST_DIRECT_TRANS( ref_matrix, ref_grid, mat_to_grid )

TEST_DIRECT_TRANS( ref_block, ref_matrix, blk_to_mat )
TEST_DIRECT_TRANS( ref_block, ref_block, blk_to_blk )
TEST_DIRECT_TRANS( ref_block, ref_grid, blk_to_grid )

TEST_DIRECT_TRANS( ref_grid, ref_matrix, grid_to_mat )
TEST_DIRECT_TRANS( ref_grid, ref_block, grid_to_blk )
TEST_DIRECT_TRANS( ref_grid, ref_grid, grid_to_grid )

BEGIN_MAIN_SUITE
	ADD_TPACK( direct_trans_mat_to_mat )
	ADD_TPACK( direct_trans_mat_to_blk )
	ADD_TPACK( direct_trans_mat_to_grid )

	ADD_TPACK( direct_trans_blk_to_mat )
	ADD_TPACK( direct_trans_blk_to_blk )
	ADD_TPACK( direct_trans_blk_to_grid )

	ADD_TPACK( direct_trans_grid_to_mat )
	ADD_TPACK( direct_trans_grid_to_blk )
	ADD_TPACK( direct_trans_grid_to_grid )
END_MAIN_SUITE

