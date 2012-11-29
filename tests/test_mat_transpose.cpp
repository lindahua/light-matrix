/**
 * @file test_mat_transpose.cpp
 *
 * Unit testing for matrix transpose
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include "multimat_supp.h"

#include <light_mat/matexpr/matrix_arith.h>
#include <light_mat/matexpr/matrix_transpose.h>
#include <string>

using namespace lmat;
using namespace lmat::test;

template<class Mat>
void fill_lin(Mat& A)
{
	const index_t m = A.nrows();
	const index_t n = A.ncolumns();

	int v = 1;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			A(i, j) = double(v++);
		}
	}
}


template<
	template<typename T1, int M1, int N1> class ClassT1,
	template<typename T2, int M2, int N2> class ClassT2,
	int M, int N
>
void test_mat_transpose()
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

	dmat = trans(smat);

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

	dmat2 = trans(smat2);

	dense_matrix<double, M, N> rmat2(m, n);
	for (index_t i = 0; i < m; ++i)
	{
		for (index_t j = 0; j < n; ++j) rmat2(i, j) = smat2(j, i);
	}

	ASSERT_MAT_EQ(m, n, dmat2, rmat2);
}


template<
	template<typename T2, int M2, int N2> class ClassT2,
	int M, int N
>
void test_xpr_transpose()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	const index_t max_dim = m > n ? m : n;

	dblock<double> dblk(LDim * max_dim, zero());

	// M x N ==> N x M

	dense_matrix<double, M, N> smat(m, n);
	for (index_t i = 0; i < m * n; ++i) smat[i] = double(i+1);

	typedef typename mat_maker<ClassT2, N, M>::mat_t mat_t;
	mat_t dmat = mat_maker<ClassT2, N, M>::get_mat(dblk.ptr_data(), n, m);

	dmat = trans(smat * 2.0 + 1.0);

	dense_matrix<double, N, M> rmat(n, m);
	for (index_t i = 0; i < m; ++i)
	{
		for (index_t j = 0; j < n; ++j) rmat(j, i) = smat(i, j) * 2.0 + 1.0;
	}

	ASSERT_MAT_EQ(n, m, dmat, rmat);

	// N x M ==> M x N

	dense_matrix<double, N, M> smat2(n, m);
	for (index_t i = 0; i < n * m; ++i) smat2[i] = double(i+1);

	typedef typename mat_maker<ClassT2, M, N>::mat_t mat2_t;
	mat2_t dmat2 = mat_maker<ClassT2, M, N>::get_mat(dblk.ptr_data(), m, n);

	dmat2 = trans(smat2 * 2.0 + 1.0);

	dense_matrix<double, M, N> rmat2(m, n);
	for (index_t i = 0; i < m; ++i)
	{
		for (index_t j = 0; j < n; ++j) rmat2(i, j) = smat2(j, i) * 2.0 + 1.0;
	}

	ASSERT_MAT_EQ(m, n, dmat2, rmat2);
}



#define TEST_MAT_TRANS( sform, dform, name ) \
	MN_CASE( mat_trans, name ) \
	{ test_mat_transpose<sform, dform, M, N>(); } \
	BEGIN_TPACK( mat_trans_##name ) \
	ADD_MN_CASE_3X3( mat_trans, name, DM, DN ) \
	END_TPACK

#define TEST_XPR_TRANS( dform, name ) \
	MN_CASE( mat_trans, name ) \
	{ test_xpr_transpose<dform, M, N>(); } \
	BEGIN_TPACK( mat_trans_##name ) \
	ADD_MN_CASE_3X3( mat_trans, name, DM, DN ) \
	END_TPACK

TEST_MAT_TRANS( ref_matrix, ref_matrix, mat_to_mat )
TEST_MAT_TRANS( ref_matrix, ref_block, mat_to_blk )
TEST_MAT_TRANS( ref_matrix, ref_grid, mat_to_grid )

TEST_MAT_TRANS( ref_block, ref_matrix, blk_to_mat )
TEST_MAT_TRANS( ref_block, ref_block, blk_to_blk )
TEST_MAT_TRANS( ref_block, ref_grid, blk_to_grid )

TEST_MAT_TRANS( ref_grid, ref_matrix, grid_to_mat )
TEST_MAT_TRANS( ref_grid, ref_block, grid_to_blk )
TEST_MAT_TRANS( ref_grid, ref_grid, grid_to_grid )

TEST_XPR_TRANS( ref_matrix, xpr_to_mat )
TEST_XPR_TRANS( ref_block, xpr_to_blk )
TEST_XPR_TRANS( ref_grid, xpr_to_grid )

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_trans_mat_to_mat )
	ADD_TPACK( mat_trans_mat_to_blk )
	ADD_TPACK( mat_trans_mat_to_grid )

	ADD_TPACK( mat_trans_blk_to_mat )
	ADD_TPACK( mat_trans_blk_to_blk )
	ADD_TPACK( mat_trans_blk_to_grid )

	ADD_TPACK( mat_trans_grid_to_mat )
	ADD_TPACK( mat_trans_grid_to_blk )
	ADD_TPACK( mat_trans_grid_to_grid )

	ADD_TPACK( mat_trans_xpr_to_mat )
	ADD_TPACK( mat_trans_xpr_to_blk )
	ADD_TPACK( mat_trans_xpr_to_grid )
END_MAIN_SUITE

