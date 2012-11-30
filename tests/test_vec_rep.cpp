/**
 * @file test_rep_vecs.cpp
 *
 * Unit testing of repeat vectors
 *
 * @author Dahua Lin
 */


#include "test_base.h"
#include "multimat_supp.h"

#include <light_mat/matexpr/matrix_arith.h>
#include <light_mat/matexpr/vector_repeat.h>

using namespace lmat;
using namespace lmat::test;


// Auxiliary facilities

template<class Tag1, class Tag2, int M, int N>
void test_repcols()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_maker<Tag1, double, M, 1>::cmat_t scol_t;
	typedef typename mat_maker<Tag2, double, M, N>::mat_t dmat_t;

	mat_maker<Tag1, double, M, 1> src(m, 1);
	mat_maker<Tag2, double, M, N> dst(m, n);

	src.fill_lin();
	scol_t scol = src.get_cmat();

	dense_matrix<double, M, N> rmat(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rmat(i, j) = scol[i];
		}
	}

	ASSERT_EQ( rep_col(scol, n).nrows(), m );
	ASSERT_EQ( rep_col(scol, n).ncolumns(), n );

	dmat_t dmat = dst.get_mat();

	if (N == 0)
	{
		dmat = rep_col(scol, n);
	}
	else
	{
		dmat = rep_col(scol, fix_int<N>());
	}

	ASSERT_MAT_EQ( m, n, dmat, rmat );

	// access-based evaluation

	typedef macc_scheme<percol_macc, scalar_kernel_t, M, N> scheme_t;
	fill(dmat, 0.0);
	scheme_t sch = scheme_t::get_default(rep_col(scol, n), dmat);
	sch.evaluate(rep_col(scol, n), dmat);
	ASSERT_MAT_EQ(m, n, dmat, rmat );
}

template<class Tag1, class Tag2, int M, int N>
void test_reprows()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_maker<Tag1, double, 1, N>::cmat_t srow_t;
	typedef typename mat_maker<Tag2, double, M, N>::mat_t dmat_t;

	mat_maker<Tag1, double, 1, N> src(1, n);
	mat_maker<Tag2, double, M, N> dst(m, n);

	src.fill_lin();
	srow_t srow = src.get_cmat();

	dense_matrix<double, M, N> rmat(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rmat(i, j) = srow[j];
		}
	}

	ASSERT_EQ( rep_row(srow, m).nrows(), m );
	ASSERT_EQ( rep_row(srow, m).ncolumns(), n );

	dmat_t dmat = dst.get_mat();
	if (M == 0)
	{
		dmat = rep_row(srow, m);
	}
	else
	{
		dmat = rep_row(srow, fix_int<M>());
	}

	ASSERT_MAT_EQ( m, n, dmat, rmat );

	// access-based evaluation

	typedef macc_scheme<percol_macc, scalar_kernel_t, M, N> scheme_t;
	fill(dmat, 0.0);
	scheme_t sch = scheme_t::get_default(rep_row(srow, m), dmat);
	sch.evaluate(rep_row(srow, m), dmat);
	ASSERT_MAT_EQ(m, n, dmat, rmat );
}


template<class Tag2, int M, int N>
void test_repcols_ex()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_maker<Tag2, double, M, N>::mat_t dmat_t;

	dense_col<double, M> scol(m);
	for (index_t i = 0; i < m; ++i) scol[i] = double(i+1);

	dense_matrix<double, M, N> rmat(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rmat(i, j) = scol[i] * 2.0 + 1.0;
		}
	}

	ASSERT_EQ( rep_col(scol, n).nrows(), m );
	ASSERT_EQ( rep_col(scol, n).ncolumns(), n );

	mat_maker<Tag2, double, M, N> dst(m, n);
	dmat_t dmat = dst.get_mat();
	if (N == 0)
	{
		dmat = rep_col(scol * 2.0 + 1.0, n);
	}
	else
	{
		dmat = rep_col(scol * 2.0 + 1.0, fix_int<N>());
	}

	ASSERT_MAT_EQ( m, n, dmat, rmat );

	// access-based evaluation

	typedef macc_scheme<percol_macc, scalar_kernel_t, M, N> scheme_t;
	fill(dmat, 0.0);
	scheme_t sch = scheme_t::get_default(rep_col(scol * 2.0 + 1.0, n), dmat);
	sch.evaluate(rep_col(scol * 2.0 + 1.0, n), dmat);
	ASSERT_MAT_EQ(m, n, dmat, rmat );
}


template<class Tag2, int M, int N>
void test_reprows_ex()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_maker<Tag2, double, M, N>::mat_t dmat_t;

	dense_row<double, N> srow(n);
	for (index_t j = 0; j < n; ++j) srow[j] = double(j+1);

	dense_matrix<double, M, N> rmat(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			rmat(i, j) = srow[j] * 2.0 + 1.0;
		}
	}

	ASSERT_EQ( rep_row(srow, m).nrows(), m );
	ASSERT_EQ( rep_row(srow, m).ncolumns(), n );

	mat_maker<Tag2, double, M, N> dst(m, n);
	dmat_t dmat = dst.get_mat();
	if (M == 0)
	{
		dmat = rep_row(srow * 2.0 + 1.0, m);
	}
	else
	{
		dmat = rep_row(srow * 2.0 + 1.0, fix_int<M>());
	}

	ASSERT_MAT_EQ( m, n, dmat, rmat );

	// access-based evaluation

	typedef macc_scheme<percol_macc, scalar_kernel_t, M, N> scheme_t;
	fill(dmat, 0.0);
	scheme_t sch = scheme_t::get_default(rep_row(srow * 2.0 + 1.0, m), dmat);
	sch.evaluate(rep_row(srow * 2.0 + 1.0, m), dmat);
	ASSERT_MAT_EQ(m, n, dmat, rmat );
}



#define TEST_REPCOLS( sform, dform, name ) \
	MN_CASE( test_repcols, name ) { \
		test_repcols<sform, dform, M, N>(); } \
	BEGIN_TPACK( test_repcols_##name ) \
		ADD_MN_CASE_3X3( test_repcols, name, DM, DN ) \
	END_TPACK

#define TEST_REPROWS( sform, dform, name ) \
	MN_CASE( test_reprows, name ) { \
		test_reprows<sform, dform, M, N>(); } \
	BEGIN_TPACK( test_reprows_##name ) \
		ADD_MN_CASE_3X3( test_reprows, name, DM, DN ) \
	END_TPACK

#define TEST_REPCOLS_EX( dform, name ) \
	MN_CASE( test_repcols, name ) { \
		test_repcols_ex<dform, M, N>(); } \
	BEGIN_TPACK( test_repcols_##name ) \
		ADD_MN_CASE_3X3( test_repcols, name, DM, DN ) \
	END_TPACK

#define TEST_REPROWS_EX( dform, name ) \
	MN_CASE( test_reprows, name ) { \
		test_reprows_ex<dform, M, N>(); } \
	BEGIN_TPACK( test_reprows_##name ) \
		ADD_MN_CASE_3X3( test_reprows, name, DM, DN ) \
	END_TPACK


TEST_REPCOLS( ref_matrix, ref_matrix, mat_to_mat )
TEST_REPCOLS( ref_matrix, ref_block,  mat_to_blk )
TEST_REPCOLS( ref_matrix, ref_grid,   mat_to_grid )
TEST_REPCOLS( ref_block,  ref_matrix, blk_to_mat )
TEST_REPCOLS( ref_block,  ref_block,  blk_to_blk )
TEST_REPCOLS( ref_block,  ref_grid,   blk_to_grid )
TEST_REPCOLS( ref_grid,   ref_matrix, grid_to_mat )
TEST_REPCOLS( ref_grid,   ref_block,  grid_to_blk )
TEST_REPCOLS( ref_grid,   ref_grid,   grid_to_grid )

TEST_REPROWS( ref_matrix, ref_matrix, mat_to_mat )
TEST_REPROWS( ref_matrix, ref_block,  mat_to_blk )
TEST_REPROWS( ref_matrix, ref_grid,   mat_to_grid )
TEST_REPROWS( ref_block,  ref_matrix, blk_to_mat )
TEST_REPROWS( ref_block,  ref_block,  blk_to_blk )
TEST_REPROWS( ref_block,  ref_grid,   blk_to_grid )
TEST_REPROWS( ref_grid,   ref_matrix, grid_to_mat )
TEST_REPROWS( ref_grid,   ref_block,  grid_to_blk )
TEST_REPROWS( ref_grid,   ref_grid,   grid_to_grid )

TEST_REPCOLS_EX( ref_matrix, xpr_to_mat )
TEST_REPCOLS_EX( ref_block,  xpr_to_blk )
TEST_REPCOLS_EX( ref_grid,   xpr_to_grid )

TEST_REPROWS_EX( ref_matrix, xpr_to_mat )
TEST_REPROWS_EX( ref_block,  xpr_to_blk )
TEST_REPROWS_EX( ref_grid,   xpr_to_grid )


BEGIN_MAIN_SUITE
	ADD_TPACK( test_repcols_mat_to_mat )
	ADD_TPACK( test_repcols_mat_to_blk )
	ADD_TPACK( test_repcols_mat_to_grid )
	ADD_TPACK( test_repcols_blk_to_mat )
	ADD_TPACK( test_repcols_blk_to_blk )
	ADD_TPACK( test_repcols_blk_to_grid )
	ADD_TPACK( test_repcols_grid_to_mat )
	ADD_TPACK( test_repcols_grid_to_blk )
	ADD_TPACK( test_repcols_grid_to_grid )

	ADD_TPACK( test_reprows_mat_to_mat )
	ADD_TPACK( test_reprows_mat_to_blk )
	ADD_TPACK( test_reprows_mat_to_grid )
	ADD_TPACK( test_reprows_blk_to_mat )
	ADD_TPACK( test_reprows_blk_to_blk )
	ADD_TPACK( test_reprows_blk_to_grid )
	ADD_TPACK( test_reprows_grid_to_mat )
	ADD_TPACK( test_reprows_grid_to_blk )
	ADD_TPACK( test_reprows_grid_to_grid )

	ADD_TPACK( test_repcols_xpr_to_mat )
	ADD_TPACK( test_repcols_xpr_to_blk )
	ADD_TPACK( test_repcols_xpr_to_grid )

	ADD_TPACK( test_reprows_xpr_to_mat )
	ADD_TPACK( test_reprows_xpr_to_blk )
	ADD_TPACK( test_reprows_xpr_to_grid )
END_MAIN_SUITE





