/**
 * @file test_repvecs.cpp
 *
 * Unit testing of repeat vectors
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include "../multimat_supp.h"

#include <light_mat/mateval/repvec_expr.h>
#include <light_mat/mateval/map_expr.h>

using namespace lmat;
using namespace lmat::test;


// Auxiliary facilities

template<class Tag1, class Tag2, int M, int N>
void test_repcols()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<Tag1, double, M, 1>::cmat_t scol_t;
	typedef typename mat_host<Tag2, double, M, N>::mat_t dmat_t;

	mat_host<Tag1, double, M, 1> src(m, 1);
	mat_host<Tag2, double, M, N> dst(m, n);

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

	ASSERT_EQ( repcol(scol, n).nrows(), m );
	ASSERT_EQ( repcol(scol, n).ncolumns(), n );

	dmat_t dmat = dst.get_mat();

	if (N == 0)
	{
		dmat = repcol(scol, n);
	}
	else
	{
		dmat = repcol(scol, dimension<N>());
	}

	ASSERT_MAT_EQ( m, n, dmat, rmat );

	zero(dmat);
	evaluate_by_map( repcol(scol, dimension<N>(n)), dmat );
	ASSERT_MAT_EQ( m, n, dmat, rmat );
}

template<class Tag1, class Tag2, int M, int N>
void test_reprows()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<Tag1, double, 1, N>::cmat_t srow_t;
	typedef typename mat_host<Tag2, double, M, N>::mat_t dmat_t;

	mat_host<Tag1, double, 1, N> src(1, n);
	mat_host<Tag2, double, M, N> dst(m, n);

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

	ASSERT_EQ( reprow(srow, m).nrows(), m );
	ASSERT_EQ( reprow(srow, m).ncolumns(), n );

	dmat_t dmat = dst.get_mat();
	if (M == 0)
	{
		dmat = reprow(srow, m);
	}
	else
	{
		dmat = reprow(srow, dimension<M>());
	}

	ASSERT_MAT_EQ( m, n, dmat, rmat );

	zero(dmat);
	evaluate_by_map( reprow(srow, dimension<M>(m)), dmat );
	ASSERT_MAT_EQ( m, n, dmat, rmat );
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


TEST_REPCOLS( cont, cont, mat_to_mat )
TEST_REPCOLS( cont, bloc, mat_to_blk )
TEST_REPCOLS( cont, grid, mat_to_grid )
TEST_REPCOLS( bloc, cont, blk_to_mat )
TEST_REPCOLS( bloc, bloc, blk_to_blk )
TEST_REPCOLS( bloc, grid, blk_to_grid )
TEST_REPCOLS( grid, cont, grid_to_mat )
TEST_REPCOLS( grid, bloc, grid_to_blk )
TEST_REPCOLS( grid, grid, grid_to_grid )

TEST_REPROWS( cont, cont, mat_to_mat )
TEST_REPROWS( cont, bloc, mat_to_blk )
TEST_REPROWS( cont, grid, mat_to_grid )
TEST_REPROWS( bloc, cont, blk_to_mat )
TEST_REPROWS( bloc, bloc, blk_to_blk )
TEST_REPROWS( bloc, grid, blk_to_grid )
TEST_REPROWS( grid, cont, grid_to_mat )
TEST_REPROWS( grid, bloc, grid_to_blk )
TEST_REPROWS( grid, grid, grid_to_grid )


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
END_MAIN_SUITE





