/**
 * @file test_dense_eval.cpp
 *
 * Unit test of evaluation of dense matrices
 *
 * @author Dahua Lin
 */


#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>

using namespace lmat;
using namespace lmat::test;

void fill_lin(dblock<double>& arr)
{
	for (index_t i = 0; i < arr.nelems(); ++i)
		arr[i] = double(i + 1);
}

MN_CASE( dense_mat_eval )
{
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 5 : N;

	dblock<double> s(m * n);
	fill_lin(s);

	dense_matrix<double, M, N> a(m, n, copy_from(s.ptr_data()));
	dense_matrix<double, M, N> r = eval(a);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_NE( r.ptr_data(), a.ptr_data() );

	ASSERT_MAT_EQ(m, n, a, r);
}

MN_CASE( ref_mat_eval )
{
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 5 : N;

	dblock<double> s(m * n);
	fill_lin(s);

	ref_matrix<double, M, N> a(s.ptr_data(), m, n);
	dense_matrix<double, M, N> r = eval(a);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_NE( r.ptr_data(), a.ptr_data() );

	ASSERT_MAT_EQ(m, n, a, r);
}

MN_CASE( ref_block_eval )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 5 : N;

	dblock<double> s(ldim * n);
	fill_lin(s);

	ref_block<double, M, N> a(s.ptr_data(), m, n, ldim);
	dense_matrix<double, M, N> r = eval(a);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_NE( r.ptr_data(), a.ptr_data() );

	ASSERT_MAT_EQ(m, n, a, r);
}

MN_CASE( ref_grid_eval )
{
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 5 : N;
	const index_t rs = 2;
	const index_t cs = 9;

	dblock<double> s(cs * n);
	fill_lin(s);

	ref_grid<double, M, N> a(s.ptr_data(), m, n, rs, cs);
	dense_matrix<double, M, N> r = eval(a);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_NE( r.ptr_data(), a.ptr_data() );

	ASSERT_MAT_EQ(m, n, a, r);
}


LTEST_INIT_AUTOSUITE

AUTO_TPACK( dense_mat_eval )
{
	ADD_MN_CASE_3X3( dense_mat_eval, 4, 5 );
}

AUTO_TPACK( ref_mat_eval )
{
	ADD_MN_CASE_3X3( ref_mat_eval, 4, 5 );
}

AUTO_TPACK( ref_block_eval )
{
	ADD_MN_CASE_3X3( ref_block_eval, 4, 5 );
}

AUTO_TPACK( ref_grid_eval )
{
	ADD_MN_CASE_3X3( ref_grid_eval, 4, 5 );
}



