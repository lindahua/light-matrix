/**
 * @file test_partial_reduce.cpp
 *
 * Unit testing of partial reduction
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matrix/partial_reduce.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 6;
const index_t DN = 7;
const index_t LDim = 8;

MN_CASE( colwise_reduce, sum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += A(i, j);
		r0[j] = s;
	}

	// test

	row_t r = sum(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, sum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += A(i, j);
		r0[i] = s;
	}

	// test

	col_t r = sum(A, rowwise());

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}

/*
MN_CASE( colwise_reduce, mean )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += A(i, j);
		r0[j] = s / double(m);
	}

	// test

	row_t r = mean(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_APPROX( n, r, r0, 1.0e-12 );
}
*/


BEGIN_TPACK( colwise_sum )
	ADD_MN_CASE_3X3( colwise_reduce, sum, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_sum )
	ADD_MN_CASE_3X3( rowwise_reduce, sum, DM, DN )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( colwise_sum )
	ADD_TPACK( rowwise_sum )
END_MAIN_SUITE


