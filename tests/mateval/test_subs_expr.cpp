/**
 * @file test_subs_expr.cpp
 *
 * Unit testing of subs_expr classes
 * 
 * @author Dahua Lin 
 */


#include "../test_base.h"
#include <light_mat/mateval/subs_expr.h>
#include <light_mat/mateval/mat_arith.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 9;
const index_t DN = 8;

TMN_CASE( test_subs, inds_expr )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	inds_expr<T, M, N> e = inds( type_<T>(), matrix_shape<M, N>(m, n) );

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n );

	dense_matrix<T> r = e;

	dense_matrix<T> r0(m, n);
	for (index_t i = 0; i < m * n; ++i) r0[i] = (T)(i);

	ASSERT_MAT_EQ( m, n, r, r0 );
}

TMN_CASE( test_subs, subs_i_expr )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	subs_i_expr<T, M, N> e = subs_i( type_<T>(), matrix_shape<M, N>(m, n) );

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n );

	dense_matrix<T> r = e;

	dense_matrix<T> r0(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) r0(i, j) = (T)(i);

	ASSERT_MAT_EQ( m, n, r, r0 );
}

TMN_CASE( test_subs, subs_j_expr )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	subs_j_expr<T, M, N> e = subs_j( type_<T>(), matrix_shape<M, N>(m, n) );

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n );

	dense_matrix<T> r = e;

	dense_matrix<T> r0(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) r0(i, j) = (T)(j);

	ASSERT_MAT_EQ( m, n, r, r0 );
}


TMN_CASE( test_subs, sub2ind )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	matrix_shape<M, N> shape(m, n);

	auto I = subs_i( type_<T>(), shape );
	auto J = subs_j( type_<T>(), shape );

	dense_matrix<T> r = I + J * T(m);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );

	dense_matrix<T> r0(m, n);
	for (index_t i = 0; i < m * n; ++i) r0[i] = (T)(i);

	ASSERT_MAT_EQ( m, n, r, r0 );
}


BEGIN_TPACK( test_inds )
	ADD_TMN_CASE_3X3( test_subs, inds_expr, double, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, inds_expr, float, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, inds_expr, int32_t, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, inds_expr, uint32_t, DM, DN )
END_TPACK

BEGIN_TPACK( test_subs_i )
	ADD_TMN_CASE_3X3( test_subs, subs_i_expr, double, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, subs_i_expr, float, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, subs_i_expr, int32_t, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, subs_i_expr, uint32_t, DM, DN )
END_TPACK

BEGIN_TPACK( test_subs_j )
	ADD_TMN_CASE_3X3( test_subs, subs_j_expr, double, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, subs_j_expr, float, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, subs_j_expr, int32_t, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, subs_j_expr, uint32_t, DM, DN )
END_TPACK

BEGIN_TPACK( test_sub2ind )
	ADD_TMN_CASE_3X3( test_subs, sub2ind, double, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, sub2ind, float, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, sub2ind, int32_t, DM, DN )
	ADD_TMN_CASE_3X3( test_subs, sub2ind, uint32_t, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( test_inds )
	ADD_TPACK( test_subs_i )
	ADD_TPACK( test_subs_j )
	ADD_TPACK( test_sub2ind )
END_MAIN_SUITE



