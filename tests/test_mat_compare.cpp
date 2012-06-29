/**
 * @file test_mat_compare.cpp
 *
 * Unit testing for matrix comparison
 *
 * @author Dahua Lin
 */


#include "test_base.h"

#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/matrix/ref_matrix_ex.h>
#include <light_mat/core/array.h>

using namespace lmat;
using namespace lmat::test;

MN_CASE( mat_equal, equal )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	darray<double> sa(m * n);
	for (index_t i = 0; i < m * n; ++i) sa[i] = double(i + 2);

	darray<double> sb(sa);

	ref_matrix<double, M, N> a(sa.ptr_begin(), m, n);
	ref_matrix<double, M, N> b(sb.ptr_begin(), m, n);

	ASSERT_TRUE( is_equal(a, b) );

	sb[m * n - 1] = double(-1);

	ASSERT_FALSE( is_equal(a, b) );
}


MN_CASE( mat_equal, equal_ex )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	const index_t ldim_a = 7;
	const index_t ldim_b = 8;

	darray<double> sa(ldim_a * n);  fill(sa, 0.0);
	darray<double> sb(ldim_b * n);  fill(sb, 0.0);

	ref_matrix_ex<double, M, N> a(sa.ptr_begin(), m, n, ldim_a);
	ref_matrix_ex<double, M, N> b(sb.ptr_begin(), m, n, ldim_b);

	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) a(i, j) = b(i, j) = double(1 + i + j * m);

	ASSERT_TRUE( is_equal(a, b) );

	b(m-1, n-1) = double(-1);

	ASSERT_FALSE( is_equal(a, b) );
}


MN_CASE( mat_approx, approx )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	darray<double> sa(m * n);
	for (index_t i = 0; i < m * n; ++i) sa[i] = double(i + 2);

	darray<double> sb(sa);

	ref_matrix<double, M, N> a(sa.ptr_begin(), m, n);
	ref_matrix<double, M, N> b(sb.ptr_begin(), m, n);

	ASSERT_TRUE( is_approx(a, b, 0.1) );

	sb[m * n - 1] = double(-1);

	ASSERT_FALSE( is_approx(a, b, 0.1) );
	ASSERT_TRUE( is_approx(a, b, 1000.0) );
}


MN_CASE( mat_approx, approx_ex )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	const index_t ldim_a = 7;
	const index_t ldim_b = 8;

	darray<double> sa(ldim_a * n);  fill(sa, 0.0);
	darray<double> sb(ldim_b * n);  fill(sb, 0.0);

	ref_matrix_ex<double, M, N> a(sa.ptr_begin(), m, n, ldim_a);
	ref_matrix_ex<double, M, N> b(sb.ptr_begin(), m, n, ldim_b);

	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) a(i, j) = b(i, j) = double(1 + i + j * m);

	ASSERT_TRUE( is_approx(a, b, 0.1) );

	b(m-1, n-1) = double(-1);

	ASSERT_FALSE( is_approx(a, b, 0.1) );
	ASSERT_TRUE( is_approx(a, b, 1000.0) );
}



BEGIN_TPACK( mat_equal )
	ADD_MN_CASE( mat_equal, equal, 0, 0 );
	ADD_MN_CASE( mat_equal, equal, 0, 1 );
	ADD_MN_CASE( mat_equal, equal, 0, 6 );
	ADD_MN_CASE( mat_equal, equal, 1, 0 );
	ADD_MN_CASE( mat_equal, equal, 1, 1 );
	ADD_MN_CASE( mat_equal, equal, 1, 6 );
	ADD_MN_CASE( mat_equal, equal, 5, 0 );
	ADD_MN_CASE( mat_equal, equal, 5, 1 );
	ADD_MN_CASE( mat_equal, equal, 5, 6 );
END_TPACK

BEGIN_TPACK( mat_equal_ex )
	ADD_MN_CASE( mat_equal, equal_ex, 0, 0 );
	ADD_MN_CASE( mat_equal, equal_ex, 0, 1 );
	ADD_MN_CASE( mat_equal, equal_ex, 0, 6 );
	ADD_MN_CASE( mat_equal, equal_ex, 1, 0 );
	ADD_MN_CASE( mat_equal, equal_ex, 1, 1 );
	ADD_MN_CASE( mat_equal, equal_ex, 1, 6 );
	ADD_MN_CASE( mat_equal, equal_ex, 5, 0 );
	ADD_MN_CASE( mat_equal, equal_ex, 5, 1 );
	ADD_MN_CASE( mat_equal, equal_ex, 5, 6 );
END_TPACK

BEGIN_TPACK( mat_approx )
	ADD_MN_CASE( mat_approx, approx, 0, 0 );
	ADD_MN_CASE( mat_approx, approx, 0, 1 );
	ADD_MN_CASE( mat_approx, approx, 0, 6 );
	ADD_MN_CASE( mat_approx, approx, 1, 0 );
	ADD_MN_CASE( mat_approx, approx, 1, 1 );
	ADD_MN_CASE( mat_approx, approx, 1, 6 );
	ADD_MN_CASE( mat_approx, approx, 5, 0 );
	ADD_MN_CASE( mat_approx, approx, 5, 1 );
	ADD_MN_CASE( mat_approx, approx, 5, 6 );
END_TPACK

BEGIN_TPACK( mat_approx_ex )
	ADD_MN_CASE( mat_approx, approx_ex, 0, 0 );
	ADD_MN_CASE( mat_approx, approx_ex, 0, 1 );
	ADD_MN_CASE( mat_approx, approx_ex, 0, 6 );
	ADD_MN_CASE( mat_approx, approx_ex, 1, 0 );
	ADD_MN_CASE( mat_approx, approx_ex, 1, 1 );
	ADD_MN_CASE( mat_approx, approx_ex, 1, 6 );
	ADD_MN_CASE( mat_approx, approx_ex, 5, 0 );
	ADD_MN_CASE( mat_approx, approx_ex, 5, 1 );
	ADD_MN_CASE( mat_approx, approx_ex, 5, 6 );
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_equal )
	ADD_TPACK( mat_equal_ex )
	ADD_TPACK( mat_approx )
	ADD_TPACK( mat_approx_ex )
END_MAIN_SUITE



