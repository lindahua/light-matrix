/*
 * @file test_mat_ebool.cpp
 *
 * Unit testing of element-wise boolean operations
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matrix/matrix_boolop.h>

using namespace lmat;
using namespace lmat::test;

const int default_m = 8;
const int default_n = 6;
const index_t LDim = 12;

MN_CASE( mat_ebool, not )
{
	typedef dense_matrix<bool, M, N> bmat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	bmat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef unary_ewise_expr<not_op, bmat_t> R_t;
	static_assert(is_same<decltype(!A), R_t>::value, "Expression type verification failed.");
#endif

	bmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = !(A[i]);

	bmat_t R = !A;
	ASSERT_TRUE( is_equal(R, R_r) );

	bmat_t R_s(m, n);
	evaluate_by_scalars(!A, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_ebool, and )
{
	typedef dense_matrix<bool, M, N> bmat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	bmat_t A(m, n);
	bmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef binary_ewise_expr<and_op, bmat_t, bmat_t> R_t;
	static_assert(is_same<decltype(A && B), R_t>::value, "Expression type verification failed.");
#endif

	bmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] && B[i]);

	bmat_t R = (A && B);
	ASSERT_TRUE( is_equal(R, R_r) );

	bmat_t R_s(m, n);
	evaluate_by_scalars(A && B, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_ebool, or )
{
	typedef dense_matrix<bool, M, N> bmat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	bmat_t A(m, n);
	bmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef binary_ewise_expr<or_op, bmat_t, bmat_t> R_t;
	static_assert(is_same<decltype(A || B), R_t>::value, "Expression type verification failed.");
#endif

	bmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] || B[i]);

	bmat_t R = (A || B);
	ASSERT_TRUE( is_equal(R, R_r) );

	bmat_t R_s(m, n);
	evaluate_by_scalars(A || B, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}



BEGIN_TPACK( mat_ebool_not )
	ADD_MN_CASE( mat_ebool, not, 0, 0 )
	ADD_MN_CASE( mat_ebool, not, 0, 1 )
	ADD_MN_CASE( mat_ebool, not, 0, default_n )
	ADD_MN_CASE( mat_ebool, not, 1, 0 )
	ADD_MN_CASE( mat_ebool, not, 1, 1 )
	ADD_MN_CASE( mat_ebool, not, 1, default_n )
	ADD_MN_CASE( mat_ebool, not, default_m, 0 )
	ADD_MN_CASE( mat_ebool, not, default_m, 1 )
	ADD_MN_CASE( mat_ebool, not, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_ebool_and )
	ADD_MN_CASE( mat_ebool, and, 0, 0 )
	ADD_MN_CASE( mat_ebool, and, 0, 1 )
	ADD_MN_CASE( mat_ebool, and, 0, default_n )
	ADD_MN_CASE( mat_ebool, and, 1, 0 )
	ADD_MN_CASE( mat_ebool, and, 1, 1 )
	ADD_MN_CASE( mat_ebool, and, 1, default_n )
	ADD_MN_CASE( mat_ebool, and, default_m, 0 )
	ADD_MN_CASE( mat_ebool, and, default_m, 1 )
	ADD_MN_CASE( mat_ebool, and, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_ebool_or )
	ADD_MN_CASE( mat_ebool, or, 0, 0 )
	ADD_MN_CASE( mat_ebool, or, 0, 1 )
	ADD_MN_CASE( mat_ebool, or, 0, default_n )
	ADD_MN_CASE( mat_ebool, or, 1, 0 )
	ADD_MN_CASE( mat_ebool, or, 1, 1 )
	ADD_MN_CASE( mat_ebool, or, 1, default_n )
	ADD_MN_CASE( mat_ebool, or, default_m, 0 )
	ADD_MN_CASE( mat_ebool, or, default_m, 1 )
	ADD_MN_CASE( mat_ebool, or, default_m, default_n )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_ebool_not )
	ADD_TPACK( mat_ebool_and )
	ADD_TPACK( mat_ebool_or )
END_MAIN_SUITE



