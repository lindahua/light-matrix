/*
 * @file test_mat_elogical.cpp
 *
 * Unit testing of element-wise logical operations
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matexpr/matrix_elogical.h>

using namespace lmat;
using namespace lmat::test;

const int DM = 8;
const int DN = 6;
const index_t LDim = 12;

MN_CASE( mat_elogical, not )
{
	typedef dense_matrix<bool, M, N> bmat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	bmat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);

	bmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = ~(A[i]);

	bmat_t R = ~A;
	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( mat_elogical, and )
{
	typedef dense_matrix<bool, M, N> bmat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	bmat_t A(m, n);
	bmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	bmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] & B[i]);

	bmat_t R = (A & B);
	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( mat_elogical, or )
{
	typedef dense_matrix<bool, M, N> bmat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	bmat_t A(m, n);
	bmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	bmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] | B[i]);

	bmat_t R = (A | B);
	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( mat_elogical, xor )
{
	typedef dense_matrix<bool, M, N> bmat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	bmat_t A(m, n);
	bmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	bmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] ^ B[i]);

	bmat_t R = (A ^ B);
	ASSERT_TRUE( is_equal(R, R_r) );
}



BEGIN_TPACK( mat_elogical_not )
	ADD_MN_CASE_3X3( mat_elogical, not, DM, DN )
END_TPACK

BEGIN_TPACK( mat_elogical_and )
	ADD_MN_CASE_3X3( mat_elogical, and, DM, DN )
END_TPACK

BEGIN_TPACK( mat_elogical_or )
	ADD_MN_CASE_3X3( mat_elogical, or, DM, DN )
END_TPACK

BEGIN_TPACK( mat_elogical_xor )
	ADD_MN_CASE_3X3( mat_elogical, xor, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_elogical_not )
	ADD_TPACK( mat_elogical_and )
	ADD_TPACK( mat_elogical_or )
	ADD_TPACK( mat_elogical_xor )
END_MAIN_SUITE



