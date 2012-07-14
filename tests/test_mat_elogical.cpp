/*
 * @file test_mat_elogical.cpp
 *
 * Unit testing of element-wise logical operations
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matrix/matrix_elogical.h>

using namespace lmat;
using namespace lmat::test;

const int default_m = 8;
const int default_n = 6;
const index_t LDim = 12;

typedef mask_t<double> dmask_t;

MN_CASE( mat_elogical, not )
{
	typedef dense_matrix<dmask_t, M, N> mmat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mmat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);

	mmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = ~(A[i]);

	mmat_t R = ~A;
	ASSERT_TRUE( is_equal(R, R_r) );

	mmat_t R_s(m, n);
	evaluate_by_scalars(~A, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_elogical, and )
{
	typedef dense_matrix<dmask_t, M, N> mmat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mmat_t A(m, n);
	mmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	mmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] & B[i]);

	mmat_t R = (A & B);
	ASSERT_TRUE( is_equal(R, R_r) );

	mmat_t R_s(m, n);
	evaluate_by_scalars(A & B, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_elogical, or )
{
	typedef dense_matrix<dmask_t, M, N> mmat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mmat_t A(m, n);
	mmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	mmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] | B[i]);

	mmat_t R = (A | B);
	ASSERT_TRUE( is_equal(R, R_r) );

	mmat_t R_s(m, n);
	evaluate_by_scalars(A | B, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_elogical, xor )
{
	typedef dense_matrix<dmask_t, M, N> mmat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mmat_t A(m, n);
	mmat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2);
	for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3);

	mmat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] ^ B[i]);

	mmat_t R = (A ^ B);
	ASSERT_TRUE( is_equal(R, R_r) );

	mmat_t R_s(m, n);
	evaluate_by_scalars(A ^ B, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}



BEGIN_TPACK( mat_elogical_not )
	ADD_MN_CASE_3X3( mat_elogical, not, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_elogical_and )
	ADD_MN_CASE_3X3( mat_elogical, and, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_elogical_or )
	ADD_MN_CASE_3X3( mat_elogical, or, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_elogical_xor )
	ADD_MN_CASE_3X3( mat_elogical, xor, default_m, default_n )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_elogical_not )
	ADD_TPACK( mat_elogical_and )
	ADD_TPACK( mat_elogical_or )
	ADD_TPACK( mat_elogical_xor )
END_MAIN_SUITE



