/**
 * @file test_rep_veval.cpp
 *
 * Test vector evaluation of repeated vectors
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/repeat_vecs_veval.h>

using namespace lmat;
using namespace lmat::test;

typedef dense_matrix<double, 0, 1> dcol_t;
typedef dense_matrix<double, 1, 0> drow_t;

#ifdef LMAT_USE_STATIC_ASSERT
static_assert(is_linear_vector_evaluator<single_vec_linear_evaluator<dcol_t>, double>::value,
		"Evaluator interface check failed");
static_assert(is_linear_vector_evaluator<rep_scalar_linear_evaluator<dcol_t>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<single_vec_percol_evaluator<dcol_t>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<rep_scalar_percol_evaluator<dcol_t>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<repcol_percol_evaluator<dcol_t>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<reprow_percol_evaluator<drow_t>, double>::value,
		"Evaluator interface check failed");
#endif


MN_CASE( linear_veval, repcol_linear )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, 1> col_t;
	typedef dense_matrix<double, M, N> mat_t;

	col_t col(m, 1);
	for (index_t i = 0; i < m; ++i) col[i] = double(i+2);

	repeat_col_expr<col_t, N> expr(col, n);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = col(i, 0);
	}

	mat_t R(m, n);
	lmat::detail::linear_mateval_internal<double, M*N>::evaluate(expr, R);

	ASSERT_TRUE( is_equal(R, R_r) );
}

MN_CASE( percol_veval, repcol_percol )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, 1> col_t;
	typedef dense_matrix<double, M, N> mat_t;

	col_t col(m, 1);
	for (index_t i = 0; i < m; ++i) col[i] = double(i+2);

	repeat_col_expr<col_t, N> expr(col, n);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = col(i, 0);
	}

	mat_t R(m, n);
	lmat::detail::percol_mateval_internal<double, M, N>::evaluate(expr, R);

	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( linear_veval, reprow_linear )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;
	typedef dense_matrix<double, M, N> mat_t;

	row_t row(1, n);
	for (index_t j = 0; j < n; ++j) row[j] = double(j+2);

	repeat_row_expr<row_t, M> expr(row, m);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = row(0, j);
	}

	mat_t R(m, n);
	lmat::detail::linear_mateval_internal<double, M*N>::evaluate(expr, R);

	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( percol_veval, reprow_percol )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;
	typedef dense_matrix<double, M, N> mat_t;

	row_t row(1, n);
	for (index_t j = 0; j < n; ++j) row[j] = double(j+2);

	repeat_row_expr<row_t, M> expr(row, m);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = row(0, j);
	}

	mat_t R(m, n);
	lmat::detail::percol_mateval_internal<double, M, N>::evaluate(expr, R);

	ASSERT_TRUE( is_equal(R, R_r) );
}




BEGIN_TPACK( repcol_linear_eval )
	ADD_MN_CASE( linear_veval, repcol_linear, 0, 0 )
	ADD_MN_CASE( linear_veval, repcol_linear, 0, 1 )
	ADD_MN_CASE( linear_veval, repcol_linear, 0, 6 )
	ADD_MN_CASE( linear_veval, repcol_linear, 1, 0 )
	ADD_MN_CASE( linear_veval, repcol_linear, 1, 1 )
	ADD_MN_CASE( linear_veval, repcol_linear, 1, 6 )
	ADD_MN_CASE( linear_veval, repcol_linear, 5, 0 )
	ADD_MN_CASE( linear_veval, repcol_linear, 5, 1 )
	ADD_MN_CASE( linear_veval, repcol_linear, 5, 6 )
END_TPACK

BEGIN_TPACK( repcol_percol_eval )
	ADD_MN_CASE( percol_veval, repcol_percol, 0, 0 )
	ADD_MN_CASE( percol_veval, repcol_percol, 0, 1 )
	ADD_MN_CASE( percol_veval, repcol_percol, 0, 6 )
	ADD_MN_CASE( percol_veval, repcol_percol, 1, 0 )
	ADD_MN_CASE( percol_veval, repcol_percol, 1, 1 )
	ADD_MN_CASE( percol_veval, repcol_percol, 1, 6 )
	ADD_MN_CASE( percol_veval, repcol_percol, 5, 0 )
	ADD_MN_CASE( percol_veval, repcol_percol, 5, 1 )
	ADD_MN_CASE( percol_veval, repcol_percol, 5, 6 )
END_TPACK

BEGIN_TPACK( reprow_linear_eval )
	ADD_MN_CASE( linear_veval, reprow_linear, 0, 0 )
	ADD_MN_CASE( linear_veval, reprow_linear, 0, 1 )
	ADD_MN_CASE( linear_veval, reprow_linear, 0, 6 )
	ADD_MN_CASE( linear_veval, reprow_linear, 1, 0 )
	ADD_MN_CASE( linear_veval, reprow_linear, 1, 1 )
	ADD_MN_CASE( linear_veval, reprow_linear, 1, 6 )
	ADD_MN_CASE( linear_veval, reprow_linear, 5, 0 )
	ADD_MN_CASE( linear_veval, reprow_linear, 5, 1 )
	ADD_MN_CASE( linear_veval, reprow_linear, 5, 6 )
END_TPACK

BEGIN_TPACK( reprow_percol_eval )
	ADD_MN_CASE( percol_veval, reprow_percol, 0, 0 )
	ADD_MN_CASE( percol_veval, reprow_percol, 0, 1 )
	ADD_MN_CASE( percol_veval, reprow_percol, 0, 6 )
	ADD_MN_CASE( percol_veval, reprow_percol, 1, 0 )
	ADD_MN_CASE( percol_veval, reprow_percol, 1, 1 )
	ADD_MN_CASE( percol_veval, reprow_percol, 1, 6 )
	ADD_MN_CASE( percol_veval, reprow_percol, 5, 0 )
	ADD_MN_CASE( percol_veval, reprow_percol, 5, 1 )
	ADD_MN_CASE( percol_veval, reprow_percol, 5, 6 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( repcol_linear_eval )
	ADD_TPACK( repcol_percol_eval )
	ADD_TPACK( reprow_linear_eval )
	ADD_TPACK( reprow_percol_eval )
END_MAIN_SUITE

