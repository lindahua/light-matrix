/**
 * @file test_rep_veval.cpp
 *
 * Test vector evaluation of repeated vectors
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/repeat_vecs_eval.h>

using namespace lmat;
using namespace lmat::test;

typedef dense_matrix<double, 0, 1> dcol_t;
typedef dense_matrix<double, 1, 0> drow_t;

#ifdef LMAT_USE_STATIC_ASSERT
static_assert(is_linear_vector_evaluator<single_vec_linear_evaluator<dcol_t, false>, double>::value,
		"Evaluator interface check failed");
static_assert(is_linear_vector_evaluator<rep_scalar_linear_evaluator<dcol_t, false>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<single_vec_percol_evaluator<dcol_t, false>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<rep_scalar_percol_evaluator<dcol_t, false>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<repcol_percol_evaluator<dcol_t, false>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<reprow_percol_evaluator<drow_t, false>, double>::value,
		"Evaluator interface check failed");
#endif


MN_CASE( repcols, eval )
{
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, 1> col_t;

	col_t col(m, 1);
	for (index_t i = 0; i < m; ++i) col[i] = double(i + 2);

	typedef repeat_col_expr<col_t, N, false> expr_t;
	expr_t expr(col, n);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );

	dense_matrix<double, M, N> R( expr );

	dense_matrix<double, M, N> R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = col(i, 0);
	}

	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( reprows, eval )
{
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;

	row_t row(1, n);
	for (index_t j = 0; j < n; ++j) row[j] = double(j + 2);

	typedef repeat_row_expr<row_t, M, false> expr_t;
	expr_t expr(row, m);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );

	dense_matrix<double, M, N> R( expr );

	dense_matrix<double, M, N> R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = row(0, j);
	}

	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( linear_veval, repcol_linear )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, 1> col_t;
	typedef dense_matrix<double, M, N> mat_t;

	col_t col(m, 1);
	for (index_t i = 0; i < m; ++i) col[i] = double(i+2);

	repeat_col_expr<col_t, N, false> expr(col, n);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = col(i, 0);
	}

	mat_t R(m, n);
	lmat::detail::ewise_linear_eval_internal<double, M*N>::evaluate(expr, R);

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

	repeat_col_expr<col_t, N, false> expr(col, n);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = col(i, 0);
	}

	mat_t R(m, n);
	lmat::detail::ewise_percol_eval_internal<double, M, N>::evaluate(expr, R);

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

	repeat_row_expr<row_t, M, false> expr(row, m);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = row(0, j);
	}

	mat_t R(m, n);
	lmat::detail::ewise_linear_eval_internal<double, M*N>::evaluate(expr, R);

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

	repeat_row_expr<row_t, M, false> expr(row, m);

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = row(0, j);
	}

	mat_t R(m, n);
	lmat::detail::ewise_percol_eval_internal<double, M, N>::evaluate(expr, R);

	ASSERT_TRUE( is_equal(R, R_r) );
}


BEGIN_TPACK( repcols_eval )
	ADD_MN_CASE( repcols, eval, 0, 0 )
	ADD_MN_CASE( repcols, eval, 0, 1 )
	ADD_MN_CASE( repcols, eval, 0, 6 )
	ADD_MN_CASE( repcols, eval, 1, 0 )
	ADD_MN_CASE( repcols, eval, 1, 1 )
	ADD_MN_CASE( repcols, eval, 1, 6 )
	ADD_MN_CASE( repcols, eval, 4, 0 )
	ADD_MN_CASE( repcols, eval, 4, 1 )
	ADD_MN_CASE( repcols, eval, 4, 6 )
END_TPACK


BEGIN_TPACK( reprows_eval )
	ADD_MN_CASE( reprows, eval, 0, 0 )
	ADD_MN_CASE( reprows, eval, 0, 1 )
	ADD_MN_CASE( reprows, eval, 0, 6 )
	ADD_MN_CASE( reprows, eval, 1, 0 )
	ADD_MN_CASE( reprows, eval, 1, 1 )
	ADD_MN_CASE( reprows, eval, 1, 6 )
	ADD_MN_CASE( reprows, eval, 4, 0 )
	ADD_MN_CASE( reprows, eval, 4, 1 )
	ADD_MN_CASE( reprows, eval, 4, 6 )
END_TPACK

BEGIN_TPACK( repcol_linear_eval )
	ADD_MN_CASE_3X3( linear_veval, repcol_linear, 5, 6 )
END_TPACK

BEGIN_TPACK( repcol_percol_eval )
	ADD_MN_CASE_3X3( percol_veval, repcol_percol, 5, 6 )
END_TPACK

BEGIN_TPACK( reprow_linear_eval )
	ADD_MN_CASE_3X3( linear_veval, reprow_linear, 5, 6 )
END_TPACK

BEGIN_TPACK( reprow_percol_eval )
	ADD_MN_CASE_3X3( percol_veval, reprow_percol, 5, 6 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( repcols_eval )
	ADD_TPACK( reprows_eval )
	ADD_TPACK( repcol_linear_eval )
	ADD_TPACK( repcol_percol_eval )
	ADD_TPACK( reprow_linear_eval )
	ADD_TPACK( reprow_percol_eval )
END_MAIN_SUITE

