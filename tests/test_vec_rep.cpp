/**
 * @file test_rep_vecs.cpp
 *
 * Unit testing of repeat vectors
 *
 * @author Dahua Lin
 */


#include "test_base.h"

#include <light_mat/matrix/const_matrix.h>
#include <light_mat/matrix/dense_matrix.h>
#include <light_mat/matrix/matrix_repeat_eval.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

typedef dense_matrix<double, 0, 1> dcol_t;
typedef dense_matrix<double, 1, 0> drow_t;

template class lmat::horizontal_repeat_expr<ref_arg_t, dcol_t, 0>;
template class lmat::horizontal_repeat_expr<ref_arg_t, dcol_t, 6>;
template class lmat::vertical_repeat_expr<ref_arg_t, drow_t, 0>;
template class lmat::vertical_repeat_expr<ref_arg_t, drow_t, 4>;


N_CASE( repcols, generic )
{
	const index_t m = N == 0 ? 4 : N;
	const index_t n = 6;

	typedef dense_matrix<double, N, 1> col_t;
	typedef horizontal_repeat_expr<ref_arg_t, col_t, 0> expr_t;

	col_t a(m, 1);

	expr_t expr = hrep(a, n);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( &(expr.arg()), &a );
}


MN_CASE( repcols, generic_fix )
{
	check_arg(N > 0, "N must be positive");

	const index_t m = M == 0 ? 4 : M;
	const index_t n = N;

	typedef dense_matrix<double, M, 1> col_t;
	typedef horizontal_repeat_expr<ref_arg_t, col_t, N> expr_t;

	col_t a(m, 1);

	expr_t expr = hrep(a, fix_int<N>());

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( &(expr.arg()), &a );
}


N_CASE( repcols, const )
{
	const index_t m = N == 0 ? 4 : N;
	const index_t n = 6;
	const double val = 12.5;

	typedef const_matrix<double, N, 1> col_t;
	typedef const_matrix<double, N, 0> expr_t;

	col_t a(m, 1, val);

	expr_t expr = hrep(a, n);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( expr.value(), val );
}


MN_CASE( repcols, const_fix )
{
	check_arg(N > 0, "N must be positive");

	const index_t m = M == 0 ? 4 : M;
	const index_t n = N;
	const double val = 12.5;

	typedef const_matrix<double, M, 1> col_t;
	typedef const_matrix<double, M, N> expr_t;

	col_t a(m, 1, val);

	expr_t expr = hrep(a, fix_int<N>());

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( expr.value(), val );
}


N_CASE( reprows, generic )
{
	const index_t m = 4;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;
	typedef vertical_repeat_expr<ref_arg_t, row_t, 0> expr_t;

	row_t a(1, n);

	expr_t expr = vrep(a, m);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( &(expr.arg()), &a );
}


MN_CASE( reprows, generic_fix )
{
	check_arg(M > 0, "M must be positive");

	const index_t m = M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;
	typedef vertical_repeat_expr<ref_arg_t, row_t, M> expr_t;

	row_t a(1, n);

	expr_t expr = vrep(a, fix_int<M>());

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( &(expr.arg()), &a );
}


N_CASE( reprows, const )
{
	const index_t m = 4;
	const index_t n = N == 0 ? 6 : N;
	const double val = 12.5;

	typedef const_matrix<double, 1, N> row_t;
	typedef const_matrix<double, 0, N> expr_t;

	row_t a(1, n, val);
	expr_t expr = vrep(a, m);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( expr.value(), val );
}

MN_CASE( reprows, const_fix )
{
	check_arg(M > 0, "M must be positive");

	const index_t m = M;
	const index_t n = N == 0 ? 6 : N;
	const double val = 12.5;

	typedef const_matrix<double, 1, N> row_t;
	typedef const_matrix<double, M, N> expr_t;

	row_t a(1, n, val);
	expr_t expr = vrep(a, fix_int<M>());

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( expr.value(), val );
}


MN_CASE( repcols, eval )
{
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, 1> col_t;

	col_t col(m, 1);
	for (index_t i = 0; i < m; ++i) col[i] = double(i + 2);

	typedef horizontal_repeat_expr<ref_arg_t, col_t, N> expr_t;
	expr_t expr = make_expr(hrep_t<N>(n), ref_arg(col));

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );

	dense_matrix<double, M, N> R( expr );

	typedef horizontal_repeat_expr<copy_arg_t, col_t, N> expr_et;
	expr_et expr_e = make_expr(hrep_t<N>(n), copy_arg(col));

	ASSERT_EQ( expr_e.nrows(), m );
	ASSERT_EQ( expr_e.ncolumns(), n );

	dense_matrix<double, M, N> Re(expr_e);

	dense_matrix<double, M, N> R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = col(i, 0);
	}

	ASSERT_TRUE( is_equal(R, R_r) );
	ASSERT_TRUE( is_equal(Re, R_r) );
}


MN_CASE( reprows, eval )
{
	const index_t m = M == 0 ? 4 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;

	row_t row(1, n);
	for (index_t j = 0; j < n; ++j) row[j] = double(j + 2);

	typedef vertical_repeat_expr<ref_arg_t, row_t, M> expr_t;
	expr_t expr = make_expr(vrep_t<M>(m), ref_arg(row));

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );

	dense_matrix<double, M, N> R( expr );

	typedef vertical_repeat_expr<copy_arg_t, row_t, M> expr_et;
	expr_et expr_e = make_expr(vrep_t<M>(m), copy_arg(row));

	ASSERT_EQ( expr_e.nrows(), m );
	ASSERT_EQ( expr_e.ncolumns(), n );

	dense_matrix<double, M, N> Re(expr_e);

	dense_matrix<double, M, N> R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = row(0, j);
	}

	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( repcols, vec_eval )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, 1> col_t;
	typedef dense_matrix<double, M, N> mat_t;

	col_t col(m, 1);
	for (index_t i = 0; i < m; ++i) col[i] = double(i+2);

	horizontal_repeat_expr<ref_arg_t, col_t, N> expr =
			make_expr(hrep_t<N>(n), ref_arg(col));

	horizontal_repeat_expr<copy_arg_t, col_t, N> expr_e =
			make_expr(hrep_t<N>(n), copy_arg(col));

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = col(i, 0);
	}

	mat_t R(m, n);

	fill(R, 0.0);
	linear_by_scalars_evaluate(expr, R);
	ASSERT_TRUE( is_equal(R, R_r) );

	fill(R, 0.0);
	linear_by_scalars_evaluate(expr_e, R);
	ASSERT_TRUE( is_equal(R, R_r) );

	fill(R, 0.0);
	percol_by_scalars_evaluate(expr, R);
	ASSERT_TRUE( is_equal(R, R_r) );

	fill(R, 0.0);
	percol_by_scalars_evaluate(expr_e, R);
	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( reprows, vec_eval )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;
	typedef dense_matrix<double, M, N> mat_t;

	row_t row(1, n);
	for (index_t j = 0; j < n; ++j) row[j] = double(j+2);

	vertical_repeat_expr<ref_arg_t, row_t, M> expr =
			make_expr(vrep_t<M>(m), ref_arg(row));

	vertical_repeat_expr<copy_arg_t, row_t, M> expr_e =
			make_expr(vrep_t<M>(m), copy_arg(row));

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) R_r(i, j) = row(0, j);
	}

	mat_t R(m, n);

	fill(R, 0.0);
	linear_by_scalars_evaluate(expr, R);
	ASSERT_TRUE( is_equal(R, R_r) );

	fill(R, 0.0);
	linear_by_scalars_evaluate(expr_e, R);
	ASSERT_TRUE( is_equal(R, R_r) );

	fill(R, 0.0);
	percol_by_scalars_evaluate(expr, R);
	ASSERT_TRUE( is_equal(R, R_r) );

	fill(R, 0.0);
	percol_by_scalars_evaluate(expr_e, R);
	ASSERT_TRUE( is_equal(R, R_r) );
}


BEGIN_TPACK( repcols_generic )
	ADD_N_CASE( repcols, generic, 0 )
	ADD_N_CASE( repcols, generic, 1 )
	ADD_N_CASE( repcols, generic, 4 )
END_TPACK

BEGIN_TPACK( repcols_generic_fix )
	ADD_MN_CASE( repcols, generic_fix, 0, 1 )
	ADD_MN_CASE( repcols, generic_fix, 1, 1 )
	ADD_MN_CASE( repcols, generic_fix, 4, 1 )
	ADD_MN_CASE( repcols, generic_fix, 0, 6 )
	ADD_MN_CASE( repcols, generic_fix, 1, 6 )
	ADD_MN_CASE( repcols, generic_fix, 4, 6 )
END_TPACK

BEGIN_TPACK( repcols_const )
	ADD_N_CASE( repcols, const, 0 )
	ADD_N_CASE( repcols, const, 1 )
	ADD_N_CASE( repcols, const, 4 )
END_TPACK

BEGIN_TPACK( repcols_const_fix )
	ADD_MN_CASE( repcols, const_fix, 0, 1 )
	ADD_MN_CASE( repcols, const_fix, 1, 1 )
	ADD_MN_CASE( repcols, const_fix, 4, 1 )
	ADD_MN_CASE( repcols, const_fix, 0, 6 )
	ADD_MN_CASE( repcols, const_fix, 1, 6 )
	ADD_MN_CASE( repcols, const_fix, 4, 6 )
END_TPACK


BEGIN_TPACK( reprows_generic )
	ADD_N_CASE( reprows, generic, 0 )
	ADD_N_CASE( reprows, generic, 1 )
	ADD_N_CASE( reprows, generic, 6 )
END_TPACK

BEGIN_TPACK( reprows_generic_fix )
	ADD_MN_CASE( reprows, generic_fix, 1, 0 )
	ADD_MN_CASE( reprows, generic_fix, 1, 1 )
	ADD_MN_CASE( reprows, generic_fix, 1, 6 )
	ADD_MN_CASE( reprows, generic_fix, 4, 0 )
	ADD_MN_CASE( reprows, generic_fix, 4, 1 )
	ADD_MN_CASE( reprows, generic_fix, 4, 6 )
END_TPACK

BEGIN_TPACK( reprows_const )
	ADD_N_CASE( reprows, const, 0 )
	ADD_N_CASE( reprows, const, 1 )
	ADD_N_CASE( reprows, const, 6 )
END_TPACK

BEGIN_TPACK( reprows_const_fix )
	ADD_MN_CASE( reprows, const_fix, 1, 0 )
	ADD_MN_CASE( reprows, const_fix, 1, 1 )
	ADD_MN_CASE( reprows, const_fix, 1, 6 )
	ADD_MN_CASE( reprows, const_fix, 4, 0 )
	ADD_MN_CASE( reprows, const_fix, 4, 1 )
	ADD_MN_CASE( reprows, const_fix, 4, 6 )
END_TPACK



BEGIN_TPACK( repcols_eval )
	ADD_MN_CASE_3X3( repcols, eval, 4, 6 )
END_TPACK

BEGIN_TPACK( reprows_eval )
	ADD_MN_CASE_3X3( reprows, eval, 4, 6 )
END_TPACK

BEGIN_TPACK( repcols_vec_eval )
	ADD_MN_CASE_3X3( repcols, vec_eval, 4, 6 )
END_TPACK

BEGIN_TPACK( reprows_vec_eval )
	ADD_MN_CASE_3X3( reprows, vec_eval, 4, 6 )
END_TPACK



BEGIN_MAIN_SUITE
	ADD_TPACK( repcols_generic )
	ADD_TPACK( repcols_generic_fix )
	ADD_TPACK( repcols_const )
	ADD_TPACK( repcols_const_fix )

	ADD_TPACK( reprows_generic )
	ADD_TPACK( reprows_generic_fix )
	ADD_TPACK( reprows_const )
	ADD_TPACK( reprows_const_fix )

	ADD_TPACK( repcols_eval )
	ADD_TPACK( reprows_eval )
	ADD_TPACK( repcols_vec_eval )
	ADD_TPACK( reprows_vec_eval )
END_MAIN_SUITE





