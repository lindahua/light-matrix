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
#include <light_mat/matrix/repeat_vecs.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

typedef dense_matrix<double, 0, 1> dcol_t;
typedef dense_matrix<double, 1, 0> drow_t;

template class lmat::repeat_col_expr<dcol_t, DynamicDim>;
template class lmat::repeat_col_expr<dcol_t, 6>;
template class lmat::repeat_row_expr<drow_t, DynamicDim>;
template class lmat::repeat_row_expr<drow_t, 4>;


#ifdef LMAT_USE_STATIC_ASSERT
static_assert(lmat::is_mat_xpr<lmat::repeat_col_expr<dcol_t> >::value, "Interface verification failed.");
static_assert(lmat::is_mat_xpr<lmat::repeat_row_expr<drow_t> >::value, "Interface verification failed.");
#endif


N_CASE( repcols, generic )
{
	const index_t m = N == 0 ? 4 : N;
	const index_t n = 6;

	typedef dense_matrix<double, N, 1> col_t;
	typedef repeat_col_expr<col_t> expr_t;

	col_t a(m, 1);

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert( is_same<decltype(repcol(a, n)), expr_t>::value, "Expression type verification failed" );
	static_assert( ct_rows<expr_t>::value == N, "CT-size verification failed." );
	static_assert( ct_cols<expr_t>::value == 0, "CT-size verification failed." );
#endif

	expr_t expr = repcol(a, n);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( &(expr.column()), &a );
}


N_CASE( reprows, generic )
{
	const index_t m = 4;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, 1, N> row_t;
	typedef repeat_row_expr<row_t> expr_t;

	row_t a(1, n);

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert( is_same<decltype(reprow(a, m)), expr_t>::value, "Expression type verification failed" );
	static_assert( ct_rows<expr_t>::value == 0, "CT-size verification failed." );
	static_assert( ct_cols<expr_t>::value == N, "CT-size verification failed." );
#endif

	expr_t expr = reprow(a, m);

	ASSERT_EQ( expr.nrows(), m );
	ASSERT_EQ( expr.ncolumns(), n );
	ASSERT_EQ( expr.nelems(), m * n);
	ASSERT_EQ( expr.size(), size_t(m * n) );

	ASSERT_EQ( &(expr.row()), &a );
}


BEGIN_TPACK( repcols_generic )
	ADD_N_CASE( repcols, generic, 0 )
	ADD_N_CASE( repcols, generic, 1 )
	ADD_N_CASE( repcols, generic, 4 )
END_TPACK

BEGIN_TPACK( reprows_generic )
	ADD_N_CASE( reprows, generic, 0 )
	ADD_N_CASE( reprows, generic, 1 )
	ADD_N_CASE( reprows, generic, 6 )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( repcols_generic )
	ADD_TPACK( reprows_generic )
END_MAIN_SUITE





