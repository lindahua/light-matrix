/**
 * @file test_mat_transpose.cpp
 *
 * Unit testing for matrix transpose
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include <light_mat/matrix/matrix_eval.h>
#include <string>

using namespace lmat;
using namespace lmat::test;

const int DM = 6;
const int DN = 7;
const int LDim = 9;

template<class Mat>
void fill_lin(Mat& A)
{
	const index_t m = A.nrows();
	const index_t n = A.ncolumns();

	int v = 1;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			A(i, j) = double(v++);
		}
	}
}

template<class SMat, class DMat>
void my_transpose(const SMat& S, DMat& D)
{
	const index_t m = S.nrows();
	const index_t n = S.ncolumns();

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			D(j, i) = S(i, j);
		}
	}
}


MN_CASE( mat_trans, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	std::string base_name =
			N == 1 ? "contcol" : (M == 1 ? "controw" : "dense");

	dense_matrix<double, M, N> S(m, n);
	fill_lin(S);

	ASSERT_STREQ( S.trans().trans_base_type_name(), base_name  );

	dense_matrix<double, N, M> T0(n, m);
	my_transpose(S, T0);

	dense_matrix<double, N, M> T = S.trans();

	ASSERT_EQ( T.nrows(), n );
	ASSERT_EQ( T.ncolumns(), m );

	ASSERT_MAT_EQ( n, m, T, T0 );

	dense_matrix<double, M, N> R = T.trans();

	ASSERT_EQ( R.nrows(), m );
	ASSERT_EQ( R.ncolumns(), n );

	ASSERT_MAT_EQ( m, n, R, S );
}


MN_CASE( mat_trans, refex )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	const char *base_name =
			N == 1 ? "contcol" : (M == 1 ? "regular_row" : "dense");

	scoped_array<double> sarr(LDim * n, -1.0);

	ref_matrix_ex<double, M, N> S(sarr.ptr_begin(), m, n, LDim);
	fill_lin(S);

	ASSERT_STREQ( S.trans().trans_base_type_name(), base_name  );

	dense_matrix<double, N, M> T0(n, m);
	my_transpose(S, T0);

	dense_matrix<double, N, M> T = S.trans();

	ASSERT_EQ( T.nrows(), n );
	ASSERT_EQ( T.ncolumns(), m );

	ASSERT_MAT_EQ( n, m, T, T0 );
}


MN_CASE( mat_trans, unary_ewise )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	std::string base_name =
			M == 1 ? "rowxpr" : (N == 1 ? "colxpr" : "generic");

	dense_matrix<double, M, N> S(m, n);
	fill_lin(S);

	ASSERT_STREQ( sqr(S).trans().trans_base_type_name(), base_name );

	dense_matrix<double, M, N> R = sqr(S);
	dense_matrix<double, N, M> T0(n, m);
	my_transpose(R, T0);

	dense_matrix<double, N, M> T = sqr(S).trans();

	ASSERT_EQ( T.nrows(), n );
	ASSERT_EQ( T.ncolumns(), m );

	ASSERT_MAT_EQ( n, m, T, T0 );

	scoped_array<double> darr(LDim * m, -1.0);
	ref_matrix_ex<double, N, M> T2(darr.ptr_begin(), n, m, LDim);

	T2 = sqr(S).trans();

	ASSERT_MAT_EQ( n, m, T2, T0 );
}


MN_CASE( mat_trans, binary_ewise )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	std::string base_name =
			M == 1 ? "rowxpr" : (N == 1 ? "colxpr" : "generic");

	dense_matrix<double, M, N> S(m, n);
	fill_lin(S);
	dense_matrix<double, M, N> S2 = sqr(S);

	ASSERT_STREQ( (S + S2).trans().trans_base_type_name(), base_name );

	dense_matrix<double, M, N> R = S + S2;
	dense_matrix<double, N, M> T0(n, m);
	my_transpose(R, T0);

	dense_matrix<double, N, M> T = (S + S2).trans();

	ASSERT_EQ( T.nrows(), n );
	ASSERT_EQ( T.ncolumns(), m );

	ASSERT_MAT_EQ( n, m, T, T0 );

	scoped_array<double> darr(LDim * m, -1.0);
	ref_matrix_ex<double, N, M> T2(darr.ptr_begin(), n, m, LDim);

	T2 = (S + S2).trans();

	ASSERT_MAT_EQ( n, m, T2, T0 );
}


MN_CASE( mat_trans, colwise_reduce )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	std::string base_name = "rowxpr";

	dense_matrix<double, M, N> S(m, n);
	fill_lin(S);
	ASSERT_STREQ( sum(S, colwise()).trans().trans_base_type_name(), base_name );

	dense_matrix<double, 1, N> R = sum(S, colwise());
	dense_matrix<double, N, 1> T0(n, 1);
	my_transpose(R, T0);

	dense_matrix<double, N, 1> T = sum(S, colwise()).trans();

	ASSERT_EQ( T.nrows(), n );
	ASSERT_EQ( T.ncolumns(), 1 );

	ASSERT_MAT_EQ( n, 1, T, T0 );
}


MN_CASE( mat_trans, rowwise_reduce )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	std::string base_name = M == 1 ? "rowxpr" : "colxpr";

	dense_matrix<double, M, N> S(m, n);
	fill_lin(S);

	ASSERT_STREQ( sum(S, rowwise()).trans().trans_base_type_name(), base_name );

	dense_matrix<double, M, 1> R = sum(S, rowwise());
	dense_matrix<double, 1, M> T0(1, m);
	my_transpose(R, T0);

	dense_matrix<double, 1, M> T = sum(S, rowwise()).trans();

	ASSERT_EQ( T.nrows(), 1 );
	ASSERT_EQ( T.ncolumns(), m );

	ASSERT_MAT_EQ( 1, m, T, T0 );
}




BEGIN_TPACK( dense_trans )
	ADD_MN_CASE_3X3( mat_trans, dense, DM, DN )
END_TPACK

BEGIN_TPACK( refex_trans )
	ADD_MN_CASE_3X3( mat_trans, refex, DM, DN )
END_TPACK

BEGIN_TPACK( unary_ewise_trans )
	ADD_MN_CASE_3X3( mat_trans, unary_ewise, DM, DN )
END_TPACK

BEGIN_TPACK( binary_ewise_trans )
	ADD_MN_CASE_3X3( mat_trans, binary_ewise, DM, DN )
END_TPACK


BEGIN_TPACK( colwise_reduce_trans )
	ADD_MN_CASE_3X3( mat_trans, colwise_reduce, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_reduce_trans )
	ADD_MN_CASE_3X3( mat_trans, rowwise_reduce, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( dense_trans )
	ADD_TPACK( refex_trans )
	ADD_TPACK( unary_ewise_trans )
	ADD_TPACK( binary_ewise_trans )
	ADD_TPACK( colwise_reduce_trans )
	ADD_TPACK( rowwise_reduce_trans )
END_MAIN_SUITE

