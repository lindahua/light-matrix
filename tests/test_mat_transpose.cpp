/**
 * @file test_mat_transpose.cpp
 *
 * Unit testing for matrix transpose
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include <light_mat/matrix/matrix_classes.h>

using namespace lmat;
using namespace lmat::test;

const int DM = 6;
const int DN = 7;


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

	dense_matrix<double, M, N> S(m, n);
	fill_lin(S);

	dense_matrix<double, N, M> T0(n, m);
	my_transpose(S, T0);

	// dense_matrix<double, N, M> T = S.trans();

	// ASSERT_EQ( T.nrows(), n );
	// ASSERT_EQ( T.ncolumns(), m );

	// ASSERT_MAT_EQ( n, m, T, T0 );
}


BEGIN_TPACK( dense_trans )
	ADD_MN_CASE( mat_trans, dense, 0, 0 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( dense_trans )
END_MAIN_SUITE

