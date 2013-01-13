/**
 * @file test_mat_zip.cpp
 *
 * @brief Unit testing for matrix zip/unzip
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/mat_zip.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 9;
const index_t DN = 8;


MN_CASE( mat_zip_pair_aa )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<int32_t, M, N> b(m, n);

	for (index_t i = 0; i < m * n; ++i)
	{
		a[i] = double(i+2);
		b[i] = 2 * i + 3;
	}

	typedef std::pair<double, int32_t> T;

	dense_matrix<T, M, N> t = zip_pair(a, b);
	ASSERT_EQ( t.nrows(), m );
	ASSERT_EQ( t.ncolumns(), n );

	dense_matrix<T, M, N> r(m, n);
	for (index_t i = 0; i < m * n; ++i)
		r[i] = std::make_pair(a[i], b[i]);

	ASSERT_MAT_EQ( m, n, t, r );

	dense_matrix<double, M, N> e0  = zip_e(t, int_<0>());
	dense_matrix<int32_t, M, N> e1 = zip_e(t, int_<1>());

	ASSERT_EQ( e0.nrows(), m );
	ASSERT_EQ( e0.ncolumns(), n );
	ASSERT_EQ( e1.nrows(), m );
	ASSERT_EQ( e1.ncolumns(), n );

	ASSERT_MAT_EQ( m, n, e0, a );
	ASSERT_MAT_EQ( m, n, e1, b );

	dense_matrix<double, M, N> u0;
	dense_matrix<double, M, N> u1;
	unzip(t, u0, u1);

	ASSERT_EQ( u0.nrows(), m );
	ASSERT_EQ( u0.ncolumns(), n );
	ASSERT_EQ( u1.nrows(), m );
	ASSERT_EQ( u1.ncolumns(), n );

	ASSERT_MAT_EQ( m, n, u0, a );
	ASSERT_MAT_EQ( m, n, u1, b );
}

MN_CASE( mat_zip_pair_av )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	int32_t bv = 13;

	for (index_t i = 0; i < m * n; ++i)
	{
		a[i] = double(i+2);
	}

	typedef std::pair<double, int32_t> T;

	dense_matrix<T, M, N> t = zip_pair_fix2(a, bv);
	ASSERT_EQ( t.nrows(), m );
	ASSERT_EQ( t.ncolumns(), n );

	dense_matrix<T, M, N> r(m, n);
	for (index_t i = 0; i < m * n; ++i)
		r[i] = std::make_pair(a[i], bv);

	ASSERT_MAT_EQ( m, n, t, r );
}


MN_CASE( mat_zip_pair_va )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	double av = 12.5;
	dense_matrix<int32_t, M, N> b(m, n);

	for (index_t i = 0; i < m * n; ++i)
	{
		b[i] = 2 * i + 3;
	}

	typedef std::pair<double, int32_t> T;

	dense_matrix<T, M, N> t = zip_pair_fix1(av, b);
	ASSERT_EQ( t.nrows(), m );
	ASSERT_EQ( t.ncolumns(), n );

	dense_matrix<T, M, N> r(m, n);
	for (index_t i = 0; i < m * n; ++i)
		r[i] = std::make_pair(av, b[i]);

	ASSERT_MAT_EQ( m, n, t, r );
}


AUTO_TPACK( mat_zip_pair_aa )
{
	ADD_MN_CASE_3X3( mat_zip_pair_aa, DM, DN )
}

AUTO_TPACK( mat_zip_pair_av )
{
	ADD_MN_CASE_3X3( mat_zip_pair_av, DM, DN )
}

AUTO_TPACK( mat_zip_pair_va )
{
	ADD_MN_CASE_3X3( mat_zip_pair_va, DM, DN )
}




