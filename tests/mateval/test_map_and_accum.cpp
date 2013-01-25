/**
 * @file test_map_and_accum.cpp
 *
 * @brief Unit testing of ewise map and accumulation
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#define DEFAULT_M_VALUE 13
#define DEFAULT_N_VALUE 9

#include "../multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/math/basic_functors.h>
#include <light_mat/mateval/ewise_eval.h>


using namespace lmat;
using namespace lmat::test;


template<int M, int N>
void test_ewise_map()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> sa(m, n);
	dense_matrix<double, M, N> sb(m, n);
	dense_matrix<double, M, N> dst(m, n, zero());
	dense_matrix<double, M, N> r1(m, n);
	dense_matrix<double, M, N> r2(m, n);

	for (index_t i = 0; i < m * n; ++i)
	{
		sa[i] = double(i + 2);
		sb[i] = double(i * i + 1);

		r1[i] = math::sqr(sa[i]);
		r2[i] = sa[i] + sb[i];
	}

	matrix_shape<M, N> shape(m, n);

	map(sqr_fun<double>())(shape, out_(dst), in_(sa));
	ASSERT_MAT_EQ( m, n, dst, r1 );

	map(add_fun<double>())(shape, out_(dst), in_(sa), in_(sb));
	ASSERT_MAT_EQ( m, n, dst, r2 );
}


MN_CASE( ewise_map )
{
	test_ewise_map<M, N>();
}



MN_CASE( ewise_map_to )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> s1(m, n);
	dense_matrix<double, M, N> s2(m, n);

	do_fill_rand( s1.ptr_data(), m * n );
	do_fill_rand( s2.ptr_data(), m * n );

	matrix_shape<M, N> shape(m, n);
	dense_matrix<double, M, N> r(m, n);

	for (index_t i = 0; i < m * n; ++i) r[i] = s1[i] + s2[i];

	map_to(a, add_fun<double>(), in_(s1), in_(s2));

	ASSERT_MAT_EQ(m, n, a, r);
}


MN_CASE( ewise_accum_to )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> s(m, n);
	dense_matrix<double, M, N> c(m, n);

	do_fill_rand( a.ptr_data(), m * n );
	do_fill_rand( s.ptr_data(), m * n );
	do_fill_rand( c.ptr_data(), m * n );

	double cv = c[0] * 2.0;

	matrix_shape<M, N> shape(m, n);
	dense_matrix<double, M, N> r(m, n);

	for (index_t i = 0; i < m * n; ++i) r[i] = a[i] + s[i];
	accum_to(a, s);

	ASSERT_MAT_EQ( m, n, a, r );

	for (index_t i = 0; i < m * n; ++i) r[i] = a[i] + cv * s[i];
	accum_to(a, cv, s);

	ASSERT_MAT_EQ( m, n, a, r );

	for (index_t i = 0; i < m * n; ++i) r[i] = a[i] + c[i] * s[i];
	accum_to(a, c, s);

	ASSERT_MAT_EQ(m, n, a, r );
}



AUTO_TPACK( ewise_map )
{
	ADD_MN_CASE_3X3( ewise_map, DM, DN )
}


AUTO_TPACK( ewise_map_to )
{
	ADD_MN_CASE_3X3( ewise_map_to, DM, DN )
}

AUTO_TPACK( ewise_accum_to )
{
	ADD_MN_CASE_3X3( ewise_accum_to, DM, DN )
}


