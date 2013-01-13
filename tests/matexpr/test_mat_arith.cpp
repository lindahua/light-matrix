/**
 * @file test_mat_arith.cpp
 *
 * Unit testing for Basic matrix arithmetics
 *
 * @author Dahua Lin
 */

#include "matfun_test_base.h"

#include <light_mat/matexpr/mat_arith.h>

using namespace lmat;
using namespace lmat::test;

typedef dense_matrix<double> dmat_t;


MN_CASE( mat_add )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// check policy

	check_policy(A + B);
	check_policy(A + c);
	check_policy(c + B);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] + B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] + c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c + B[i];

	// default evaluation

	mat_t AB = A + B;
	ASSERT_EQ( AB.nrows(), m );
	ASSERT_EQ( AB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AB, AB_r);

	mat_t AC = A + c;
	ASSERT_EQ( AC.nrows(), m );
	ASSERT_EQ( AC.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AC, AC_r );

	mat_t CB = c + B;
	ASSERT_EQ( CB.nrows(), m );
	ASSERT_EQ( CB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, CB, CB_r );

}


MN_CASE( mat_add_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] + B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] + c;

	// default evaluation

	mat_t AB(A);
	AB += B;
	ASSERT_MAT_EQ( m, n, AB, AB_r );

	mat_t AC(A);
	AC += c;
	ASSERT_MAT_EQ( m, n, AC, AC_r );
}



MN_CASE( mat_sub )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// check policy

	check_policy(A - B);
	check_policy(A - c);
	check_policy(c - B);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] - B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] - c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c - B[i];

	// default evaluation

	mat_t AB = A - B;
	ASSERT_EQ( AB.nrows(), m );
	ASSERT_EQ( AB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AB, AB_r);

	mat_t AC = A - c;
	ASSERT_EQ( AC.nrows(), m );
	ASSERT_EQ( AC.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AC, AC_r );

	mat_t CB = c - B;
	ASSERT_EQ( CB.nrows(), m );
	ASSERT_EQ( CB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, CB, CB_r );

}

MN_CASE( mat_sub_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] - B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] - c;

	// default evaluation

	mat_t AB(A);
	AB -= B;
	ASSERT_MAT_EQ( m, n, AB, AB_r );

	mat_t AC(A);
	AC -= c;
	ASSERT_MAT_EQ( m, n, AC, AC_r );
}


MN_CASE( mat_mul )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// check policy

	check_policy(A * B);
	check_policy(A * c);
	check_policy(c * B);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] * B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] * c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c * B[i];

	// default evaluation

	mat_t AB = A * B;
	ASSERT_EQ( AB.nrows(), m );
	ASSERT_EQ( AB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AB, AB_r);

	mat_t AC = A * c;
	ASSERT_EQ( AC.nrows(), m );
	ASSERT_EQ( AC.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AC, AC_r );

	mat_t CB = c * B;
	ASSERT_EQ( CB.nrows(), m );
	ASSERT_EQ( CB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, CB, CB_r );
}

MN_CASE( mat_mul_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] * B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] * c;

	// default evaluation

	mat_t AB(A);
	AB *= B;
	ASSERT_MAT_EQ( m, n, AB, AB_r );

	mat_t AC(A);
	AC *= c;
	ASSERT_MAT_EQ( m, n, AC, AC_r );
}


MN_CASE( mat_div )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 4.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(1 << (i % 5));

	// check policy

	check_policy(A / B);
	check_policy(A / c);
	check_policy(c / B);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] / B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] / c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c / B[i];

	// default evaluation

	mat_t AB = A / B;
	ASSERT_EQ( AB.nrows(), m );
	ASSERT_EQ( AB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AB, AB_r);

	mat_t AC = A / c;
	ASSERT_EQ( AC.nrows(), m );
	ASSERT_EQ( AC.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, AC, AC_r );

	mat_t CB = c / B;
	ASSERT_EQ( CB.nrows(), m );
	ASSERT_EQ( CB.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, CB, CB_r );
}


MN_CASE( mat_div_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 4.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(1 << (i % 5));

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] / B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] / c;

	// default evaluation

	mat_t AB(A);
	AB /= B;
	ASSERT_MAT_EQ( m, n, AB, AB_r );

	mat_t AC(A);
	AC /= c;
	ASSERT_MAT_EQ( m, n, AC, AC_r );
}


MN_CASE( mat_neg )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// check policy

	check_policy(-A);

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = -A[i];

	// default evaluation

	mat_t R = -A;
	ASSERT_EQ( R.nrows(), m );
	ASSERT_EQ( R.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, R, R_r  );
}


// test packs for arithmetic operators

AUTO_TPACK( mat_add )
{
	ADD_MN_CASE_3X3( mat_add, DM, DN )
}

AUTO_TPACK( mat_add_ip )
{
	ADD_MN_CASE_3X3( mat_add_ip, DM, DN )
}

AUTO_TPACK( mat_sub )
{
	ADD_MN_CASE_3X3( mat_sub, DM, DN )
}

AUTO_TPACK( mat_sub_ip )
{
	ADD_MN_CASE_3X3( mat_sub_ip, DM, DN )
}

AUTO_TPACK( mat_mul )
{
	ADD_MN_CASE_3X3( mat_mul, DM, DN )
}

AUTO_TPACK( mat_mul_ip )
{
	ADD_MN_CASE_3X3( mat_mul_ip, DM, DN )
}

AUTO_TPACK( mat_div )
{
	ADD_MN_CASE_3X3( mat_div, DM, DN )
}

AUTO_TPACK( mat_div_ip )
{
	ADD_MN_CASE_3X3( mat_div_ip, DM, DN )
}

AUTO_TPACK( mat_neg )
{
	ADD_MN_CASE_3X3( mat_neg, DM, DN )
}


// test packs for other matrix functions

DEFINE_MATFUN_TESTS3a( fma, -2.0, 2.0, 1.0e-15 )

DEFINE_MATFUN_TESTS2( min, -2.0, 2.0, -2.0, 2.0, 1.0e-16 )
DEFINE_MATFUN_TESTS2( max, -2.0, 2.0, -2.0, 2.0, 1.0e-16 )
DEFINE_MATFUN_TESTS3a( clamp, -2.0, 2.0, 1.0e-15 )

// simple powers

DEFINE_MATFUN_TESTS1( abs, -2.0, 2.0, 1.0e-16 )
DEFINE_MATFUN_TESTS1( sqr, -2.0, 2.0, 1.0e-16 )
DEFINE_MATFUN_TESTS1( cube, -2.0, 2.0, 1.0e-15 )

DEFINE_MATFUN_TESTS1( rcp, -2.0, 2.0, 1.0e-15 )
DEFINE_MATFUN_TESTS1( sqrt, -2.0, 2.0, 1.0e-15 )
DEFINE_MATFUN_TESTS1( rsqrt, -2.0, 2.0, 1.0e-15 )

// rounding

DEFINE_MATFUN_TESTS1( floor, -10.0, 10.0, 1.0e-16 )
DEFINE_MATFUN_TESTS1( ceil,  -10.0, 10.0, 1.0e-16 )
DEFINE_MATFUN_TESTS1( round, -10.0, 10.0, 1.0e-16 )
DEFINE_MATFUN_TESTS1( trunc, -10.0, 10.0, 1.0e-16 )

