/**
 * @file test_cpd_ewise.cpp
 *
 * Unit testing of compound ewise expressions
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/matexpr/mat_arith.h>
#include <light_mat/matexpr/mat_emath.h>
#include <light_mat/matexpr/mat_cast.h>
#include <light_mat/matexpr/repvec_expr.h>

#include <light_mat/math/math_constants.h>
#include <cstdlib>

using namespace lmat;
using namespace lmat::test;


// auxiliary functions

template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& X, double a, double b)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = a + (double(std::rand()) / RAND_MAX) * (b - a);
	}
}

// test cases

MN_CASE( ewise_sqdist )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_matrix<double, M, N> A(m, n);
	dense_matrix<double, M, N> B(m, n);

	fill_ran(A, -5.0, 5.0);
	fill_ran(B, -5.0, 5.0);
	const double tol = 1.0e-12;

	dense_matrix<double, M, N> Y = sqr(A - B);

	dense_matrix<double, M, N> R(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			R(i, j) = math::sqr(A(i,j) - B(i,j));
		}
	}

	ASSERT_EQ( Y.nrows(), m );
	ASSERT_EQ( Y.ncolumns(), n );
	ASSERT_MAT_APPROX( m, n, Y, R, tol );
}


MN_CASE( ewise_normpdf )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_matrix<double, M, N> X(m, n);

	fill_ran(X, -4.0, 6.0);
	const double tol = 1.0e-12;

	const double mu = 1.2;
	const double sigma2 = 2.5;

	double hrpi = lmat::math::consts<double>::rcp_two_rpi();
	const double c = math::sqrt(hrpi / sigma2);

	dense_matrix<double, M, N> Y = c * exp( - sqr(X - mu) / (2.0 * sigma2) );

	dense_matrix<double, M, N> R(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			R(i,j) = c * math::exp( - math::sqr(X(i,j) - mu) / (2.0 * sigma2) );
		}
	}

	ASSERT_EQ( Y.nrows(), m );
	ASSERT_EQ( Y.ncolumns(), n );
	ASSERT_MAT_APPROX( m, n, Y, R, tol );
}


MN_CASE( axpy_cast )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	const double a = 1.5;
	dense_matrix<double, M, N> X(m, n);
	dense_matrix<double, M, N> Y(m, n);

	fill_ran(X, -5.0, 5.0);
	fill_ran(Y, -5.0, 5.0);
	const float tol = 1.0e-5f;

	dense_matrix<float, M, N> Z = to_f32(a * X + Y);
	dense_matrix<float, M, N> R(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			R(i, j) = float(a * X(i,j) + Y(i,j));
		}
	}

	ASSERT_EQ( Z.nrows(), m );
	ASSERT_EQ( Z.ncolumns(), n );
	ASSERT_MAT_APPROX( m, n, Z, R, tol );
}

MN_CASE( pwise_sqdist )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_col<double, M> vx(m);
	dense_row<double, N> vy(n);

	fill_ran(vx, -5.0, 5.0);
	fill_ran(vy, -5.0, 5.0);
	const double tol = 1.0e-12;

	dense_matrix<double, M, N> Y =
			  repcol(eval(sqr(vx)), n)
			+ reprow(eval(sqr(vy)), m)
			- 2.0 * repcol(vx, n) * reprow(vy, m);

	dense_matrix<double, M, N> R(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			R(i,j) = math::sqr(vx[i] - vy[j]);
		}
	}

	ASSERT_EQ( Y.nrows(), m );
	ASSERT_EQ( Y.ncolumns(), n );
	ASSERT_MAT_APPROX( m, n, Y, R, tol );
}


LTEST_INIT_AUTOSUITE

AUTO_TPACK( ewise_sqdist )
{
	ADD_MN_CASE_3X3( ewise_sqdist, 5, 6 )
}

AUTO_TPACK( ewise_normpdf )
{
	ADD_MN_CASE_3X3( ewise_normpdf, 5, 6 )
}

AUTO_TPACK( axpy_cast )
{
	ADD_MN_CASE_3X3( axpy_cast, 5, 6 )
}

AUTO_TPACK( pwise_sqdist )
{
	ADD_MN_CASE_3X3( pwise_sqdist, 5, 6 )
}



