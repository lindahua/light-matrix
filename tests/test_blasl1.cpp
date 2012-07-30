/**
 * @file test_blasl1.cpp
 *
 * Test BLAS Level 1
 *
 * @author Dahua Lin
 */


#include "linalg_test_base.h"

#include <light_mat/matrix/matrix_arith.h>
#include <light_mat/matrix/matrix_ewise_eval.h>

#include <light_mat/linalg/matrix_blasl1.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;


template<typename T> struct tolvalue;

template<> struct tolvalue<float>
{
	static float get() { return 1.0e-5f; }
};


template<> struct tolvalue<double>
{
	static double get() { return 1.0e-12; }
};


const int DM = 6;
const int DN = 7;
const index_t LDim = 8;


// asum

TMN_CASE( asum, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> X(m, n);
	fill_rand(X, T(-1), T(1));

	T tol = tolvalue<T>::get();

	double s = 0.0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			s += std::abs(X(i,j));
		}
	}

	T r0 = T(s);
	T r = lmat::blas::asum(X);

	ASSERT_APPROX( r, r0, tol );
}


TMN_CASE( asum, refex )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	scoped_array<T> Xbuf(LDim * n, 0.0);
	ref_matrix_ex<T, M, N> X(Xbuf.ptr_begin(), m, n, LDim);

	fill_rand(X, T(-1), T(1));

	T tol = tolvalue<T>::get();

	double s = 0.0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			s += std::abs(X(i,j));
		}
	}

	T r0 = T(s);
	T r = lmat::blas::asum(X);

	ASSERT_APPROX( r, r0, tol );
}


TMN_CASE( asum, mexpr )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> A(m, n);
	dense_matrix<T, M, N> B(m, n);

	fill_rand(A, T(-1), T(1));
	fill_rand(B, T(-1), T(1));
	dense_matrix<T, M, N> X = A + B;

	T tol = tolvalue<T>::get();

	double s = 0.0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			s += std::abs(X(i,j));
		}
	}

	T r0 = T(s);
	T r = lmat::blas::asum(A + B);

	ASSERT_APPROX( r, r0, tol );
}


// axpy

TMN_CASE( axpy, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> X(m, n);
	fill_rand(X, T(-1), T(1));

	dense_matrix<T, M, N> Y(m, n);
	fill_rand(Y, T(-1), T(-1));

	dense_matrix<T, M, N> R0(m, n);

	T a = T(1.8);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			R0(i,j) = a * X(i,j) + Y(i,j);
		}
	}

	lmat::blas::axpy(a, X, Y);

	T tol = tolvalue<T>::get();
	ASSERT_MAT_APPROX(m, n, Y, R0, tol);
}


TMN_CASE( axpy, refex )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	scoped_array<T> Xbuf(LDim * n);
	ref_matrix_ex<T, M, N> X(Xbuf.ptr_begin(), m, n, LDim);
	fill_rand(X, T(-1), T(1));

	dense_matrix<T, M, N> Y(m, n);
	fill_rand(Y, T(-1), T(-1));

	dense_matrix<T, M, N> R0(m, n);

	T a = T(1.8);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			R0(i,j) = a * X(i,j) + Y(i,j);
		}
	}

	lmat::blas::axpy(a, X, Y);

	T tol = tolvalue<T>::get();
	ASSERT_MAT_APPROX(m, n, Y, R0, tol);
}


TMN_CASE( axpy, mexpr )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> A(m, n);
	dense_matrix<T, M, N> B(m, n);

	fill_rand(A, T(-1), T(1));
	fill_rand(B, T(-1), T(1));

	dense_matrix<T, M, N> X = A + B;

	dense_matrix<T, M, N> Y(m, n);
	fill_rand(Y, T(-1), T(-1));

	dense_matrix<T, M, N> R0(m, n);

	T a = T(1.8);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			R0(i,j) = a * X(i,j) + Y(i,j);
		}
	}

	lmat::blas::axpy(a, A + B, Y);

	T tol = tolvalue<T>::get();
	ASSERT_MAT_APPROX(m, n, Y, R0, tol);
}


// dot

TMN_CASE( dot, dense_ex )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> X(m, n);

	scoped_array<T> Ybuf(LDim * n);
	ref_matrix_ex<T, M, N> Y(Ybuf.ptr_begin(), m, n, LDim);

	fill_rand(X, T(-1), T(1));
	fill_rand(Y, T(-1), T(1));

	T tol = tolvalue<T>::get();

	double s = 0.0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			s += X(i,j) * Y(i,j);
		}
	}

	T r0 = T(s);
	T r = lmat::blas::dot(X, Y);

	ASSERT_APPROX( r, r0, tol );
}


// nrm2

TMN_CASE( nrm2, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> X(m, n);
	fill_rand(X, T(-1), T(1));

	T tol = tolvalue<T>::get();

	double s = 0.0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			s += X(i,j) * X(i,j);
		}
	}

	T r0 = (T)std::sqrt(s);
	T r = lmat::blas::nrm2(X);

	ASSERT_APPROX( r, r0, tol );
}


// rot

TMN_CASE( rot, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> X(m, n);
	dense_matrix<T, M, N> Y(m, n);

	fill_rand(X, T(-1), T(1));
	fill_rand(Y, T(-1), T(1));

	T tol = tolvalue<T>::get();

	dense_matrix<T, M, N> Xr(m, n);
	dense_matrix<T, M, N> Yr(m, n);

	T c = T(2.5);
	T s = T(3.2);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			Xr(i,j) = c * X(i,j) + s * Y(i,j);
			Yr(i,j) = c * Y(i,j) - s * X(i,j);
		}
	}

	lmat::blas::rot(X, Y, c, s);

	ASSERT_MAT_APPROX( m, n, X, Xr, tol );
	ASSERT_MAT_APPROX( m, n, Y, Yr, tol );
}



BEGIN_TPACK( asum_dense )
	ADD_TMN_CASE_3X3( asum, dense, float,  DM, DN )
	ADD_TMN_CASE_3X3( asum, dense, double, DM, DN )
END_TPACK

BEGIN_TPACK( asum_refex )
	ADD_TMN_CASE_3X3( asum, refex, float,  DM, DN )
	ADD_TMN_CASE_3X3( asum, refex, double, DM, DN )
END_TPACK

BEGIN_TPACK( asum_mexpr )
	ADD_TMN_CASE_3X3( asum, mexpr, float,  DM, DN )
	ADD_TMN_CASE_3X3( asum, mexpr, double, DM, DN )
END_TPACK


BEGIN_TPACK( axpy_dense )
	ADD_TMN_CASE_3X3( axpy, dense, float,  DM, DN )
	ADD_TMN_CASE_3X3( axpy, dense, double, DM, DN )
END_TPACK

BEGIN_TPACK( axpy_refex )
	ADD_TMN_CASE_3X3( axpy, refex, float,  DM, DN )
	ADD_TMN_CASE_3X3( axpy, refex, double, DM, DN )
END_TPACK

BEGIN_TPACK( axpy_mexpr )
	ADD_TMN_CASE_3X3( axpy, mexpr, float,  DM, DN )
	ADD_TMN_CASE_3X3( axpy, mexpr, double, DM, DN )
END_TPACK


BEGIN_TPACK( dot_dense_ex )
	ADD_TMN_CASE_3X3( dot, dense_ex, float,  DM, DN )
	ADD_TMN_CASE_3X3( dot, dense_ex, double, DM, DN )
END_TPACK

BEGIN_TPACK( nrm2_dense )
	ADD_TMN_CASE_3X3( nrm2, dense, float,  DM, DN )
	ADD_TMN_CASE_3X3( nrm2, dense, double, DM, DN )
END_TPACK

BEGIN_TPACK( rot_dense )
	ADD_TMN_CASE_3X3( rot, dense, float,  DM, DN )
	ADD_TMN_CASE_3X3( rot, dense, double, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( asum_dense )
	ADD_TPACK( asum_refex )
	ADD_TPACK( asum_mexpr )

	ADD_TPACK( axpy_dense )
	ADD_TPACK( axpy_refex )
	ADD_TPACK( axpy_mexpr )

	ADD_TPACK( dot_dense_ex )

	ADD_TPACK( nrm2_dense )

	ADD_TPACK( rot_dense )
END_MAIN_SUITE






