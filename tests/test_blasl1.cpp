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


BEGIN_MAIN_SUITE
	ADD_TPACK( asum_dense )
	ADD_TPACK( asum_refex )
	ADD_TPACK( asum_mexpr )
END_MAIN_SUITE






