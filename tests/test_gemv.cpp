/*
 * @file test_gemv.cpp
 *
 * Test of generic matrix-vector multiplication
 *
 * @author Dahua Lin
 */


#include "linalg_test_base.h"

#include <light_mat/matrix/matrix_arith.h>
#include <light_mat/matrix/matrix_ewise_eval.h>

#include <light_mat/linalg/matrix_blasl2.h>

using namespace lmat;
using namespace lmat::test;


template<typename T> struct tolvalue;

template<> struct tolvalue<float>
{
	static float get() { return 2.0e-5f; }
};


template<> struct tolvalue<double>
{
	static double get() { return 1.0e-12; }
};


const int DM = 6;
const int DN = 7;
const index_t LDim = 8;


TMN_CASE( gemv_n, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> A(m, n);
	dense_matrix<T, M, 1> x(m, 1);
	dense_matrix<T, N, 1> y(n, 1, zeros<T>());
	dense_matrix<T, N, 1> r(n, 1, zeros<T>());

	fill_rand(A, T(-1), T(1));
	fill_rand(x, T(-1), T(1));

	gemv(A, x, y);

	naive_mtimes(A, x, r);

	T tol = tolvalue<T>::get();
	ASSERT_MAT_APPROX(n, 1, y, r, tol);
}


BEGIN_TPACK( gemv_n_dense )
	ADD_TMN_CASE( gemv_n, dense, double, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( gemv_n_dense )
END_MAIN_SUITE
